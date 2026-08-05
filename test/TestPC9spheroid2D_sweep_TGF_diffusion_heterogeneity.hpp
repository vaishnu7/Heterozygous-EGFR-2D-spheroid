/*
 * TestPC9spheroid2D_sweep_TGF_diffusion_heterogeneity.hpp
 *
 * Created on: 24 Oct, 2025
 * Author: Vaishnudebi Dutta
 *
 */

#ifndef TESTPC9SPHEROID2D_HPP_
#define TESTPC9SPHEROID2D_HPP_

#include <cxxtest/TestSuite.h>                          // Unit testing framework
#include "CellBasedSimulationArchiver.hpp"              // Save/load simulation checkpoints
#include "AbstractCellBasedTestSuite.hpp"               // Base class for cell-based tests
#include "CheckpointArchiveTypes.hpp"                  // Serialization support for checkpoints
#include "SmartPointers.hpp"                            // Shared pointer utilities

// Core simulation classes
#include "HoneycombMeshGenerator.hpp"                   // Generate 2D honeycomb mesh geometry
#include "OffLatticeSimulation.hpp"                     // Off-lattice cell-based simulation engine
#include "MeshBasedCellPopulation.hpp"                  // Container for cells on mesh with interactions
#include "VoronoiDataWriter.hpp"                        // Output Voronoi tessellation data
#include "DifferentiatedCellProliferativeType.hpp"      // Non-dividing cell type marker
#include "TransitCellProliferativeType.hpp"             // Dividing cell type marker
#include "StemCellProliferativeType.hpp"
#include "FakePetscSetup.hpp"                           // Mock parallel solver setup for testing

// Cell cycle and mutation classes
#include "SimpleOxygenBasedCellCycleModel.hpp"          // Oxygen-dependent cell cycle model
#include "WildTypeCellMutationState.hpp"                // Normal/healthy cell state marker
#include "ApcOneHitCellMutationState.hpp"               // Hypoxic cell state marker
#include "ApcTwoHitCellMutationState.hpp"               // Necrotic/dead cell state marker

// Force classes
#include "GeneralisedLinearSpringForce.hpp"             // Spring forces between neighboring cells (Meineke model)

// Contact inhibition
#include "SimpleTargetAreaModifier.hpp"                 // Cell growth limiting by target area

// Utility classes
#include <boost/lexical_cast.hpp>                       // String conversion utility
#include "CellAgesWriter.hpp"                           // Output cell age data
#include "CellLabelWriter.hpp"                          // Output cell labels
#include "CellMutationStatesWriter.hpp"                 // Output mutation state data
#include "CellAncestorWriter.hpp"                       // Output cell lineage/ancestry
#include "OutputFileHandler.hpp"                        // File I/O management
#include "CellDataWriter.hpp"                             // Custom cell data output
#include "CellData.hpp"                                 // Generic cell property storage
#include <iomanip>                                      // Output formatting
#include <iostream>                                     // Console I/O
#include <cmath>                                        // Math functions

// Custom header files
#include "HypoxiaSignalingModifier.hpp"                 // HIF-1a and TGF-a signaling under hypoxia
#include "EGFRSignalingModifier.hpp"                    // EGFR-TGF-a signaling pathway
#include "SpheroidOxygenGradientModifier.hpp"           // Radial oxygen gradient in spheroid
#include "HypoxiaCellStateModifier.hpp"                 // Hypoxic cell state changes
#include "HIF1AlphaWriter.hpp"                          // Output HIF-1a data
#include "TGFAlphaWriter.hpp"                           // Output TGF-a data
#include "EGFRActivationWriter.hpp"                     // Output EGFR activation data

// ============================================================================
// EGFRMutationType  (added for per-cell heterogeneity)
// ============================================================================
enum class EGFRMutationType
{
    EXON19_DEL,  // dE746-A750 — PC9 canonical mutation
    L858R        // Leu858Arg  — activation loop mutation
};

double GetEGFRBaseRescueFraction(EGFRMutationType mutationType)
{
    switch (mutationType)
    {
        case EGFRMutationType::EXON19_DEL: return 0.5;
        case EGFRMutationType::L858R:      return 0.35;
        default:                           return 0.5;
    }
}

std::string GetMutationLabel(EGFRMutationType mutationType)
{
    switch (mutationType)
    {
        case EGFRMutationType::EXON19_DEL: return "Exon19del_dE746-A750";
        case EGFRMutationType::L858R:      return "L858R";
        default:                           return "Unknown";
    }
}

static double SampleBeta(double alpha, double beta)
{
    double mean = alpha / (alpha + beta);
    double variance = (alpha * beta) /
        ((alpha + beta) * (alpha + beta) * (alpha + beta + 1.0));
    double sd = std::sqrt(variance);

    double u = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();
    double sample = mean + sd * u;

    if (sample < 0.01) sample = 0.01;
    if (sample > 0.99) sample = 0.99;
    return sample;
}

class TestPC9Spheroid2D : public AbstractCellBasedTestSuite
{
public:
    void RunSpheroidSimulation(double tgf_diffusion_radius, int sweep_index, unsigned replicate_index)
    {
        std::cout << "\n=== Sweep " << sweep_index << ", Replicate " << replicate_index
                  << ": TGF-alpha diffusion radius = " << tgf_diffusion_radius 
                  << " ===" << std::endl;

        // [USER SETTING] — EGFR MUTATION TYPE and BETA DISTRIBUTION SHAPE PARAMETERS
        // (added for per-cell heterogeneity, matching TestPC9Spheroid2D_Heterogeneous_Resumable.hpp)
        EGFRMutationType egfr_mutation = EGFRMutationType::EXON19_DEL;
        const double beta_alpha = 1.0;
        const double beta_beta  = 9.0;
        double base_rescue_fraction = GetEGFRBaseRescueFraction(egfr_mutation);
        const double g1_min         = 54.0;
        const double g1_extra_hours = 54.0;
        
        // Simulation parameters
        double dt = 0.01; // 36 sec
        double end_time = 240.0;  // 7days 
        double sampling_timestep = 4.0/dt;  // Sample every hour
        
        // Spheroid parameters 
        double spheroid_radius = 8.0;
        unsigned mesh_size = 30;
        
        RandomNumberGenerator::Instance()->Reseed(replicate_index);
        
        // Generate CIRCULAR honeycomb mesh
        HoneycombMeshGenerator generator(mesh_size, mesh_size, 0);
        boost::shared_ptr<MutableMesh<2,2>> p_mesh_shared = generator.GetCircularMesh(spheroid_radius);
        MutableMesh<2,2>* p_mesh = p_mesh_shared.get();
        
        // ============================================================
        // ADD SPATIAL NOISE: jitter each node position slightly
        // This mimics the irregular cell packing seen in real spheroids.
        // A perturbation of ~0.1 cell diameters is small enough to not
        // break the mesh topology but large enough to break symmetry.
        // ============================================================
        double jitter_magnitude = 0.1;  // fraction of cell diameter (~1.0 in Chaste units)

        for (unsigned i = 0; i < p_mesh->GetNumNodes(); i++)
        {
            c_vector<double, 2>& location = p_mesh->GetNode(i)->rGetModifiableLocation();
            
            // Random displacement in x: uniform in [-jitter_magnitude, +jitter_magnitude]
            double dx = jitter_magnitude * (2.0 * RandomNumberGenerator::Instance()->ranf() - 1.0);
            double dy = jitter_magnitude * (2.0 * RandomNumberGenerator::Instance()->ranf() - 1.0);
            
            location[0] += dx;
            location[1] += dy;
        }
        // Get all node indices
        std::vector<unsigned> all_indices = generator.GetCellLocationIndices();
        
        // Spheroid center at origin
        c_vector<double, 2> spheroid_center = zero_vector<double>(2);
        
        // Filter nodes to create circular spheroid
        std::vector<unsigned> location_indices;
        for (unsigned i = 0; i < all_indices.size(); i++)
        {
            unsigned node_index = all_indices[i];
            c_vector<double, 2> node_location = p_mesh->GetNode(node_index)->rGetLocation();
            double radius = norm_2(node_location - spheroid_center);
            
            if (radius <= spheroid_radius)
            {
                location_indices.push_back(node_index);
            }
        }
        
        std::cout << "Created circular mesh with " << location_indices.size() 
                  << " cells within radius" << std::endl;
        
        // Cell properties
        boost::shared_ptr<AbstractCellProperty> p_wildtype = 
            CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();
        boost::shared_ptr<AbstractCellProperty> p_hypoxic_state = 
            CellPropertyRegistry::Instance()->Get<ApcOneHitCellMutationState>();
        boost::shared_ptr<AbstractCellProperty> p_necrotic_state = 
            CellPropertyRegistry::Instance()->Get<ApcTwoHitCellMutationState>();
        //boost::shared_ptr<AbstractCellProperty> p_transit_type = 
           // CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>();
        boost::shared_ptr<AbstractCellProperty> p_stem_type = 
        CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>();
        //boost::shared_ptr<AbstractCellProperty> p_diff_type = 
            //CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
        
        // Create cells
        std::vector<CellPtr> cells;
        
        // Register cell data items for output
        CellPropertyRegistry::Instance()->Clear();

        // Tracking for per-cell EGFR mut:wt heterogeneity (added)
        double sum_mut_frac_init = 0.0;
        double sum_sq_mut_frac_init = 0.0;
        double min_mut_frac_init = 1.0;
        double max_mut_frac_init = 0.0;
        
        for (unsigned i = 0; i < location_indices.size(); i++)
        {
            // Sample this cell's EGFR mut:wt fraction from the Beta distribution (added)
            double cell_mut_frac = SampleBeta(beta_alpha, beta_beta);
            double cell_wt_frac  = 1.0 - cell_mut_frac;
            double cell_g1 = g1_min + g1_extra_hours * cell_wt_frac;
            double cell_rescue = base_rescue_fraction * cell_mut_frac;
            double cell_total_cycle = cell_g1 + 16.0 + 8.0 + 2.0;

            // Create oxygen-sensitive cell cycle model
            SimpleOxygenBasedCellCycleModel* p_cycle_model = 
                new SimpleOxygenBasedCellCycleModel();
            p_cycle_model->SetDimension(2);
            
            // Set oxygen thresholds 
            p_cycle_model->SetQuiescentConcentration(0.08);
            p_cycle_model->SetHypoxicConcentration(0.05);
            p_cycle_model->SetCriticalHypoxicDuration(27.0);
            
            // All cells start as viable cells
            CellPtr p_cell = CellPtr(new Cell(p_wildtype, p_cycle_model));
            p_cell->SetCellProliferativeType(p_stem_type);
            
            // Initialize cell cycle model first
            //p_cell->InitialiseCellCycleModel();
            
            p_cycle_model->SetStemCellG1Duration(cell_g1);   // G1 (heterogeneous, was 54.0)
            //p_cycle_model->SetTransitCellG1Duration(24.0);  
            p_cycle_model->SetSDuration(16.0);            // S
            p_cycle_model->SetG2Duration(8.0);            // G2
            p_cycle_model->SetMDuration(2.0);             // M

            // Initialize cell cycle model first
            p_cell->InitialiseCellCycleModel();

            // Randomize birth time to avoid synchronization
            double birth_time = cell_total_cycle * RandomNumberGenerator::Instance()->ranf();
            p_cycle_model->SetBirthTime(-birth_time);
                                
            // Initial oxygen (will be updated by modifier immediately)
            p_cell->GetCellData()->SetItem("oxygen", 1.0);
            p_cell->GetCellData()->SetItem("cell_state", 1.0);
            p_cell->GetCellData()->SetItem("time_became_hypoxic", SimulationTime::Instance()->GetTime());
            p_cell->GetCellData()->SetItem("HIF1Alpha", 0.0);
            p_cell->GetCellData()->SetItem("TGFAlpha", 0.0);
            p_cell->GetCellData()->SetItem("EGFR_activation", 0.0);

            // Per-cell EGFR mut:wt heterogeneity data (added)
            p_cell->GetCellData()->SetItem("egfr_mut_fraction", cell_mut_frac);
            p_cell->GetCellData()->SetItem("egfr_wt_fraction",  cell_wt_frac);
            p_cell->GetCellData()->SetItem("egfr_rescue",        cell_rescue);

            cells.push_back(p_cell);

            sum_mut_frac_init    += cell_mut_frac;
            sum_sq_mut_frac_init += cell_mut_frac * cell_mut_frac;
            if (cell_mut_frac < min_mut_frac_init) min_mut_frac_init = cell_mut_frac;
            if (cell_mut_frac > max_mut_frac_init) max_mut_frac_init = cell_mut_frac;
        }
        
        std::cout << "Created " << cells.size() << " PC9 cells" << std::endl;

        // Report initial mut_frac distribution (added)
        unsigned n_cells_init = location_indices.size();
        double obs_mean = sum_mut_frac_init / n_cells_init;
        double obs_var  = (sum_sq_mut_frac_init / n_cells_init) - obs_mean * obs_mean;
        double obs_sd   = std::sqrt(std::max(0.0, obs_var));

        std::cout << "  Initial mut_frac distribution: mean = " << obs_mean
                  << "  sd = " << obs_sd
                  << "  min = " << min_mut_frac_init
                  << "  max = " << max_mut_frac_init << std::endl;
        
        // Create cell population (no ghost nodes for spheroid)
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);
        
        // Add cell writers - these generate .dat files
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellAncestorWriter>();
        //cell_population.AddCellWriter<CellDataWriter>();  // Writes oxygen, HIF1Alpha, TGFAlpha, EGFR_activation, etc.
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.AddPopulationWriter<VoronoiDataWriter>();
        //MAKE_PTR(SpheroidRadiusWriter, p_radius_writer);
        //cell_population.AddPopulationWriter(p_radius_writer);

        // Create output directory for this run
        //std::string output_dir = "PC9_Spheroid2D_TGFDiff_Try2_" + boost::lexical_cast<std::string>(sweep_index);
        
        // Create output directory: sweep index + replicate index
        std::string output_dir = "PC9_Hetero_TGFDiffusion_sweep_" 
            + boost::lexical_cast<std::string>(sweep_index) 
            + "_rep_" 
            + boost::lexical_cast<std::string>(replicate_index);

        boost::shared_ptr<CellDataWriter<2,2>> p_cell_data_writer = boost::shared_ptr<CellDataWriter<2,2>>(new CellDataWriter<2,2>());
        p_cell_data_writer->SetOutputDirectory(output_dir);
        cell_population.AddCellWriter(p_cell_data_writer);

        // Create simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_dir);
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(sampling_timestep);
        simulator.SetEndTime(end_time);
        
        std::cout << "Simulation parameters:" << std::endl;
        std::cout << "  Duration: " << end_time << " hours" << std::endl;
        std::cout << "  Output directory: " << output_dir << std::endl;
        
        // Add spring force (tumor cells)
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetMeinekeSpringStiffness(15.0);       // Stiffer springs for spheroid cells
        p_force->SetCutOffLength(1.5);                  // Interact up to 1.5 cell diameters
        simulator.AddForce(p_force);
        
        // Literature: oxygen penetrates ~100-200 um = ~5-10 cell diameters
        MAKE_PTR_ARGS(SpheroidOxygenGradientModifier, p_oxygen_modifier, 
             (spheroid_radius,  // Initial spheroid radius (12 units)
              5.0,              // Oxygen penetration depth (~100 um)
              0.08,             // Max O2 at surface (normoxic)
              0.001,             // Min O2 in necrotic core
              9.0));           // Hypoxia onset at radius ~10 units (~200 um)
        simulator.AddSimulationModifier(p_oxygen_modifier);
        
        // Add HIF-1a and TGF-a signaling
        // SWEPT PARAMETER: TGF-alpha diffusion radius
        MAKE_PTR_ARGS(HypoxiaSignalingModifier, p_hypoxia_signaling, 
            (0.05,                  // Hypoxia threshold
                1.0,                // Max HIF-1a concentration
                0.01,                // HIF-1a degradation rate (half life is 0.1)
                0.004,               // TGF-a production rate (assumption)
                tgf_diffusion_radius,// TGF-a diffusion radius (SWEPT PARAMETER)
                0.001,               // TGF-a degradation rate (assumption)
                1.0,                  // Max TGF-a concentration
                0.005));              // dt
        simulator.AddSimulationModifier(p_hypoxia_signaling);

        // Add EGFR-TGF-a signaling (reads HIF-1a and TGF-a from above)
        MAKE_PTR_ARGS(EGFRSignalingModifier, p_egfr_signaling, 
            (0.001,                   // TGF-a activation threshold
                1.0,                // WT EGFR sensitivity
                1.0));                // Mutant EGFR basal activity
        simulator.AddSimulationModifier(p_egfr_signaling);

        // Add hypoxic cell state modifier (marks necrotic cells without removing them)
        MAKE_PTR_ARGS(HypoxicCellStateModifier, p_hypoxic_modifier, 
            (0.001,                // Severe hypoxia threshold
             27.0));              // Time to necrosis under severe hypoxia
        simulator.AddSimulationModifier(p_hypoxic_modifier);
        
        // Run simulation
        std::cout << "Starting simulation...\n" << std::endl;
        
        simulator.Solve();
        
        // Final statistics
        unsigned total_cells = 0;
        unsigned proliferative_cells = 0;
        unsigned quiescent_cells = 0;
        unsigned dead_cells = 0;

        // Per-mut_frac final statistics (added)
        double final_sum_mut_frac    = 0.0;
        double final_sum_sq_mut_frac = 0.0;
        double final_min_mut_frac    = 1.0;
        double final_max_mut_frac    = 0.0;

        double normoxic_sum_mut = 0.0;
        double hypoxic_sum_mut  = 0.0;
        double necrotic_sum_mut = 0.0;
        unsigned n_normoxic = 0, n_hypoxic = 0, n_necrotic = 0;
        
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            total_cells++;
            double oxygen = cell_iter->GetCellData()->GetItem("oxygen");
            double mut_frac = cell_iter->GetCellData()->GetItem("egfr_mut_fraction");

            final_sum_mut_frac    += mut_frac;
            final_sum_sq_mut_frac += mut_frac * mut_frac;
            if (mut_frac < final_min_mut_frac) final_min_mut_frac = mut_frac;
            if (mut_frac > final_max_mut_frac) final_max_mut_frac = mut_frac;
            
            if (oxygen > 0.05)
            {
                proliferative_cells++;
                normoxic_sum_mut += mut_frac;
                n_normoxic++;
            }
            else if (oxygen > 0.01)
            {
                quiescent_cells++;
                hypoxic_sum_mut += mut_frac;
                n_hypoxic++;
            }
            else
            {
                dead_cells++;
                necrotic_sum_mut += mut_frac;
                n_necrotic++;
            }
        }

        double f_mean = (total_cells > 0) ? final_sum_mut_frac / total_cells : 0.0;
        double f_var  = (total_cells > 0)
            ? (final_sum_sq_mut_frac / total_cells) - f_mean * f_mean : 0.0;
        double f_sd   = std::sqrt(std::max(0.0, f_var));

        std::cout << "\n--- Sweep " << sweep_index << ", Replicate " << replicate_index << " Results ---" << std::endl;
        std::cout << "\n=== Final Results for TGF Diffusion Radius " << tgf_diffusion_radius << " ===" << std::endl;
        std::cout << "Total cells: " << total_cells << std::endl;
        std::cout << "Normoxic (O2>0.05): " << proliferative_cells 
                  << " (" << 100.0*proliferative_cells/total_cells << "%)"
                  << "  mean mut_frac = " << (n_normoxic > 0 ? normoxic_sum_mut/n_normoxic : 0.0) << std::endl;
        std::cout << "Hypoxic (0.01<O2<0.05): " << quiescent_cells 
                  << " (" << 100.0*quiescent_cells/total_cells << "%)"
                  << "  mean mut_frac = " << (n_hypoxic > 0 ? hypoxic_sum_mut/n_hypoxic : 0.0) << std::endl;
        std::cout << "Severely hypoxic/necrotic (O2<0.01): " << dead_cells 
                  << " (" << 100.0*dead_cells/total_cells << "%)"
                  << "  mean mut_frac = " << (n_necrotic > 0 ? necrotic_sum_mut/n_necrotic : 0.0) << std::endl;
        std::cout << "Final mut_frac distribution: mean = " << f_mean << "  sd = " << f_sd
                  << "  min = " << final_min_mut_frac << "  max = " << final_max_mut_frac << std::endl;
        
        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
    }

    void TestPC9SpheroidCrossSection()
    {
        std::cout << "\n=== PC9 Spheroid 2D Cross-Section Parameter Sweep (with replicates) ===\n" << std::endl;
        
        std::vector<double> tgf_diffusion_radii = {0.25, 1.0, 2.0, 4.0};

        unsigned start_rep = 1;
        unsigned num_reps = 5;

        // --- RESUME SETTINGS ---
        unsigned resume_sweep_index = 1;   // 1-based sweep index to resume from (tgf_diffusion_radius = 2.0)
        unsigned resume_rep = 1;           // replicate to resume from within that sweep

        std::cout << "Resuming from Sweep " << resume_sweep_index 
                << ", Replicate " << resume_rep << std::endl;
        std::cout << "Total runs: " << tgf_diffusion_radii.size() * num_reps << std::endl;

        for (unsigned i = 0; i < tgf_diffusion_radii.size(); i++)
        {
            unsigned sweep_index = i + 1;  // 1-based

            // Skip sweeps before the resume point
            if (sweep_index < resume_sweep_index) continue;

            for (unsigned rep = start_rep; rep < start_rep + num_reps; rep++)
            {
                // Within the resume sweep, skip replicates before the resume rep
                if (sweep_index == resume_sweep_index && rep < resume_rep) continue;

                RunSpheroidSimulation(tgf_diffusion_radii[i], sweep_index, rep);
            }
        }
        
        std::cout << "\n=== All runs completed! ===" << std::endl;
        std::cout << "4 TGF production rates x 5 replicates = 20 simulations" << std::endl;
    }
};

#endif /* TESTPC9SPHEROID2D_HPP_ */
