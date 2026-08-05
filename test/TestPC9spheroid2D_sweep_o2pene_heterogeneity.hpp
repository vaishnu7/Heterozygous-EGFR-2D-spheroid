/*
 * TestPC9spheroid2D_sweep_o2pene_heterogeneity.hpp
 *
 * Created on: 24 Oct, 2025
 * Author: Vaishnudebi Dutta
 * Modified: Sweep on oxygen penetration depth
 * FIXED: 16 Feb, 2026
 * HETEROGENEITY ADDED: 29 Jun, 2026 (aligned to TestPC9Spheroid2D_Heterogeneous_Resumable.hpp)
 * MESH/CELL MISMATCH BUG FIXED: 29 Jun, 2026 (same fix applied to the radius sweep)
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
#include <algorithm>                                    // std::min / std::max for clamping

// Custom header files
#include "HypoxiaSignalingModifier.hpp"                 // HIF-1a and TGF-a signaling under hypoxia
#include "EGFRSignalingModifier.hpp"                    // EGFR-TGF-a signaling pathway
#include "SpheroidOxygenGradientModifier.hpp"           // Radial oxygen gradient in spheroid
#include "HypoxiaCellStateModifier.hpp"                 // Hypoxic cell state changes
#include "HIF1AlphaWriter.hpp"                          // Output HIF-1a data
#include "TGFAlphaWriter.hpp"                           // Output TGF-a data
#include "EGFRActivationWriter.hpp"                     // Output EGFR activation data


// ============================================================================
// EGFRMutationType
// ============================================================================
enum class EGFRMutationType
{
    EXON19_DEL,  
    L858R        
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

// ============================================================================
// SampleBeta
// ============================================================================
// Samples a Beta(alpha, beta) random variate using Chaste's own
// StandardNormalRandomDeviate(), via a normal approximation:
//   mean = alpha/(alpha+beta)
//   variance = alpha*beta / ((alpha+beta)^2 * (alpha+beta+1))
// clamped to [0.01, 0.99] to avoid degenerate 0/1 fractions.
// ============================================================================
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
    void RunSpheroidSimulation(double spheroid_radius, double oxygen_penetration_depth, int sweep_index, unsigned replicate_index)
    {
        std::cout << "\n=== Sweep " << sweep_index << ", Replicate " << replicate_index
                  << ": Spheroid radius = " << spheroid_radius 
                  << " cell units (~" << spheroid_radius * 15 << " um), O2 penetration = " 
                  << oxygen_penetration_depth << " cell units (~" << oxygen_penetration_depth * 15 
                  << " um) ===" << std::endl;

        // --------------------------------------------------------------
        // Heterogeneity parameters (EGFR mut:wt ratio, per cell)
        // Fixed across the whole sweep: EXON19_DEL mutation type,
        // Beta(5,5) shape parameters (matches the radius-sweep script).
        // --------------------------------------------------------------
        EGFRMutationType egfr_mutation = EGFRMutationType::EXON19_DEL;
        const double beta_alpha = 5.0;
        const double beta_beta  = 5.0;
        double base_rescue_fraction = GetEGFRBaseRescueFraction(egfr_mutation);
        const double g1_min         = 54.0;
        const double g1_extra_hours = 54.0;
        
        // Simulation parameters
        double dt = 0.01; // 36 sec
        double end_time = 240.0;  // 7 days 
        double sampling_timestep = 4.0/dt;  // Sample every hour
        
        // Spheroid parameters 
        unsigned mesh_size = 30;
        
        // Different seed per replicate
        RandomNumberGenerator::Instance()->Reseed(replicate_index);
        
        // Generate CIRCULAR honeycomb mesh
        HoneycombMeshGenerator generator(mesh_size, mesh_size, 0);
        boost::shared_ptr<MutableMesh<2,2>> p_mesh_shared = generator.GetCircularMesh(spheroid_radius);
        MutableMesh<2,2>* p_mesh = p_mesh_shared.get();

        // location_indices must cover every node GetCircularMesh produced,
        // one cell per node, with no extras and no omissions.
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Add spatial jitter AFTER fixing location_indices, applied to all
        // nodes in the (already circular) mesh -- jitter only perturbs
        // positions, it must never change which nodes exist.
        double jitter_magnitude = 0.1;
        for (unsigned i = 0; i < p_mesh->GetNumNodes(); i++)
        {
            c_vector<double, 2>& location = p_mesh->GetNode(i)->rGetModifiableLocation();
            double dx = jitter_magnitude * (2.0 * RandomNumberGenerator::Instance()->ranf() - 1.0);
            double dy = jitter_magnitude * (2.0 * RandomNumberGenerator::Instance()->ranf() - 1.0);
            location[0] += dx;
            location[1] += dy;
        }
        
        // Spheroid center at origin (kept for any downstream radius-based
        // calculations, e.g. logging/diagnostics -- no longer used to
        // re-filter nodes).
        c_vector<double, 2> spheroid_center = zero_vector<double>(2);
        
        std::cout << "Created circular mesh with " << location_indices.size() 
                  << " cells within radius" << std::endl;
        
        // Cell properties
        boost::shared_ptr<AbstractCellProperty> p_wildtype = 
            CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();
        boost::shared_ptr<AbstractCellProperty> p_hypoxic_state = 
            CellPropertyRegistry::Instance()->Get<ApcOneHitCellMutationState>();
        boost::shared_ptr<AbstractCellProperty> p_necrotic_state = 
            CellPropertyRegistry::Instance()->Get<ApcTwoHitCellMutationState>();
        boost::shared_ptr<AbstractCellProperty> p_stem_type = 
        CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>();
        
        // Create cells
        std::vector<CellPtr> cells;
        
        // Register cell data items for output
        CellPropertyRegistry::Instance()->Clear();

        double sum_mut_frac_init = 0.0;
        double sum_sq_mut_frac_init = 0.0;
        double min_mut_frac_init = 1.0;
        double max_mut_frac_init = 0.0;
        
        for (unsigned i = 0; i < location_indices.size(); i++)
        {
            // ----------------------------------------------------------
            // HETEROGENEITY: per-cell EGFR mut:wt fraction
            //
            // Each cell is an independent agent with its own EGFR allelic
            // composition, sampled from a Beta(alpha, beta) distribution.
            // Biological basis: PC9 tumors harbour pre-existing EGFR-low
            // subclones (Alsaed et al. 2025); within-subline variability
            // is epigenetic/stochastic (Camp et al. 2021).
            // ----------------------------------------------------------
            double cell_mut_frac = SampleBeta(beta_alpha, beta_beta);
            double cell_wt_frac  = 1.0 - cell_mut_frac;

            // Per-cell G1: longer in WT-enriched cells (slower cycling)
            double cell_g1 = g1_min + g1_extra_hours * cell_wt_frac;

            // Per-cell rescue fraction: scales linearly with mut content.
            // Full mutant -> full base rescue; full WT -> zero rescue.
            double cell_rescue = base_rescue_fraction * cell_mut_frac;

            double cell_total_cycle = cell_g1 + 16.0 + 8.0 + 2.0; // G1 + S + G2 + M

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
            
            // Per-cell G1 duration is the key heterogeneous cell-cycle parameter
            p_cycle_model->SetStemCellG1Duration(cell_g1);   // G1 (per-cell)
            p_cycle_model->SetSDuration(16.0);            // S
            p_cycle_model->SetG2Duration(8.0);            // G2
            p_cycle_model->SetMDuration(2.0);             // M

            // Initialize cell cycle model first
            p_cell->InitialiseCellCycleModel();

            // Randomize birth time across this cell's own (heterogeneous)
            // total cycle length so cells start distributed across all phases
            double birth_time = cell_total_cycle * RandomNumberGenerator::Instance()->ranf();
            p_cycle_model->SetBirthTime(-birth_time);
                                
            // Initial oxygen (will be updated by modifier immediately)
            p_cell->GetCellData()->SetItem("oxygen", 1.0);
            p_cell->GetCellData()->SetItem("cell_state", 1.0);
            p_cell->GetCellData()->SetItem("time_became_hypoxic", SimulationTime::Instance()->GetTime());
            p_cell->GetCellData()->SetItem("hif1alpha", 0.0);
            p_cell->GetCellData()->SetItem("tgfalpha", 0.0);
            p_cell->GetCellData()->SetItem("egfr_activation", 0.0);

            // ----------------------------------------------------------
            // HETEROGENEITY: store per-cell EGFR fractions in CellData,
            // read by HypoxiaSignalingModifier / EGFRSignalingModifier.
            // ----------------------------------------------------------
            p_cell->GetCellData()->SetItem("egfr_mut_fraction", cell_mut_frac);
            p_cell->GetCellData()->SetItem("egfr_wt_fraction",  cell_wt_frac);
            p_cell->GetCellData()->SetItem("egfr_rescue",        cell_rescue);

            cells.push_back(p_cell);

            sum_mut_frac_init    += cell_mut_frac;
            sum_sq_mut_frac_init += cell_mut_frac * cell_mut_frac;
            if (cell_mut_frac < min_mut_frac_init) min_mut_frac_init = cell_mut_frac;
            if (cell_mut_frac > max_mut_frac_init) max_mut_frac_init = cell_mut_frac;
        }
        
        std::cout << "Created " << cells.size() << " PC9 cells (heterogeneous EGFR mut:wt)" << std::endl;

        unsigned n_cells_init = location_indices.size();
        double init_mean = sum_mut_frac_init / n_cells_init;
        double init_var  = (sum_sq_mut_frac_init / n_cells_init) - init_mean * init_mean;
        double init_sd   = std::sqrt(std::max(0.0, init_var));

        std::cout << "  Initial mut_frac distribution: mean = " << init_mean
                  << "  sd = " << init_sd
                  << "  min = " << min_mut_frac_init
                  << "  max = " << max_mut_frac_init << std::endl;
        
        // Create cell population (no ghost nodes for spheroid)
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);
        
        // Add cell writers - these generate .dat files
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellAncestorWriter>();
        cell_population.AddCellWriter<HIF1AlphaWriter>();
        cell_population.AddCellWriter<TGFAlphaWriter>();
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        // Create output directory for this run
        // (renamed to flag this as the heterogeneous-EGFR version of the
        // O2-penetration-depth sweep)
        std::string output_dir = "PC9_Spheroid2D_O2Penetration_Heterogeneous_EGFR_"
            + GetMutationLabel(egfr_mutation)
            + "_O2pen" + boost::lexical_cast<std::string>(oxygen_penetration_depth)
            + "_rep_" + boost::lexical_cast<std::string>(replicate_index);
        
        // Create CellDataWriter explicitly with output directory
        boost::shared_ptr<CellDataWriter<2,2>> p_cell_data_writer(new CellDataWriter<2,2>());
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
        
        // SWEPT PARAMETER: Oxygen penetration depth
        // Literature: oxygen penetrates ~100-200 um = ~5-10 cell diameters
        MAKE_PTR_ARGS(SpheroidOxygenGradientModifier, p_oxygen_modifier, 
             (spheroid_radius,              // Fixed spheroid radius (8 units)
              oxygen_penetration_depth,     // SWEPT: oxygen penetration depth
              0.08,                         // Max O2 at surface (normoxic)
              0.001,                        // Min O2 in necrotic core
              9.0));                        // Hypoxia onset at radius ~10 units (~200 um)
        simulator.AddSimulationModifier(p_oxygen_modifier);
        
        // Add HIF-1a and TGF-a signaling
        MAKE_PTR_ARGS(HypoxiaSignalingModifier, p_hypoxia_signaling, 
            (0.05,                  // Hypoxia threshold
                1.0,                // Max HIF-1a concentration
                0.01,                // HIF-1a degradation rate (half life is 0.1)
                0.004,                // TGF-a production rate (assumption)
                0.5,                // TGF-a diffusion radius (assumption)
                0.001,               // TGF-a degradation rate (assumption)
                1.0,                  // Max TGF-a concentration
                dt,                   // dt
                0.05,                 // HIF basal synthesis rate (confirmed)
                0.08,                 // Max O2-dependent HIF degradation (confirmed)
                0.003,                // O2-independent HIF degradation rate (confirmed)
                0.0,                  // UNVERIFIED — copied as-is from reference script
                0.0));                // UNVERIFIED — copied as-is from reference script
        simulator.AddSimulationModifier(p_hypoxia_signaling);

        // Add hypoxic cell state modifier (marks necrotic cells without removing them)
        MAKE_PTR_ARGS(HypoxicCellStateModifier, p_hypoxic_modifier, 
            (0.001,                // Severe hypoxia threshold
             27.0));              // Time to necrosis under severe hypoxia
        simulator.AddSimulationModifier(p_hypoxic_modifier);

        // Add EGFR-TGF-a signaling (reads HIF-1a and TGF-a from above, and
        // each cell's own "egfr_mut_fraction" / "egfr_wt_fraction" /
        // "egfr_rescue" from CellData)
        MAKE_PTR_ARGS(EGFRSignalingModifier, p_egfr_signaling, 
            (0.001,                   // TGF-a activation threshold
                1.0,                // WT EGFR sensitivity
                1.0,                // Mutant EGFR basal activity
                0.3,                // UNVERIFIED — copied as-is from reference script
                2,                  // UNVERIFIED — copied as-is from reference script
                0.25,               // UNVERIFIED — copied as-is from reference script
                dt,                 // dt (inferred by position/context only)
                0.0));              // UNVERIFIED — copied as-is from reference script
        simulator.AddSimulationModifier(p_egfr_signaling);
        
        // Run simulation
        std::cout << "Starting simulation...\n" << std::endl;
        
        simulator.Solve();
        
        // Final statistics
        unsigned total_cells = 0;
        unsigned proliferative_cells = 0;
        unsigned quiescent_cells = 0;
        unsigned dead_cells = 0;

        double final_sum_mut_frac = 0.0;
        double final_sum_sq_mut_frac = 0.0;
        double final_min_mut_frac = 1.0;
        double final_max_mut_frac = 0.0;

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

        double final_mean = (total_cells > 0) ? final_sum_mut_frac / total_cells : 0.0;
        double final_var  = (total_cells > 0)
            ? (final_sum_sq_mut_frac / total_cells) - final_mean * final_mean : 0.0;
        double final_sd   = std::sqrt(std::max(0.0, final_var));
        
        std::cout << "\n--- Sweep " << sweep_index << ", Replicate " << replicate_index 
                  << " (O2 Penetration " << oxygen_penetration_depth << ") ---" << std::endl;
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
        std::cout << "Final mut_frac distribution: mean = " << final_mean
                  << "  sd = " << final_sd
                  << "  min = " << final_min_mut_frac
                  << "  max = " << final_max_mut_frac << std::endl;
        
        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
    }

    void TestPC9SpheroidCrossSection()
    {
        std::cout << "\n=== PC9 Spheroid 2D Cross-Section: O2 Penetration Parameter Sweep (Heterogeneous EGFR) ===\n" << std::endl;
        
        double spheroid_radius = 8.0;
        
        std::vector<double> o2_penetration_depths = {2.0, 4.0, 6.0, 8.0};
        
        unsigned start_rep = 1;
        unsigned num_reps = 5;

        // --- RESUME SETTINGS (index-based, unchanged from original sweep) ---
        unsigned resume_sweep_index = 1;   // 1-based sweep index to resume from
        unsigned resume_rep = 1;           // replicate to resume from within that sweep

        std::cout << "Fixed spheroid radius: " << spheroid_radius << " cell units (~" 
                << spheroid_radius * 15 << " um)" << std::endl;
        std::cout << "Resuming from Sweep " << resume_sweep_index 
                << ", Replicate " << resume_rep << std::endl;
        std::cout << "Total runs: " << o2_penetration_depths.size() * num_reps << std::endl;
        std::cout << std::endl;

        for (unsigned i = 0; i < o2_penetration_depths.size(); i++)
        {
            unsigned sweep_index = i + 1;  // 1-based

            // Skip sweeps before the resume point
            if (sweep_index < resume_sweep_index) continue;

            for (unsigned rep = start_rep; rep < start_rep + num_reps; rep++)
            {
                // Within the resume sweep, skip replicates before the resume rep
                if (sweep_index == resume_sweep_index && rep < resume_rep) continue;

                RunSpheroidSimulation(spheroid_radius, o2_penetration_depths[i], sweep_index, rep);
            }
        }
        
        std::cout << "\n=== Parameter sweep completed! ===" << std::endl;
        std::cout << "4 O2 penetration depths x 5 replicates = 20 simulations" << std::endl;
    }
};

#endif /* TESTPC9SPHEROID2D_HPP_ */
