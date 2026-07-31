/*
 * TestPC9Spheroid2D_TGFaProdDiffDegSweep_Resumable.hpp
 *
 * 2D cross-section model of cancer tumor spheroid with hypoxic zones.
 * JOINT parameter sweep on TGF-alpha PRODUCTION RATE x TGF-alpha DIFFUSION
 * RADIUS x TGF-alpha DEGRADATION RATE (i.e. a 3D grid of every production
 * rate against every diffusion radius against every degradation rate),
 * with replicates at each grid point.
 *
 * Models radial oxygen gradient from spheroid surface to hypoxic/necrotic core.
 * Incorporates HIF and TGF-a signaling under hypoxia, affecting proliferation.
 * Based on off-lattice cell-based modeling framework in Chaste.
 *
 * Per-cell intratumoral heterogeneity in the EGFR mut:wt ratio (Beta-sampled),
 * carried over from TestPC9Spheroid2D_Heterogeneous_Resumable.hpp.
 *
 * RESUME LOGIC: same pattern as TestPC9Spheroid2D_Heterogeneous_Resumable.hpp —
 *   - a ".run_complete" marker file is written inside each run's own output
 *     directory only after that run finishes successfully
 *   - before running, check for that marker; if present, skip it
 *
 * Original single-parameter sweep author: Vaishnudebi Dutta (24 Oct 2025)
 * Heterogeneous resumable reference:      19 Jun 2026
 * Extended to joint production x diffusion sweep + heterogeneity: 22 Jul 2026
 * Extended to joint production x diffusion x degradation sweep:   23 Jul 2026
 *
 * NOTE: Modifier signatures/values, cell-data key casing, mesh jitter
 * ordering, and cell writers below are all matched to the reference script
 * TestPC9Spheroid2D_Heterogeneous_Resumable.hpp so that only the TGF-a
 * production rate, diffusion radius, and degradation rate differ from that
 * reference (as the three swept axes) — everything else is unchanged from it.
 */

#ifndef TESTPC9SPHEROID2D_TGFAPRODDIFFDEGSWEEP_RESUMABLE_HPP_
#define TESTPC9SPHEROID2D_TGFAPRODDIFFDEGSWEEP_RESUMABLE_HPP_

#include <cxxtest/TestSuite.h>
#include "CellBasedSimulationArchiver.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"

// Core simulation classes
#include "HoneycombMeshGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "VoronoiDataWriter.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "FakePetscSetup.hpp"

// Cell cycle and mutation classes
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"

// Force classes
#include "GeneralisedLinearSpringForce.hpp"

// Utility classes
#include <boost/lexical_cast.hpp>
#include "CellAgesWriter.hpp"
#include "OutputFileHandler.hpp"
#include "CellDataWriter.hpp"
#include "CellData.hpp"
#include "HIF1AlphaWriter.hpp"
#include "TGFAlphaWriter.hpp"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>

// Custom header files
#include "HypoxiaSignalingModifier.hpp"
#include "EGFRSignalingModifier.hpp"
#include "SpheroidOxygenGradientModifier.hpp"
#include "HypoxiaCellStateModifier.hpp"

// ============================================================================
// EGFRMutationType — carried over from the heterogeneous reference script
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

class TestPC9Spheroid2D_TGFaProdDiffDegSweep_Resumable : public AbstractCellBasedTestSuite
{
private:

    /**
     * Path to the completion marker file for a given run's output directory.
     * Lives inside that run's own output folder so it can never be confused
     * with another run's marker.
     */
    std::string GetCompletionMarkerPath(const std::string& rOutputDir)
    {
        std::string output_dir = OutputFileHandler::GetChasteTestOutputDirectory()
            + rOutputDir + "/";
        return output_dir + ".run_complete";
    }

    /** True if the given run's output directory already has a completion marker. */
    bool IsRunComplete(const std::string& rOutputDir)
    {
        std::ifstream marker(GetCompletionMarkerPath(rOutputDir));
        return marker.good();
    }

    /** Write the completion marker once a run finishes successfully. */
    void MarkRunComplete(const std::string& rOutputDir)
    {
        std::ofstream marker(GetCompletionMarkerPath(rOutputDir));
        marker << "completed" << std::endl;
        marker.close();
    }

public:

    /**
     * Runs a single simulation at one (production rate, diffusion radius,
     * degradation rate, replicate) grid point, with per-cell EGFR mut:wt
     * heterogeneity.
     *
     * @param egfrMutation        EGFR mutation subtype (sets base rescue fraction)
     * @param betaAlpha           Beta distribution alpha shape param for mut_frac
     * @param betaBeta            Beta distribution beta shape param for mut_frac
     * @param tgfProductionRate   TGF-a production rate under hypoxia (SWEPT: axis 1)
     * @param tgfDiffusionRadius  TGF-a diffusion radius (SWEPT: axis 2)
     * @param tgfDegradationRate  TGF-a degradation rate (SWEPT: axis 3)
     * @param prodSweepIndex      1-based index into the production-rate list
     * @param diffSweepIndex      1-based index into the diffusion-radius list
     * @param degSweepIndex       1-based index into the degradation-rate list
     * @param replicateIndex      1-based replicate number within this grid point
     */
    void RunSpheroidSimulation(EGFRMutationType egfrMutation, double betaAlpha, double betaBeta,
                                double tgfProductionRate, double tgfDiffusionRadius,
                                double tgfDegradationRate,
                                unsigned prodSweepIndex, unsigned diffSweepIndex,
                                unsigned degSweepIndex, unsigned replicateIndex)
    {
        std::string output_dir = "PC9_ProdDiffDegSweep_Hetero_"
            + GetMutationLabel(egfrMutation)
            + "_BetaA" + boost::lexical_cast<std::string>(static_cast<int>(betaAlpha))
            + "B" + boost::lexical_cast<std::string>(static_cast<int>(betaBeta))
            + "_Prod" + boost::lexical_cast<std::string>(prodSweepIndex)
            + "_Diff" + boost::lexical_cast<std::string>(diffSweepIndex)
            + "_Deg" + boost::lexical_cast<std::string>(degSweepIndex)
            + "_rep" + boost::lexical_cast<std::string>(replicateIndex);

        // ---------------------------------------------------------------
        // RESUME LOGIC: skip this run if it already completed successfully
        // ---------------------------------------------------------------
        if (IsRunComplete(output_dir))
        {
            std::cout << "\n=== Prod " << prodSweepIndex
                      << " (rate=" << tgfProductionRate << ")"
                      << ", Diff " << diffSweepIndex
                      << " (radius=" << tgfDiffusionRadius << ")"
                      << ", Deg " << degSweepIndex
                      << " (rate=" << tgfDegradationRate << ")"
                      << ", Replicate " << replicateIndex
                      << ": " << output_dir
                      << " — ALREADY COMPLETED, skipping. ===" << std::endl;
            return;
        }

        std::cout << "\n=== Prod " << prodSweepIndex
                  << " (TGF-a production rate = " << tgfProductionRate << ")"
                  << ", Diff " << diffSweepIndex
                  << " (TGF-a diffusion radius = " << tgfDiffusionRadius << ")"
                  << ", Deg " << degSweepIndex
                  << " (TGF-a degradation rate = " << tgfDegradationRate << ")"
                  << ", Replicate " << replicateIndex
                  << ": " << output_dir << " ===" << std::endl;

        double base_rescue_fraction = GetEGFRBaseRescueFraction(egfrMutation);
        const double g1_min         = 54.0;
        const double g1_extra_hours = 54.0;

        // Simulation parameters — matched to the heterogeneous reference script
        double dt                = 0.005;
        double end_time          = 240.0;
        double sampling_timestep = 1.0 / dt;

        // Spheroid parameters
        double spheroid_radius  = 8.0;
        unsigned mesh_size      = 30;
        double jitter_magnitude = 0.1;

        // Different seed per replicate (kept independent of prod/diff/deg
        // indices so the same replicate uses the same initial cell-cycle
        // phasing, mesh jitter, and per-cell mut_frac draws across the
        // whole grid, isolating the effect of the three swept parameters)
        RandomNumberGenerator::Instance()->Reseed(replicateIndex);

        // Generate CIRCULAR honeycomb mesh
        HoneycombMeshGenerator generator(mesh_size, mesh_size, 0);
        boost::shared_ptr<MutableMesh<2,2>> p_mesh_shared = generator.GetCircularMesh(spheroid_radius);
        MutableMesh<2,2>* p_mesh = p_mesh_shared.get();

        // Add spatial jitter to node positions (BEFORE filtering to a
        // circle, matching the heterogeneous reference script's ordering)
        for (unsigned i = 0; i < p_mesh->GetNumNodes(); i++)
        {
            c_vector<double, 2>& location = p_mesh->GetNode(i)->rGetModifiableLocation();
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
        boost::shared_ptr<AbstractCellProperty> p_stem_type =
            CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>();

        // Create cells
        std::vector<CellPtr> cells;

        double sum_mut_frac_init    = 0.0;
        double sum_sq_mut_frac_init = 0.0;
        double min_mut_frac_init    = 1.0;
        double max_mut_frac_init    = 0.0;

        for (unsigned i = 0; i < location_indices.size(); i++)
        {
            double cell_mut_frac    = SampleBeta(betaAlpha, betaBeta);
            double cell_wt_frac     = 1.0 - cell_mut_frac;
            double cell_g1          = g1_min + g1_extra_hours * cell_wt_frac;
            double cell_rescue      = base_rescue_fraction * cell_mut_frac;
            double cell_total_cycle = cell_g1 + 16.0 + 8.0 + 2.0;

            SimpleOxygenBasedCellCycleModel* p_cycle_model =
                new SimpleOxygenBasedCellCycleModel();
            p_cycle_model->SetDimension(2);

            p_cycle_model->SetQuiescentConcentration(0.08);
            p_cycle_model->SetHypoxicConcentration(0.05);
            p_cycle_model->SetCriticalHypoxicDuration(27.0);

            CellPtr p_cell = CellPtr(new Cell(p_wildtype, p_cycle_model));
            p_cell->SetCellProliferativeType(p_stem_type);

            p_cycle_model->SetStemCellG1Duration(cell_g1);
            p_cycle_model->SetSDuration(16.0);
            p_cycle_model->SetG2Duration(8.0);
            p_cycle_model->SetMDuration(2.0);

            p_cell->InitialiseCellCycleModel();

            double birth_time = cell_total_cycle * RandomNumberGenerator::Instance()->ranf();
            p_cycle_model->SetBirthTime(-birth_time);

            p_cell->GetCellData()->SetItem("oxygen", 1.0);
            p_cell->GetCellData()->SetItem("cell_state", 1.0);
            p_cell->GetCellData()->SetItem("time_became_hypoxic", SimulationTime::Instance()->GetTime());
            p_cell->GetCellData()->SetItem("hif1alpha", 0.0);
            p_cell->GetCellData()->SetItem("tgfalpha", 0.0);
            p_cell->GetCellData()->SetItem("egfr_activation", 0.0);

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

        unsigned n_cells_init = location_indices.size();
        double obs_mean = sum_mut_frac_init / n_cells_init;
        double obs_var  = (sum_sq_mut_frac_init / n_cells_init) - obs_mean * obs_mean;
        double obs_sd   = std::sqrt(std::max(0.0, obs_var));

        std::cout << "  Initial mut_frac distribution: mean = " << obs_mean
                  << "  sd = " << obs_sd
                  << "  min = " << min_mut_frac_init
                  << "  max = " << max_mut_frac_init << std::endl;

        // Create cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<HIF1AlphaWriter>();
        cell_population.AddCellWriter<TGFAlphaWriter>();
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        boost::shared_ptr<CellDataWriter<2,2>> p_cell_data_writer(new CellDataWriter<2,2>());
        p_cell_data_writer->SetOutputDirectory(output_dir);
        cell_population.AddCellWriter(p_cell_data_writer);

        // Create simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_dir);
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(sampling_timestep);
        simulator.SetEndTime(end_time);

        std::cout << "\nSimulation setup:" << std::endl;
        std::cout << "  Duration: " << end_time << " hours" << std::endl;
        std::cout << "  Output directory: " << output_dir << std::endl;

        // Add forces and modifiers
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetMeinekeSpringStiffness(15.0);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(SpheroidOxygenGradientModifier, p_oxygen_modifier,
            (spheroid_radius, 5.0, 0.08, 0.001, 9.0));
        simulator.AddSimulationModifier(p_oxygen_modifier);

        // HIF-1a / TGF-a signaling — production rate, diffusion radius, and
        // degradation rate are the three swept axes here; all other args
        // match the heterogeneous reference script exactly.
        MAKE_PTR_ARGS(HypoxiaSignalingModifier, p_hypoxia_signaling,
            (0.05,                  // Hypoxia threshold
             1.0,                   // Max HIF-1a concentration
             0.01,                  // HIF-1a degradation rate
             tgfProductionRate,     // TGF-a production rate   (SWEPT: axis 1)
             tgfDiffusionRadius,    // TGF-a diffusion radius  (SWEPT: axis 2)
             tgfDegradationRate,    // TGF-a degradation rate  (SWEPT: axis 3)
             1.0,                   // Max TGF-a concentration
             dt,                    // dt
             0.05, 0.08, 0.003, 0.0, 0.0));  // fixed, carried over from reference script
        simulator.AddSimulationModifier(p_hypoxia_signaling);

        MAKE_PTR_ARGS(HypoxicCellStateModifier, p_hypoxic_modifier,
            (0.001, 27.0));
        simulator.AddSimulationModifier(p_hypoxic_modifier);

        MAKE_PTR_ARGS(EGFRSignalingModifier, p_egfr_signaling,
            (0.001, 1.0, 1.0, 0.3, 2, 0.25, dt, 0.0));
        simulator.AddSimulationModifier(p_egfr_signaling);

        // Run simulation
        std::cout << "Starting simulation...\n" << std::endl;
        simulator.Solve();

        // Final statistics
        unsigned total_cells    = 0;
        unsigned normoxic_cells = 0;
        unsigned hypoxic_cells  = 0;
        unsigned dead_cells     = 0;

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
            double oxygen   = cell_iter->GetCellData()->GetItem("oxygen");
            double mut_frac = cell_iter->GetCellData()->GetItem("egfr_mut_fraction");

            final_sum_mut_frac    += mut_frac;
            final_sum_sq_mut_frac += mut_frac * mut_frac;
            if (mut_frac < final_min_mut_frac) final_min_mut_frac = mut_frac;
            if (mut_frac > final_max_mut_frac) final_max_mut_frac = mut_frac;

            if (oxygen > 0.05)
            {
                normoxic_cells++;
                normoxic_sum_mut += mut_frac;
                n_normoxic++;
            }
            else if (oxygen > 0.01)
            {
                hypoxic_cells++;
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

        std::cout << "\n--- Prod " << prodSweepIndex << " (rate=" << tgfProductionRate << ")"
                  << ", Diff " << diffSweepIndex << " (radius=" << tgfDiffusionRadius << ")"
                  << ", Deg " << degSweepIndex << " (rate=" << tgfDegradationRate << ")"
                  << ", Replicate " << replicateIndex << " ---" << std::endl;
        std::cout << "Total cells: " << total_cells << std::endl;
        std::cout << "Normoxic (O2>0.05): " << normoxic_cells
                  << " (" << 100.0*normoxic_cells/total_cells << "%)"
                  << "  mean mut_frac = " << (n_normoxic > 0 ? normoxic_sum_mut/n_normoxic : 0.0) << std::endl;
        std::cout << "Hypoxic (0.01<O2<0.05): " << hypoxic_cells
                  << " (" << 100.0*hypoxic_cells/total_cells << "%)"
                  << "  mean mut_frac = " << (n_hypoxic > 0 ? hypoxic_sum_mut/n_hypoxic : 0.0) << std::endl;
        std::cout << "Necrotic (O2<0.01): " << dead_cells
                  << " (" << 100.0*dead_cells/total_cells << "%)"
                  << "  mean mut_frac = " << (n_necrotic > 0 ? necrotic_sum_mut/n_necrotic : 0.0) << std::endl;
        std::cout << "Final mut_frac distribution: mean = " << f_mean << "  sd = " << f_sd
                  << "  min = " << final_min_mut_frac << "  max = " << final_max_mut_frac << std::endl;

        // ---------------------------------------------------------------
        // RESUME LOGIC: write completion marker only after full success
        // ---------------------------------------------------------------
        MarkRunComplete(output_dir);
        std::cout << "Completion marker written for " << output_dir << std::endl;

        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
    }

    void TestPC9SpheroidTGFaProdDiffDegSweepHeterogeneous()
    {
        std::cout << "\n=== Spheroid 2D — Joint TGF-a Production x Diffusion x Degradation Sweep, "
                  << "Per-Cell EGFR Heterogeneity (Resumable) ===\n" << std::endl;

        // [USER SETTING 1] — EGFR MUTATION TYPE
        EGFRMutationType egfr_mutation = EGFRMutationType::EXON19_DEL;

        // [USER SETTING 2] — BETA DISTRIBUTION SHAPE PARAMETERS
        const double beta_alpha = 5.0;
        const double beta_beta  = 5.0;

        // [USER SETTING 3] — TGF-a PRODUCTION RATES TO SWEEP
        std::vector<double> tgf_production_rates = {0.001, 0.002, 0.008, 0.016};

        // [USER SETTING 4] — TGF-a DIFFUSION RADII TO SWEEP
        std::vector<double> tgf_diffusion_radii = {0.25, 1.0, 2.0, 4.0};

        // [USER SETTING 5] — TGF-a DEGRADATION RATES TO SWEEP
        std::vector<double> tgf_degradation_rates = {0.001, 0.005, 0.01, 0.05};

        // [USER SETTING 6] — REPLICATES PER GRID POINT
        unsigned start_rep = 1;
        unsigned num_reps  = 5;

        unsigned total_grid_points = tgf_production_rates.size()
            * tgf_diffusion_radii.size()
            * tgf_degradation_rates.size();
        unsigned total_runs        = total_grid_points * num_reps;

        std::cout << tgf_production_rates.size() << " production rates x "
                  << tgf_diffusion_radii.size() << " diffusion radii x "
                  << tgf_degradation_rates.size() << " degradation rates = "
                  << total_grid_points << " grid points" << std::endl;
        std::cout << total_grid_points << " grid points x " << num_reps
                  << " replicates = " << total_runs << " total simulations" << std::endl;

        // Nested sweep: production rate (outer) x diffusion radius (middle)
        // x degradation rate (inner) x replicate (innermost). Each run
        // checks its own completion marker before starting, so this loop
        // can simply be re-run in full after an interruption.
        for (unsigned p = 0; p < tgf_production_rates.size(); p++)
        {
            unsigned prod_sweep_index = p + 1;  // 1-based

            for (unsigned d = 0; d < tgf_diffusion_radii.size(); d++)
            {
                unsigned diff_sweep_index = d + 1;  // 1-based

                for (unsigned g = 0; g < tgf_degradation_rates.size(); g++)
                {
                    unsigned deg_sweep_index = g + 1;  // 1-based

                    for (unsigned rep = start_rep; rep < start_rep + num_reps; rep++)
                    {
                        RunSpheroidSimulation(egfr_mutation, beta_alpha, beta_beta,
                                              tgf_production_rates[p], tgf_diffusion_radii[d],
                                              tgf_degradation_rates[g],
                                              prod_sweep_index, diff_sweep_index,
                                              deg_sweep_index, rep);
                    }
                }
            }
        }

        std::cout << "\n=== Joint production x diffusion x degradation sweep (heterogeneous) completed! ===" << std::endl;
        std::cout << total_grid_points << " grid points x " << num_reps
                  << " replicates = " << total_runs << " simulations" << std::endl;
    }
};

#endif /* TESTPC9SPHEROID2D_TGFAPRODDIFFDEGSWEEP_RESUMABLE_HPP_ */