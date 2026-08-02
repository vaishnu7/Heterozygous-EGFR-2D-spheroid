/*
 * TestPC9Spheroid2D_Heterogeneous_Resumable.hpp
 *
 * 2D cross-section model of cancer tumor spheroid with hypoxic zones.
 * Models radial oxygen gradient from spheroid surface to hypoxic/necrotic core.
 * Incorporates HIF and TGF-a signaling under hypoxia, affecting proliferation.
 * Based on off-lattice cell-based modeling framework in Chaste.
 *
 * Per-cell intratumoral heterogeneity in the EGFR mut:wt ratio (Beta-sampled).
 *
 * RESUME LOGIC: same pattern as TestPC9Spheroid2D.hpp —
 *   - a ".run_complete" marker file is written inside each replicate's own
 *     output directory only after that replicate finishes successfully
 *   - before running a replicate, check for that marker; if present, skip it
 *
 * Original author:  Vaishnudebi Dutta
 * Created:           24 Oct 2025
 * Edited:         19 Jun 2026
 */

#ifndef TESTPC9SPHEROID2D_MULTI_EXPERI_HETEROGENOUS_RESUMABLE_HPP_
#define TESTPC9SPHEROID2D_MULTI_EXPERI_HETEROGENOUS_RESUMABLE_HPP_

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
// EGFRMutationType
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

class TestPC9spheroid2D_multi_experi_heterogenous_resumable : public AbstractCellBasedTestSuite
{
private:

    /**
     * Returns the path to the completion marker file for a given replicate.
     * The marker is written inside the run's own output directory so it lives
     * alongside the simulation results and is never confused with another run.
     */
    std::string GetCompletionMarkerPath(const std::string& rOutputDir)
    {
        std::string output_dir = OutputFileHandler::GetChasteTestOutputDirectory()
            + rOutputDir + "/";
        return output_dir + ".run_complete";
    }

    /** Returns true if the given replicate already finished successfully. */
    bool IsRunComplete(const std::string& rOutputDir)
    {
        std::ifstream marker(GetCompletionMarkerPath(rOutputDir));
        return marker.good();
    }

    /** Writes the completion marker after a successful run. */
    void MarkRunComplete(const std::string& rOutputDir)
    {
        std::ofstream marker(GetCompletionMarkerPath(rOutputDir));
        marker << "completed" << std::endl;
        marker.close();
    }

public:
    void RunSpheroidSimulation(EGFRMutationType egfrMutation, double betaAlpha, double betaBeta,
                                unsigned replicateIndex)
    {
        std::string output_dir = "PC9_Hetero_"
            + GetMutationLabel(egfrMutation)
            + "_BetaA" + boost::lexical_cast<std::string>(static_cast<int>(betaAlpha))
            + "B" + boost::lexical_cast<std::string>(static_cast<int>(betaBeta))
            + "_rep" + boost::lexical_cast<std::string>(replicateIndex);

        // ---------------------------------------------------------------
        // RESUME LOGIC: skip this run if it already completed successfully
        // ---------------------------------------------------------------
        if (IsRunComplete(output_dir))
        {
            std::cout << "\n=== Replicate " << replicateIndex
                      << ": " << output_dir
                      << " — ALREADY COMPLETED, skipping. ===" << std::endl;
            return;
        }

        std::cout << "\n=== Replicate " << replicateIndex
                  << ": " << output_dir << " ===" << std::endl;

        double base_rescue_fraction = GetEGFRBaseRescueFraction(egfrMutation);
        const double g1_min         = 54.0;
        const double g1_extra_hours = 54.0;

        // Simulation parameters
        double dt                = 0.005;
        double end_time          = 2.0; // please change this value to something reasonable like 240 or 450 as per your need
        double sampling_timestep = 1.0 / dt;

        // Spheroid parameters
        double spheroid_radius = 8.0;
        unsigned mesh_size     = 30;
        double jitter_magnitude = 0.1;

        // Different seed per replicate
        RandomNumberGenerator::Instance()->Reseed(replicateIndex);

        // Generate CIRCULAR honeycomb mesh
        HoneycombMeshGenerator generator(mesh_size, mesh_size, 0);
        boost::shared_ptr<MutableMesh<2,2>> p_mesh_shared = generator.GetCircularMesh(spheroid_radius);
        MutableMesh<2,2>* p_mesh = p_mesh_shared.get();

        // Add spatial jitter to node positions
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

        double sum_mut_frac_init = 0.0;
        double sum_sq_mut_frac_init = 0.0;
        double min_mut_frac_init = 1.0;
        double max_mut_frac_init = 0.0;

        for (unsigned i = 0; i < location_indices.size(); i++)
        {
            double cell_mut_frac = SampleBeta(betaAlpha, betaBeta);
            double cell_wt_frac  = 1.0 - cell_mut_frac;
            double cell_g1 = g1_min + g1_extra_hours * cell_wt_frac;
            double cell_rescue = base_rescue_fraction * cell_mut_frac;
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

        MAKE_PTR_ARGS(HypoxiaSignalingModifier, p_hypoxia_signaling,
            (0.05, 1.0, 0.01, 0.004, 0.5, 0.001, 1.0, dt, 0.05, 0.08, 0.003,
             0.0, 0.0));
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

        std::cout << "\n--- Replicate " << replicateIndex << " ---" << std::endl;
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
        std::cout << "Completion marker written for Replicate " << replicateIndex << std::endl;

        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
    }

    void TestPC9SpheroidHeterogeneousMultipleRuns()
    {
        std::cout << "\n=== Spheroid 2D — Per-Cell EGFR Heterogeneity (Multi-Run, Resumable) ===\n" << std::endl;

        // [USER SETTING 1] — EGFR MUTATION TYPE
        EGFRMutationType egfr_mutation = EGFRMutationType::L858R;

        // [USER SETTING 2] — BETA DISTRIBUTION SHAPE PARAMETERS
        const double beta_alpha = 9.0;
        const double beta_beta  = 1.0;

        unsigned start_rep = 1;
        unsigned num_reps  = 5;

        for (unsigned rep = start_rep; rep < start_rep + num_reps; rep++)
        {
            RunSpheroidSimulation(egfr_mutation, beta_alpha, beta_beta, rep);
        }
    }
};

#endif /* TESTPC9SPHEROID2D_MULTI_EXPERI_HETEROGENOUS_RESUMABLE_HPP_ */
