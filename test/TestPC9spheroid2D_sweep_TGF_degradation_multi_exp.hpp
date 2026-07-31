/*
 * TestPC9Spheroid2D.hpp
 *
 * 2D cross-section model of PC9 or general NSCLC tumor spheroid with hypoxic zones
 * Parameter sweep on TGF-alpha degradation rate
 * Models radial oxygen gradient from spheroid surface to hypoxic/necrotic core
 * Incorporates HIF-1a and TGF-a signaling under hypoxia, affecting proliferation
 * Based on off-lattice cell-based modeling framework in Chaste
 *
 * Created on: 24 Oct, 2025
 * Author: Vaishnudebi Dutta
 * FIXED: 16 Feb, 2026
 *
 */

#ifndef TESTPC9SPHEROID2D_HPP_
#define TESTPC9SPHEROID2D_HPP_

#include <cxxtest/TestSuite.h>
#include "CellBasedSimulationArchiver.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"

#include "HoneycombMeshGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "VoronoiDataWriter.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "FakePetscSetup.hpp"

#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"

#include "GeneralisedLinearSpringForce.hpp"
#include "SimpleTargetAreaModifier.hpp"

#include <boost/lexical_cast.hpp>
#include "CellAgesWriter.hpp"
#include "CellLabelWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "OutputFileHandler.hpp"
#include "CellDataWriter.hpp"
#include "CellData.hpp"
#include <iomanip>
#include <iostream>
#include <cmath>
#include <set>

#include "HypoxiaSignalingModifier.hpp"
#include "EGFRSignalingModifier.hpp"
#include "SpheroidOxygenGradientModifier.hpp"
#include "HypoxiaCellStateModifier.hpp"
#include "HIF1AlphaWriter.hpp"
#include "TGFAlphaWriter.hpp"
#include "EGFRActivationWriter.hpp"

class TestPC9Spheroid2D : public AbstractCellBasedTestSuite
{
public:
    void RunSpheroidSimulation(double tgf_degradation_rate, int sweep_index, unsigned replicate_index)
    {
        std::cout << "\n=== Sweep " << sweep_index << ", Replicate " << replicate_index
                  << ": TGF-alpha degradation rate = " << tgf_degradation_rate 
                  << " ===" << std::endl;
        
        // Simulation parameters
        double dt = 0.01; // 36 sec
        double end_time = 240.0;  // 7 days
        double sampling_timestep = 4.0/dt;  // Sample every hour
        
        // Spheroid parameters 
        double spheroid_radius = 8.0;
        unsigned mesh_size = 30;
        
        // Different seed per replicate
        RandomNumberGenerator::Instance()->Reseed(replicate_index);
        
        // Generate CIRCULAR honeycomb mesh
        HoneycombMeshGenerator generator(mesh_size, mesh_size, 0);
        boost::shared_ptr<MutableMesh<2,2>> p_mesh_shared = generator.GetCircularMesh(spheroid_radius);
        MutableMesh<2,2>* p_mesh = p_mesh_shared.get();
        
        std::vector<unsigned> all_indices = generator.GetCellLocationIndices();
        
        c_vector<double, 2> spheroid_center = zero_vector<double>(2);
        
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

        std::set<unsigned> location_set(location_indices.begin(), location_indices.end());
        double jitter_magnitude = 0.1;
        for (unsigned i = 0; i < p_mesh->GetNumNodes(); i++)
        {
            if (location_set.count(i) > 0)
            {
                c_vector<double, 2>& location = p_mesh->GetNode(i)->rGetModifiableLocation();
                double dx = jitter_magnitude * (2.0 * RandomNumberGenerator::Instance()->ranf() - 1.0);
                double dy = jitter_magnitude * (2.0 * RandomNumberGenerator::Instance()->ranf() - 1.0);
                location[0] += dx;
                location[1] += dy;
            }
        }
        
        std::cout << "Created circular mesh with " << location_indices.size() 
                  << " cells within radius" << std::endl;
        
        boost::shared_ptr<AbstractCellProperty> p_wildtype = 
            CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();
        boost::shared_ptr<AbstractCellProperty> p_hypoxic_state = 
            CellPropertyRegistry::Instance()->Get<ApcOneHitCellMutationState>();
        boost::shared_ptr<AbstractCellProperty> p_necrotic_state = 
            CellPropertyRegistry::Instance()->Get<ApcTwoHitCellMutationState>();
        boost::shared_ptr<AbstractCellProperty> p_stem_type = 
            CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>();
        
        std::vector<CellPtr> cells;
        
        CellPropertyRegistry::Instance()->Clear();
        
        for (unsigned i = 0; i < location_indices.size(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_cycle_model = 
                new SimpleOxygenBasedCellCycleModel();
            p_cycle_model->SetDimension(2);
            
            p_cycle_model->SetQuiescentConcentration(0.08);
            p_cycle_model->SetHypoxicConcentration(0.05);
            p_cycle_model->SetCriticalHypoxicDuration(27.0);
            
            CellPtr p_cell = CellPtr(new Cell(p_wildtype, p_cycle_model));
            p_cell->SetCellProliferativeType(p_stem_type);
            
            p_cycle_model->SetStemCellG1Duration(54.0);
            p_cycle_model->SetSDuration(16.0);
            p_cycle_model->SetG2Duration(8.0);
            p_cycle_model->SetMDuration(2.0);

            p_cell->InitialiseCellCycleModel();

            double birth_time = 80.0 * RandomNumberGenerator::Instance()->ranf();
            p_cycle_model->SetBirthTime(-birth_time);
                                
            p_cell->GetCellData()->SetItem("oxygen", 1.0);
            p_cell->GetCellData()->SetItem("cell_state", 1.0);
            p_cell->GetCellData()->SetItem("time_became_hypoxic", SimulationTime::Instance()->GetTime());
            p_cell->GetCellData()->SetItem("HIF1Alpha", 0.0);
            p_cell->GetCellData()->SetItem("TGFAlpha", 0.0);
            p_cell->GetCellData()->SetItem("EGFR_activation", 0.0);

            cells.push_back(p_cell);
        }
        
        std::cout << "Created " << cells.size() << " PC9 cells" << std::endl;
        
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);
        
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellAncestorWriter>();
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        std::string output_dir = "PC9_Spheroid2D_TGFDeg_MultiExp_CorrectedO2_EGFR5050_" 
            + boost::lexical_cast<std::string>(sweep_index)
            + "_rep_"
            + boost::lexical_cast<std::string>(replicate_index);
        
        boost::shared_ptr<CellDataWriter<2,2>> p_cell_data_writer(new CellDataWriter<2,2>());
        p_cell_data_writer->SetOutputDirectory(output_dir);
        cell_population.AddCellWriter(p_cell_data_writer);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_dir);
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(sampling_timestep);
        simulator.SetEndTime(end_time);
        
        std::cout << "Simulation parameters:" << std::endl;
        std::cout << "  Duration: " << end_time << " hours" << std::endl;
        std::cout << "  Output directory: " << output_dir << std::endl;
        
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetMeinekeSpringStiffness(15.0);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);
        
        MAKE_PTR_ARGS(SpheroidOxygenGradientModifier, p_oxygen_modifier, 
             (spheroid_radius,
              5.0,
              0.08,
              0.001,
              9.0));
        simulator.AddSimulationModifier(p_oxygen_modifier);
        
        // Add HIF-1a and TGF-a signaling
        // SWEPT PARAMETER: TGF-alpha degradation rate
        MAKE_PTR_ARGS(HypoxiaSignalingModifier, p_hypoxia_signaling, 
            (0.05,                  // Hypoxia threshold
             1.0,                   // Max HIF-1a concentration
             0.01,                  // HIF-1a degradation rate
             0.001,                 // TGF-a production rate (fixed from production sweep)
             0.5,                   // TGF-a diffusion radius
             tgf_degradation_rate,  // TGF-a degradation rate (SWEPT PARAMETER)
             1.0,                   // Max TGF-a concentration
             dt,                    // dt
             0.05,                  // HIF basal synthesis rate
             0.08,                  // Max O2-dependent HIF degradation
             0.003));               // O2-independent HIF degradation rate
        simulator.AddSimulationModifier(p_hypoxia_signaling);

        MAKE_PTR_ARGS(EGFRSignalingModifier, p_egfr_signaling, 
            (0.001,
             1.0,
             1.0));
        simulator.AddSimulationModifier(p_egfr_signaling);

        MAKE_PTR_ARGS(HypoxicCellStateModifier, p_hypoxic_modifier, 
            (0.001,
             27.0));
        simulator.AddSimulationModifier(p_hypoxic_modifier);
        
        std::cout << "Starting simulation...\n" << std::endl;
        
        simulator.Solve();
        
        // Final statistics
        unsigned total_cells = 0;
        unsigned proliferative_cells = 0;
        unsigned quiescent_cells = 0;
        unsigned dead_cells = 0;
        
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            total_cells++;
            double oxygen = cell_iter->GetCellData()->GetItem("oxygen");
            
            if (oxygen > 0.05)
                proliferative_cells++;
            else if (oxygen > 0.01)
                quiescent_cells++;
            else
                dead_cells++;
        }
        
        std::cout << "\n--- Sweep " << sweep_index << ", Replicate " << replicate_index 
                  << " (TGF Degradation Rate " << tgf_degradation_rate << ") ---" << std::endl;
        std::cout << "Total cells: " << total_cells << std::endl;
        std::cout << "Normoxic (O2>0.05): " << proliferative_cells 
                  << " (" << 100.0*proliferative_cells/total_cells << "%)" << std::endl;
        std::cout << "Hypoxic (0.01<O2<0.05): " << quiescent_cells 
                  << " (" << 100.0*quiescent_cells/total_cells << "%)" << std::endl;
        std::cout << "Severely hypoxic/necrotic (O2<0.01): " << dead_cells 
                  << " (" << 100.0*dead_cells/total_cells << "%)" << std::endl;
        
        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
    }

    void TestPC9SpheroidCrossSection()
    {
        std::cout << "\n=== PC9 Spheroid 2D Cross-Section Parameter Sweep: TGF-alpha Degradation Rate ===\n" << std::endl;
        
        // Sweep range covers: fast local decay → slow persistent signaling
        // 0.001 = very slow (TGF-a persists, broad field effect)
        // 0.005 = moderate-slow
        // 0.01  = moderate (baseline from production sweep)
        // 0.05  = fast (TGF-a degrades quickly, short-range paracrine)
        std::vector<double> tgf_degradation_rates = {0.001, 0.005, 0.01, 0.05};
        
        unsigned start_rep = 1;
        unsigned num_reps = 5;

        // --- RESUME SETTINGS ---
        unsigned resume_sweep_index = 1;   // 1-based sweep index to resume from
        unsigned resume_rep = 1;           // replicate to resume from within that sweep

        std::cout << "Resuming from Sweep " << resume_sweep_index 
                << ", Replicate " << resume_rep << std::endl;
        std::cout << "Total runs: " << tgf_degradation_rates.size() * num_reps << std::endl;

        for (unsigned i = 0; i < tgf_degradation_rates.size(); i++)
        {
            unsigned sweep_index = i + 1;  // 1-based
            
            if (sweep_index < resume_sweep_index) continue;
            
            for (unsigned rep = start_rep; rep < start_rep + num_reps; rep++)
            {
                if (sweep_index == resume_sweep_index && rep < resume_rep) continue;
                
                RunSpheroidSimulation(tgf_degradation_rates[i], sweep_index, rep);
            }
        }
        
        std::cout << "\n=== Parameter sweep completed! ===" << std::endl;
        std::cout << "4 TGF degradation rates x 5 replicates = 20 simulations" << std::endl;
    }
};

#endif /* TESTPC9SPHEROID2D_HPP_ */