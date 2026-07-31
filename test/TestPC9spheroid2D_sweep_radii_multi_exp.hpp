/*
 * TestPC9Spheroid2D.hpp
 *
 * 2D cross-section model of PC9 or general NSCLC tumor spheroid with hypoxic zones
 * Parameter sweep on initial spheroid radius (4 different values)
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
//#include "SpheroidRadiusWriter.hpp"                     // Output spheroid radius data
//#include "SpheroidBoundaryModifier.hpp"                 // Enforce spheroid boundary constraints

class TestPC9Spheroid2D : public AbstractCellBasedTestSuite
{
public:
    void RunSpheroidSimulation(double spheroid_radius, int sweep_index, unsigned replicate_index)
    {
        std::cout << "\n=== Sweep " << sweep_index << ", Replicate " << replicate_index
                  << ": Spheroid radius = " << spheroid_radius 
                  << " cell units (~" << spheroid_radius * 15 << " um) ===" << std::endl;
        
        // Simulation parameters
        double dt = 0.01; // 36 sec
        double end_time = 240.0;  // 7days 
        double sampling_timestep = 4.0/dt;  // Sample every hour
        
        // Spheroid parameters 
        unsigned mesh_size = 30;
        
        // Different seed per replicate
        RandomNumberGenerator::Instance()->Reseed(replicate_index);
        
        // Generate CIRCULAR honeycomb mesh
        HoneycombMeshGenerator generator(mesh_size, mesh_size, 0);
        boost::shared_ptr<MutableMesh<2,2>> p_mesh_shared = generator.GetCircularMesh(spheroid_radius);
        MutableMesh<2,2>* p_mesh = p_mesh_shared.get();
        
        // Get all node indices BEFORE jitter (these are the nodes in the circular mesh)
        std::vector<unsigned> all_indices = generator.GetCellLocationIndices();
        
        // Spheroid center at origin
        c_vector<double, 2> spheroid_center = zero_vector<double>(2);
        
        // Filter nodes to create circular spheroid BEFORE jitter
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

        // Apply spatial jitter AFTER filtering — only to nodes that will have cells
        // Build a set for fast lookup
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
        
        for (unsigned i = 0; i < location_indices.size(); i++)
        {
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
            
            p_cycle_model->SetStemCellG1Duration(54.0);   // G1
            p_cycle_model->SetSDuration(16.0);            // S
            p_cycle_model->SetG2Duration(8.0);            // G2
            p_cycle_model->SetMDuration(2.0);             // M

            // Initialize cell cycle model first
            p_cell->InitialiseCellCycleModel();

            // Randomize birth time to avoid synchronization
            double birth_time = 80.0 * RandomNumberGenerator::Instance()->ranf();
            p_cycle_model->SetBirthTime(-birth_time);
                                
            // Initial oxygen (will be updated by modifier immediately)
            p_cell->GetCellData()->SetItem("oxygen", 1.0);
            p_cell->GetCellData()->SetItem("cell_state", 1.0);
            p_cell->GetCellData()->SetItem("time_became_hypoxic", SimulationTime::Instance()->GetTime());
            p_cell->GetCellData()->SetItem("HIF1Alpha", 0.0);
            p_cell->GetCellData()->SetItem("TGFAlpha", 0.0);
            p_cell->GetCellData()->SetItem("EGFR_activation", 0.0);

            cells.push_back(p_cell);
        }
        
        std::cout << "Created " << cells.size() << " PC9 cells" << std::endl;
        
        // Create cell population (no ghost nodes for spheroid)
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);
        
        // Add cell writers - these generate .dat files
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellAncestorWriter>();
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.AddPopulationWriter<VoronoiDataWriter>();
        //MAKE_PTR(SpheroidRadiusWriter, p_radius_writer);
        //cell_population.AddPopulationWriter(p_radius_writer);

        // Create output directory for this run
        std::string output_dir = "PC9_Spheroid2D_Radius_" 
            + boost::lexical_cast<std::string>(spheroid_radius)
            + "_rep_"
            + boost::lexical_cast<std::string>(replicate_index);
        
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
        
        // Literature: oxygen penetrates ~100-200 um = ~5-10 cell diameters
        MAKE_PTR_ARGS(SpheroidOxygenGradientModifier, p_oxygen_modifier, 
             (spheroid_radius,  // Initial spheroid radius
              5.0,              // Oxygen penetration depth (~100 um)
              0.08,             // Max O2 at surface (normoxic)
              0.001,             // Min O2 in necrotic core
              9.0));           // Hypoxia onset at radius ~10 units (~200 um)
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
        
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            total_cells++;
            double oxygen = cell_iter->GetCellData()->GetItem("oxygen");
            
            if (oxygen > 0.08)
                proliferative_cells++;
            else if (oxygen > 0.02)
                quiescent_cells++;
            else
                dead_cells++;
        }
        
        std::cout << "\n--- Sweep " << sweep_index << ", Replicate " << replicate_index 
                  << " (Radius " << spheroid_radius << ") ---" << std::endl;
        std::cout << "Total cells: " << total_cells << std::endl;
        std::cout << "Normoxic (O2>0.08): " << proliferative_cells 
                  << " (" << 100.0*proliferative_cells/total_cells << "%)" << std::endl;
        std::cout << "Hypoxic (0.02<O2<0.08): " << quiescent_cells 
                  << " (" << 100.0*quiescent_cells/total_cells << "%)" << std::endl;
        std::cout << "Severely hypoxic/necrotic (O2<0.02): " << dead_cells 
                  << " (" << 100.0*dead_cells/total_cells << "%)" << std::endl;
        
        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
    }

    void TestPC9SpheroidCrossSection()
    {
        std::cout << "\n=== PC9 Spheroid 2D Cross-Section: Radius Parameter Sweep ===\n" << std::endl;
        
        // Define spheroid radii to sweep
        std::vector<double> spheroid_radii = {4.0, 6.0, 10.0, 12.0};
        
        // Replicate settings
        unsigned start_rep = 1;
        unsigned num_reps = 5;

        std::cout << "Total runs: " << spheroid_radii.size() * num_reps << std::endl;

        // Run simulation for each parameter x replicate
        for (unsigned i = 0; i < spheroid_radii.size(); i++)
        {
            for (unsigned rep = start_rep; rep < start_rep + num_reps; rep++)
            {
                RunSpheroidSimulation(spheroid_radii[i], i + 1, rep);
            }
        }
        
        std::cout << "\n=== Parameter sweep completed! ===" << std::endl;
        std::cout << "4 radii x 5 replicates = 20 simulations" << std::endl;
    }
};

#endif /* TESTPC9SPHEROID2D_HPP_ */