/*
 * CellDataWriter.cpp
 * Author: Vaishnudebi Dutta
 * Date: 24 Oct, 2025
 */

#include "CellDataWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "SimulationTime.hpp"
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <filesystem>

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellDataWriter<ELEMENT_DIM, SPACE_DIM>::CellDataWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("cell_data"),
      mOutputDirectory(""),
      mStoredFrictionForce(0.0)
{
    this->mVtkCellDataName = "CellData";
    this->mOutputScalarData = true;
    this->mOutputVectorData = false;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellDataWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(
    CellPtr pCell,
    AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    try
    {
        // Get simulation time information
        double current_time = SimulationTime::Instance()->GetTime();
        unsigned time_step = SimulationTime::Instance()->GetTimeStepsElapsed();

        // Get cell location and index
        unsigned cell_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
        c_vector<double, SPACE_DIM> cell_position = pCellPopulation->GetLocationOfCellCentre(pCell);
        
        // Calculate distance from spheroid center
        double distance_from_center = 0.0;
        for (unsigned i = 0; i < SPACE_DIM; i++)
        {
            distance_from_center += cell_position[i] * cell_position[i];
        }
        distance_from_center = sqrt(distance_from_center);

        // Get cell cycle and age information
        double cell_age = pCell->GetAge();
        double cell_birth_time = pCell->GetBirthTime();

        // Helper lambda to safely get cell data with default value
        auto GetCellDataItem = [pCell](const std::string& key, double default_val) -> double {
            try
            {
                return pCell->GetCellData()->GetItem(key);
            }
            catch (const Exception&)
            {
                return default_val;
            }
        };

        // Retrieve all cell data
        double oxygen = GetCellDataItem("oxygen", 0.0);
        double cell_state = GetCellDataItem("cell_state", 1.0);
        double time_became_hypoxic = GetCellDataItem("time_became_hypoxic", current_time);
        double hif1alpha = GetCellDataItem("hif1alpha", 0.0);
        double tgf_alpha = GetCellDataItem("tgfalpha", 0.0);
        double egfr_activation = GetCellDataItem("egfr_activation", 0.0);
        double distance_to_boundary = GetCellDataItem("distance_to_boundary", 0.0);

        // Get mutation and proliferative type
        std::string mutation_state = pCell->GetMutationState()->GetIdentifier();
        std::string proliferative_type = pCell->GetCellProliferativeType()->GetIdentifier();

        // Get output directory - use environment variable set by Chaste
        const char* chaste_output = std::getenv("CHASTE_TEST_OUTPUT");
        std::string base_dir = (chaste_output != nullptr) ? std::string(chaste_output) : ".";
        
        // Find the most recently modified directory (current simulation)
        std::string output_dir = base_dir;
        try
        {
            auto latest_time = std::filesystem::file_time_type::min();
            if (!mOutputDirectory.empty())
            {
                output_dir = base_dir + "/" + mOutputDirectory;
            }
            else
            {
                for (const auto& entry : std::filesystem::directory_iterator(base_dir))
                {
                    if (entry.is_directory() && entry.path().filename().string() != ".git")
                    {
                        auto mod_time = std::filesystem::last_write_time(entry);
                        if (mod_time > latest_time)
                        {
                            latest_time = mod_time;
                            output_dir = entry.path().string();
                        }
                    }
                }
            }
        }
        catch (const std::exception&)
        {
            // If anything fails, use base directory
        }

        std::string cell_data_dir = output_dir + "/cellData/";

        // Create directory
        try
        {
            if (system(("mkdir -p " + cell_data_dir).c_str()) != 0)
            {
                std::cerr << "ERROR: Could not create directory: " << cell_data_dir << std::endl;
                return;
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "ERROR: Exception creating directory: " << e.what() << std::endl;
            return;
        }

        // Define filename based on time step
        std::stringstream filename;
        filename << "cell_data_time_" << time_step << ".dat";
        std::string filepath = cell_data_dir + filename.str();

        // Check if file already has content (header written)
        bool file_has_content = false;
        std::ifstream check_file(filepath);
        if (check_file.good())
        {
            check_file.seekg(0, std::ios::end);
            file_has_content = (check_file.tellg() > 0);
        }
        check_file.close();

        // Open file in append mode
        std::ofstream cell_data_file(filepath, std::ios::app);
        if (!cell_data_file.is_open())
        {
            std::cerr << "ERROR at timestep " << time_step << ", cell " << cell_index 
                      << ": Could not open file: " << filepath << std::endl;
            return;
        }

        // Write header on first write to this file
        if (!file_has_content && cell_index == 0)
        {
            cell_data_file << "Time\tTimeStep\tCellIndex\tDistanceFromCenter\t"
                           << "PositionX\tPositionY";
            
            if (SPACE_DIM == 3)
                cell_data_file << "\tPositionZ";
            
            cell_data_file << "\tCellAge\tBirthTime\tOxygen\tCellState\tTimeHypoxic\t"
                           << "HIF1Alpha\tTGFAlpha\tEGFRActivation\tDistanceToBoundary\t"
                           << "MutationState\tProliferativeType\n";
        }

        // Write data row for this cell
        cell_data_file << std::fixed << std::setprecision(6)
            << current_time << "\t" << time_step << "\t" << cell_index << "\t"
            << distance_from_center << "\t"
            << cell_position[0] << "\t" << cell_position[1];
        
        if (SPACE_DIM == 3)
            cell_data_file << "\t" << cell_position[2];
        
        cell_data_file << "\t" << cell_age << "\t" << cell_birth_time << "\t"
            << oxygen << "\t" << cell_state << "\t"
            << (current_time - time_became_hypoxic) << "\t"
            << hif1alpha << "\t" << tgf_alpha << "\t"
            << egfr_activation << "\t" << distance_to_boundary << "\t"
            << mutation_state << "\t" << proliferative_type << "\n";

        cell_data_file.close();
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL ERROR in VisitCell: " << e.what() << std::endl;
    }
}

// Explicit template instantiation
template class CellDataWriter<1,1>;
template class CellDataWriter<2,2>;
template class CellDataWriter<3,3>;
