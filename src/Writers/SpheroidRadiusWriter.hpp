// File: SpheroidRadiusWriter.hpp
#ifndef SPHEROIDRADIUSWRITER_HPP_
#define SPHEROIDRADIUSWRITER_HPP_

#include "AbstractCellPopulationWriter.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include <cmath>

class SpheroidRadiusWriter : public AbstractCellPopulationWriter<2, 2>
{
private:
    double CalculateSpheroidRadius(AbstractCellPopulation<2>* pCellPopulation,
                                   double& meanRadius)
    {
        c_vector<double, 2> center = zero_vector<double>(2);
        double max_radius = 0.0;
        meanRadius = 0.0;
        unsigned num_cells = 0;

        for (AbstractCellPopulation<2>::Iterator cell_iter = pCellPopulation->Begin();
             cell_iter != pCellPopulation->End();
             ++cell_iter)
        {
            c_vector<double, 2> location = pCellPopulation->GetLocationOfCellCentre(*cell_iter);
            double radius = norm_2(location - center);
            
            max_radius = std::max(max_radius, radius);
            meanRadius += radius;
            num_cells++;
        }

        if (num_cells > 0)
        {
            meanRadius /= num_cells;
        }
        
        return max_radius;
    }

public:
    SpheroidRadiusWriter() : AbstractCellPopulationWriter<2, 2>("spheroid_radius.txt")
    {
    }

    void Visit(MeshBasedCellPopulation<2>* pCellPopulation)
    {
        double mean_radius = 0.0;
        double max_radius = CalculateSpheroidRadius(pCellPopulation, mean_radius);
        
        double time = SimulationTime::Instance()->GetTime();
        *this->mpOutStream << time << "\t" << max_radius << "\t" << mean_radius << std::endl;
    }

    void Visit(NodeBasedCellPopulation<2>* pCellPopulation)
    {
        double mean_radius = 0.0;
        double max_radius = CalculateSpheroidRadius(pCellPopulation, mean_radius);
        
        double time = SimulationTime::Instance()->GetTime();
        *this->mpOutStream << time << "\t" << max_radius << "\t" << mean_radius << std::endl;
    }

    void Visit(VertexBasedCellPopulation<2>* pCellPopulation)
    {
        double mean_radius = 0.0;
        double max_radius = CalculateSpheroidRadius(pCellPopulation, mean_radius);
        
        double time = SimulationTime::Instance()->GetTime();
        *this->mpOutStream << time << "\t" << max_radius << "\t" << mean_radius << std::endl;
    }

    void Visit(ImmersedBoundaryCellPopulation<2>* pCellPopulation)
    {
        double mean_radius = 0.0;
        double max_radius = CalculateSpheroidRadius(pCellPopulation, mean_radius);
        
        double time = SimulationTime::Instance()->GetTime();
        *this->mpOutStream << time << "\t" << max_radius << "\t" << mean_radius << std::endl;
    }

    void Visit(PottsBasedCellPopulation<2>* pCellPopulation)
    {
        double mean_radius = 0.0;
        double max_radius = CalculateSpheroidRadius(pCellPopulation, mean_radius);
        
        double time = SimulationTime::Instance()->GetTime();
        *this->mpOutStream << time << "\t" << max_radius << "\t" << mean_radius << std::endl;
    }

    void Visit(CaBasedCellPopulation<2>* pCellPopulation)
    {
        double mean_radius = 0.0;
        double max_radius = CalculateSpheroidRadius(pCellPopulation, mean_radius);
        
        double time = SimulationTime::Instance()->GetTime();
        *this->mpOutStream << time << "\t" << max_radius << "\t" << mean_radius << std::endl;
    }
};

#endif /* SPHEROIDRADIUSWRITER_HPP_ */