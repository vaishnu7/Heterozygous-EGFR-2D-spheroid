/*
* Created on: 21 Oct 2025
* Last modified: 21 Oct 2025
* 		Author: Vaishnudebi Dutta
* Modified: Added EGFR ratio shift in hypoxia (Step 4)
*/

#include "CellData.hpp"
#include "ApcTwoHitCellMutationState.hpp"

// Oxygen gradient modifier for spheroid cross-section
class SpheroidOxygenGradientModifier : public AbstractCellBasedSimulationModifier<2>
{
private:
    double mInitialSpheroidRadius;     // Initial radius at t=0
    double mOxygenPenetrationDepth;    // How deep oxygen penetrates (~100-200 μm)
    double mMaxOxygenConcentration;    // O2 at surface
    double mMinOxygenConcentration;    // O2 in deep core
    double mHypoxiaOnsetRadius;        // Radius at which hypoxia begins (~10 cell units)
    c_vector<double, 2> mSpheroidCenter;
    
public:
    SpheroidOxygenGradientModifier(double initialRadius = 12.0,
                                   double penetrationDepth = 5.0,
                                   double maxOxygen = 0.11,
                                   double minOxygen = 0.01,
                                   double hypoxiaOnset = 13)
        : AbstractCellBasedSimulationModifier<2>(),
          mInitialSpheroidRadius(initialRadius),
          mOxygenPenetrationDepth(penetrationDepth),
          mMaxOxygenConcentration(maxOxygen),
          mMinOxygenConcentration(minOxygen),
          mHypoxiaOnsetRadius(hypoxiaOnset)
    {
        mSpheroidCenter(0) = 0.0;
        mSpheroidCenter(1) = 0.0;
    }
    
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<2,2>& rCellPopulation)
    {
        // Calculate actual current maximum radius
        double actual_radius = 0.0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            c_vector<double,2> location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
            double radius = norm_2(location - mSpheroidCenter);
            if (radius > actual_radius)
            {
                actual_radius = radius;
            }
        }
        
        // Use actual radius for oxygen calculations (accounts for real growth)
        double current_spheroid_radius = actual_radius;
        
        // Determine if hypoxia should occur based on spheroid size
        // Small spheroids start fully normoxic, hypoxia develops as they grow

          bool has_hypoxic_region = (current_spheroid_radius >= mHypoxiaOnsetRadius);

          // Calculate target core oxygen based on spheroid size (gradual development)
          double target_core_oxygen;
          if (!has_hypoxic_region)
          {
              // Small spheroid: core stays highly oxygenated
              target_core_oxygen = mMaxOxygenConcentration; 
          }
          else
          {
              // Large spheroid: core oxygen decreases gradually as it grows
              double size_beyond_onset = current_spheroid_radius - mHypoxiaOnsetRadius;
              double max_growth = 5.0;  // Core fully hypoxic after growing 5 more units

              double hypoxia_factor = std::min(1.0, size_beyond_onset / max_growth);

              // Gradually transition core oxygen from max to min
              target_core_oxygen = mMaxOxygenConcentration - (mMaxOxygenConcentration - mMinOxygenConcentration) * hypoxia_factor;
          }

          // Update oxygen for each cell based on radial distance from center
          for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
               cell_iter != rCellPopulation.End();
               ++cell_iter)
          {
              c_vector<double,2> location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

              // Calculate radial distance from spheroid center
              double radius = norm_2(location - mSpheroidCenter);

              // Distance from surface (key for oxygen diffusion)
              double distance_from_surface = current_spheroid_radius - radius;

              double oxygen;

              if (distance_from_surface <= 0.0)
              {
                  // At or beyond surface - always max oxygen
                  oxygen = mMaxOxygenConcentration;
              }
              //else if (distance_from_surface >= mOxygenPenetrationDepth)
              //{
                  // Deep core beyond oxygen penetration
                  //oxygen = target_core_oxygen;  // Use gradual core oxygen
              //}
              else
              {
                  // Transition zone: exponential decay from surface to core
                  double decay_rate = 1.0;
                  double normalized_depth = distance_from_surface / mOxygenPenetrationDepth;

                  // Exponential decay: high O2 at surface → target core O2
                  oxygen = target_core_oxygen +
                          (mMaxOxygenConcentration - target_core_oxygen) *
                          exp(-decay_rate * normalized_depth);
               }
            
            // Set oxygen in cell data
            cell_iter->GetCellData()->SetItem("oxygen", oxygen);
            
            // Track cell state for visualization
            if (oxygen > 0.05)
            {
                cell_iter->GetCellData()->SetItem("cell_state", 1.0); // Normoxic/proliferative
            }
            else if (oxygen > 0.01)
            {
                cell_iter->GetCellData()->SetItem("cell_state", 2.0); // Mild hypoxia/quiescent
            }
            //else if (oxygen > 0.015)
            //{
                //cell_iter->GetCellData()->SetItem("cell_state", 3.0); // Severe hypoxia
            //}
            else
            {
                cell_iter->GetCellData()->SetItem("cell_state", 3.0); // Necrotic
            }
        }
    }
    
    virtual void SetupSolve(AbstractCellPopulation<2,2>& rCellPopulation, 
                           std::string outputDirectory)
    {
        UpdateAtEndOfTimeStep(rCellPopulation);
    }
    
    virtual void OutputSimulationModifierParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<InitialSpheroidRadius>" << mInitialSpheroidRadius 
                     << "</InitialSpheroidRadius>\n";
        *rParamsFile << "\t\t\t<OxygenPenetrationDepth>" << mOxygenPenetrationDepth 
                     << "</OxygenPenetrationDepth>\n";
        *rParamsFile << "\t\t\t<HypoxiaOnsetRadius>" << mHypoxiaOnsetRadius 
                     << "</HypoxiaOnsetRadius>\n";
        
        AbstractCellBasedSimulationModifier<2>::OutputSimulationModifierParameters(rParamsFile);
    }
};