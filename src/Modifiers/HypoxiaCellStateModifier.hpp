/*
* Created on: 21 Oct 2025
* Last modified: 21 Oct 2025
* 		Author: Vaishnudebi Dutta
* Modified: Added EGFR ratio shift in hypoxia (Step 4)
*/

#include "CellData.hpp"                                 // Generic cell property storage
#include "WildTypeCellMutationState.hpp"                // Normal/healthy cell state marker
#include "ApcOneHitCellMutationState.hpp"               // Hypoxic cell state marker
#include "ApcTwoHitCellMutationState.hpp"               // Necrotic/dead cell state marker
#include "SimpleOxygenBasedCellCycleModel.hpp"          // Oxygen-dependent cell cycle model

// Hypoxic cell state modifier - marks cells as necrotic without removing them
class HypoxicCellStateModifier : public AbstractCellBasedSimulationModifier<2>
{
private:
    double mHypoxicThreshold;
    double mHypoxicDuration;
    
public:
    HypoxicCellStateModifier(double hypoxicThreshold = 0.02,
                            double hypoxicDuration = 48.0)
        : AbstractCellBasedSimulationModifier<2>(),
        mHypoxicThreshold(hypoxicThreshold),
        mHypoxicDuration(hypoxicDuration)
    {
    }
    
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<2,2>& rCellPopulation)
    {
        // Get mutation states
        boost::shared_ptr<AbstractCellProperty> p_wildtype = 
            CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();
        boost::shared_ptr<AbstractCellProperty> p_hypoxic = 
            CellPropertyRegistry::Instance()->Get<ApcOneHitCellMutationState>();
        boost::shared_ptr<AbstractCellProperty> p_necrotic = 
            CellPropertyRegistry::Instance()->Get<ApcTwoHitCellMutationState>();
        
        for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
            cell_iter != rCellPopulation.End();
            ++cell_iter)
        {
            double oxygen = cell_iter->GetCellData()->GetItem("oxygen");
            
            // Check if cell has been severely hypoxic long enough to become necrotic
            if (oxygen < mHypoxicThreshold)
            {
                double time_since_hypoxic = SimulationTime::Instance()->GetTime() - 
                                        cell_iter->GetCellData()->GetItem("time_became_hypoxic");
                
                if (time_since_hypoxic > mHypoxicDuration)
                {
                    // Mark as necrotic (don't kill, just change state)
                    cell_iter->SetMutationState(p_necrotic);
                    //cell_iter->GetCellData()->SetItem("cell_state", 4.0); // Necrotic state
                }
                else
                {
                    // Mark as hypoxic but still viable
                    cell_iter->SetMutationState(p_hypoxic);
                }
            }
            else if (oxygen >= 0.02 && oxygen < 0.08)
            {
                // Hypoxic but not severely
                cell_iter->SetMutationState(p_hypoxic);
                // Reset timer if oxygen improved
                cell_iter->GetCellData()->SetItem("time_became_hypoxic", 
                                                SimulationTime::Instance()->GetTime());
            }
            else
            {
                // Normoxic - back to wildtype
                cell_iter->SetMutationState(p_wildtype);
                // Reset timer
                cell_iter->GetCellData()->SetItem("time_became_hypoxic", 
                                                SimulationTime::Instance()->GetTime());
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
        *rParamsFile << "\t\t\t<HypoxicThreshold>" << mHypoxicThreshold 
                    << "</HypoxicThreshold>\n";
        *rParamsFile << "\t\t\t<HypoxicDuration>" << mHypoxicDuration 
                    << "</HypoxicDuration>\n";
        
        AbstractCellBasedSimulationModifier<2>::OutputSimulationModifierParameters(rParamsFile);
    }
};