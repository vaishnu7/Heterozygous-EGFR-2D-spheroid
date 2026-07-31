/*
* Created on: 21 Oct 2025
* Last modified: 16 Jan 2026
* Author: Vaishnudebi Dutta
* Cleaned version: Removed EGFR ratio shift complexity, keeping constant 90:10 ratio
*/

#include "CellData.hpp"
#include "ApcTwoHitCellMutationState.hpp"

// ============================================================================
// Simulates HIF-1a stabilisation in hypoxic regions and TGF-a secretion
class HypoxiaSignalingModifier : public AbstractCellBasedSimulationModifier<2>
{
private:
    double mHypoxiaThreshold;           // O2 threshold for HIF-1a activation (~0.08)
    double mMaxHIF1AlphaProduction;     // Max HIF-1a concentration
    double mHIF1AlphaDegradationRate;   // Degradation rate in normoxic conditions
    double mTGFAlphaProductionRate;     // How much TGF-a per unit of HIF-1a
    double mTGFAlphaDiffusionRadius;    // How far TGF-a diffuses to neighbors (~2-3 cell units)
    double mTGFAlphaDegradationRate;    // TGF-a degradation/internalization
    double mMaxTGFAlphaConcentration;   // Max TGF-a concentration 
    double mDt;                         // Simulation timestep
    double mhifBasalSynthesisRate;      // HIF basal synthesis rate
    double mk_deg_O2_max;               // Max O2 dependent degradation rate of HIF
    double mk_deg_const;                // O2 independent degradation rate of HIF
    double MUTANT_EGFR_FRACTION;  //  mutant EGFR
    double WT_EGFR_FRACTION;      //  WT EGFR
    
public:
    HypoxiaSignalingModifier(double hypoxiaThreshold = 0.08,
                          double maxHIF1Alpha = 1.0,
                          double hif1AlphaDegradation = 0.1,
                          double tgfAlphaProductionRate = 0.3,
                          double tgfAlphaDiffusionRadius = 2.5,
                          double tgfAlphaDegradation = 0.05,
                          double maxTGFAlphaConcentration = 1.0,
                          double dt = 0.005,
                          double hifBasalSynthesisRate = 0.05,
                          double k_deg_O2_max = 0.08,
                          double k_deg_const = 0.003,
                          double mutant_egfr_fraction = 0.9,
                          double wt_egfr_fraction = 0.1)

        : AbstractCellBasedSimulationModifier<2>(),
          mHypoxiaThreshold(hypoxiaThreshold),
          mMaxHIF1AlphaProduction(maxHIF1Alpha),
          mHIF1AlphaDegradationRate(hif1AlphaDegradation),
          mTGFAlphaProductionRate(tgfAlphaProductionRate),
          mTGFAlphaDiffusionRadius(tgfAlphaDiffusionRadius),
          mTGFAlphaDegradationRate(tgfAlphaDegradation),
          mMaxTGFAlphaConcentration(maxTGFAlphaConcentration),
          mDt(dt),
          mhifBasalSynthesisRate(hifBasalSynthesisRate),
          mk_deg_O2_max(k_deg_O2_max),
          mk_deg_const(k_deg_const),
            MUTANT_EGFR_FRACTION(mutant_egfr_fraction),
            WT_EGFR_FRACTION(wt_egfr_fraction)
    {
    }
    
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<2,2>& rCellPopulation)
    {
        // Step 1: Update HIF-1a based on oxygen levels
        for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            double oxygen = cell_iter->GetCellData()->GetItem("oxygen");
            double current_hif1alpha = 0.0;
            
            if (cell_iter->GetCellData()->HasItem("hif1alpha"))
            {
                current_hif1alpha = cell_iter->GetCellData()->GetItem("hif1alpha");
            }
            
            double new_hif1alpha = current_hif1alpha;

            if (oxygen >= mHypoxiaThreshold)
            {
                // Normoxic: HIF-1a degrades (VHL-mediated degradation)
                new_hif1alpha = current_hif1alpha * (1.0 - mHIF1AlphaDegradationRate);
            }
            else if (oxygen > 0.01 && oxygen < mHypoxiaThreshold)
            {
                // Hypoxic range (0.01-0.05): HIF-1a accumulation
                // d[HIF-1a]/dt = k_0 - (k_deg_O2 * [O2]/(K_m + [O2]) + k_deg_const) * [HIF-1a]
                
                //double k_0 = mhifBasalSynthesisRate;                      // Basal synthesis rate
                double k_deg_O2_max = mk_deg_O2_max;             // Max oxygen-dependent degradation
                double K_m = 0.03;                      // Michaelis constant (% O2)
                double k_deg_const = mk_deg_const;             // Oxygen-independent degradation
                
                // Synthesis rate
                double synthesis_rate = mhifBasalSynthesisRate;
                
                // Oxygen-dependent degradation rate 
                double k_deg_O2 = k_deg_O2_max * (oxygen / (K_m + oxygen));
                
                // Total degradation
                double degradation_rate = (k_deg_O2 + k_deg_const) * current_hif1alpha;
                
                // Update: new = current + (synthesis - degradation) * dt
                new_hif1alpha = current_hif1alpha + (synthesis_rate - degradation_rate) * mDt;
                new_hif1alpha = std::max(0.0, new_hif1alpha);
            }
            else  
            {
                // Severe hypoxia (=< 1%): cells too quiescent, degradation dominates
                new_hif1alpha = current_hif1alpha * (1.0 - mHIF1AlphaDegradationRate);
                new_hif1alpha = std::max(0.0, new_hif1alpha);
            }
            
            cell_iter->GetCellData()->SetItem("hif1alpha", new_hif1alpha);
        }
        
        // Step 2: TGF-a production driven by HIF-1a
        for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {   
            bool cell_is_necrotic = cell_iter->GetMutationState()->IsType<ApcTwoHitCellMutationState>();

            // Skip TGF-alpha production for necrotic cells
            if (cell_is_necrotic)
            {
                cell_iter->GetCellData()->SetItem("tgfalpha", 0.0);
                continue;
            }
            
            double hif1alpha = cell_iter->GetCellData()->GetItem("hif1alpha");
            double current_tgfalpha = 0.0;
            
            if (cell_iter->GetCellData()->HasItem("tgfalpha"))
            {
                current_tgfalpha = cell_iter->GetCellData()->GetItem("tgfalpha");
            }
            
            // TGF-a production is proportional to HIF-1a levels
            double tgfalpha_production = hif1alpha * mTGFAlphaProductionRate;
            double new_tgfalpha = current_tgfalpha + tgfalpha_production -
                                (current_tgfalpha * mTGFAlphaDegradationRate);

            // Cap TGF-a at maximum concentration 
            new_tgfalpha = std::max(0.0, std::min(mMaxTGFAlphaConcentration, new_tgfalpha));

            cell_iter->GetCellData()->SetItem("tgfalpha", new_tgfalpha);
        }
        
        // Step 3: TGF-a diffusion to neighboring cells
        DiffuseTGFAlpha(rCellPopulation);
    }
    
    void DiffuseTGFAlpha(AbstractCellPopulation<2,2>& rCellPopulation)
    {
        auto p_mesh_pop = dynamic_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);
        MutableMesh<2,2>& mesh = p_mesh_pop->rGetMesh();

        // Build node_index → CellPtr map (O(N), done once per call)
        std::map<unsigned, CellPtr> node_to_cell;
        for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
            cell_iter != rCellPopulation.End();
            ++cell_iter)
        {
            unsigned node_idx = p_mesh_pop->GetLocationIndexUsingCell(*cell_iter);
            node_to_cell[node_idx] = *cell_iter;
        }

        // Accumulate transfers into a separate map (node_index → transfer amount)
        // Separating accumulation from application prevents a cell that receives
        // transfer in one pass from incorrectly re-propagating it in the same step.
        std::map<unsigned, double> tgfalpha_transfer;

        for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
            cell_iter != rCellPopulation.End();
            ++cell_iter)
        {
            double cell_tgfalpha = cell_iter->GetCellData()->GetItem("tgfalpha");
            if (cell_tgfalpha < 1e-6) continue;

            c_vector<double, 2> cell_location =
                rCellPopulation.GetLocationOfCellCentre(*cell_iter);
            unsigned node_index = p_mesh_pop->GetLocationIndexUsingCell(*cell_iter);

            // Collect unique neighbour nodes from mesh elements (~6 in honeycomb)
            const std::set<unsigned>& element_indices =
                mesh.GetNode(node_index)->rGetContainingElementIndices();

            std::set<unsigned> neighbour_nodes;
            for (unsigned elem_idx : element_indices)
            {
                Element<2,2>* p_elem = mesh.GetElement(elem_idx);
                for (unsigned i = 0; i < p_elem->GetNumNodes(); ++i)
                {
                    unsigned nb_node = p_elem->GetNodeGlobalIndex(i);
                    if (nb_node != node_index)
                    {
                        neighbour_nodes.insert(nb_node);
                    }
                }
            }

            // O(k) lookup — directly look up each neighbour in the map
            for (unsigned nb_node : neighbour_nodes)
            {
                auto it = node_to_cell.find(nb_node);
                if (it == node_to_cell.end()) continue;  // boundary/ghost node

                c_vector<double, 2> nb_location =
                    rCellPopulation.GetLocationOfCellCentre(it->second);
                double distance = norm_2(nb_location - cell_location);

                if (distance > 0.0 && distance <= mTGFAlphaDiffusionRadius)
                {
                    double transfer = cell_tgfalpha
                                    * (1.0 - distance / mTGFAlphaDiffusionRadius)
                                    * 0.15;
                    tgfalpha_transfer[nb_node] += transfer;
                }
            }
        }

        // Apply accumulated transfers (one O(N) pass at most)
        for (auto& kv : tgfalpha_transfer)
        {
            auto it = node_to_cell.find(kv.first);
            if (it == node_to_cell.end()) continue;

            CellPtr p_cell = it->second;
            double new_amount = p_cell->GetCellData()->GetItem("tgfalpha") + kv.second;
            p_cell->GetCellData()->SetItem(
                "tgfalpha",
                std::max(0.0, std::min(mMaxTGFAlphaConcentration, new_amount)));
        }
    }
     
    virtual void SetupSolve(AbstractCellPopulation<2,2>& rCellPopulation, 
                           std::string outputDirectory)
    {
        // Initialize HIF-1a, TGF-a, and EGFR fractions in all cells
        for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            if (!cell_iter->GetCellData()->HasItem("hif1alpha"))
            {
                cell_iter->GetCellData()->SetItem("hif1alpha", 0.0);
            }
            if (!cell_iter->GetCellData()->HasItem("tgfalpha"))
            {
                cell_iter->GetCellData()->SetItem("tgfalpha", 0.0);
            }
            // Initialize EGFR fractions to constant 90:10 ratio
            //if (!cell_iter->GetCellData()->HasItem("egfr_mut_fraction"))
            //{
                //cell_iter->GetCellData()->SetItem("egfr_mut_fraction", MUTANT_EGFR_FRACTION);
            //}
            //if (!cell_iter->GetCellData()->HasItem("egfr_wt_fraction"))
            //{
                //cell_iter->GetCellData()->SetItem("egfr_wt_fraction", WT_EGFR_FRACTION);
            //}
        }
    }
    
    virtual void OutputSimulationModifierParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<HypoxiaThreshold>" << mHypoxiaThreshold 
                     << "</HypoxiaThreshold>\n";
        *rParamsFile << "\t\t\t<MaxHIF1Alpha>" << mMaxHIF1AlphaProduction 
                     << "</MaxHIF1Alpha>\n";
        *rParamsFile << "\t\t\t<TGFAlphaProductionRate>" << mTGFAlphaProductionRate 
                     << "</TGFAlphaProductionRate>\n";
        *rParamsFile << "\t\t\t<TGFAlphaDiffusionRadius>" << mTGFAlphaDiffusionRadius 
                     << "</TGFAlphaDiffusionRadius>\n";
        *rParamsFile << "\t\t\t<MutantEGFRFraction>" << MUTANT_EGFR_FRACTION 
                     << "</MutantEGFRFraction>\n";
        *rParamsFile << "\t\t\t<WTEGFRFraction>" << WT_EGFR_FRACTION 
                     << "</WTEGFRFraction>\n";
        *rParamsFile << "\t\t\t<MaxTGFAlphaConcentration>" << mMaxTGFAlphaConcentration
                     << "</MaxTGFAlphaConcentration>\n";
        
        AbstractCellBasedSimulationModifier<2>::OutputSimulationModifierParameters(rParamsFile);
    }
};