/*
* Created on: 21 Oct 2025
* Last modified: 01 Jun 2026
* 		Author: Vaishnudebi Dutta
*/
#include "CellData.hpp"
#include "ApcTwoHitCellMutationState.hpp"

class EGFRSignalingModifier : public AbstractCellBasedSimulationModifier<2>
{
private:
    double mTGFAlphaActivationThreshold;  // Min TGF-alpha to activate WT EGFR
    double mWildtypeEGFRSensitivity;      // WT EGFR sensitivity to TGF-alpha
    double mMutantEGFRBasalActivity;      // Mutant EGFR constitutive activity level
    double mDownstreamK;                  // Hill K for downstream signal
    int    mDownstreamN;                  // Hill n for downstream signal
    double mReceptorDecayRate;  // 0.25/h
    double mDt;                 // simulation timestep, passed from constructor
    double mMaxOxygenRescueFraction;

public:
    EGFRSignalingModifier(double tgfAlphaThreshold        = 0.1,
                          double wildtypeEGFRSensitivity   = 1.0,
                          double mutantEGFRBasalActivity   = 0.8,
                          double downstreamK               = 0.3,
                          int    downstreamN               = 2,
                          double receptorDecayRate           = 0.25,
                          double dt                        = 0.005,
                          double maxOxygenRescueFraction   = 0.5)
        : AbstractCellBasedSimulationModifier<2>(),
          mTGFAlphaActivationThreshold(tgfAlphaThreshold),
          mWildtypeEGFRSensitivity(wildtypeEGFRSensitivity),
          mMutantEGFRBasalActivity(mutantEGFRBasalActivity),
          mDownstreamK(downstreamK),
          mDownstreamN(downstreamN),
          mReceptorDecayRate(receptorDecayRate),
          mDt(dt),  // default value
          mMaxOxygenRescueFraction(maxOxygenRescueFraction)
    {
    }

    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<2,2>& rCellPopulation)
    {
        for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            bool cell_is_necrotic = cell_iter->GetMutationState()->IsType<ApcTwoHitCellMutationState>();

            if (cell_is_necrotic)
            {
 
                double necrosis_onset = cell_iter->GetCellData()->GetItem("necrosis_onset_time");
                double now = SimulationTime::Instance()->GetTime();
                if (necrosis_onset < 0.0)
                {
                    cell_iter->GetCellData()->SetItem("necrosis_onset_time", now);
                }

                double decay = std::exp(-mReceptorDecayRate * mDt);

                double new_mut = cell_iter->GetCellData()->GetItem("mut_egfr_activation") * decay;
                double new_wt  = cell_iter->GetCellData()->GetItem("wt_egfr_activation")  * decay;

                cell_iter->GetCellData()->SetItem("mut_egfr_activation", new_mut);
                cell_iter->GetCellData()->SetItem("wt_egfr_activation",  new_wt);

                double wt_frac  = cell_iter->GetCellData()->GetItem("egfr_wt_fraction");
                double mut_frac = cell_iter->GetCellData()->GetItem("egfr_mut_fraction");
                double new_blended = (wt_frac * new_wt) + (mut_frac * new_mut);
                cell_iter->GetCellData()->SetItem("egfr_activation", new_blended);

                double egfr_act_n  = std::pow(new_blended, mDownstreamN);
                double K_n         = std::pow(mDownstreamK,  mDownstreamN);
                double downstream  = egfr_act_n / (K_n + egfr_act_n);
                cell_iter->GetCellData()->SetItem("downstream_signal", downstream);

                double oxygen = cell_iter->GetCellData()->GetItem("oxygen");
                cell_iter->GetCellData()->SetItem("effective_oxygen", oxygen); 

                continue;
            }

            double tgfalpha        = cell_iter->GetCellData()->GetItem("tgfalpha");
            double egfr_wt_fraction  = cell_iter->GetCellData()->GetItem("egfr_wt_fraction");
            double egfr_mut_fraction = cell_iter->GetCellData()->GetItem("egfr_mut_fraction");

            // ===== WT EGFR Response (ligand-dependent) =====
            double wt_egfr_activation = 0.0;

            if (tgfalpha > mTGFAlphaActivationThreshold)
            {
                // Sigmoidal dose-response; requires TGF-alpha above threshold
                wt_egfr_activation = (tgfalpha * mWildtypeEGFRSensitivity) /
                                     (mTGFAlphaActivationThreshold + tgfalpha * mWildtypeEGFRSensitivity);
                wt_egfr_activation = std::min(1.0, wt_egfr_activation);
            }

            // ===== Mutant EGFR Response (constitutive) =====
            // Both EXON19_DEL and L858R are constitutively active. The mutation type only
            // affects mMaxOxygenRescueFraction, set in the constructor.
            double mut_egfr_activation = std::min(1.0, mMutantEGFRBasalActivity);

            // ===== Blended EGFR Activation =====
            double blended_egfr_activation = (egfr_wt_fraction  * wt_egfr_activation) +
                                             (egfr_mut_fraction * mut_egfr_activation);
            blended_egfr_activation = std::min(1.0, blended_egfr_activation);

            cell_iter->GetCellData()->SetItem("egfr_activation",    blended_egfr_activation);
            cell_iter->GetCellData()->SetItem("wt_egfr_activation",  wt_egfr_activation);
            cell_iter->GetCellData()->SetItem("mut_egfr_activation",  mut_egfr_activation);

            // ================================================================
            // DOWNSTREAM PATHWAY SIGNAL
            //
            // Abstract Hill-function representation of the EGFR-triggered
            // intracellular cascade (PI3K/AKT, RAS/MAPK). Switch-like
            // behaviour (n=2) reflects the threshold nature of cyclin D1
            // induction. Output is normalised to [0,1].
            // ================================================================
            double egfr_act_n       = std::pow(blended_egfr_activation, mDownstreamN);
            double K_n              = std::pow(mDownstreamK,             mDownstreamN);
            double downstream_signal = egfr_act_n / (K_n + egfr_act_n);
            downstream_signal        = std::min(1.0, downstream_signal);

            cell_iter->GetCellData()->SetItem("downstream_signal", downstream_signal);

            // ================================================================
            // EFFECTIVE OXYGEN — coupling downstream_signal to the cell cycle
            // ================================================================
            double oxygen = cell_iter->GetCellData()->GetItem("oxygen");
            double effective_oxygen = oxygen;  // default: no rescue

            if (oxygen >= 0.01 && oxygen < 0.05)
            {
                double gap    = 0.05 - oxygen;  // always > 0 in this band
                double rescue = gap
                                * mMaxOxygenRescueFraction
                                * downstream_signal;

                effective_oxygen = oxygen + rescue;
            }
            // Outside 0.01–0.05: effective_oxygen stays = oxygen

            cell_iter->GetCellData()->SetItem("effective_oxygen", effective_oxygen);
        }
    }

    virtual void SetupSolve(AbstractCellPopulation<2,2>& rCellPopulation,
                            std::string outputDirectory)
    {
        // Initialise EGFR fractions and activation in all cells
        for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            if (!cell_iter->GetCellData()->HasItem("egfr_wt_fraction"))
            {
                cell_iter->GetCellData()->SetItem("egfr_wt_fraction",  0.1);
            }
            if (!cell_iter->GetCellData()->HasItem("egfr_mut_fraction"))
            {
                cell_iter->GetCellData()->SetItem("egfr_mut_fraction", 0.9);
            }
            if (!cell_iter->GetCellData()->HasItem("egfr_activation"))
            {
                cell_iter->GetCellData()->SetItem("egfr_activation",    0.0);
            }
            if (!cell_iter->GetCellData()->HasItem("wt_egfr_activation"))
            {
                cell_iter->GetCellData()->SetItem("wt_egfr_activation", 0.0);
            }
            if (!cell_iter->GetCellData()->HasItem("mut_egfr_activation"))
            {
                cell_iter->GetCellData()->SetItem("mut_egfr_activation", mMutantEGFRBasalActivity);
            }
            if (!cell_iter->GetCellData()->HasItem("effective_oxygen"))
            {
                cell_iter->GetCellData()->SetItem("effective_oxygen", 1.0);
            }
            if (!cell_iter->GetCellData()->HasItem("downstream_signal"))
            {
                cell_iter->GetCellData()->SetItem("downstream_signal", 0.0);
            }
            if (!cell_iter->GetCellData()->HasItem("necrosis_onset_time"))
            {
                cell_iter->GetCellData()->SetItem("necrosis_onset_time", -1.0);
            }
        }
        UpdateAtEndOfTimeStep(rCellPopulation);
    }

    virtual void OutputSimulationModifierParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<TGFAlphaActivationThreshold>" << mTGFAlphaActivationThreshold
                     << "</TGFAlphaActivationThreshold>\n";
        *rParamsFile << "\t\t\t<WildtypeEGFRSensitivity>" << mWildtypeEGFRSensitivity
                     << "</WildtypeEGFRSensitivity>\n";
        *rParamsFile << "\t\t\t<MutantEGFRBasalActivity>" << mMutantEGFRBasalActivity
                     << "</MutantEGFRBasalActivity>\n";
        *rParamsFile << "\t\t\t<DownstreamK>" << mDownstreamK << "</DownstreamK>\n";
        *rParamsFile << "\t\t\t<DownstreamN>" << mDownstreamN << "</DownstreamN>\n";
        *rParamsFile << "\t\t\t<MaxOxygenRescueFraction>" << mMaxOxygenRescueFraction
                     << "</MaxOxygenRescueFraction>\n";

        AbstractCellBasedSimulationModifier<2>::OutputSimulationModifierParameters(rParamsFile);
    }
};
