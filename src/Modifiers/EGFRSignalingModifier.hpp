/*
* Created on: 21 Oct 2025
* Last modified: 01 Jun 2026
* 		Author: Vaishnudebi Dutta
* Modified: Blended WT and mutant EGFR response
* Modified: effective_oxygen rewritten to couple downstream_signal
*           to the cell cycle via SimpleOxygenBasedCellCycleModel.
* Modified: EGFRMutationType enum moved to test file. Modifier now
*           accepts mMaxOxygenRescueFraction directly as a constructor
*           argument. The test file owns the mutation type decision
*           and passes the appropriate fraction here.
*/
#include "CellData.hpp"
#include "ApcTwoHitCellMutationState.hpp"

// ============================================================================
// EGFR-TGF-alpha Signaling Modifier (Blended WT + Mutant)
// ============================================================================
// PC9 cells: 90% mutant EGFR (constitutive) + 10% WT EGFR (ligand-dependent)
//
// WT EGFR: requires TGF-alpha to activate
// Mutant EGFR: constitutively active (independent of TGF-alpha)
//
// mMaxOxygenRescueFraction is passed in from the test file, where the
// EGFRMutationType enum lives. Values:
//   EXON19_DEL → 0.5    (αC-helix stabilised, more stress-resilient)
//   L858R      → 0.35   (activation loop mutation, less stress-resilient)
// ============================================================================
class EGFRSignalingModifier : public AbstractCellBasedSimulationModifier<2>
{
private:
    double mTGFAlphaActivationThreshold;  // Min TGF-alpha to activate WT EGFR
    double mWildtypeEGFRSensitivity;      // WT EGFR sensitivity to TGF-alpha
    double mMutantEGFRBasalActivity;      // Mutant EGFR constitutive activity level
    double mDownstreamK;                  // Hill K for downstream signal
    int    mDownstreamN;                  // Hill n for downstream signal
    double mReceptorDecayRate;  // 0.25/h → t½ ≈ 2.8h
    double mDt;                 // simulation timestep, passed from constructor
    // =========================================================================
    // mMaxOxygenRescueFraction
    //
    // Passed in from the test file via EGFRMutationType enum lookup.
    //
    // Biological meaning:
    //   The maximum fraction of the mild-hypoxia O2 gap (0.05 - oxygen) that
    //   EGFR downstream signalling (PI3K/AKT, RAS/MAPK → cyclin D1) can
    //   compensate for, at full signal strength. Both WT and mutant EGFR
    //   contribute via downstream_signal (which already encodes the blended
    //   WT + mutant activation weighted by their fractions). The mutation type
    //   sets the ceiling because exon 19 del maintains a more structurally
    //   stable active conformation under hypoxic stress than L858R.
    // =========================================================================
    double mMaxOxygenRescueFraction;

public:
    EGFRSignalingModifier(double tgfAlphaThreshold        = 0.1,
                          double wildtypeEGFRSensitivity   = 1.0,
                          double mutantEGFRBasalActivity   = 0.8,
                          double downstreamK               = 0.3,
                          int    downstreamN               = 2,
                          double receptorDecayRate           = 0.25,
                          double maxOxygenRescueFraction   = 0.5)
        : AbstractCellBasedSimulationModifier<2>(),
          mTGFAlphaActivationThreshold(tgfAlphaThreshold),
          mWildtypeEGFRSensitivity(wildtypeEGFRSensitivity),
          mMutantEGFRBasalActivity(mutantEGFRBasalActivity),
          mDownstreamK(downstreamK),
          mDownstreamN(downstreamN),
          mReceptorDecayRate(receptorDecayRate),
          mDt(0.01),  // default value; should be overridden by test via constructor arg
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
                // Necrotic cells have no EGFR activity
                //cell_iter->GetCellData()->SetItem("egfr_activation",     0.0);
                //cell_iter->GetCellData()->SetItem("wt_egfr_activation",  0.0);
                //cell_iter->GetCellData()->SetItem("mut_egfr_activation",  0.0);
                //cell_iter->GetCellData()->SetItem("effective_oxygen",    0.0);
                //cell_iter->GetCellData()->SetItem("downstream_signal",   0.0);
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
                cell_iter->GetCellData()->SetItem("effective_oxygen", oxygen); // no rescue

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
            // Both EXON19_DEL and L858R are constitutively active — they do
            // not require TGF-alpha. mMutantEGFRBasalActivity is the same
            // for both (1.0 when passed from test). The mutation type only
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
            //
            // SimpleOxygenBasedCellCycleModel reads "effective_oxygen" in
            // preference to "oxygen" (in both UpdateCellCyclePhase and
            // UpdateHypoxicDuration) when the key is present. By raising
            // effective_oxygen within the mild hypoxia band, we reduce the
            // rate at which G1 is lengthened each timestep and delay
            // hypoxic death onset — both direct cell cycle effects.
            //
            // Rescue is active only in mild hypoxia: 0.01 <= oxygen < 0.05
            //
            // Formula:
            //   effective_oxygen = oxygen
            //                    + (0.05 - oxygen)          [gap to hypoxic upper bound]
            //                    * mMaxOxygenRescueFraction  [mutation type ceiling]
            //                    * downstream_signal         [total EGFR signal strength]
            //
            // Note: downstream_signal is derived from blended_egfr_activation,
            // which already encodes both WT (TGF-alpha dependent) and mutant
            // (constitutive) contributions weighted by their fractions.
            // egfr_mut_fraction is NOT multiplied separately — that would
            // double-penalise the WT contribution and incorrectly suppress
            // rescue in cells where WT EGFR is active via TGF-alpha.
            //
            // Properties:
            //   - In normoxia or above hypoxic threshold (oxygen >= 0.05):
            //     no rescue. Cell cycle governed by real oxygen alone.
            //   - In mild hypoxia (0.01 <= oxygen < 0.05): rescue proportional
            //     to total EGFR signal strength, ceilinged by mutation type.
            //   - In severe hypoxia (oxygen < 0.01): rescue disabled.
            //     EGFR signalling cannot overcome critical O2 deficit.
            //   - effective_oxygen is naturally bounded below 0.05 because
            //     gap > 0 and mMaxOxygenRescueFraction < 1.
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
