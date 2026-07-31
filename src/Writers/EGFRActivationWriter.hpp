/*
* Created on: 21 Oct 2025
* Last modified: 21 Oct 2025
* 		Author: Vaishnudebi Dutta
*/

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class EGFRActivationWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
{
public:
    EGFRActivationWriter()
        : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("egfr_activation.dat")
    {
        this->mVtkCellDataName = "EGFRActivation";
    }

    double GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
    {
        if (pCell->GetCellData()->HasItem("egfr_activation"))
        {
            return pCell->GetCellData()->GetItem("egfr_activation");
        }
        return 0.0;
    }

    void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
    {
        double val = GetCellDataForVtkOutput(pCell, pCellPopulation);
        *this->mpOutStream << val << " ";
    }
};