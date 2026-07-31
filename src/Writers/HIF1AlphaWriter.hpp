/*
* Created on: 21 Oct 2025
* Last modified: 21 Oct 2025
* 		Author: Vaishnudebi Dutta
*/

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class HIF1AlphaWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
{
public:
    HIF1AlphaWriter()
        : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("hif1alpha.dat")
    {
        this->mVtkCellDataName = "HIF1Alpha";
    }

    double GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
    {
        if (pCell->GetCellData()->HasItem("hif1alpha"))
        {
            return pCell->GetCellData()->GetItem("hif1alpha");
        }
        return 0.0;
    }

    void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
    {
        double val = GetCellDataForVtkOutput(pCell, pCellPopulation);
        *this->mpOutStream << val << " ";
    }
};