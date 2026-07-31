/*
* Created on: 21 Oct 2025
* Last modified: 21 Oct 2025
* 		Author: Vaishnudebi Dutta
*/

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class TGFAlphaWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
{
public:
    TGFAlphaWriter()
        : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("tgfalpha.dat")
    {
        this->mVtkCellDataName = "TGFAlpha";
    }

    double GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
    {
        if (pCell->GetCellData()->HasItem("tgfalpha"))
        {
            return pCell->GetCellData()->GetItem("tgfalpha");
        }
        return 0.0;
    }

    void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
    {
        double val = GetCellDataForVtkOutput(pCell, pCellPopulation);
        *this->mpOutStream << val << " ";
    }
};