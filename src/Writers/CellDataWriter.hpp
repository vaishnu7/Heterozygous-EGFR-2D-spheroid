/*
 * CellDataWriter.hpp
 *
 * Custom cell data writer for PC9 spheroid simulation
 * Outputs oxygen gradient, hypoxia signaling, EGFR activation, and cell state data
 * Outputs force components and cell properties
 *
 * Author: Vaishnudebi Dutta
 * Date: 24 Oct, 2025
 */

#ifndef CELL_DATA_WRITER_HPP_
#define CELL_DATA_WRITER_HPP_

#include "AbstractCellWriter.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "OutputFileHandler.hpp"

/**
 * A cell writer class that outputs custom cell data for spheroid simulations.
 * Writes oxygen concentration, HIF-1a, TGF-a, EGFR activation, cell state,
 * and force information to structured output files.
 *
 * Template parameters:
 * ELEMENT_DIM - Mesh dimension (typically 2 for cross-section)
 * SPACE_DIM - Physical space dimension (typically 2)
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellDataWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
{
private:

    /** Output directory path for cell data files */
    std::string mOutputDirectory;

    /** Storage for computed friction force values */
    double mStoredFrictionForce;

    friend class boost::serialization::access;
    
    template<class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
        archive & mOutputDirectory;
        archive & mStoredFrictionForce;
    }

public:

    /**
     * Default constructor.
     * Sets up output directory and VTK naming conventions
     */
    CellDataWriter();
    
    /**
     * Visit each cell and write its data to file
     * Outputs spatial, temporal, and force information
     * 
     * @param pCell pointer to the cell
     * @param pCellPopulation pointer to the cell population
     */
    virtual void VisitCell(CellPtr pCell, 
                          AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation) override;
    
    /**
     * Set the output directory path for cell data files
     * @param rOutputDirectory Path to output directory
     */
    void SetOutputDirectory(const std::string& rOutputDirectory)
    {
        mOutputDirectory = rOutputDirectory;
    }
    
    /**
     * Get the current output directory path
     * @return Output directory path
     */
    std::string GetOutputDirectory() const
    {
        return mOutputDirectory;
    }

    /**
     * Set stored friction force value
     * @param storedFrictionForce Friction force value
     */
    void SetStoredFrictionForce(double storedFrictionForce)
    {
        mStoredFrictionForce = storedFrictionForce;
    }

    /**
     * Get stored friction force value
     * @return Friction force value
     */
    double GetStoredFrictionForce() const
    {
        return mStoredFrictionForce;
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellDataWriter)

#endif /* CELL_DATA_WRITER_HPP_ */