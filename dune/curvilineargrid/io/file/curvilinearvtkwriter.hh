/*******************************************************************
 * Curvilinear VTK writer
 * 
 * author: Aleksejs Fomins
 * date: 01.09.2014 - created
 * 
 * description:
 * Generates VTK, VTU and PVTU files to visualise curvilinear elements.
 * Curvilinear elements added to the writer are automatically discretized using
 * fixed interval linear elements, which are consequently written to a file
 * 
 *******************************************************************/


#ifndef DUNE_CURVILINEARVTKWRITER_HH
#define DUNE_CURVILINEARVTKWRITER_HH

#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <stdio.h>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>

#include <dune/curvilineargeometry/interpolation/lagrangeinterpolator.hh>
#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>

#include <dune/curvilineargrid/common/constant.hh>
#include <dune/curvilineargrid/common/loggingmessage.hh>
#include <dune/curvilineargrid/io/file/curvilinearvtkformat.hh>
#include <dune/curvilineargrid/io/file/curvilinearvtkentitysubset.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilineargridstorage.hh>



namespace Dune
{

namespace CurvGrid {


// [TODO] VTKWriter consumes a lot of RAM, since it stores the entire subsampled mesh before writing it to file.
//				The problem is that for each appended element, several places in the output file change (vertices, elements, fields)
//				Solution: Write a temporary file for each dataset during insert. During write combine all temporary files into output file and delete temporary files
// [TODO] Output files are gigantic.
//				Problem: Using "ascii" format for very large data. Vertices shared by nearby elements written twice
//				Solution: Implement format="binary" output format, as defined in VTK Legacy File Formats Reference.
//				Solution: Somehow find which vertices are shared.
//					Note: This is hard because we do not want to store all data in the RAM
//					Note: Benefit only for vertices - can't in general use for shared edges and faces, since they sometimes need to be written twice with different tags
// [TODO] Implement abstract SCALAR CELL DATA formalism.
//				Rewrite PhysicalTag, StructuralType, Rank, ProcessBind through this formalism.
//				Do not store them explicitly, but add them as fields through the VTKGridWriter
// [TODO] For complex fields we use phase animation - rotate all fields by the same time-dependent phase factor.
//				Currently our FEM code splits the phase animation into 360 steps, and writes 360 real-valued VTK files
//				Geometry is rewritten 360 times, although exactly the same, and 360 fields for each data point are exhaustively defined by 2 numbers
//				Solution: Write only one file. Store two real fields for each complex field (Re and Im).
//				Question: How to animate a simple linear function of 2 fields in paraview: (Cory Quammen) suggests that this can be done using "Programmable Filter" -> "Calculator"

  template<class GridType>
  class CurvilinearVTKWriter
  {

/** \brief This class takes curved elements, samples them on a grid, connects points into mesh of
 * small straight-sided triangles and edges, ands writes them to the .vtk file
 *
 */
  public:

	  typedef CurvilinearVTKWriter<GridType> This;

	  typedef typename GridType::ctype  ctype;
      static const int dimension = GridType::dimension;
      static const bool isCached = GridType::is_cached;

      // Codimensions of entity types for better code readability
      static const int   VERTEX_CODIM   = 3;
      static const int   EDGE_CODIM     = 2;
      static const int   FACE_CODIM     = 1;
      static const int   ELEMENT_CODIM  = 0;

      typedef FieldVector< ctype, dimension >             GlobalCoordinate;
      typedef std::vector<GlobalCoordinate>               GlobalCoordinateVec;
      typedef std::vector<int>                            IndexVector;
      typedef std::vector<int>                            TagVector;
      typedef std::vector<IndexVector >                   SubEntityIndexVector;
      typedef std::map<std::vector<int>, int>             LocalCoordinate2GlobalIdMap;
      typedef std::vector< LocalCoordinate2GlobalIdMap >  Coord2GlobalMapVector;
      typedef typename LocalCoordinate2GlobalIdMap::iterator  LC2GIMapIter;

      typedef std::map<std::string, unsigned int>         FieldNameMap;
      typedef std::pair<std::string, unsigned int>        FieldNamePair;
      typedef std::map<int, ctype>                        FieldScalarMap;
      typedef std::pair<int, ctype>                       FieldScalarPair;
      typedef std::map<int, GlobalCoordinate>             FieldCoordMap;
      typedef std::pair<int, GlobalCoordinate>            FieldCoordPair;

      typedef typename FieldNameMap::iterator             FieldNameMapIter;
      typedef typename FieldNameMap::const_iterator       FieldNameMapConstIter;
      typedef typename FieldScalarMap::iterator           FieldScalarMapIter;
      typedef typename FieldScalarMap::const_iterator     FieldScalarMapConstIter;
      typedef typename FieldCoordMap::iterator            FieldCoordMapIter;
      typedef typename FieldCoordMap::const_iterator      FieldCoordMapConstIter;

      typedef std::vector< IndexVector >                  ElemGridEnumerate;

	  typedef std::vector< std::size_t >                  SizeTVec;
      typedef std::vector< std::string >                  StrVec;

  public:

	CurvilinearVTKWriter (int rank, int size, bool withPeriodicBind = false) :
		nInterval_(0),         // Defined later
		rank_(rank),
		size_(size)
	{

	}


    CurvilinearVTKWriter (MPIHelper &mpihelper)
    	: This(mpihelper.rank(), mpihelper.size())
    {}




    /** \brief Takes curvilinear element, discretizes it into linear element, adds linear elements to the writer
     *
     *  \param[in]  geomtype                    Geometry Type of the element
     *  \param[in]  nodeSet                 Coordinates of the element's interpolation points
     *  \param[in]  tagSet                  Set of 0 to 3 tags, in priority order. [Physical Tag, Structural Type, ProcessRank]
     *  \param[in]  elementOrder                   Curvilinear interpolation order of the element
     *  \param[in]  nDiscretizationPoint          Number of discretization points-per-edge
     *  \param[in]  interpolate                    If set to false, instead of interpolating, the provided interpolation points will be used as discretization points. Otherwise, a sub-sampling grid over the element will be interpolated
     *  \param[in]  explode	                    there will be gaps between all elements, obtained by scaling them away from the center of mass
     *  \param[in]  writeEdgeData                  If true, discretization points are connected using linear edges forming triangular mesh, and added to writer
     *  \param[in]  writeTriangleData              If true, discretization points are connected using linear triangles forming triangular mesh, and added to writer
     *
     *    \note Discretization is performed using a regular grid over the reference element, for example,
     *    below the stars represent the triangular discretization points with nDiscretizationPoint = 5
     *
     *  *
     *  **
     *  ***
     *  ****
     *  *****
     */
    template<int mydim, int subdim>
    void addCurvilinearElement(
            const Dune::GeometryType & geomtype,
            const GlobalCoordinateVec & nodeSet,
            const TagVector & tagSet,
            int elementOrder,
            int nDiscretizationPoint,
            bool interpolate,
            bool explode)
    {
    	assert((mydim > 0) && (mydim <= 3));  // Forbid writing vertices and hypergeometric entities
    	assert((subdim > 0) && (subdim <= 3));  // Forbid subrefinement vertices and hypergeometric entities
    	assert(subdim <= mydim);  // Subentity refinement dimension can not be larger than the entity dimension

        int   thisElmPhysTag        = (tagSet.size() > 0 ) ? tagSet[0] : 0;
        int   thisElmPartitionType  = (tagSet.size() > 1 ) ? tagSet[1] : 0;
        int   thisElmProcessRank    = (tagSet.size() > 2 ) ? tagSet[2] : 0;

        // Treat boundary segments in a special way
        std::string pname;
        bool magnify;

        // [TODO] Boundary types do not belong here
        // - pass magnification as a parameter to addCurvilinearElement
        // - implement partitionNameExtended in constant.hh
        if (thisElmPartitionType == PERIODIC_GHOST_PARTITION_TYPE) {
        	magnify = false;
        	pname = "PeriodicGhost";
        } else if (thisElmPartitionType == BOUNDARY_SEGMENT_PARTITION_TYPE) {
        	magnify = true;
        	pname = "DomainBoundarySegment";
        } else if (thisElmPartitionType == PERIODIC_BOUNDARY_PARTITION_TYPE) {
            	magnify = true;
            	pname = "PeriodicBoundarySegment";
        } else if (thisElmPartitionType == INTERIOR_BOUNDARY_SEGMENT_PARTITION_TYPE) {
            	magnify = false;
            	pname = "InteriorBoundarySegment";
        } else if (thisElmPartitionType == PERIODIC_GHOST_BIND_EDGE_TYPE) {
            	magnify = false;
            	pname = "PeriodicBindEdge";
        } else {
        	magnify = false;
        	pname = Dune::PartitionName(static_cast<Dune::PartitionType> (thisElmPartitionType));
        }

        std::stringstream log_message;
        log_message << "VTK_WRITER: Adding a curvilinear element Type=" << CurvilinearGeometryHelper::geometryName(geomtype);
        log_message << " Order="               << elementOrder;
        log_message << " PhysicalTag="         << thisElmPhysTag;
        log_message << " StructuralType="      << pname;
        log_message << " ProcessRank="         << thisElmProcessRank;
        log_message << " nDiscretization="     << nDiscretizationPoint;
        log_message << " useInterpolation="    << interpolate;
        log_message << " explodeElements="     << explode;
        log_message << " subrefinement_dim="      << subdim;
        //log_message << " vertices=" << VectorHelper::vector2string(nodeSet);
        LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, log_message.str());

        addCurvilinearSimplex<mydim, subdim>(geomtype, nodeSet, tagSet, elementOrder, nDiscretizationPoint, explode, magnify, interpolate);
    }


    // Initialize field name. Useful, if master process does not have any entities of a given type, but still needs to write a correct .pvtu file
    template <class VTKVectorFunction>
    void initScalarField(const VTKVectorFunction & vtkfunction) {
    	LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: Initializing scalar field " + vtkfunction.name());

    	typedef typename VTKVectorFunction::Domain  Domain;
    	typedef typename VTKVectorFunction::Range   Range;
    	const int mydimension = VTKVectorFunction::mydimension;

        // The index of the field associated with this field name
    	int fieldIndex = updateFieldIndex(scalarFieldName2Index_, vtkFieldScalar_, vtkfunction.name());
    }

    // Initialize field name. Useful, if master process does not have any entities of a given type, but still needs to write a correct .pvtu file
    template <class VTKVectorFunction>
    void initVectorField(const VTKVectorFunction & vtkfunction) {
    	LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: Initializing vector field " + vtkfunction.name());

    	typedef typename VTKVectorFunction::Domain  Domain;
    	typedef typename VTKVectorFunction::Range   Range;
    	const int mydimension = VTKVectorFunction::mydimension;

        // The index of the field associated with this field name
    	int fieldIndex = updateFieldIndex(vectorFieldName2Index_, vtkFieldVector_, vtkfunction.name());
    }



    template <class VTKVectorFunction>
    void addScalarField(const VTKVectorFunction & vtkfunction) {
    	LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: Adding element field " + vtkfunction.name());

    	typedef typename VTKVectorFunction::Domain  Domain;
    	typedef typename VTKVectorFunction::Range   Range;
    	const int mydimension = VTKVectorFunction::mydimension;

        // The index of the field associated with this field name
    	int fieldIndex = updateFieldIndex(scalarFieldName2Index_, vtkFieldScalar_, vtkfunction.name());

    	// Loop over all points
    	for (LC2GIMapIter iter = parameter2Index_.begin(); iter != parameter2Index_.end(); iter++)
    	{
    		std::vector<int> integerCoordinate = (*iter).first;
    		unsigned int     vertexIndex       = (*iter).second;

    		// Compute local coordinate
    		Domain local;
    		for (int i = 0; i < mydimension; i++)  { local[i] = (ctype((*iter).first[i])) / nInterval_; }

    		// Evaluate field
    		Range field = vtkfunction.evaluate(local);

    		//std::cout << "VTKWriter sampling local coordinate " << local << " coordinate index " << vertexIndex << " field=" << field << std::endl;

    		// Append coordinate index and field
    		vtkFieldScalar_[fieldIndex].insert(FieldScalarPair(vertexIndex, field));
    	}
    }

    template <class VTKVectorFunction>
    void addVectorField(const VTKVectorFunction & vtkfunction) {
    	LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: Adding element field " + vtkfunction.name());

    	typedef typename VTKVectorFunction::Domain  Domain;
    	typedef typename VTKVectorFunction::Range   Range;
    	const int mydimension = VTKVectorFunction::mydimension;

        // The index of the field associated with this field name
    	int fieldIndex = updateFieldIndex(vectorFieldName2Index_, vtkFieldVector_, vtkfunction.name());

    	// Loop over all points
    	for (LC2GIMapIter iter = parameter2Index_.begin(); iter != parameter2Index_.end(); iter++)
    	{
    		std::vector<int> integerCoordinate = (*iter).first;
    		unsigned int     vertexIndex       = (*iter).second;

    		// Compute local coordinate
    		Domain local;
    		for (int i = 0; i < mydimension; i++)  { local[i] = (ctype((*iter).first[i])) / nInterval_; }

    		// Evaluate field
    		Range field = vtkfunction.evaluate(local);

    		//std::cout << "VTKWriter sampling local coordinate " << local << " coordinate index " << vertexIndex << " field=" << field << std::endl;

    		// Append coordinate index and field
    		vtkFieldVector_[fieldIndex].insert(FieldCoordPair(vertexIndex, field));
    	}
    }


    // Writes serial VTK file
    // Takes a vector of vertices, a vector of edges and a vector of triangles, writes them all to a .VTK file
    void writeVTK(std::string filename, std::string format) const
    {
        std::ofstream vtkFile(filename.c_str());

        std::vector<std::size_t> nEntity
        {
        	vtkCodimVertexIndex_[ELEMENT_CODIM].size(),
        	vtkCodimVertexIndex_[FACE_CODIM].size(),
        	vtkCodimVertexIndex_[EDGE_CODIM].size(),
        	vtkPoint_.size()
        };

        unsigned int nEntityTot =   nEntity[EDGE_CODIM] +   nEntity[FACE_CODIM] +   nEntity[ELEMENT_CODIM];
        unsigned int nCellTot   = 3*nEntity[EDGE_CODIM] + 4*nEntity[FACE_CODIM] + 5*nEntity[ELEMENT_CODIM];

        // Check if the datasize is consistent for writing
        writerSelfCheck(nEntity);

        vtkFile << "# vtk DataFile Version 2.0" << std::endl;
        vtkFile << "CurvilinearGmshReader test output" << std::endl;
        vtkFile << format << std::endl;  // e.g. ASCII
        vtkFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

        // Write all points
        vtkFile << "POINTS " << vtkPoint_.size() << " double" << std::endl;
        VTKFormat::writeCoordinateArray(vtkFile, vtkPoint_, "\n", format);

        // Write all elements
        vtkFile << std::endl;
        vtkFile << "CELLS " << nEntityTot << " " << nCellTot << std::endl;
        for (int iCodim = EDGE_CODIM; iCodim >= ELEMENT_CODIM; iCodim--)  // Write edges first, elements last
        {
        	IndexVector codimVertexIndex;
            for (auto const & vertexIndex : vtkCodimVertexIndex_[iCodim])     {
            	codimVertexIndex.push_back(vertexIndex.size());
				for (auto const & thisIndex : vertexIndex) { codimVertexIndex.push_back(thisIndex); }
            }
            VTKFormat::writeArray(vtkFile, codimVertexIndex, " ", format);
        }


        // Write edge and triangle cell types
        vtkFile << std::endl;
        vtkFile << "CELL_TYPES " << nEntityTot << std::endl;
    	VTKFormat::writeArray(vtkFile, nEntity[EDGE_CODIM], "3", "\n", format);
    	VTKFormat::writeArray(vtkFile, nEntity[FACE_CODIM], "5", "\n", format);
    	VTKFormat::writeArray(vtkFile, nEntity[ELEMENT_CODIM], "10", "\n", format);

        // If are defined fields associated with vertices, then write them too
        if ((vectorFieldName2Index_.size() > 0) || (scalarFieldName2Index_.size() > 0))
        {
        	vtkFile << std::endl;
        	vtkFile << "POINT_DATA " << vtkPoint_.size() << std::endl;

        	// Write all scalar fields
        	for (const auto & iterName : scalarFieldName2Index_) {
        		std::string fieldName = iterName.first;
        		int fieldIndex  = iterName.second;

        		vtkFile << "SCALARS " << fieldName << " FLOAT" << std::endl;
        		{
        			std::vector<ctype> scalarField;
                    for (unsigned int i = 0; i < vtkPoint_.size(); i++ ) {
                    	// If there is no field defined for this vertex, just print a zero vector for consistency
                    	FieldScalarMapConstIter iterField = vtkFieldScalar_[fieldIndex].find(i);
                    	if (iterField == vtkFieldScalar_[fieldIndex].end())  {
                    		scalarField.push_back(0.0);
                    	} else {
                    		scalarField.push_back((*iterField).second);
                    	}
                    }
                    VTKFormat::writeArray(vtkFile, scalarField, "\n", format);
        		}
        	}

        	// Write all vector fields
        	for (const auto & iterName : vectorFieldName2Index_) {
        		std::string fieldName = iterName.first;
        		int fieldIndex  = iterName.second;

        		vtkFile << "VECTORS " << fieldName << " FLOAT" << std::endl;
        		{
        			std::vector<GlobalCoordinate> vectorField;
                    for (unsigned int i = 0; i < vtkPoint_.size(); i++ ) {
                    	// If there is no field defined for this vertex, just print a zero vector for consistency
                    	FieldCoordMapConstIter iterField = vtkFieldVector_[fieldIndex].find(i);
                    	if (iterField == vtkFieldVector_[fieldIndex].end())  {
                    		vectorField.push_back(GlobalCoordinate(0.0));
                    	} else {
                    		vectorField.push_back((*iterField).second);
                    	}
                    }
                    VTKFormat::writeCoordinateArray(vtkFile, vectorField, "\n", format);
        		}
        	}
        }

        vtkFile << std::endl;
        vtkFile << "CELL_DATA " << nEntityTot << std::endl;

        // Write edge and triangle Structural type
        vtkFile << "SCALARS physicalTag FLOAT" << std::endl;
        vtkFile << "LOOKUP_TABLE default" << std::endl;
        VTKFormat::writeArray(vtkFile, vtkCodimPhysicalTag_[EDGE_CODIM], "\n", format);
        VTKFormat::writeArray(vtkFile, vtkCodimPhysicalTag_[FACE_CODIM], "\n", format);
        VTKFormat::writeArray(vtkFile, vtkCodimPhysicalTag_[ELEMENT_CODIM], "\n", format);
        vtkFile << std::endl;

        // Write edge and triangle physicalTags
        vtkFile << "SCALARS structuralType FLOAT" << std::endl;
        vtkFile << "LOOKUP_TABLE default" << std::endl;
        VTKFormat::writeArray(vtkFile, vtkCodimStructuralType_[EDGE_CODIM], "\n", format);
        VTKFormat::writeArray(vtkFile, vtkCodimStructuralType_[FACE_CODIM], "\n", format);
        VTKFormat::writeArray(vtkFile, vtkCodimStructuralType_[ELEMENT_CODIM], "\n", format);
        vtkFile << std::endl;

        // Write edge and triangle provider process ranks
        vtkFile << "SCALARS processRank FLOAT" << std::endl;
        vtkFile << "LOOKUP_TABLE default" << std::endl;
        VTKFormat::writeArray(vtkFile, vtkCodimProcessRank_[EDGE_CODIM], "\n", format);
        VTKFormat::writeArray(vtkFile, vtkCodimProcessRank_[FACE_CODIM], "\n", format);
        VTKFormat::writeArray(vtkFile, vtkCodimProcessRank_[ELEMENT_CODIM], "\n", format);
        vtkFile << std::endl;

        // Empty line at the end of file
        vtkFile << std::endl;
        vtkFile.close();
    }


    // Writes a PVTU parallel file (no data in this file)
    // Requite path and filename as separate strings, because
    void writePVTU(std::string path, std::string filenameBody, int size, std::string format) const
    {
    	std::string filename = path + filenameBody + ".pvtu";
    	std::ofstream pvtuFile(filename.c_str());
    	bool selfClose = true;  // In PVTU file most XML closes are self-closing


        // Write header
        // *****************************************************
        pvtuFile << "<?xml version=\"" << VTK_XML_VERSION << "\"?>" << std::endl;
        {
        	VTKFormat::XMLClause clauseVTKFile(pvtuFile, "VTKFile", StrVec{"type", "version", "byte_order"}, StrVec{"P"+VTK_GRID_TYPE, VTK_VTU_VERSION, VTK_BYTE_ORDER});
        	{
        		VTKFormat::XMLClause clauseVTKGrid(pvtuFile, "P" + VTK_GRID_TYPE, StrVec{"GhostLevel"}, StrVec{"0"});

                // Write PointData - scalar and vector fields sampled over all vertices
                // *****************************************************
                bool haveScalar = scalarFieldName2Index_.size() > 0;
                bool haveVector = vectorFieldName2Index_.size() > 0;
                if (haveScalar || haveVector)
                {
                	// Generate and write PointData preamble
                    std::stringstream pointDataNames;
                    if (haveScalar)  { pointDataNames << "Scalars=\"" << namesPile(scalarFieldName2Index_) << "\"";  if (haveVector) { pointDataNames << " "; } }
                    if (haveVector)  { pointDataNames << "Vectors=\"" << namesPile(vectorFieldName2Index_) << "\"";                                             }

                    {
                        VTKFormat::XMLClause clausePPointData(pvtuFile, "PPointData", pointDataNames.str());

                    	// Write point data arrays
                    	for (const auto & iterName : scalarFieldName2Index_) {
                			VTKFormat::XMLClauseSingle(pvtuFile, "DataArray",
                					StrVec{"type", "Name", "NumberOfComponents", "format"},
        							StrVec{"Float32", iterName.first, "1", format}, true
                			);
                    	}

                    	for (const auto & iterName : vectorFieldName2Index_) {
                			VTKFormat::XMLClauseSingle(pvtuFile, "DataArray",
                					StrVec{"type", "Name", "NumberOfComponents", "format"},
        							StrVec{"Float32", iterName.first, "3", format}, true
                			);
                    	}
                    }
                }


                // Write scalars associated with entities
                //*****************************************************
                {
                    VTKFormat::XMLClause clausePCellData(pvtuFile, "PCellData", StrVec{"Scalars"}, StrVec{"physicalTag"});

                    VTKFormat::XMLClauseSingle(pvtuFile, "PDataArray",
                    		StrVec{"type", "Name", "NumberOfComponents"},
    						StrVec{"Float32","physicalTag", "1"}, true);

                    VTKFormat::XMLClauseSingle(pvtuFile, "PDataArray",
                    		StrVec{"type", "Name", "NumberOfComponents"},
    						StrVec{"Float32","structuralType", "1"}, true);

                    VTKFormat::XMLClauseSingle(pvtuFile, "PDataArray",
                    		StrVec{"type", "Name", "NumberOfComponents"},
    						StrVec{"Float32","processRank", "1"}, true);
                }


                // Write coordinates of vertices
                // *****************************************************
                {
                	VTKFormat::XMLClause clausePPoints(pvtuFile, "PPoints");
                    VTKFormat::XMLClauseSingle(pvtuFile, "PDataArray", StrVec{"type", "NumberOfComponents"}, StrVec{"Float32", "3"}, true);
                }

                // Write element information
                // *****************************************************
                {
                	VTKFormat::XMLClause clausePCells(pvtuFile, "PCells");
                    VTKFormat::XMLClauseSingle(pvtuFile, "PDataArray", StrVec{"type", "Name"}, StrVec{"Int32", "connectivity"}, true);
                    VTKFormat::XMLClauseSingle(pvtuFile, "PDataArray", StrVec{"type", "Name"}, StrVec{"Int32", "offsets"}, true);
                    VTKFormat::XMLClauseSingle(pvtuFile, "PDataArray", StrVec{"type", "Name"}, StrVec{"UInt8", "types"}, true);
                }

                // Write all .vtu data files file
                // *****************************************************
                for (int iProc = 0; iProc < size; iProc++ )
                {
                    std::string vtuFilename = filenameBody + "_process_" + std::to_string(iProc) + ".vtu";
                    VTKFormat::XMLClauseSingle(pvtuFile, "Piece", StrVec{"Source"}, StrVec{vtuFilename}, true);
                }
        	}
        }
        pvtuFile.close();
    }


    // Writes serial VTU file
    void writeVTU(std::string filename, std::string format) const
    {
    	assert((format == "ascii") || (format == "binary"));

    	std::ofstream vtkFile(filename.c_str());

        std::vector<std::size_t> nEntity
        {
        	vtkCodimVertexIndex_[ELEMENT_CODIM].size(),
        	vtkCodimVertexIndex_[FACE_CODIM].size(),
        	vtkCodimVertexIndex_[EDGE_CODIM].size(),
        	vtkPoint_.size()
        };

        std::size_t nCells = nEntity[EDGE_CODIM] + nEntity[FACE_CODIM] + nEntity[ELEMENT_CODIM];

        // Check if the datasize is consistent for writing
        writerSelfCheck(nEntity);

        // Compute offsets which is a general way to determine the number of vertices per element
        std::vector<unsigned int> offsets;
        int tmp_offset = 0;
        for (unsigned int i = 0; i < nEntity[EDGE_CODIM]; i++)     { tmp_offset += 2;  offsets.push_back(tmp_offset); }
        for (unsigned int i = 0; i < nEntity[FACE_CODIM]; i++)     { tmp_offset += 3;  offsets.push_back(tmp_offset); }
        for (unsigned int i = 0; i < nEntity[ELEMENT_CODIM]; i++)  { tmp_offset += 4;  offsets.push_back(tmp_offset); }

        // Write header
        // *****************************************************

        vtkFile << "<?xml version=\"" << VTK_XML_VERSION << "\"?>" << std::endl;
        {
        	VTKFormat::XMLClause clauseVTKFile(vtkFile, "VTKFile", StrVec{"type", "version", "byte_order"}, StrVec{VTK_GRID_TYPE, VTK_VTU_VERSION, VTK_BYTE_ORDER});
        	{
        		VTKFormat::XMLClause clauseVTKGrid(vtkFile, VTK_GRID_TYPE);
        		{
        			VTKFormat::XMLClause clauseVTKPiece(vtkFile, "Piece", StrVec{"NumberOfPoints", "NumberOfCells"}, SizeTVec{nEntity[VERTEX_CODIM], nCells});
                    // Write PointData - scalar and vector fields sampled over all vertices
                    // *****************************************************
                    bool haveScalar = scalarFieldName2Index_.size() > 0;
                    bool haveVector = vectorFieldName2Index_.size() > 0;
                    if (haveScalar || haveVector)
                    {
                    	// Generate and write PointData preamble
                        std::stringstream pointDataNames;
                        if (haveScalar)  { pointDataNames << "Scalars=\"" << namesPile(scalarFieldName2Index_) << "\"";  if (haveVector) { pointDataNames << " "; } }
                        if (haveVector)  { pointDataNames << "Vectors=\"" << namesPile(vectorFieldName2Index_) << "\"";                                             }

                        { // Write point data arrays
                        	VTKFormat::XMLClause clauseVTKPointData(vtkFile, "PointData", pointDataNames.str());

                        	for (const auto & iterName : scalarFieldName2Index_)
                        	{
                        		std::string fieldName = iterName.first;
                        		int fieldIndex  = iterName.second;

                        		{
                        			VTKFormat::XMLClause clauseVTKScalarField(vtkFile, "DataArray",
                        					StrVec{"type", "Name", "NumberOfComponents", "format"},
											StrVec{"Float32", fieldName, "1", format}
                        			);

                        			std::vector<ctype> scalarField;
                                    for (unsigned int i = 0; i < vtkPoint_.size(); i++ ) {
                                    	// If there is no field defined for this vertex, just print a zero vector for consistency
                                    	FieldScalarMapConstIter iterField = vtkFieldScalar_[fieldIndex].find(i);
                                    	if (iterField == vtkFieldScalar_[fieldIndex].end())  {
                                    		scalarField.push_back(0.0);
                                    	} else {
                                    		scalarField.push_back((*iterField).second);
                                    	}
                                    }
                                    VTKFormat::writeArray(vtkFile, scalarField, " ", format);
                        		}
                        	}

                        	for (const auto & iterName : vectorFieldName2Index_)
                        	{
                        		std::string fieldName = iterName.first;
                        		int fieldIndex  = iterName.second;

                        		{
                        			VTKFormat::XMLClause clauseVTKVectorField(vtkFile, "DataArray",
                        					StrVec{"type", "Name", "NumberOfComponents", "format"},
											StrVec{"Float32", fieldName, "3", format}
                        			);

                        			std::vector<GlobalCoordinate> vectorField;
                                    for (unsigned int i = 0; i < vtkPoint_.size(); i++ ) {
                                    	// If there is no field defined for this vertex, just print a zero vector for consistency
                                    	FieldCoordMapConstIter iterField = vtkFieldVector_[fieldIndex].find(i);
                                    	if (iterField == vtkFieldVector_[fieldIndex].end())  {
                                    		vectorField.push_back(GlobalCoordinate(0.0));
                                    	} else {
                                    		vectorField.push_back((*iterField).second);
                                    	}
                                    }
                                    VTKFormat::writeCoordinateArray(vtkFile, vectorField, "\n", format);
                        		}
                        	}
                        }
                    }


                    // Write edge and triangle physicalTags and structural types
                    // *****************************************************
                    {
                    	VTKFormat::XMLClause clauseCellData(vtkFile, "CellData", StrVec{"Scalars"}, StrVec{"physicalTag"});
                        { // Write edge and triangle physicalTags
                			VTKFormat::XMLClause clauseTag(vtkFile, "DataArray",
                					StrVec{"type", "Name", "NumberOfComponents", "format"},
									StrVec{"Float32", "physicalTag", "1", format}
                			);

                			VTKFormat::writeArray(vtkFile, vtkCodimPhysicalTag_[EDGE_CODIM], " ", format);
                			VTKFormat::writeArray(vtkFile, vtkCodimPhysicalTag_[FACE_CODIM], " ", format);
                			VTKFormat::writeArray(vtkFile, vtkCodimPhysicalTag_[ELEMENT_CODIM], " ", format);
                        }

                        { // Write edge and triangle structural type
                			VTKFormat::XMLClause clauseTag(vtkFile, "DataArray",
                					StrVec{"type", "Name", "NumberOfComponents", "format"},
									StrVec{"Float32", "structuralType", "1", format}
                			);

                			VTKFormat::writeArray(vtkFile, vtkCodimStructuralType_[EDGE_CODIM], " ", format);
                			VTKFormat::writeArray(vtkFile, vtkCodimStructuralType_[FACE_CODIM], " ", format);
                			VTKFormat::writeArray(vtkFile, vtkCodimStructuralType_[ELEMENT_CODIM], " ", format);
                        }

                        { // Write edge and triangle provider process ranks
                			VTKFormat::XMLClause clauseTag(vtkFile, "DataArray",
                					StrVec{"type", "Name", "NumberOfComponents", "format"},
									StrVec{"Float32", "processRank", "1", format}
                			);

                			VTKFormat::writeArray(vtkFile, vtkCodimProcessRank_[EDGE_CODIM], " ", format);
                			VTKFormat::writeArray(vtkFile, vtkCodimProcessRank_[FACE_CODIM], " ", format);
                			VTKFormat::writeArray(vtkFile, vtkCodimProcessRank_[ELEMENT_CODIM], " ", format);
                        }
                    }

                    // Write coordinates of vertices
                    // *****************************************************
                    {
                    	VTKFormat::XMLClause xmlPoints(vtkFile, "Points");
                        { // Write all points
                			VTKFormat::XMLClause clauseTag(vtkFile, "DataArray",
                					StrVec{"type", "NumberOfComponents", "format"},
									StrVec{"Float32", "3", format}
                			);

                			VTKFormat::writeCoordinateArray(vtkFile, vtkPoint_, "\n", format);
                        }
                    }

                    // Write element information
                    // *****************************************************
                    {
                    	VTKFormat::XMLClause xmlCells(vtkFile, "Cells");
                        {
                			VTKFormat::XMLClause clauseTag(vtkFile, "DataArray",
                					StrVec{"type", "Name", "format"},
									StrVec{"Int32", "connectivity", format}
                			);

                            for (int iCodim = EDGE_CODIM; iCodim >= ELEMENT_CODIM; iCodim--)  // Write edges first, elements last
                            {
                            	for (const auto & entityVertexIndex : vtkCodimVertexIndex_[iCodim]) {
                            		VTKFormat::writeArray(vtkFile, entityVertexIndex, " ", format);
                            	}
                            }
                        }

                        {
                			VTKFormat::XMLClause clauseTag(vtkFile, "DataArray",
                					StrVec{"type", "Name", "format"},
									StrVec{"Int32", "offsets", format}
                			);

                			VTKFormat::writeArray(vtkFile, offsets, " ", format);
                        }

                        {   // Element types
                			VTKFormat::XMLClause clauseTag(vtkFile, "DataArray",
                					StrVec{"type", "Name", "format"},
									StrVec{"UInt8", "types", format}
                			);

                        	VTKFormat::writeArray(vtkFile, nEntity[EDGE_CODIM], "3", " ", format);
                        	VTKFormat::writeArray(vtkFile, nEntity[FACE_CODIM], "5", " ", format);
                        	VTKFormat::writeArray(vtkFile, nEntity[ELEMENT_CODIM], "10", " ", format);
                        }
                    }
        		}
        	}
        }
        vtkFile.close();
    }


    // Writes a VTU file on all processes and a PVTU on Master Process
    void writeParallelVTU(std::string path, std::string filenameBody, std::string format) const
    {
        // Write a PVTU file on master process
        if (rank_ == 0) { writePVTU(path, filenameBody, size_, format); }

        // Write a VTU file on all processes
        std::string vtuFileName = path + filenameBody  + "_process_" + std::to_string(rank_) + ".vtu";
        writeVTU(vtuFileName, format);
    }


    // Checks if the storage arrays are consistent before writing to file
    void writerSelfCheck(std::vector<std::size_t> nEntity) const
    {
        for (int iCodim = ELEMENT_CODIM; iCodim <= EDGE_CODIM; iCodim++)
        {
            assert(vtkCodimPhysicalTag_[iCodim].size()    == nEntity[iCodim]);
            assert(vtkCodimStructuralType_[iCodim].size() == nEntity[iCodim]);
            assert(vtkCodimProcessRank_[iCodim].size()    == nEntity[iCodim]);
        }
    }

  protected:

    // ***********************************************
    // Auxiliary Methods
    // ***********************************************

    /** \brief Adds a set of linear discretization entities which will be explicitly written to the file
     *  Note: Only add tags if they are provided by the user
     * */
    template<int codim, int subcodim>
    void refineEntity(
    		const ElemGridEnumerate & simplexEnumerateReduced,
			int nInterval,
			const LocalCoordinate2GlobalIdMap & parametricToIndex,
			const std::vector<int> & tagSet)
    {
    	LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKWriter: Computing and writing refinement-edges" );
    	SubEntityIndexVector thisEntitySubset = VTKEntitySubset::refineEntitySubset<codim, subcodim>(simplexEnumerateReduced, nInterval, parametricToIndex);

    	for (const auto & subEntity : thisEntitySubset) {
            vtkCodimVertexIndex_[subcodim].push_back(subEntity);
            vtkCodimPhysicalTag_[subcodim].push_back(    (tagSet.size() > 0) ? tagSet[0] : 0);
            vtkCodimStructuralType_[subcodim].push_back( (tagSet.size() > 1) ? tagSet[1] : 0);
            vtkCodimProcessRank_[subcodim].push_back(    (tagSet.size() > 2) ? tagSet[2] : 0);
    	}
    }

    /** \brief Calculates the centre of mass of a vector of points (equal-weighted) */
    GlobalCoordinate vectorCentreOfMass( const GlobalCoordinateVec & cornerVector) const
    {
        GlobalCoordinate rez;
        for (unsigned int i = 0; i < cornerVector.size(); i++) { rez += cornerVector[i]; }
        rez /= cornerVector.size();
        return rez;
    }


    /** \brief Checks if the tetrahedral discretization point corresponds to tetrahedral boundary
     *
     *  \param[in]  p                               tetrahedral discretization point key
     *  \param[in]  nInterval                     Number of discretization intervals per edge (number of discretization points-per-edge minus 1)
     *
     */
    bool onTetrahedronBoundary(const std::vector<int> & p, const int nInterval) const
    {
          return ((p[0] == 0) || (p[1] == 0) || (p[2] == 0) || (p[0] + p[1] + p[2] == nInterval));
    }


    // Pass CoM, because it is more optimal to only calculate it once for all samples over the element
    void shrinkMagnify(GlobalCoordinate & point, const GlobalCoordinate & CoM, bool explode, bool magnify) const
    {
        ctype shrinkMagnitude       = explode ? 0.2 : 0.0;  // 0.0 - no shrinking, 0.99 - very small element (must be < 1)
        ctype boundaryMagnification = magnify ? 1.2 : 1.0;  // Expand all domain boundary surfaces a little bit so that they do not interlay with element surfaces

        for (int d = 0; d < dimension; d++)  {
        	point[d] = (point[d] + (CoM[d] - point[d]) * shrinkMagnitude) * boundaryMagnification;
        }
    }


    template <int mydim, int subdim>
    void addCurvilinearSimplex(
            const Dune::GeometryType & geomtype,
            const GlobalCoordinateVec & nodeSet,
            const std::vector<int> & tagSet,
            int elementOrder,
            int nDiscretizationPoint,
            bool shrink,
            bool magnify,
            bool interpolate)
    {
    	const int codim = dimension - mydim;
    	const int subcodim = dimension - subdim;
    	typedef FieldVector< ctype, mydim >      LocalVector;

        LagrangeInterpolator<ctype, mydim, dimension> elementInterpolator(geomtype, nodeSet, elementOrder);




        // *******************************************************************************
        // Step 1. Find coordinates of the corners and the (linear) center of mass
        // *******************************************************************************
        int nDofThis    = nodeSet.size();
        int nCornerThis = elementInterpolator.nCorner();


        LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKWriter: Calculating CoM" );
        GlobalCoordinateVec cornerVector;
        for (int i = 0; i < nCornerThis; i++) { cornerVector.push_back(elementInterpolator.corner(i)); }
        GlobalCoordinate CoM = vectorCentreOfMass(cornerVector);


        // *******************************************************************************
        // Step 2: Construct sampling grid over the element. Sample the points, store them
        //   and map to their insertion coordinate
        // For a Generic triangle the sampling looks like this
        //
        //  (i,j)=(n,0) -> *
        //                 **
        //                 ***
        //                 ****
        //  (i,j)=(0,0) -> ****** <- (i,j) = (0,n)
        //
        // Step 3: Store the index of vertices used for future use
        // *******************************************************************************
        nInterval_ = interpolate ? nDiscretizationPoint - 1 : elementOrder;
        parameter2Index_.clear();


        LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKWriter: Calculating Enumerators" );
        ElemGridEnumerate  simplexEnumerate        = CurvilinearGeometryHelper::simplexGridEnumerate<mydim>(nInterval_);
        ElemGridEnumerate  simplexEnumerateReduced = CurvilinearGeometryHelper::simplexGridEnumerate<mydim>(nInterval_ - 1);
        std::vector< LocalVector > simplexLocalGrid = CurvilinearGeometryHelper::simplexGridCoordinateSet<ctype, mydim>(simplexEnumerate, nInterval_);

        LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKWriter: Computing and inserting refinement vertices" );
        for (unsigned int i = 0; i < simplexEnumerate.size(); i++)
        {
            // Find if this vertex is internal or boundary
            // bool isBoundaryPoint = (mydim == 3) ? onTetrahedronBoundary(simplexEnumerate[i], nInterval_) : true;

            // If we interpolate, then all points will be taken from new sample grid
            // Otherwise we take the intrinsic interpolation point grid which has the same shape
            GlobalCoordinate tmpPoint = interpolate ? elementInterpolator.global(simplexLocalGrid[i]) : nodeSet[i];

            //Optionally, transform the element. Do not try to multipy element by 1.0 if it is not going to be transformed
            if (shrink || magnify) { shrinkMagnify(tmpPoint, CoM, shrink, magnify); }

            // Add coordinates to the coordinates array
            vtkPoint_.push_back(tmpPoint);

            // Add point to the point map
            parameter2Index_[simplexEnumerate[i]] = vtkPoint_.size() - 1;
        }

        // *******************************************************************************
        // Step 4: Write entities that discretize this element to VTK
        // Current paradigm is to either discretize all inserted entities with edges, making sort of a net
        // Or to discretize each inserted entity with smaller entities of the same codim
        // *******************************************************************************
        refineEntity<codim, subcodim>(simplexEnumerateReduced, nInterval_, parameter2Index_, tagSet);

        //if (writeCodim[EDGE_CODIM])  { refineEntity<codim, EDGE_CODIM>(simplexEnumerateReduced, nInterval_, parameter2Index_, tagSet);  }
        //if (writeCodim[FACE_CODIM])     { refineEntity<codim, FACE_CODIM>(simplexEnumerateReduced, nInterval, parametricToIndex, tagSet, periodicBindIndex);  }
        //if (writeCodim[ELEMENT_CODIM])  { refineEntity<codim, ELEMENT_CODIM>(simplexEnumerateReduced, nInterval, parametricToIndex, tagSet, periodicBindIndex);  }
    }


    template <class Map, class Storage>
    unsigned int updateFieldIndex(Map & namemap, Storage & storage, std::string name) {
    	FieldNameMapIter iter = namemap.find(name);

    	if (iter != namemap.end() )  {       // If field exists, find its index
    		return (*iter).second;
    	} else {                                      // Otherwise, make new index and add field name to them map
    		unsigned int rez = storage.size();
    		namemap.insert(FieldNamePair(name, rez));
    		storage.push_back(typename Storage::value_type());
    		return rez;
    	}
    }


    // Stacks up all names from the map in a single string, separated by commas, with no extra spaces
    std::string namesPile(const FieldNameMap & map) const {
    	std::stringstream rez;

    	int i = 0;
    	for (FieldNameMapConstIter iterName = map.begin(); iterName != map.end(); ++iterName) {
    		if (i > 0) { rez << ","; }
    		rez << (*iterName).first;
    		i++;
     	}

    	return rez.str();
    }




  private:
    // MPI
    int rank_;
    int size_;

    // Temporary storage of current entity before it is written to storage arrays
    int nInterval_;
    LocalCoordinate2GlobalIdMap parameter2Index_;   // Stores the set of vertex indices used when inserting the last entity

    GlobalCoordinateVec vtkPoint_;

    // Discretization entity storage
    std::vector<IndexVector> vtkCodimVertexIndex_[4];   // Vertices that discretize the entity of a given codimension
    TagVector vtkCodimPhysicalTag_[4];                  // Physical tag of the entity
    TagVector vtkCodimStructuralType_[4];               // Structural type of the entity
    TagVector vtkCodimProcessRank_[4];                  // Process rank of the entity

    // Field storage
    FieldNameMap                 scalarFieldName2Index_;
    FieldNameMap                 vectorFieldName2Index_;
    std::vector<FieldScalarMap>  vtkFieldScalar_;   // For each field maps vertex index to field
    std::vector<FieldCoordMap>   vtkFieldVector_;   // For each field maps vertex index to field



  };
} // namespace CurvGrid

}  // namespace Dune

#endif /** DUNE_CURVILINEARVTKWRITER_HH **/
