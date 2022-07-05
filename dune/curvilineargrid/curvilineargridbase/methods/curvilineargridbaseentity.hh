#ifndef CURVILINEAR_GRID_BASE_ENTITY
#define CURVILINEAR_GRID_BASE_ENTITY

#include <dune/common/exceptions.hh>

#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>


namespace Dune {

namespace CurvGrid {


template<typename GridBaseType>
class CurvilinearGridBaseEntity {

public:

	typedef typename GridBaseType::GridStorageType  GridStorageType;

    typedef typename GridBaseType::ctype  ctype;
	static const int   dimension  = GridBaseType::dimension;
	static const bool is_cached = GridBaseType::is_cached;

    typedef typename GridStorageType::GlobalIndexType           GlobalIndexType;
    typedef typename GridStorageType::LocalIndexType            LocalIndexType;
	typedef typename GridStorageType::InternalIndexType	        InternalIndexType;
	typedef typename GridStorageType::StructuralType            StructuralType;
	typedef typename GridStorageType::PhysicalTagType           PhysicalTagType;
	typedef typename GridStorageType::InterpolatoryOrderType    InterpolatoryOrderType;
	typedef typename GridStorageType::IdType                 	IdType;

    typedef typename GridStorageType::GlobalCoordinate       GlobalCoordinate;
    typedef typename GridStorageType::VertexStorage          VertexStorage;
    typedef typename GridStorageType::EdgeStorage            EdgeStorage;
    typedef typename GridStorageType::FaceStorage            FaceStorage;
    typedef typename GridStorageType::EntityStorage          EntityStorage;

    // Codimensions of entity types for better code readability
    static const int   VERTEX_CODIM   = GridStorageType::VERTEX_CODIM;
    static const int   EDGE_CODIM     = GridStorageType::EDGE_CODIM;
    static const int   FACE_CODIM     = GridStorageType::FACE_CODIM;
    static const int   ELEMENT_CODIM  = GridStorageType::ELEMENT_CODIM;

    typedef Dune::ReferenceElement< ctype, dimension > ReferenceElement3d;
//    typedef Dune::ReferenceElements<ctype, dimension>  ReferenceElements3d;
    typedef Dune::Geo::template ReferenceElements<ctype, dimension>  ReferenceElements3d;



public:

    // Defines the Elementary geometry
    template <int codim>
    struct Codim
    {
    	typedef typename std::conditional<is_cached,
    	  CachedCurvilinearGeometry<ctype, dimension - codim, dimension>,
    	  CurvilinearGeometry<ctype, dimension - codim, dimension>
    	>::type  EntityGeometry;
    };




public:

	CurvilinearGridBaseEntity(GridBaseType & gridbase) :
		gridstorage_(gridbase.gridstorage())
	{

	}


	/** \brief Get the level of the entity in the refinement structure
	 *   NOTE: Grid Refinement not implemented!
	 * */
    int level(int codim, LocalIndexType localIndex)  const { return 0; }


    /** \brief Get the GeometryType of entities on this process  */
    Dune::GeometryType geometryType(int codim, LocalIndexType localIndex) const
    {
    	switch(codim)
    	{
    	case VERTEX_CODIM  : return Dune::GeometryTypes::vertex;                         break;
    	case EDGE_CODIM    : return Dune::GeometryTypes::line;                           break;
    	case FACE_CODIM    : return gridstorage_.face_[localIndex].geometryType;         break;
    	case ELEMENT_CODIM : return gridstorage_.element_[localIndex].geometryType;      break;
    	}
    }


    /** \brief Get physical tag based on codimension  */
    PhysicalTagType physicalTag(int codim, LocalIndexType localIndex) const
    {
    	switch(codim)
    	{
    	case VERTEX_CODIM  : return 0;                                              break;
    	case EDGE_CODIM    : return 0;                                              break;
    	case FACE_CODIM    : return gridstorage_.face_[localIndex].physicalTag;     break;
    	case ELEMENT_CODIM : return gridstorage_.element_[localIndex].physicalTag;  break;
    	}
    }

    /** \brief Get Partition Type of an entity */
    PartitionType partitionType(int codim, LocalIndexType localIndex) const
    {
    	switch(codim)
    	{
    	case VERTEX_CODIM  : assert(localIndex < gridstorage_.point_.size());    return gridstorage_.point_[localIndex].ptype;    break;
    	case EDGE_CODIM    : assert(localIndex < gridstorage_.edge_.size());     return gridstorage_.edge_[localIndex].ptype;     break;
    	case FACE_CODIM    : assert(localIndex < gridstorage_.face_.size());     return gridstorage_.face_[localIndex].ptype;     break;
    	case ELEMENT_CODIM : assert(localIndex < gridstorage_.element_.size());  return gridstorage_.element_[localIndex].ptype;  break;
    	default :
    	{
    		LoggingMessage::template write<LOG_MSG_PERSISTENT>( __FILE__, __LINE__, "CurvilinearGridBase: unexpected codim " + std::to_string(codim));
    		DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected codimension");
    		break;
    	}
    	}
    }


    /** \brief Get Interpolatory order of the entity */
    InterpolatoryOrderType interpolationOrder(int codim, LocalIndexType localIndex) const
    {
    	if (codim >= VERTEX_CODIM)  { return 0; }

    	LocalIndexType localAssocElementIndex;

    	switch (codim)
    	{
    	case ELEMENT_CODIM : localAssocElementIndex = localIndex;  break;
    	case FACE_CODIM    : localAssocElementIndex = gridstorage_.face_[localIndex].element1Index;  break;
    	case EDGE_CODIM    : localAssocElementIndex = gridstorage_.edge_[localIndex].elementIndex;  break;
    	}

    	return gridstorage_.element_[localAssocElementIndex].interpOrder;
    }


    /** Storage data related to this entity, except of explicit vertex coordinates
     *  \param[in] codim                 codimension of the entity
     *  \param[in] localIndex            local edge index
     *
     *  Note: This data is stored only partially - vertex indices are extracted from associated element
     * */
    EntityStorage data(int codim, LocalIndexType localIndex) const
    {
    	switch (codim)
    	{
    	case VERTEX_CODIM  : return pointData(localIndex);     break;
    	case EDGE_CODIM    : return edgeData(localIndex);      break;
    	case FACE_CODIM    : return faceData(localIndex);      break;
    	case ELEMENT_CODIM : return elementData(localIndex);   break;
    	default:
    	{
    		DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Requested unexpected codim for entity data ");
    	} break;

    	}
    }

    /** \brief Retrieves the geometry class of the entity. This procedure is expensive, especially for cached geometries */
    template<int codim>
    typename Codim<codim>::EntityGeometry
    geometry(LocalIndexType localIndex) const
    {
    	//std::cout << "attempting to create a codim " << codim << " entity with localIndex=" << localIndex << std::endl;
    	return geometryConstructor<codim>(data(codim, localIndex));
    }



    /** \brief Finds global index using local index and codimension of entity.
     *   Returns false if requested local index does not correspond to an entity on this process  */
    bool findGlobalIndex(int codim, LocalIndexType localIndex, GlobalIndexType & globalIndex) const
    {
    	if (localIndex < 0)  { DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Received negative index");  }

    	switch(codim)
    	{
    	case VERTEX_CODIM  : if (localIndex >= gridstorage_.point_.size())    { return false; }  else { globalIndex = gridstorage_.point_[localIndex].globalIndex; }    break;
    	case EDGE_CODIM    : if (localIndex >= gridstorage_.edge_.size())     { return false; }  else { globalIndex = gridstorage_.edge_[localIndex].globalIndex; }     break;
    	case FACE_CODIM    : if (localIndex >= gridstorage_.face_.size())     { return false; }  else { globalIndex = gridstorage_.face_[localIndex].globalIndex; }     break;
    	case ELEMENT_CODIM : if (localIndex >= gridstorage_.element_.size())  { return false; }  else { globalIndex = gridstorage_.element_[localIndex].globalIndex; }  break;
    	}
    	return true;
    }


    /** \brief Finds local index using global index and codimension of entity.
     *   Returns false if requested global index does not correspond to an entity on this process  */
    bool findLocalIndex(int codim, GlobalIndexType globalIndex, LocalIndexType & localIndex) const
    {
    	if (globalIndex < 0)  { DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Received negative index");  }

    	const auto tmpIter = gridstorage_.entityIndexMap_[codim].find(globalIndex);
    	if (tmpIter != gridstorage_.entityIndexMap_[codim].end())  { localIndex = (*tmpIter).second;  return true; }
    	else  { return false; }
    }


    /** Generates a vector of local coordinates of corners of requested entity
     * The order of corners corresponds to the one provided by the ref.subEntity()
     *   */
    std::vector<LocalIndexType> cornerLocalIndex(int codim, LocalIndexType entityLocalIndex) const
	{
    	std::vector<LocalIndexType> rez;

    	if (codim == VERTEX_CODIM)  { rez.push_back(entityLocalIndex);  return rez; }

    	LocalIndexType    localElementIndex;
    	InternalIndexType subentityIndex;

    	switch(codim)
    	{
    	case ELEMENT_CODIM : localElementIndex = entityLocalIndex;  break;
    	case FACE_CODIM : {
    		localElementIndex = gridstorage_.face_[entityLocalIndex].element1Index;
    		subentityIndex  = gridstorage_.face_[entityLocalIndex].element1SubentityIndex;

    	} break;
    	case EDGE_CODIM : {
    		localElementIndex = gridstorage_.edge_[entityLocalIndex].elementIndex;
    		subentityIndex  = gridstorage_.edge_[entityLocalIndex].subentityIndex;
    	} break;
    	}


    	const EntityStorage & thisElem = gridstorage_.element_[localElementIndex];
        std::vector<LocalIndexType> elementCornerLocalIndexSet =
        		CurvilinearGeometryHelper::entityVertexCornerSubset<ctype, 3>(thisElem.geometryType, thisElem.vertexIndexSet, thisElem.interpOrder);

        // If we are interested in corners of the element, we are interested in all corners
        if (codim == ELEMENT_CODIM)  { return elementCornerLocalIndexSet; }

        // Otherwise, we need to calculate the subset of corners wrt selected subentity

        Dune::GeometryType elemGT=Dune::GeometryTypes::simplex(dimension);
        const auto & refElem = ReferenceElements3d::general(elemGT);
        int thisEntityCornerNumber = refElem.size(0, codim, dimension);

        for (int i = 0; i < thisEntityCornerNumber; i++)  {
        	//std::cout << "-- Attempting ref.subentity() using gt.dim=" << elemGT.dim() << " cdim=" << cdim << " codim="<<codim << " subcodim="<< VERTEX_CODIM << " subIndex=" << subentityIndex << " subsubindex=" << i << std::endl;

        	InternalIndexType thisCornerSubIndex = refElem.subEntity(subentityIndex, codim, i, VERTEX_CODIM);
        	rez.push_back(elementCornerLocalIndexSet[thisCornerSubIndex]);
        }

        return rez;
	}


    /** Returns the local index of a subentity of a given entity
     *  \param[in] entityIndex              local index of the entity
     *  \param[in] codim                    codimension of the entity
     *  \param[in] subcodim                 codimension of the subentity
     *  \param[in] subentityInternalIndex   subentity internal index wrt entity
     *
     *  \note subcodim > codim required
     *
     *  Algorithm:
     *  1) Find the parent tetrahedron index from (entityIndex, codim)
     *  2) Find the requested subentity internal index wrt parent tetrahedron
     *  3) Access tetrahedral subentity index storage to return stored answer
     *
     *  Note: To comply with DUNE regulations, if index of a vertex are requested,
     *  that index will be over the corners of the element as if it was a linear element.
     *  Curvilinear interpolation vertices will not be accessible by this method
     *
     * */
    LocalIndexType subentityLocalIndex (LocalIndexType entityIndex, int codim, int subcodim, InternalIndexType subentityInternalIndex) const
    {
    	// Stage 1) Find this subentity as an element subentity
    	// **************************************************************************
    	Dune::GeometryType tetrahedronGeometry=Dune::GeometryTypes::simplex(3);
    	const auto & thisRefElement = ReferenceElements3d::general(tetrahedronGeometry);

    	if (subcodim == codim)  { return entityIndex; }  // In this case return itself as own subentity
    	if (subcodim < codim) {                          // Wrong by definition
    		LoggingMessage::template write<LOG_MSG_PERSISTENT>( __FILE__, __LINE__, "CurvilinearGridBase: subentityIndex(): Unexpected codim-subcodim pair = (" + std::to_string(codim) + "," + std::to_string(subcodim) + ")");
    		DUNE_THROW(Dune::IOError, "CurvilinearGridBase: subentityIndex(): Unexpected codim-subcodim pair");
    	}

    	LocalIndexType elementLocalIndex;
    	InternalIndexType elementSubentityInternalIndex1;

    	switch (codim)
    	{
    		case ELEMENT_CODIM :
    		{
    			elementLocalIndex = entityIndex;
    			elementSubentityInternalIndex1 = subentityInternalIndex;
    		} break;
    		case FACE_CODIM :
    		{
    			InternalIndexType elementSubentityInternalIndex2 = gridstorage_.face_[entityIndex].element1SubentityIndex;

    			elementLocalIndex = gridstorage_.face_[entityIndex].element1Index;
    			elementSubentityInternalIndex1 = thisRefElement.subEntity(elementSubentityInternalIndex2, codim, subentityInternalIndex, subcodim);
    		} break;
    		case EDGE_CODIM :
    		{
    			InternalIndexType elementSubentityInternalIndex2 = gridstorage_.edge_[entityIndex].subentityIndex;

    			elementLocalIndex = gridstorage_.edge_[entityIndex].elementIndex;
    			elementSubentityInternalIndex1 = thisRefElement.subEntity(elementSubentityInternalIndex2, codim, subentityInternalIndex, subcodim);
    		} break;
    	}


    	// Stage 2) Find index of the element subentity
    	// **************************************************************************

    	int rez;

    	switch (subcodim)
    	{
    		// Face
    		case FACE_CODIM   :  rez = gridstorage_.elementSubentityCodim1_[elementLocalIndex][elementSubentityInternalIndex1];  break;
    		// Edge
    		case EDGE_CODIM   :  rez = gridstorage_.elementSubentityCodim2_[elementLocalIndex][elementSubentityInternalIndex1];  break;
    		// Corner
    		case VERTEX_CODIM :
    		{
    			InterpolatoryOrderType interpolationOrder = gridstorage_.element_[elementLocalIndex].interpOrder;
    			InternalIndexType cornerInternalIndex = CurvilinearGeometryHelper::cornerIndex(tetrahedronGeometry, interpolationOrder, elementSubentityInternalIndex1);
    			rez = gridstorage_.element_[elementLocalIndex].vertexIndexSet[cornerInternalIndex];
    		} break;
    	}

    	//std::cout << "Requested subentity index of entity=" << entityIndex << " codim=" << codim << " subcodim=" << subcodim << " internalIndex=" << subentityInternalIndex << " rez=" << rez << std::endl;
    	return rez;
    }


    /** \brief Returns the globalId of the entity  */
    IdType globalId(int codim, LocalIndexType localIndex) const
    {
    	GlobalIndexType globalIndex;
    	if (! findGlobalIndex(codim, localIndex, globalIndex))  { DUNE_THROW(Dune::IOError, "CurvilinearGridBase: requested local index does not point to an entity");  }
    	return IdType(globalIndex);
    }


    /** \brief Returns the local index of the associated globalId of the entity  */
    LocalIndexType localIndexFromGlobalId(int codim, IdType globalId) const
    {
    	GlobalIndexType globalIndex = globalId.globalindex_;
    	LocalIndexType localIndex;
    	if (! findLocalIndex(codim, globalIndex, localIndex))  { DUNE_THROW(Dune::IOError, "CurvilinearGridBase: requested globalID does not have an associated local index on this process");  }
    	return localIndex;
    }


protected:

    /* ***************************************************************************
     * Section: Auxiliary Methods
     * ***************************************************************************/

    EntityStorage pointData(LocalIndexType localIndex) const
    {
    	assert(localIndex < gridstorage_.point_.size());
    	const VertexStorage & thisPointData =  gridstorage_.point_[localIndex];

        EntityStorage thisPoint;
        thisPoint.geometryType=Dune::GeometryTypes::vertex;
        thisPoint.globalIndex  = thisPointData.globalIndex;
        thisPoint.ptype        = thisPointData.ptype;
        thisPoint.interpOrder  = 0;                    // Note: Points do not have an interpolation order
        thisPoint.physicalTag  = -1;                   // Note: Points do not have a physical tag
        thisPoint.vertexIndexSet.push_back(localIndex);   // Note: Point has only one vertex, and its index is the point index

        return thisPoint;
    }

    EntityStorage edgeData(LocalIndexType localIndex) const
    {
    	assert(localIndex < gridstorage_.edge_.size());
    	const EdgeStorage & thisEdgeData =    gridstorage_.edge_[localIndex];
        const EntityStorage & assocElement =  gridstorage_.element_[thisEdgeData.elementIndex];

        EntityStorage thisEdge;
        thisEdge.geometryType=Dune::GeometryTypes::line;
        thisEdge.globalIndex  = thisEdgeData.globalIndex;
        thisEdge.ptype        = thisEdgeData.ptype;
        thisEdge.interpOrder  = assocElement.interpOrder;
        thisEdge.physicalTag  = -1;        // Note: Edges do not have a physical tag

        // Get the internal element vertex indices associated with this face as a subentity
        std::vector<InternalIndexType> subentityVertexIndices =
            CurvilinearGeometryHelper::subentityInternalCoordinateSet<ctype, dimension>(assocElement.geometryType, thisEdge.interpOrder, 2, thisEdgeData.subentityIndex);

        // Calculate the localIndices's of vertices of this face by extracting them from the element vertex Ids
        for(unsigned int i = 0; i < subentityVertexIndices.size(); i++) { thisEdge.vertexIndexSet.push_back(assocElement.vertexIndexSet[subentityVertexIndices[i]]); }

        return thisEdge;
    }

    EntityStorage faceData(LocalIndexType localIndex) const
    {
    	assert(localIndex < gridstorage_.face_.size());
    	const FaceStorage & thisFaceData = gridstorage_.face_[localIndex];
    	assert(thisFaceData.element1Index < gridstorage_.element_.size());
        const EntityStorage & assocElement = gridstorage_.element_[thisFaceData.element1Index];

        EntityStorage thisFace;
        thisFace.geometryType=Dune::GeometryTypes::triangle;
        thisFace.globalIndex  = thisFaceData.globalIndex;
        thisFace.ptype        = thisFaceData.ptype;
        //if (thisFaceData.boundaryType == GridStorageType::FaceBoundaryType::DomainBoundary)  { thisFace.ptype = BOUNDARY_SEGMENT_PARTITION_TYPE; }

        thisFace.interpOrder  = assocElement.interpOrder;
        thisFace.physicalTag  = thisFaceData.physicalTag;

        // Get the internal element vertex indices associated with this face as a subentity
        std::vector<InternalIndexType> subentityVertexIndices =
            CurvilinearGeometryHelper::subentityInternalCoordinateSet<ctype, dimension>(assocElement.geometryType, thisFace.interpOrder, 1, thisFaceData.element1SubentityIndex);

        // Calculate the localIndices's of vertices of this face by extracting them from the element vertex Ids
        for(unsigned int i = 0; i < subentityVertexIndices.size(); i++) { thisFace.vertexIndexSet.push_back(assocElement.vertexIndexSet[subentityVertexIndices[i]]); }

        return thisFace;
    }

    EntityStorage elementData(LocalIndexType localIndex) const {
    	assert(localIndex < gridstorage_.element_.size());
    	return gridstorage_.element_[localIndex];
    }


    // Get curved geometry of an entity
    // TODO: assert mydim == element geometry type dim
    template<int codim>
    typename Codim<codim>::EntityGeometry
    geometryConstructor(EntityStorage thisData) const
    {
    	assert(thisData.geometryType.dim() == dimension - codim);

        std::vector<GlobalCoordinate> entityVertices;
        for (unsigned int i = 0; i < thisData.vertexIndexSet.size(); i++) {
        	LocalIndexType thisIndex = thisData.vertexIndexSet[i];
        	entityVertices.push_back(gridstorage_.point_[thisIndex].coord);
        }

        return typename Codim<codim>::EntityGeometry (thisData.geometryType, entityVertices, thisData.interpOrder);
    }





private:
	    GridStorageType & gridstorage_;


};


}

}


#endif //CURVILINEAR_GRID_BASE_ENTITY
