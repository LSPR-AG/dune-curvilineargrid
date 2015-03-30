/***************************************************************************
                          tetmesh.h  -  description
                             -------------------
    begin                : Mon Dec 15 2003
    copyright            : (C) 2003 by Roman Geus
    email                : roman.geus@psi.ch
    edited by            : Hua Guo, Sep 2010
***************************************************************************/

/***************************************************************************
                          curvilineargridbase.hh
                             -------------------
    begin                : Tue Nov 25 2014
    copyright            : (C) 2014 by Aleksejs Fomins, LSPR AG
    description          : Upgraded the mesh to curvilinear grid
***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DUNE_CURVILINEARGRIDCONSTRUCTOR_HH
#define DUNE_CURVILINEARGRIDCONSTRUCTOR_HH

#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <iostream>
#include <assert.h>

#include <fstream>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parallel/mpicollectivecommunication.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/curvilineargeometry/interpolation/curvilinearelementinterpolator.hh>
#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>
#include <dune/curvilineargeometry/curvilineargeometry.hh>

#include <dune/curvilineargrid/common/loggingmessage.hh>
#include <dune/curvilineargrid/common/vectorhelper.hh>

#include <dune/curvilineargrid/curvilineargridbase/curvilineargridstorage.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilineargridbase.hh>

#include <dune/curvilineargrid/curvilineargridbase/curvilinearghostconstructor.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilinearpostconstructor.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilinearoctreenode.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilinearlooseoctree.hh>


namespace Dune {




// Forwards-declatation of the base class
// **********************************************
template <class ct, int cdim, bool isCached>
class CurvilinearGridBase;


// Constructor class
// **********************************************
template <class ct, int cdim, bool isCached>
class CurvilinearGridConstructor {
public:

    /* public types */
    typedef Dune::CurvilinearGridStorage<ct, cdim, isCached>    GridStorageType;
    typedef Dune::CurvilinearGridBase<ct, cdim, isCached>       GridBaseType;

    typedef Dune::CurvilinearGhostConstructor<ct,cdim, isCached> GridGhostConstructor;
    typedef Dune::CurvilinearPostConstructor<ct,cdim, isCached>  GridPostConstructor;


    typedef typename GridStorageType::GlobalIndexType           GlobalIndexType;
    typedef typename GridStorageType::LocalIndexType            LocalIndexType;
	typedef typename GridStorageType::InternalIndexType	        InternalIndexType;
	typedef typename GridStorageType::StructuralType            StructuralType;
	typedef typename GridStorageType::PhysicalTagType           PhysicalTagType;
	typedef typename GridStorageType::InterpolatoryOrderType    InterpolatoryOrderType;

    typedef typename GridStorageType::Vertex                    Vertex;
    typedef typename GridStorageType::VertexStorage             VertexStorage;
    typedef typename GridStorageType::EdgeStorage               EdgeStorage;
    typedef typename GridStorageType::FaceStorage               FaceStorage;
    typedef typename GridStorageType::EntityStorage             EntityStorage;

    typedef typename GridStorageType::EdgeKey                   EdgeKey;
    typedef typename GridStorageType::FaceKey                   FaceKey;

    typedef std::map<EdgeKey, LocalIndexType>                   EdgeKey2EdgeIndexMap;
    typedef std::map<FaceKey, LocalIndexType>                   FaceKey2FaceIndexMap;
    typedef typename EdgeKey2EdgeIndexMap::iterator             EdgeMapIterator;
    typedef typename FaceKey2FaceIndexMap::iterator             FaceMapIterator;

    typedef typename GridStorageType::Global2LocalMap           Global2LocalMap;
    typedef typename GridStorageType::Global2LocalIterator      Global2LocalIterator;
    typedef typename GridStorageType::Local2LocalMap            Local2LocalMap;
    typedef typename GridStorageType::Local2LocalIterator       Local2LocalIterator;

    typedef typename GridStorageType::LocalIndexSet             LocalIndexSet;
    typedef typename GridStorageType::IndexSetIterator          IndexSetIterator;


	typedef typename GridStorageType::NodeType                  NodeType;
    typedef typename GridStorageType::CurvilinearLooseOctree    CurvilinearLooseOctree;

    // Codimensions of entity types for better code readability
    static const int   VERTEX_CODIM   = GridStorageType::VERTEX_CODIM;
    static const int   EDGE_CODIM     = GridStorageType::EDGE_CODIM;
    static const int   FACE_CODIM     = GridStorageType::FACE_CODIM;
    static const int   ELEMENT_CODIM  = GridStorageType::ELEMENT_CODIM;

    // Face Structural Type
    static const unsigned int DomainBoundaryType   = GridStorageType::PartitionType::DomainBoundary;
    static const unsigned int ProcessBoundaryType  = GridStorageType::PartitionType::ProcessBoundary;
    static const unsigned int InternalType         = GridStorageType::PartitionType::Internal;
    static const unsigned int GhostType            = GridStorageType::PartitionType::Ghost;

    // Logging Message Typedefs
    static const unsigned int LOG_PHASE_DEV       = Dune::LoggingMessage::Phase::DEVELOPMENT_PHASE;
    static const unsigned int LOG_CATEGORY_DEBUG  = Dune::LoggingMessage::Category::DEBUG;
    static const unsigned int LOG_CATEGORY_ERROR  = Dune::LoggingMessage::Category::ERROR;



public: /* public methods */

    /** Parallel constructor - USE THIS CONSTRUCTOR*/
    CurvilinearGridConstructor(
    		GridStorageType & gridstorage,
    		GridBaseType & gridbase,
    		MPIHelper &mpihelper,
    		LoggingMessage &loggingmessage) :
        gridstorage_(gridstorage),
        gridbase_(gridbase),
        mpihelper_(mpihelper),
        loggingmessage_(loggingmessage)
    {
        rank_ = mpihelper_.rank();
        size_ = mpihelper_.size();

        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "Initialized CurvilinearGridConstructor");
    }


public:

    /* ***************************************************************************
     * Section: Loading the mesh
     * ***************************************************************************/


    /** \brief Add a new vertex to the mesh
     * \param[in] p                coordinate of this vertex
     * \param[in] globalIndex      global index of this vertex
     * \param[in] structtype       (optional) partition type of this vertex. User SHOULD use default value
     * */
    void insertVertex(Vertex p, GlobalIndexType globalIndex, StructuralType structtype = InternalType)
    {
        VertexStorage point;
        point.coord = p;
        point.globalIndex = globalIndex;
        point.structuralType = structtype;

        gridstorage_.entityIndexMap_[VERTEX_CODIM][globalIndex] = gridstorage_.point_.size();
        gridstorage_.point_.push_back(point);

        std::stringstream log_stream;
        log_stream << "CurvilinearGridConstructor: Inserted vertex LocalIndex=" << gridstorage_.point_.size()-1 << " GlobalIndex=" << globalIndex;
        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, log_stream.str());
    }

    /** \brief Insert an element into the mesh
     * \param[in] gt               geometry type of the element (hopefully tetrahedron)
     * \param[in] vertexIndexSet   local indices of interpolatory vertices. Local with respect to the order the vertices were inserted
     * \param[in] order            interpolatory order of the element
     * \param[in] physicalTag      physical tag of the element
     *
     * Note: Even though we pass the globalId as a parameter from GMSH, it is a globalIndex for the set Elements+BoundaryFaces,
     * therefore obtaining globalIndex for elements from it is not possible, it will have to be communicated
     *
     * */
    void insertElement(
    	Dune::GeometryType gt,
    	const std::vector<LocalIndexType> & vertexIndexSet,
    	InterpolatoryOrderType order,
    	PhysicalTagType physicalTag)
    {
        if (!gt.isTetrahedron() || (vertexIndexSet.size() != Dune::CurvilinearGeometryHelper::dofPerOrder(gt, order)))  {
        	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_ERROR>( __FILE__, __LINE__, "CurvilinearGridConstructor: insertElement() unexpected element type or number of interpolatory points");
            DUNE_THROW(Dune::IOError, "CurvilinearGrid: insertElement() unexpected element type or number of interpolatory points");
        }

        EntityStorage thisElement;

        thisElement.geometryType = gt;
        thisElement.globalIndex = 0;        // At this stage globalIndex is not known yet
        thisElement.structuralType = GridStorageType::PartitionType::Internal;
        thisElement.interpOrder = order;
        thisElement.physicalTag = physicalTag;
        thisElement.vertexIndexSet = vertexIndexSet;

        LocalIndexType thisLocalIndex = gridstorage_.element_.size();
        gridstorage_.element_.push_back(thisElement);

        std::stringstream log_stream;
        log_stream << "CurvilinearGridConstructor: Inserted Element Type=" << Dune::CurvilinearGeometryHelper::geometryName(gt);
        log_stream << " LocalIndex=" << thisLocalIndex;
        log_stream << " Order=" << order;
        log_stream << " PhysicalTag=" << physicalTag;
        log_stream << " VertexIndices=(" << Dune::VectorHelper::vector2string(vertexIndexSet) << ")";
        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, log_stream.str());
    }

    /** Insert a boundary segment into the mesh
     *
     *     Note: It is expected that all faces - domain an process boundaries - are inserted by the factory before finalising
     *     Note: Only domain boundary faces have initial globalId given by GMSH. Therefore, we ignore it, and generate our own
     *     globalId for all faces at a later stage.
     *
     *  \param[in] gt                       geometry type of the face (should be a triangle)
     *  \param[in] associatedElementIndex   local index of the element this face is associated to
     *  \param[in] vertexIndexSet           local indices of the interpolatory vertices of this face
     *  \param[in] order                    interpolatory order of the face
     *  \param[in] physicalTag              physical tag of the element (material property)
     *
     * */

    void insertBoundarySegment(
    	Dune::GeometryType gt,
    	LocalIndexType associatedElementIndex,
    	const std::vector<LocalIndexType> & vertexIndexSet,
    	InterpolatoryOrderType order,
    	PhysicalTagType physicalTag)
    {
        if (!gt.isTriangle() || (vertexIndexSet.size() != Dune::CurvilinearGeometryHelper::dofPerOrder(gt, order)))  {
        	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_ERROR>( __FILE__, __LINE__, "CurvilinearGridConstructor: insertBoundarySegment() unexpected number of interpolatory points");
            DUNE_THROW(Dune::IOError, "CurvilinearGridConstructor: insertBoundarySegment() unexpected number of interpolatory points");
        }


        // Get corners of this face
        // **********************************************************************************
        std::vector<LocalIndexType> faceCorners = Dune::CurvilinearGeometryHelper::entityVertexCornerSubset<ct, 2>(gt, vertexIndexSet, order);
        FaceKey thisFaceKey;
        thisFaceKey.node0 = gridstorage_.point_[faceCorners[0]].globalIndex;
        thisFaceKey.node1 = gridstorage_.point_[faceCorners[1]].globalIndex;
        thisFaceKey.node2 = gridstorage_.point_[faceCorners[2]].globalIndex;

        // Sort in ascending order
        thisFaceKey.sort();


        // Take associated element, get all its corners, get all keys, compare to face key
        // **********************************************************************************
        std::vector<LocalIndexType> elementCorners = Dune::CurvilinearGeometryHelper::entityVertexCornerSubset<ct, 3>(
        		gridstorage_.element_[associatedElementIndex].geometryType,
        		gridstorage_.element_[associatedElementIndex].vertexIndexSet,
        		gridstorage_.element_[associatedElementIndex].interpOrder
        );



        // Search for the face among subentities of the element
        // **********************************************************************************
        int j = 0;
        int nFacePerTetrahedron = 4;
        bool found_face = false;

        while (!found_face)
        {
            if (j == nFacePerTetrahedron)  {
            	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_ERROR>( __FILE__, __LINE__, "CurvilinearGridConstructor: insertBoundarySegment() did not find the face in the associated element");
                DUNE_THROW(Dune::IOError, "CurvilinearGrid: insertBoundarySegment() did not find the face in the associated element");
            }


            // Get internal indices of the corners of this face wrt its associated element
            Dune::GeometryType assocElementGeometryType = gridstorage_.element_[associatedElementIndex].geometryType;
            InternalIndexType node0SubIndex = Dune::ReferenceElements<ct, cdim>::general(assocElementGeometryType).subEntity(j, FACE_CODIM, 0, VERTEX_CODIM);
            InternalIndexType node1SubIndex = Dune::ReferenceElements<ct, cdim>::general(assocElementGeometryType).subEntity(j, FACE_CODIM, 1, VERTEX_CODIM);
            InternalIndexType node2SubIndex = Dune::ReferenceElements<ct, cdim>::general(assocElementGeometryType).subEntity(j, FACE_CODIM, 2, VERTEX_CODIM);

            // Define (key = sorted localIndices of corners)
            FaceKey thisKey;
            thisKey.node0 = gridstorage_.point_[elementCorners[node0SubIndex]].globalIndex;
            thisKey.node1 = gridstorage_.point_[elementCorners[node1SubIndex]].globalIndex;
            thisKey.node2 = gridstorage_.point_[elementCorners[node2SubIndex]].globalIndex;


            // Sort in ascending order
            thisKey.sort();

            // By comparison find internalIndex of this face
            if (thisKey == thisFaceKey)
            {
                found_face = true;

                // Store face in a Map (key -> faceIndex) for constructor purposes
                // Also create a domain boundary index for future indexing
                LocalIndexType localFaceIndex = gridstorage_.face_.size();
                LocalIndexType localFaceDBIndex = gridstorage_.boundarySegmentIndexMap_.size();
                domainBoundaryFaceKey2LocalIndexMap_[thisKey] = localFaceIndex;
                gridstorage_.boundarySegmentIndexMap_[localFaceIndex] = localFaceDBIndex;

                // Store Vector (faceId -> associated element)
                FaceStorage thisFaceAsSubentity;
                thisFaceAsSubentity.geometryType = gt;
                thisFaceAsSubentity.globalIndex  = 0;                  // At this stage the globalId is not known yet
                thisFaceAsSubentity.structuralType = DomainBoundaryType;
                thisFaceAsSubentity.element1Index = associatedElementIndex;
                thisFaceAsSubentity.element2Index = -1;              // Boundary Segments do not have a 2nd neighbor
                thisFaceAsSubentity.element1SubentityIndex = j;
                thisFaceAsSubentity.physicalTag = physicalTag;    // Here physical tag is very important as it need not match the tag of the element

                gridstorage_.face_.push_back(thisFaceAsSubentity);


                std::stringstream log_stream;
                log_stream << "CurvilinearGridConstructor: Inserted BoundarySegment Type=" << Dune::CurvilinearGeometryHelper::geometryName(gt);
                log_stream << " LocalIndex=" << gridstorage_.face_.size()-1;
                log_stream << " Order=" << order;
                log_stream << " PhysicalTag=" << physicalTag;
                log_stream << " AssociatedElementIndex=" << associatedElementIndex;
                log_stream << " InternalSubentityIndex=" << j;
                loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, log_stream.str());
            }

            j++;
        }
    }


    /* ***************************************************************************
     * Section: Constructing the mesh
     * ***************************************************************************/

    /** Calls the subroutines for transforming the inserted data into a functional mesh.
     *  Note: It is expected, that all necessary data (vertices, elements and boundary segments) have been added before this function is called.
     * */

    // 0) Compute process bounding box - coordinates of a box - compact, but sufficient to fit in the elements of this process
    // 1) Generate all edges
    // 2) Generate all faces - process boundary and internal (as domain boundary already inserted)
    // 3) Generate process boundary corner set - from process boundary faces
    // 4) Generate process boundary edge set  - from process boundary faces
    // 5) Generate Global Indices
    //        - Communicate neighbor ranks for all process boundary corners
    //      - Compute neighbor ranks for all edges from set intersection. Mark as owned if your rank is top
    //      - Compute neighbor ranks for all faces from set intersection. Mark as owned if your rank is top
    //      - [Nothing changes here] Enumerate owned
    //      - Communicate non-owned. No such thing as complex vertex at this stage - just communicate to all other processes who own it
    // 6) Generate Ghost Elements
    // 7) Construct OCTree

    void generateMesh() {

        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Initializing mesh");

        // Construct missing parts of the mesh
        // ************************************************************
        gridstorage_.nEntityTotal_[EDGE_CODIM] = 0;  // Will be updated later
        gridstorage_.nEntityTotal_[FACE_CODIM] = 0;  // Will be updated later

        generateEdges();
        generateFaces();
        markBoundaryVertexStructuralType();
        markBoundaryEdgeStructuralType();

        if (size_ > 1)
        {
#if HAVE_MPI
        	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Assembling Parallel Grid");

        	// Parallel case
            generateGlobalIndices();
            if (gridstorage_.withGhostElements_)
            {
            	GridGhostConstructor ghostConstructor(gridstorage_, mpihelper_, loggingmessage_);
            	ghostConstructor.generate();
            }
#endif
        }
        else
        {
        	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Assembling Serial Grid");
            // Serial case:
            // * Boundary Neighbors not necessary, since all boundaries are domain boundaries
            // * No ghost elements, even if requested by user
            // * Fake globalIndex by making it equal to localIndex

            gridstorage_.withGhostElements_ = false;
            for (int i = 0; i < gridstorage_.edge_.size();    i++)  { gridstorage_.edge_[i].globalIndex = i;     gridstorage_.entityIndexMap_[EDGE_CODIM][i] = i; }
            for (int i = 0; i < gridstorage_.face_.size();    i++)  { gridstorage_.face_[i].globalIndex = i;     gridstorage_.entityIndexMap_[FACE_CODIM][i] = i; }
            for (int i = 0; i < gridstorage_.element_.size(); i++)  { gridstorage_.element_[i].globalIndex = i;  gridstorage_.entityIndexMap_[ELEMENT_CODIM][i] = i;  gridstorage_.entityInternalIndexSet_[ELEMENT_CODIM].insert(i); }

            for (FaceMapIterator faceIter = internalFaceKey2LocalIndexMap_.begin(); faceIter != internalFaceKey2LocalIndexMap_.end(); faceIter++)              { LocalIndexType localIndex = (*faceIter).second;        gridstorage_.entityInternalIndexSet_[FACE_CODIM].insert(localIndex); }
            for (FaceMapIterator faceIter = domainBoundaryFaceKey2LocalIndexMap_.begin(); faceIter != domainBoundaryFaceKey2LocalIndexMap_.end(); faceIter++)  { LocalIndexType localIndex = (*faceIter).second;        gridstorage_.entityDomainBoundaryIndexSet_[FACE_CODIM].insert(localIndex); }

            gridstorage_.nEntityTotal_[EDGE_CODIM] = gridstorage_.edge_.size();
            gridstorage_.nEntityTotal_[FACE_CODIM] = gridstorage_.face_.size();
        }

        // Deletes all temporary memory
        // ************************************************************
        edgeKey2LocalIndexMap_.clear();
        internalFaceKey2LocalIndexMap_.clear();
        domainBoundaryFaceKey2LocalIndexMap_.clear();
        processBoundaryFaceKey2LocalIndexMap_.clear();


        // Create sets that will be used for iteration over the map
        // Create maps of all entity subsets which can be communicated over
        // Find neighbor ranks for all entities that can be communicated over
        // ************************************************************
        GridPostConstructor postConstructor(gridstorage_, gridbase_, mpihelper_, loggingmessage_);
        postConstructor.generateCornerIndex();
        postConstructor.generateIteratorSets();



        // The PB-PB communication interface is available by default. The below procedures enable the communication interfaces
        // involving ghost entities, and require existence of ghost entities
        if ((size_ > 1)&& gridstorage_.withGhostElements_)
        {
#if HAVE_MPI
            postConstructor.generateCommunicationMaps();
            postConstructor.communicateCommunicationEntityNeighborRanks();
#endif
        }


        // Construct OCTree
        // ************************************************************
        computeProcessBoundingBox();
        //constructOctree();
    }




protected:



    /* ***************************************************************************
     * Section: Constructing the mesh
     * ***************************************************************************/

    /** Compute center and extent (halved) of the bounding box of the whole mesh.
     *
     * FIXME: Current implementation underestimates bounding box due to curvilinear effects.
     * Suggested correction:
     * 1) Min and max over all interpolation points [oh, this is already done]
     * 2) Correction factor for Excess curvature
     * 2.1) For example, compute average element linear size, and take half of that, and enlarge extent by it
     *
     * */
    void computeProcessBoundingBox()
    {
        Vertex min = gridstorage_.point_[0].coord;
        Vertex max = min;

        for (int i = 1; i < gridstorage_.point_.size(); i ++) {
            min[0] = std::min(min[0], gridstorage_.point_[i].coord[0]);
            min[1] = std::min(min[1], gridstorage_.point_[i].coord[1]);
            min[2] = std::min(min[2], gridstorage_.point_[i].coord[2]);

            max[0] = std::max(max[0], gridstorage_.point_[i].coord[0]);
            max[1] = std::max(max[1], gridstorage_.point_[i].coord[1]);
            max[2] = std::max(max[2], gridstorage_.point_[i].coord[2]);
        }
        gridstorage_.boundingBoxCenter_ = min + max;  gridstorage_.boundingBoxCenter_ *= 0.5;
        gridstorage_.boundingBoxExtent_ = max - min;  gridstorage_.boundingBoxExtent_ *= 0.5;
    }


    /** \brief Generates all edges
     *
     * Algorithm:
     * 1) Loop over all elements
     * 2) For each element, construct EdgeKeys of all its edges, based on the corners of the element
     * 3) Add each EdgeKey to the map. If the edge did not exist before, give it a local index
     * 4) Note for each edge the local index of the parent that created it
     * 5) Mark edge local index as a subentity of containing element
     *
     * \note Orientation of elements ensured via Dune::ReferenceElement::subEntity()
     * [TODO]  Use more generic functions when extending to any mesh other than tetrahedral
     *
     * */
    void generateEdges()
    {
    	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Started constructing edges");

    	// Init the subentity index vector
    	int nEdgePerTetrahedron = 6;
    	int nElem = gridstorage_.element_.size();
    	gridstorage_.elementSubentityCodim2_ = std::vector<std::vector<int> > (nElem, std::vector<int>(nEdgePerTetrahedron));

        // Loop over all elements and their edges
        for (int iElem = 0; iElem < nElem; iElem++)
        {
        	EntityStorage & thisElem = gridstorage_.element_[iElem];
            std::vector<LocalIndexType> elementCornerLocalIndexSet = Dune::CurvilinearGeometryHelper::entityVertexCornerSubset<ct, 3>(thisElem.geometryType, thisElem.vertexIndexSet, thisElem.interpOrder);

            for (int iEdge = 0; iEdge < nEdgePerTetrahedron; iEdge++)
            {

                // Get internal indices of the corners of this face wrt its associated element
                InternalIndexType node0SubIndex = Dune::ReferenceElements<ct, cdim>::general(thisElem.geometryType).subEntity(iEdge, EDGE_CODIM, 0, VERTEX_CODIM);
                InternalIndexType node1SubIndex = Dune::ReferenceElements<ct, cdim>::general(thisElem.geometryType).subEntity(iEdge, EDGE_CODIM, 1, VERTEX_CODIM);

                // Define (key = sorted globalIndices of corners)
                EdgeKey thisKey;
                thisKey.node0 = gridstorage_.point_[elementCornerLocalIndexSet[node0SubIndex]].globalIndex;
                thisKey.node1 = gridstorage_.point_[elementCornerLocalIndexSet[node1SubIndex]].globalIndex;

                // Sort in ascending order
                thisKey.sort();
                //std::cout << "process_" << rank_ << "GenerateEdgeKey=(" << thisKey.node0 << "," << thisKey.node1 << ")" << std::endl;


                EdgeMapIterator edgeIter = edgeKey2LocalIndexMap_.find(thisKey);

                // If this edge has not been added already, add it to the map
                // Find local index of this edge and note it as this element subentity
                if (edgeIter == edgeKey2LocalIndexMap_.end())
                {
                    // Store map (key -> edgeIndex)
                	LocalIndexType localEdgeIndex = gridstorage_.edge_.size();
                    edgeKey2LocalIndexMap_[thisKey] = localEdgeIndex;

                    // Store vector (edgeId -> elemId + edgeElemIndex)
                    // Note: Edges do not have physical tag at all so we do not even store it
                    EdgeStorage thisEdge;
                    thisEdge.globalIndex = 0;        // GlobalId for edge determined later using global communication
                    thisEdge.structuralType = GridStorageType::PartitionType::Internal;   // For process boundaries will be redefined later
                    thisEdge.elementIndex = iElem;
                    thisEdge.subentityIndex = iEdge;

                    // Log output
                    std::stringstream log_stream;
                    log_stream << "CurvilinearGridConstructor: Added Edge";
                    log_stream << " LocalIndex=" << localEdgeIndex;
                    log_stream << " AssociatedElementIndex=" << iElem;
                    log_stream << " InternalSubentityIndex=" << iEdge;
                    loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, log_stream.str());

                    gridstorage_.edge_.push_back(thisEdge);

                    gridstorage_.elementSubentityCodim2_[iElem][iEdge] = localEdgeIndex;
                } else {
                	gridstorage_.elementSubentityCodim2_[iElem][iEdge] = (*edgeIter).second;
                }
            }
        }
    	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Finished constructing edges");
    }


    /** Generates Internal and ProcessBoundary Faces. (!!!) Assumes that all Domain Boundary Faces have been added.
     *
     * Algorithm:
     * 1) Loop over all elements
     * 2) For each element, construct FaceKeys of all its edges, based on the corners of the element
     * 3) Add each FaceKey to a temporary map, unless it is already in the DomainBoundary map. In the latter case, assign its local index as the containing element subentity
     * 4) For each FaceKey in the temporary map, check if it has been added once or twice, this determines if it is a process boundary or an internal face.
     * 5) Add each face to the total face map and to the corresponding structure map (internal/processboundary). Mark face local index as a subentity of containing element
     * 6) Note for each face the local index of the parent that created it
     * 7) For each process boundary face, create map from localFaceIndex to dummy creation index
     * 8) Resize the future neighbor rank array with the number of PB faces
     *
     * \note Orientation of elements ensured via Dune::ReferenceElement::subEntity()
     * [TODO]  If assume mesh with non-uniform p-refinement, it may make sense to point the face to the element which has the higher refinement
     *
     * */
    void generateFaces()
    {
    	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Started generating faces");

    	// Init the subentity index vector
        int nFacePerTetrahedron = 4;
    	int nElem = gridstorage_.element_.size();
    	gridstorage_.elementSubentityCodim1_ = std::vector<std::vector<LocalIndexType> > (nElem, std::vector<LocalIndexType>(nFacePerTetrahedron));

        typedef std::map<FaceKey, std::vector<int>> tmpFace2InfoMap;
        typedef typename tmpFace2InfoMap::iterator  tmpMapIterator;
        tmpFace2InfoMap tmpFaceMap;

        // Loop over all elements and their faces
        for (int iElem = 0; iElem < nElem; iElem++)
        {
        	EntityStorage & thisElem = gridstorage_.element_[iElem];
            std::vector<int> elementCornerLocalIndexSet = Dune::CurvilinearGeometryHelper::entityVertexCornerSubset<ct, 3>(thisElem.geometryType, thisElem.vertexIndexSet, thisElem.interpOrder);

            // Store info for all faces except of domain boundaries
            // Store it in a map, not to store internal faces twice
            for (int iFace = 0; iFace < nFacePerTetrahedron; iFace++)
            {

                // Get internal indices of the corners of this face wrt its associated element
                InternalIndexType node0SubIndex = Dune::ReferenceElements<ct, cdim>::general(thisElem.geometryType).subEntity(iFace, FACE_CODIM, 0, VERTEX_CODIM);
                InternalIndexType node1SubIndex = Dune::ReferenceElements<ct, cdim>::general(thisElem.geometryType).subEntity(iFace, FACE_CODIM, 1, VERTEX_CODIM);
                InternalIndexType node2SubIndex = Dune::ReferenceElements<ct, cdim>::general(thisElem.geometryType).subEntity(iFace, FACE_CODIM, 2, VERTEX_CODIM);

                // Define (key = sorted globalIndices of corners)
                FaceKey thisKey;
                thisKey.node0 = gridstorage_.point_[elementCornerLocalIndexSet[node0SubIndex]].globalIndex;
                thisKey.node1 = gridstorage_.point_[elementCornerLocalIndexSet[node1SubIndex]].globalIndex;
                thisKey.node2 = gridstorage_.point_[elementCornerLocalIndexSet[node2SubIndex]].globalIndex;

                // Sort in ascending order
                thisKey.sort();


                FaceMapIterator faceIter = domainBoundaryFaceKey2LocalIndexMap_.find(thisKey);

                // Mark this face for creation if it is not an already existing Domain Boundary
                // Otherwise note its local index
                if (faceIter == domainBoundaryFaceKey2LocalIndexMap_.end())
                {
                    std::vector<int> connectedFaceInfo;
                    tmpMapIterator iter = tmpFaceMap.find(thisKey);

                    if (iter != tmpFaceMap.end()) { connectedFaceInfo = std::vector<int> ( (*iter).second ); }

                    connectedFaceInfo.push_back(iElem);
                    connectedFaceInfo.push_back(iFace);

                    tmpFaceMap[thisKey] = connectedFaceInfo;


                    std::stringstream log_stream;
                    log_stream << "CurvilinearGridConstructor: Adding FaceKey=(" << thisKey.node0 << ", " << thisKey.node1 << ", " << thisKey.node2 << ") attached to total of " << connectedFaceInfo.size() / 2 << " elements";
                    loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, log_stream.str());
                } else {
                	LocalIndexType localFaceIndex = (*faceIter).second;
                	gridstorage_.elementSubentityCodim1_[iElem][iFace] = localFaceIndex;
                }


            }
        }

        // Add internal and process boundary faces to the mesh
        for (tmpMapIterator iter = tmpFaceMap.begin();  iter != tmpFaceMap.end();  iter++)
        {
            FaceStorage thisFace;
            LocalIndexType localFaceIndex = gridstorage_.face_.size();
            std::vector<int> connectedFaceInfo = (*iter).second;

            LocalIndexType    thisAssociatedElementIndex = connectedFaceInfo[0];
            InternalIndexType thisFaceSubentityIndex = connectedFaceInfo[1];

            // Store the face local index as element subentity
            gridstorage_.elementSubentityCodim1_[thisAssociatedElementIndex][thisFaceSubentityIndex] = localFaceIndex;

            // Recover parental information for this face
            thisFace.globalIndex = 0;       // GlobalId is defined at a later stage
            thisFace.element1Index = thisAssociatedElementIndex;
            thisFace.element1SubentityIndex = thisFaceSubentityIndex;
            thisFace.physicalTag = -1;    // At the moment physicalTag of an internal face is not defined as it could be inbetween two different elements

            // Find the geometry type of the face from its parent face
            Dune::GeometryType parentGeometry = gridstorage_.element_[thisFace.element1Index].geometryType;
            thisFace.geometryType = Dune::ReferenceElements<ct, cdim>::general(parentGeometry).type(thisFace.element1SubentityIndex, 1);


            std::stringstream log_stream;
            log_stream << "CurvilinearGridConstructor: Added Face";
            log_stream << " LocalIndex=" << localFaceIndex;
            log_stream << " AssociatedElementIndex=" << thisFace.element1Index;
            log_stream << " InternalSubentityIndex=" << thisFace.element1SubentityIndex;

            // Store internal, domain and process boundaries separately for faster iterators
            if (connectedFaceInfo.size() == 2)
            {
                thisFace.structuralType = ProcessBoundaryType;
                processBoundaryFaceKey2LocalIndexMap_[(*iter).first] = localFaceIndex;  // Store Map (key -> faceIndex)
                thisFace.element2Index = 0;                                             // Eventually this will be the Ghost Element Index

                // Add this face to the process boundary map
                LocalIndexType thisFaceLocalPBIndex = gridstorage_.processBoundaryIndexMap_[FACE_CODIM].size();
                gridstorage_.processBoundaryIndexMap_[FACE_CODIM][localFaceIndex] = thisFaceLocalPBIndex;

                log_stream << " StructuralType=processBoundary";
            }
            else
            {
            	// In this case 2nd neighbouring element needs to be mapped
                LocalIndexType    thisAssociatedElement2Index = connectedFaceInfo[2];
                InternalIndexType thisFaceSubentityIndex2 = connectedFaceInfo[3];
                gridstorage_.elementSubentityCodim1_[thisAssociatedElement2Index][thisFaceSubentityIndex2] = localFaceIndex;
                log_stream << " AssociatedElement2Index=" << thisAssociatedElement2Index;
                log_stream << " InternalSubentityIndex2=" << thisFaceSubentityIndex2;

                // Add this face to the internal map
                thisFace.structuralType = InternalType;
                internalFaceKey2LocalIndexMap_[(*iter).first] = localFaceIndex;    // Store Map (key -> faceIndex)
                thisFace.element2Index = thisAssociatedElement2Index;              // This is the 2nd neighbor of this internal face
                log_stream << " StructuralType=internal";
            }

            // Add face to the mesh
            loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, log_stream.str());
            gridstorage_.face_.push_back(thisFace);

            // Update neighbor index storage size
            int nPBFace = gridstorage_.processBoundaryIndexMap_[FACE_CODIM].size();
            gridstorage_.PB2PBNeighborRank_[FACE_CODIM].resize(nPBFace);
        }

        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Finished generating faces");
    }




    /** Mark correct structural type for all boundary vertices, by making it equal to the
     * structural type of the containing face. When sharing several structural types, the
     * priority order is ProcessBoundary > DomainBoundary > Internal
     *
     * Then, give all process boundary corners a unique index by storing them in the
     * processBoundaryIndexMap_[VERTEX_CODIM]. These are later used for communication
     *
     * [TODO] To make work with arbitrary geometries, must replace numbers, as well as make FaceKeys more abstract
     *
     * */
    void markBoundaryVertexStructuralType()
    {
    	// Mark domain boundary vertices
    	// ********************************************************
    	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Started marking domain boundary vertices");
        for (FaceMapIterator faceIter = domainBoundaryFaceKey2LocalIndexMap_.begin(); faceIter != domainBoundaryFaceKey2LocalIndexMap_.end(); faceIter++)
        {
        	EntityStorage thisFace = gridbase_.entityData(1, (*faceIter).second);

            for (int i = 0; i < thisFace.vertexIndexSet.size(); i++)
            {
            	LocalIndexType thisVertexLocalIndex = thisFace.vertexIndexSet[i];
            	gridstorage_.point_[thisVertexLocalIndex].structuralType = GridStorageType::PartitionType::DomainBoundary;
            }
        }

    	// Mark process boundary vertices - overwrite domain boundary vertices where necessary
        // ********************************************************
        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Started marking process boundary vertices");
        for (FaceMapIterator faceIter = processBoundaryFaceKey2LocalIndexMap_.begin(); faceIter != processBoundaryFaceKey2LocalIndexMap_.end(); faceIter++)
        {
        	EntityStorage thisFace = gridbase_.entityData(1, (*faceIter).second);

            for (int i = 0; i < thisFace.vertexIndexSet.size(); i++)
            {
            	LocalIndexType thisVertexLocalIndex = thisFace.vertexIndexSet[i];
            	gridstorage_.point_[thisVertexLocalIndex].structuralType = GridStorageType::PartitionType::ProcessBoundary;
            }
        }

        // Construct the set of process boundary corners - corners necessary to make process boundary faces on this process
        // ********************************************************
        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Started generating BoundaryCorneers");
        for (FaceMapIterator faceIter = processBoundaryFaceKey2LocalIndexMap_.begin(); faceIter != processBoundaryFaceKey2LocalIndexMap_.end(); faceIter++)
        {
            // Get global indices of the associated vertices from the map
            FaceKey thisFaceKey = (*faceIter).first;

            GlobalIndexType thisVertexKey[3] = {thisFaceKey.node0, thisFaceKey.node1, thisFaceKey.node2};

            for (int i = 0; i < 3; i++)
            {
            	LocalIndexType thisCornerLocalIndex = gridstorage_.entityIndexMap_[VERTEX_CODIM][thisVertexKey[i]];

            	// For each vertex, if it has not yet been added to the map, add it, mapping to a new entry
                // in the processBoundaryCornerNeighborRank_
                if (gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM].find(thisCornerLocalIndex) ==
                	gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM].end())
                {
                	// DO NOT combine these two lines into one - then the map returns the size after insertion which is +1 to expected
                    LocalIndexType processBoundaryCornerIndex = gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM].size();
                	gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM][thisCornerLocalIndex] = processBoundaryCornerIndex;

                    loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: -- Adding boundary corner GlobalIndex=" + std::to_string(thisVertexKey[i]));
                }
            }
        }

        // Resize the neighbor rank vector such that it can store for each process boundary corner
        gridstorage_.PB2PBNeighborRank_[VERTEX_CODIM].resize(gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM].size());

        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Finished marking vertices and corners");
    }


    /** Mark correct structural type for all boundary edges, by making it equal to the
     * structural type of the containing face. When sharing several structural types, the
     * priority order is ProcessBoundary > DomainBoundary > Internal
     *
     * Then, give all process boundary edges a unique index by storing them in the
     * processBoundaryIndexMap_[EDGE_CODIM]. These are later used for communication
     *
     * [TODO] To make work with arbitrary geometries, must replace numbers, as well as make FaceKeys more abstract
     *
     * */
    void markBoundaryEdgeStructuralType()
    {
        // Construct the set of EdgeKeys corresponding to edges of processBoundaries
        // ********************************************************
    	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Started marking domain boundary edges");
        for (FaceMapIterator faceIter = domainBoundaryFaceKey2LocalIndexMap_.begin(); faceIter != domainBoundaryFaceKey2LocalIndexMap_.end(); faceIter++)
        {
            // Get global indices of the associated vertices from the map
            FaceKey thisFaceKey = (*faceIter).first;
            EdgeKey thisEdgeKey[3];

            thisEdgeKey[0].node0 = thisFaceKey.node0;  thisEdgeKey[0].node1 = thisFaceKey.node1;
            thisEdgeKey[1].node0 = thisFaceKey.node0;  thisEdgeKey[1].node1 = thisFaceKey.node2;
            thisEdgeKey[2].node0 = thisFaceKey.node1;  thisEdgeKey[2].node1 = thisFaceKey.node2;

            for (int i = 0; i < 3; i++)
            {
            	LocalIndexType thisEdgeLocalIndex = edgeKey2LocalIndexMap_[thisEdgeKey[i]];
            	gridstorage_.edge_[thisEdgeLocalIndex].structuralType = GridStorageType::PartitionType::DomainBoundary;

                std::stringstream log_stream;
                log_stream << "CurvilinearGridConstructor: -- From face index= " << (*faceIter).second << " marking domain boundary edge index=" << thisEdgeLocalIndex << " EdgeKey= (" << thisEdgeKey[i].node0 << ", " << thisEdgeKey[i].node1 << ")";
                loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, log_stream.str());
            }
        }

        // Construct the set of EdgeKeys corresponding to edges of processBoundaries
        // ********************************************************
        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Started marking process boundary edges");
        for (FaceMapIterator faceIter = processBoundaryFaceKey2LocalIndexMap_.begin(); faceIter != processBoundaryFaceKey2LocalIndexMap_.end(); faceIter++)
        {
            // Get info of this face wrt associated element
            FaceKey thisFaceKey = (*faceIter).first;
            LocalIndexType thisFaceLocalIndex = (*faceIter).second;
            LocalIndexType assocElementIndex = gridstorage_.face_[thisFaceLocalIndex].element1Index;
            InternalIndexType thisFaceSubIndex = gridstorage_.face_[thisFaceLocalIndex].element1SubentityIndex;
            Dune::GeometryType assocElemGT = gridstorage_.element_[assocElementIndex].geometryType;

            // Get global indices of the associated vertices from the map
            EdgeKey thisEdgeKey[3];
            thisEdgeKey[0].node0 = thisFaceKey.node0;  thisEdgeKey[0].node1 = thisFaceKey.node1;
            thisEdgeKey[1].node0 = thisFaceKey.node0;  thisEdgeKey[1].node1 = thisFaceKey.node2;
            thisEdgeKey[2].node0 = thisFaceKey.node1;  thisEdgeKey[2].node1 = thisFaceKey.node2;

            for (int i = 0; i < 3; i++)
            {
            	LocalIndexType thisEdgeLocalIndex = edgeKey2LocalIndexMap_[thisEdgeKey[i]];

            	// Mark this edge as subentity of associated boundary internal element
            	// It is strongly suspected that this operation is unnecessary, because it is already done at generateEdges()
            	//InternalIndexType thisEdgeSubIndex = Dune::ReferenceElements<ct, cdim>::general(assocElemGT).subEntity(thisFaceSubIndex, FACE_CODIM, i, EDGE_CODIM);
            	//gridstorage_.elementSubentityCodim2_[assocElementIndex][thisEdgeSubIndex] = thisEdgeLocalIndex;

                // For each vertex, if it has not yet been added to the map, add it, mapping to a new entry
                // in the processBoundaryCornerNeighborRank_
                if (gridstorage_.processBoundaryIndexMap_[EDGE_CODIM].find(thisEdgeLocalIndex) ==
                	gridstorage_.processBoundaryIndexMap_[EDGE_CODIM].end())
                {
                	// Change the structural type of this edge to ProcessBoundary
                	gridstorage_.edge_[thisEdgeLocalIndex].structuralType = GridStorageType::PartitionType::ProcessBoundary;

                	// Insert this edge into the local2local map
                	// DO NOT combine these two lines into one - then the map returns the size after insertion which is +1 to expected
                    LocalIndexType processBoundaryEdgeIndex = gridstorage_.processBoundaryIndexMap_[EDGE_CODIM].size();
                	gridstorage_.processBoundaryIndexMap_[EDGE_CODIM][thisEdgeLocalIndex] = processBoundaryEdgeIndex;

                    std::stringstream log_stream;
                    log_stream << "CurvilinearGridConstructor: -- From face index= " << thisFaceLocalIndex << " marking process boundary edge index=" << thisEdgeLocalIndex << " EdgeKey= (" << thisEdgeKey[i].node0 << ", " << thisEdgeKey[i].node1 << ")";
                    loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, log_stream.str());
                }
            }
        }

        // Resize the neighbor rank vector such that it can store for each process boundary corner
        gridstorage_.PB2PBNeighborRank_[EDGE_CODIM].resize(gridstorage_.processBoundaryIndexMap_[EDGE_CODIM].size());
        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Finished marking edges");
    }


    /** Generates Global Indices for Edges, Faces and Elements
     *
     * Algorithm:
     * 1) Communicate neighbour ranks associated with each process boundary corner
     * 2) Compute (provisional) neighbour ranks of PB edges and faces by intersection of ranks of associated PB corners
     * 2.1) Sometimes, an entity does not exist on a neighbouring process, even though all associated PB corners are present
     *      This only happens if the (provisional) number of neighbors is larger than 1 (complicated PB entity),
     *      because each PB entity must have at least 1 neighbour.
     * 2.2) For each complicated PB entity (edge/face), communicate EdgeKeys and FaceKeys to all provisional neighbours
     * 2.3) For each received key, reply to sender if such entity exists on this process or not
     * 2.4) Remove neighbour ranks mapping to non-existing entities
     *
     * 3) Find ownership of each edge and face. A shared entity is owned by the process with lowest rank
     * 4) Communicate number of edges and faces owned by each process to all
     * 5) Locally enumerate all edges, faces and elements owned by this process. That is, to assign them a global index
     * 5.1) Global index for edges starts at nVertexTotal+nEdgesOwnedBeforeMe.
     * 5.2) Global index for faces starts at nVertexTotal+nEdgeTotal+nFacesOwnedBeforeMe.
     * 5.3) Global index for elements starts at nVertexTotal+nEdgesTotal+nFacesTotal+nElementsOwnedBeforeMe. Note that each process owns all its elements since they are not shared.
     * 6) Communicate missing edge and face globalIndices
     * 6.1) By analysing entity neighbours, each process can compute how many how many global indices it needs to send and to receive to each other process
     * 6.2) Each process sends to each neighbour the shared entity global indices enumerated by this process and receives those enumerated by the neighbour process
     * 7) Fill in Global2Local maps. They are required for user functionality and for construction of GhostElements
     *
     * [TODO] Communication of corner neighbour ranks via allgather very inefficient. Try to find better algorithm
     * [TODO] MinRank-Ownership paradigm non-uniform. If ever becomes bottleneck, replace by XORRank-Ownership
     *
     * */
    void generateGlobalIndices()
    {
        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Constructing Global Indices");
        Dune::CollectiveCommunication<MPI_Comm> collective_comm = mpihelper_.getCollectiveCommunication();


        // 1) Communicate process ranks associated with each process boundary corner
        // 2) Compute neighbour ranks of PB edges and faces by intersection of ranks of associated PB corners
        //    Then eliminate non-existing entities generated this way
        // *************************************************************************
        globalCommunicateCornerNeighborRank();
        globalComputeEdgeNeighborRanks();
        globalComputeFaceNeighborRank();


        // 3) Get edges and faces on this process that are not owned by this process
        // *************************************************************************
        LocalIndexSet edgeNonOwned;  // LocalIndex of edge not owned by this process
        LocalIndexSet faceNonOwned;  // LocalIndex of face not owned by this process

        // Edges
        for (Local2LocalIterator edgeIter  = gridstorage_.processBoundaryIndexMap_[EDGE_CODIM].begin();
        		                 edgeIter != gridstorage_.processBoundaryIndexMap_[EDGE_CODIM].end(); edgeIter++ )
        {
        	LocalIndexType thisEdgeLocalIndex = (*edgeIter).first;
        	LocalIndexType thisEdgePBIndex    = (*edgeIter).second;

            int edgeOwnerCandidateRank = gridstorage_.PB2PBNeighborRank_[EDGE_CODIM][thisEdgePBIndex][0];
            if (edgeOwnerCandidateRank < rank_) { edgeNonOwned.insert(thisEdgeLocalIndex); }
        }

        // Faces
        for (Local2LocalIterator faceIter = gridstorage_.processBoundaryIndexMap_[FACE_CODIM].begin();
        		                 faceIter != gridstorage_.processBoundaryIndexMap_[FACE_CODIM].end(); faceIter++ )
        {
        	LocalIndexType thisFaceLocalIndex = (*faceIter).first;
        	LocalIndexType thisFacePBIndex = (*faceIter).second;
            int faceOwnerCandidateRank = gridstorage_.PB2PBNeighborRank_[FACE_CODIM][thisFacePBIndex][0];
            if (faceOwnerCandidateRank < rank_) { faceNonOwned.insert(thisFaceLocalIndex); }
        }

        int nEdgeOwned = gridstorage_.edge_.size() - edgeNonOwned.size();
        int nFaceOwned = gridstorage_.face_.size() - faceNonOwned.size();
        int elementsOwned = gridstorage_.element_.size();


        // 4) Communicate number of edges and faces owned by each process to all
        // *************************************************************************
        std::vector<int> edgesOnProcess(size_);      // owned edges [rank]
        std::vector<int> facesOnProcess(size_);      // owned faces [rank]
        std::vector<int> elementsOnProcess(size_);  // owned elements [rank]

        collective_comm.allgather (&nEdgeOwned, 1, reinterpret_cast<int*> (edgesOnProcess.data()));
        collective_comm.allgather (&nFaceOwned, 1, reinterpret_cast<int*> (facesOnProcess.data()));
        collective_comm.allgather (&elementsOwned, 1, reinterpret_cast<int*> (elementsOnProcess.data()));

        int edgesBeforeMe = 0;      // Sum(edgesOwned : rank < thisRank)
        int facesBeforeMe = 0;      // Sum(facesOwned : rank < thisRank)
        int elementsBeforeMe = 0;   // Sum(elementsOwned : rank < thisRank)

        for (int iProc = 0; iProc < size_; iProc++)
        {
        	gridstorage_.nEntityTotal_[EDGE_CODIM] += edgesOnProcess[iProc];
        	gridstorage_.nEntityTotal_[FACE_CODIM] += facesOnProcess[iProc];

            if (iProc < rank_)
            {
                edgesBeforeMe += edgesOnProcess[iProc];
                facesBeforeMe += facesOnProcess[iProc];
                elementsBeforeMe += elementsOnProcess[iProc];
            }
        }


        // 5) Enumerate all edges, faces and elements that you own
        // *************************************************************************

        GlobalIndexType iEdgeGlobalId = edgesBeforeMe;
        GlobalIndexType iFaceGlobalId = facesBeforeMe;

        // Enumerating elements is simply shifting the local index, since all elements on this process are owned by it
        for (int i = 0; i < gridstorage_.element_.size(); i++)  { gridstorage_.element_[i].globalIndex = elementsBeforeMe + i; }

        // Faces that are not shared with other processes are automatically owned by this process
        for (FaceMapIterator faceIter = internalFaceKey2LocalIndexMap_.begin(); faceIter != internalFaceKey2LocalIndexMap_.end(); faceIter++)              { LocalIndexType localIndex = (*faceIter).second;        gridstorage_.face_[localIndex].globalIndex = iFaceGlobalId++; }
        for (FaceMapIterator faceIter = domainBoundaryFaceKey2LocalIndexMap_.begin(); faceIter != domainBoundaryFaceKey2LocalIndexMap_.end(); faceIter++)  { LocalIndexType localIndex = (*faceIter).second;        gridstorage_.face_[localIndex].globalIndex = iFaceGlobalId++; }

        // This process owns this face if it is in the processInternalMap and if it is not in the non-owned face map
        for (FaceMapIterator faceIter = processBoundaryFaceKey2LocalIndexMap_.begin(); faceIter != processBoundaryFaceKey2LocalIndexMap_.end(); faceIter++)
        {
        	LocalIndexType thisFaceLocalIndex = (*faceIter).second;
            if (faceNonOwned.find(thisFaceLocalIndex) == faceNonOwned.end())  {
            	gridstorage_.face_[thisFaceLocalIndex].globalIndex = iFaceGlobalId++;
            }
        }

        // This process owns this edge if it is in the edgemap and if it is not in the non-owned edge map
        for (LocalIndexType iEdge = 0; iEdge < gridstorage_.edge_.size(); iEdge++)
        {
            if (edgeNonOwned.find(iEdge) == edgeNonOwned.end())  { gridstorage_.edge_[iEdge].globalIndex = iEdgeGlobalId++; }
            else { loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: do not own edge localIndex=" + std::to_string(iEdge) + " of type=" + gridstorage_.PartitonTypeName[gridstorage_.edge_[iEdge].structuralType] ); }
        }


        // 6) Communicate missing edge and face global indices to their corresponding neighbors. Receive them and assign.
        // *************************************************************************

        globalDistributeMissingEdgeGlobalIndex();
        globalDistributeMissingFaceGlobalIndex();


        // 7) Fill in Global2Local maps
        // *************************************************************************
        for (LocalIndexType iEdge = 0; iEdge < gridstorage_.edge_.size(); iEdge++)    { gridstorage_.entityIndexMap_[EDGE_CODIM][gridstorage_.edge_[iEdge].globalIndex] = iEdge; }

        for (LocalIndexType iFace = 0; iFace < gridstorage_.face_.size(); iFace++)    { gridstorage_.entityIndexMap_[FACE_CODIM][gridstorage_.face_[iFace].globalIndex] = iFace; }

        for (LocalIndexType iElem = 0; iElem < gridstorage_.element_.size(); iElem++) {	gridstorage_.entityIndexMap_[ELEMENT_CODIM][gridstorage_.element_[iElem].globalIndex] = iElem; }

        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Finished Constructing Global Indices");

    }


    /** Construct OCTree for locating tetrahedrons in mesh
     *
     * Algorithm:
     * 0) Uses OctreeNodes, structures that point to the original grid and contain the localIndex to the associated element.
     * 1) For each element of the mesh, create an OCTreeNode and add it to the octree
     * 1.1) When OctreeNode is created, it computes the rectangular box in which the element fits.
     * 1.2) When OctreeNode is added, it is recursively sinks into a tree starting with the root Octant.
     * 1.3) The root Octant contains the whole mesh of this process. Each octant contains 8 children by spitting itself into 8 smaller boxes.
     * 1.4) The element is fit into a given child octant if its bounding box intersects with that of the child. Thus the element may and will likely be contained inside several octants
     *
     * [FIXME] Revisit OCTree. See how the depth is regulated. See whether all candidates for this element from all octants are returned upon request
     *
     * [TODO]  Use standard logging message
     * [TODO]  Original octree has diagnostics output under #if 0, can append at later stage
     * [TODO]  Replace OCTree pointer by just an instance
     *
     * */
    void constructOctree() {
    	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Started OCTree construction");

        // bounding box of whole mesh
        Vertex center = gridstorage_.boundingBoxCenter_;
        Vertex extent = gridstorage_.boundingBoxExtent_;

        // octree length is the largest component of extent
        double length = extent[0];
        if (extent[1] > length)  { length = extent[1]; }
        if (extent[2] > length)  { length = extent[2]; }

        // construct LooseOctree with large max depth
        gridstorage_.octree_ = new CurvilinearLooseOctree(center, length, 100, mpihelper_, loggingmessage_);

        // loop over all tets and insert them in the octree
        for (int iElem = 0; iElem < gridstorage_.element_.size(); iElem ++)
        {
            NodeType* thisNode = new NodeType(gridbase_, iElem);
            gridstorage_.octree_->addNode(thisNode);
        }


        int maxDepth, nOctant, nNode;
        double avgNodeDepth;
        gridstorage_.octree_->statistics(maxDepth, avgNodeDepth, nOctant, nNode);

        std::stringstream outputString;
        outputString << "CurvilinearGridConstructor: Constructed OCTree MaxDepth=" << maxDepth;
        outputString << ", #octants=" << nOctant;
        outputString << ", #nodes=" << nNode;
        outputString << ", avg. node depth=" << avgNodeDepth;
        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, outputString.str());
    }


    /* ***************************************************************************
     * SubSection: Auxiliary methods of generateGlobalIndices()
     * ***************************************************************************/

    /** Communicate the process ranks of neighbor processes for all process boundary vertices
     *
     * Algorithm:
     * 1) collective_comm.max() - find the maximal number of process boundary corners per process
     * 2) Loop over maximal number of process boundary corners per process
     * 2.1) collective_comm.allgather() - communicate a global index of your process boundary corner to all other processes
     * 2.2) If all process boundary corners of this process have already been communicated, create and communicate fake indices (negative)
     *      This is necessary to keep the protocol going until the last process communicates all its corners
     * 2.3) From received vertices, select ones that are on this process, and mark the sender as the neighbor
     *
     * [TODO] Algorithm possibly inefficient. No ideas how to improve at the moment
     * * Every process boundary corners is communicated to all processes, but can be used only by few
     * * All processes have to wait until the process with the largest number of corners finishes communicating, since they could receive sth from it
     *
     * [TODO] Currently uses negative LocalIndex values to pass fake vertex info. If we want a uint impl, need to have different def. for fake vertex
     *
     *
     * */

    void globalCommunicateCornerNeighborRank ()
    {
    	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Started communicating corner process boundary neighbors");

        // 1) collective_comm.max() - find the maximal number of process boundary corners per process
        // ********************************************************

        Dune::CollectiveCommunication<MPI_Comm> collective_comm = mpihelper_.getCollectiveCommunication();

        // Reserve memory for saving ranks associated to process boundary corners
        LocalIndexType thisProcessBoundarySize = gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM].size();

        LocalIndexType maxProcessBoundarySize = collective_comm.max(thisProcessBoundarySize);

        // 2) collective_comm.allgather() - communicate global index of your process boundary corner to all other processes
        // Repeat this process until every process has communicated all its corners
        // If you run out of corners, communicate fake corners
        // ********************************************************
        Local2LocalIterator procCornerIter = gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM].begin();

        for (LocalIndexType iCorner = 0; iCorner < maxProcessBoundarySize; iCorner++)
        {
            // If all process boundary corners have been sent start sending fake corners
        	GlobalIndexType thisCornerGlobalIndex = -1;
        	if (iCorner < thisProcessBoundarySize)
        	{
        		LocalIndexType thisCornerLocalIndex = (*(procCornerIter++)).first;
        		thisCornerGlobalIndex = gridstorage_.point_[thisCornerLocalIndex].globalIndex;
        	}

            // Communicate global indices
            std::vector<LocalIndexType> processBoundaryGlobalIndexSet (size_);
            collective_comm.allgather(&thisCornerGlobalIndex, 1, reinterpret_cast<LocalIndexType*> (processBoundaryGlobalIndexSet.data()) );

            // Loop over corners sent by other processes. If this corner present, note its sender rank
            for (int iProc = 0; iProc < size_; iProc++)
            {
                // Only consider non-fake corners sent by other processes
                if ((iProc != rank_) && (processBoundaryGlobalIndexSet[iProc] >= 0))
                {
                    // Attempt to find this corner global id among process boundary corners of this process
                	GlobalIndexType thisCornerGlobalIndex = processBoundaryGlobalIndexSet[iProc];
                	Global2LocalIterator tmpIter = gridstorage_.entityIndexMap_[VERTEX_CODIM].find(thisCornerGlobalIndex);

                    // If this corner is present, note its sender process
                    if (tmpIter != gridstorage_.entityIndexMap_[VERTEX_CODIM].end()) {
                    	LocalIndexType thisCornerLocalIndex = (*tmpIter).second;
                    	LocalIndexType thisCornerLocalPBIndex = gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM][thisCornerLocalIndex];
                    	gridstorage_.PB2PBNeighborRank_[VERTEX_CODIM][thisCornerLocalPBIndex].push_back(iProc);
                    }
                }
            }
        }

        // 3) Sort all neighbor rank sets, to accelerate set intersection algorithm in future
        // ********************************************************
        for (int i = 0; i < gridstorage_.PB2PBNeighborRank_[VERTEX_CODIM].size(); i++)
        {
            std::sort(gridstorage_.PB2PBNeighborRank_[VERTEX_CODIM][i].begin(), gridstorage_.PB2PBNeighborRank_[VERTEX_CODIM][i].end());
        }


        // Testing output
        std::stringstream log_stream;
        log_stream << "CurvilinearGridConstructor: -- Process boundary corner";
        for (Local2LocalIterator cornerIter = gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM].begin(); cornerIter != gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM].end(); cornerIter++)
        {
        	log_stream << " GlobalIndex=" << (*cornerIter).first;
        	log_stream << " has Neighbors=(" << Dune::VectorHelper::vector2string(gridstorage_.PB2PBNeighborRank_[VERTEX_CODIM][(*cornerIter).second]) << ")";
        }
        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, log_stream.str());

        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Finished corner process boundary neighbors");
    }


    /** Compute the process ranks of neighbor processes for all process boundary edges
     *
     * Algorithm:
     * 1) Loop over all edges in the edge map
     * 1.1) For each corner in the EdgeKey get associated neighbor ranks from provided vertex neighbor ranks
     * 1.2) Perform intersection on the two sets
     * 1.3) Following edge map write that intersection to the output array
     *
     * */
    void globalComputeEdgeNeighborRanks()
    {
    	MPI_Comm comm = Dune::MPIHelper::getCommunicator();
    	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Started computing edge process boundary neighbors");


        // For each process stores the set of edge indices local to this process, which are shared between more than two processes in total
    	typedef std::pair<LocalIndexType, EdgeKey> TmpEdgeData;

        std::vector<std::vector<TmpEdgeData > > neighborProcessComplicatedEdgePBLocalIndex(size_);

    	// 1) Compute neighbor ranks for each process boundary edge by intersecting neighbor ranks of its corners
    	// *************************************************************************************************************
        for (Local2LocalIterator edgeIter = gridstorage_.processBoundaryIndexMap_[EDGE_CODIM].begin();
        		                 edgeIter != gridstorage_.processBoundaryIndexMap_[EDGE_CODIM].end(); edgeIter++ )
        {
            LocalIndexType thisEdgeLocalIndex = (*edgeIter).first;
            LocalIndexType thisPBEdgeLocalIndex = (*edgeIter).second;

            // Get corners of the edge
            std::vector<LocalIndexType> thisCornerLocalIndices = gridbase_.entityCornerLocalIndex(EDGE_CODIM, thisEdgeLocalIndex);

            EdgeKey thisEdgeKey;
            thisEdgeKey.node0 = gridstorage_.point_[thisCornerLocalIndices[0]].globalIndex;
            thisEdgeKey.node1 = gridstorage_.point_[thisCornerLocalIndices[1]].globalIndex;
            thisEdgeKey.sort();

            // Get neighbor processes associated with each corner
            LocalIndexType thisVertexPBIndex0 = gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM][thisCornerLocalIndices[0]];
            LocalIndexType thisVertexPBIndex1 = gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM][thisCornerLocalIndices[1]];

            std::vector<int> corner0neighborset = gridstorage_.PB2PBNeighborRank_[VERTEX_CODIM][thisVertexPBIndex0];
            std::vector<int> corner1neighborset = gridstorage_.PB2PBNeighborRank_[VERTEX_CODIM][thisVertexPBIndex1];

            // Find neighbors common to both edge corners
            std::vector<int> edgeneighborset = Dune::VectorHelper::sortedSetIntersection(corner0neighborset, corner1neighborset);

            // Debug info
            std::stringstream log_stream;
            log_stream << "CurvilinearGridConstructor: --";
            log_stream << "localEdgeIndex=" << thisEdgeLocalIndex;
            log_stream << " EdgeKey=(" << thisEdgeKey.node0 << "," << thisEdgeKey.node1 << ")";
            //log_stream << "Neighbors[0]=(" << Dune::VectorHelper::vector2string(corner0neighborset) << ")";
            //log_stream << " Neighbors[1]=(" << Dune::VectorHelper::vector2string(corner1neighborset) << ")";
            log_stream << " Intersection=" << Dune::VectorHelper::vector2string(edgeneighborset);
            loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, log_stream.str());


            int nEdgeNeighbor = edgeneighborset.size();
            if (nEdgeNeighbor < 1) {
            	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_ERROR>( __FILE__, __LINE__, "CurvilinearGridConstructor: Found no neighbor processes to an edge ");
            	DUNE_THROW(Dune::IOError, "CurvilinearGrid: Found no neighbor processes to an edge ");
            }
            else if (nEdgeNeighbor > 1)
            {
            	// Add a complicated edge for further verification
            	// Store only after verification
                TmpEdgeData thisPBEdgeData(thisPBEdgeLocalIndex, thisEdgeKey);
                for (int iEdge = 0; iEdge < edgeneighborset.size(); iEdge++)  {
                	neighborProcessComplicatedEdgePBLocalIndex[edgeneighborset[iEdge]].push_back(thisPBEdgeData);
                }
            } else
            {
                // Store the edge neighbor rank set
                edgeneighborset.swap(gridstorage_.PB2PBNeighborRank_[EDGE_CODIM][thisPBEdgeLocalIndex]);
            }
        }


        // 2) Communicate to each process the number of complicated edges shared with it
        // *************************************************************************************************************
        std::vector<int> processNComplicatedEdgeRequested(size_);
        std::vector<int> processNComplicatedEdgeToSend(size_);
        for (int iProc = 0; iProc < size_; iProc++)  { processNComplicatedEdgeRequested[iProc] = neighborProcessComplicatedEdgePBLocalIndex[iProc].size(); }
        MPI_Alltoall(processNComplicatedEdgeRequested.data(), 1, MPI_INT, reinterpret_cast<int*>(processNComplicatedEdgeToSend.data()), 1, MPI_INT, comm);

        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Total complicated edges per process =(" + Dune::VectorHelper::vector2string(processNComplicatedEdgeRequested) + ")");


        // 3) Communicate to each process the shared complicated edge EdgeKeys
        // *************************************************************************************************************
        int thisCommSize = 0;
        std::vector<int> processEdgeKeyRequested, sdispls;
        std::vector<int> processEdgeKeyToSend, rdispls;

        for (int iProc = 0; iProc < size_; iProc++)
        {
        	for (int iEdge = 0; iEdge < processNComplicatedEdgeRequested[iProc]; iEdge++)
        	{
        		EdgeKey thisEdgeKey = neighborProcessComplicatedEdgePBLocalIndex[iProc][iEdge].second;
        		processEdgeKeyRequested.push_back(thisEdgeKey.node0);
        		processEdgeKeyRequested.push_back(thisEdgeKey.node1);
        	}

        	processNComplicatedEdgeRequested[iProc] *= 2;
        	processNComplicatedEdgeToSend[iProc]    *= 2;
        	thisCommSize += processNComplicatedEdgeToSend[iProc];

            sdispls.push_back((iProc == 0) ? 0 : sdispls[iProc-1] + processNComplicatedEdgeRequested[iProc-1] );
            rdispls.push_back((iProc == 0) ? 0 : rdispls[iProc-1] + processNComplicatedEdgeToSend[iProc-1] );
        }

        processEdgeKeyToSend.resize(thisCommSize);
        MPI_Alltoallv (processEdgeKeyRequested.data(), processNComplicatedEdgeRequested.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(processEdgeKeyToSend.data()), processNComplicatedEdgeToSend.data(), rdispls.data(), MPI_INT, comm );
        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Communicated complicated edge EdgeKeys");


        // 4) Communicate to each process whether requested edges exist on this process
        // *************************************************************************************************************

        // Note: now we communicate 1 int for every edge key requested, so send and recv switch places and are divided by 2

        thisCommSize = 0;
        std::vector<int> processEdgeExistToSend;      rdispls.clear();
        std::vector<int> processEdgeExistRequested;   sdispls.clear();

        int iData = 0;
        for (int iProc = 0; iProc < size_; iProc++)
        {
        	processNComplicatedEdgeRequested[iProc] /= 2;
        	processNComplicatedEdgeToSend[iProc] /= 2;
        	thisCommSize += processNComplicatedEdgeRequested[iProc];

        	for (int iEdge = 0; iEdge < processNComplicatedEdgeToSend[iProc]; iEdge++)
        	{
        		EdgeKey thisEdgeKey;
        		thisEdgeKey.node0 = processEdgeKeyToSend[iData++];
        		thisEdgeKey.node1 = processEdgeKeyToSend[iData++];

        		bool isReal = (edgeKey2LocalIndexMap_.find(thisEdgeKey) != edgeKey2LocalIndexMap_.end());
        		processEdgeExistToSend.push_back( isReal ? 1 : 0 );
        	}

            sdispls.push_back((iProc == 0) ? 0 : sdispls[iProc-1] + processNComplicatedEdgeToSend[iProc-1] );
            rdispls.push_back((iProc == 0) ? 0 : rdispls[iProc-1] + processNComplicatedEdgeRequested[iProc-1] );
        }

        processEdgeExistRequested.resize(thisCommSize, 0);
        MPI_Alltoallv (processEdgeExistToSend.data(), processNComplicatedEdgeToSend.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(processEdgeExistRequested.data()), processNComplicatedEdgeRequested.data(), rdispls.data(), MPI_INT, comm );
        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Communicated if requested EdgeKeys correspond to real edges");


        // 5) Fill in correct neighbors for complicated edges
        // *************************************************************************************************************
        iData = 0;
        for (int iProc = 0; iProc < size_; iProc++)
        {
        	for (int iEdge = 0; iEdge < processNComplicatedEdgeRequested[iProc]; iEdge++)
        	{
        		bool isReal = (processEdgeExistRequested[iData++] == 1);
        		LocalIndexType thisEdgePBLocalIndex = neighborProcessComplicatedEdgePBLocalIndex[iProc][iEdge].first;
        		if (isReal)  { gridstorage_.PB2PBNeighborRank_[EDGE_CODIM][thisEdgePBLocalIndex].push_back(iProc); }

        		std::stringstream log_stream;
        		log_stream << " complicated edge PBIndex=" << thisEdgePBLocalIndex << " marked as real=" << isReal << " by process " << iProc;
        		loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, log_stream.str());
        	}
        }
        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Finished computing edge process boundary neighbors");


        // 6) Sort all edge neighbor rank sets
        // *************************************************************************************************************
        for (int iEdge = 0; iEdge < gridstorage_.PB2PBNeighborRank_[EDGE_CODIM].size(); iEdge++)
        {
        	std::sort(gridstorage_.PB2PBNeighborRank_[EDGE_CODIM][iEdge].begin(), gridstorage_.PB2PBNeighborRank_[EDGE_CODIM][iEdge].end());
        }


    }


    /** Compute the process ranks of neighbor processes for all process boundary faces
     *
     *  Algorithm:
     *  1) Loop over all process boundary faces in the face map
     *  1.1) For each face corner in the FaceKey get associated neighbor ranks from provided vertex neighbor ranks
     *  1.2) Perform intersection on the three sets
     *  1.3) Ideally the intersection should result in one single rank, which is this face's neighbor. Otherwise throw error
     *
     * */
    void globalComputeFaceNeighborRank()
    {
    	MPI_Comm comm = Dune::MPIHelper::getCommunicator();
    	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Started computing face process boundary neighbors");


        // For each process stores the set of face indices local to this process, which are shared between more than two processes in total
    	typedef std::pair<LocalIndexType, FaceKey> TmpFaceData;
        std::vector<std::vector<TmpFaceData > > neighborProcessComplicatedFaceLocalIndex(size_);

        for (FaceMapIterator faceIter = processBoundaryFaceKey2LocalIndexMap_.begin(); faceIter != processBoundaryFaceKey2LocalIndexMap_.end(); faceIter++ )
        {
            // Get corners of the face
            FaceKey thisFaceKey = (*faceIter).first;
            LocalIndexType thisFaceLocalIndex = (*faceIter).second;

            // Get neighbor processes associated with each corner
            LocalIndexType thisVertexLocalIndex0 = gridstorage_.entityIndexMap_[VERTEX_CODIM][thisFaceKey.node0];
            LocalIndexType thisVertexLocalIndex1 = gridstorage_.entityIndexMap_[VERTEX_CODIM][thisFaceKey.node1];
            LocalIndexType thisVertexLocalIndex2 = gridstorage_.entityIndexMap_[VERTEX_CODIM][thisFaceKey.node2];

            LocalIndexType thisVertexPBIndex0 = gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM][thisVertexLocalIndex0];
            LocalIndexType thisVertexPBIndex1 = gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM][thisVertexLocalIndex1];
            LocalIndexType thisVertexPBIndex2 = gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM][thisVertexLocalIndex2];

            std::vector<int> corner0neighborset = gridstorage_.PB2PBNeighborRank_[VERTEX_CODIM][thisVertexPBIndex0];
            std::vector<int> corner1neighborset = gridstorage_.PB2PBNeighborRank_[VERTEX_CODIM][thisVertexPBIndex1];
            std::vector<int> corner2neighborset = gridstorage_.PB2PBNeighborRank_[VERTEX_CODIM][thisVertexPBIndex2];

            // Find neighbors common to all 3 face corners. Need to intersect sets twice
            std::vector<int> faceneighborset;
            faceneighborset = Dune::VectorHelper::sortedSetIntersection(corner0neighborset, corner1neighborset);
            faceneighborset = Dune::VectorHelper::sortedSetIntersection(faceneighborset,    corner2neighborset);

            std::stringstream log_stream;
            log_stream << "CurvilinearGridConstructor: --";
            log_stream << "localFaceIndex=" << thisFaceLocalIndex;
            log_stream << " FaceKey=(" << thisFaceKey.node0 << "," << thisFaceKey.node1 << "," << thisFaceKey.node2 << ")";
            //log_stream << "Neighbors[0]=(" << Dune::VectorHelper::vector2string(corner0neighborset) << ")";
            //log_stream << " Neighbors[1]=(" << Dune::VectorHelper::vector2string(corner1neighborset) << ")";
            //log_stream << " Neighbors[2]=(" << Dune::VectorHelper::vector2string(corner2neighborset) << ")";
            log_stream << " Intersection=(" << Dune::VectorHelper::vector2string(faceneighborset) << ")";
            loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, log_stream.str());

            int nFaceNeighbor = faceneighborset.size();

            if (nFaceNeighbor < 1) {
            	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_ERROR>( __FILE__, __LINE__, "CurvilinearGridConstructor: Found face with neighbor nProcess=" + std::to_string(nFaceNeighbor));
            	DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected number of neighbor processes to a face");
            }
            else if (nFaceNeighbor > 1)
            {
              	// Add a complicated face for further verification
              	// Store only after verification
            	TmpFaceData thisPBFaceData(thisFaceLocalIndex, thisFaceKey);
                for (int iFace = 0; iFace < nFaceNeighbor; iFace++)
                {
                  	neighborProcessComplicatedFaceLocalIndex[faceneighborset[iFace]].push_back(thisPBFaceData);
                }
            } else
            {
                // Store the face neighbor rank. Face is only allowed to have exactly one neighbor
            	LocalIndexType thisFaceLocalPBIndex = gridstorage_.processBoundaryIndexMap_[FACE_CODIM][thisFaceLocalIndex];
            	gridstorage_.PB2PBNeighborRank_[FACE_CODIM][thisFaceLocalPBIndex].push_back(faceneighborset[0]);
            }
        }


        // 2) Communicate to each process the number of complicated faces shared with it
        // *************************************************************************************************************
        std::vector<int> processNComplicatedFaceRequested(size_);
        std::vector<int> processNComplicatedFaceToSend(size_);
        for (int iProc = 0; iProc < size_; iProc++)  { processNComplicatedFaceRequested[iProc] = neighborProcessComplicatedFaceLocalIndex[iProc].size(); }
        MPI_Alltoall(processNComplicatedFaceRequested.data(), 1, MPI_INT, reinterpret_cast<int*>(processNComplicatedFaceToSend.data()), 1, MPI_INT, comm);
        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Complicated faces per process sent=( " + Dune::VectorHelper::vector2string(processNComplicatedFaceRequested) + ") received =(" + Dune::VectorHelper::vector2string(processNComplicatedFaceToSend) + ")");


        // 3) Communicate to each process the shared complicated face FaceKeys
        // *************************************************************************************************************
        int thisCommSize = 0;
        std::vector<int> processFaceKeyRequested, sdispls;
        std::vector<int> processFaceKeyToSend, rdispls;

        for (int iProc = 0; iProc < size_; iProc++)
        {
        	for (int iFace = 0; iFace < processNComplicatedFaceRequested[iProc]; iFace++)
        	{
        		FaceKey thisFaceKey = neighborProcessComplicatedFaceLocalIndex[iProc][iFace].second;
        		processFaceKeyRequested.push_back(thisFaceKey.node0);
        		processFaceKeyRequested.push_back(thisFaceKey.node1);
        		processFaceKeyRequested.push_back(thisFaceKey.node2);
        	}

        	processNComplicatedFaceRequested[iProc] *= 3;
        	processNComplicatedFaceToSend[iProc]    *= 3;
        	thisCommSize += processNComplicatedFaceToSend[iProc];

            sdispls.push_back((iProc == 0) ? 0 : sdispls[iProc-1] + processNComplicatedFaceRequested[iProc-1] );
            rdispls.push_back((iProc == 0) ? 0 : rdispls[iProc-1] + processNComplicatedFaceToSend[iProc-1] );
        }

        processFaceKeyToSend.resize(thisCommSize);
        MPI_Alltoallv (processFaceKeyRequested.data(), processNComplicatedFaceRequested.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(processFaceKeyToSend.data()), processNComplicatedFaceToSend.data(), rdispls.data(), MPI_INT, comm );
        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Communicated complicated face FaceKeys");

        //std::cout << "process_" << rank_ << "stage 3) sendcounts=" << Dune::VectorHelper::vector2string(processNComplicatedFaceRequested) << " recvcounts=" << Dune::VectorHelper::vector2string(processNComplicatedFaceToSend) <<" send=" << Dune::VectorHelper::vector2string(processFaceKeyRequested) << " recv=" << Dune::VectorHelper::vector2string(processFaceKeyToSend) << std::endl;


        // 4) Communicate to each process whether requested faces exist on this process
        // *************************************************************************************************************

        // Note: now we communicate 1 int for every face key requested, so send and recv switch places and are divided by 3

        std::vector<int> processFaceExistToSend;                                          rdispls.clear();
        std::vector<int> processFaceExistRequested(processFaceKeyRequested.size() / 3);   sdispls.clear();

        int iData = 0;
        for (int iProc = 0; iProc < size_; iProc++)
        {
        	processNComplicatedFaceRequested[iProc] /= 3;
        	processNComplicatedFaceToSend[iProc] /= 3;

        	for (int iFace = 0; iFace < processNComplicatedFaceToSend[iProc]; iFace++)
        	{
        		FaceKey thisFaceKey;
        		thisFaceKey.node0 = processFaceKeyToSend[iData++];
        		thisFaceKey.node1 = processFaceKeyToSend[iData++];
        		thisFaceKey.node2 = processFaceKeyToSend[iData++];

        		bool isReal = (processBoundaryFaceKey2LocalIndexMap_.find(thisFaceKey) != processBoundaryFaceKey2LocalIndexMap_.end());
        		processFaceExistToSend.push_back( isReal ? 1 : 0 );
        	}

            sdispls.push_back((iProc == 0) ? 0 : sdispls[iProc-1] + processNComplicatedFaceToSend[iProc-1] );
            rdispls.push_back((iProc == 0) ? 0 : rdispls[iProc-1] + processNComplicatedFaceRequested[iProc-1] );
        }
        MPI_Alltoallv (processFaceExistToSend.data(), processNComplicatedFaceToSend.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(processFaceExistRequested.data()), processNComplicatedFaceRequested.data(), rdispls.data(), MPI_INT, comm );
        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Communicated correspondence of requested FaceKeys correspond to real faces");

        std::cout << "process_" << rank_ << "stage 4) sendcounts=" << Dune::VectorHelper::vector2string(processNComplicatedFaceToSend) << " recvcounts=" << Dune::VectorHelper::vector2string(processNComplicatedFaceRequested) <<" send=" << Dune::VectorHelper::vector2string(processFaceExistToSend) << " recv=" << Dune::VectorHelper::vector2string(processFaceExistRequested) << std::endl;


        // 5) Fill in correct neighbors for complicated faces
        // *************************************************************************************************************
        iData = 0;
        for (int iProc = 0; iProc < size_; iProc++)
        {
        	for (int iFace = 0; iFace < processNComplicatedFaceRequested[iProc]; iFace++)
        	{
        		bool isReal = (processFaceExistRequested[iData++] == 1);
        		LocalIndexType thisFaceLocalIndex = neighborProcessComplicatedFaceLocalIndex[iProc][iFace].first;
        		FaceKey thisFaceKey = neighborProcessComplicatedFaceLocalIndex[iProc][iFace].second;

        		std::stringstream log_stream;
        		log_stream << " complicated face LocalIndex=" << thisFaceLocalIndex;
        		log_stream << " FaceKey=(" << thisFaceKey.node0 << "," << thisFaceKey.node1 << "," << thisFaceKey.node2 << ")";
        		log_stream << " marked as real=" << isReal << " by process " << iProc;
        		loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, log_stream.str());

        		if (isReal)
        		{
        			LocalIndexType thisFaceLocalPBIndex = gridstorage_.processBoundaryIndexMap_[FACE_CODIM][thisFaceLocalIndex];
        			int nNeighborAlready = gridstorage_.PB2PBNeighborRank_[FACE_CODIM][thisFaceLocalPBIndex].size();

        			// If the face neighbor has already been assigned, this face has more than 1 real neighbor process, which is impossible
        			if (nNeighborAlready != 0)
        			{
                    	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_ERROR>( __FILE__, __LINE__, "CurvilinearGridConstructor: Found face with neighbor more than two even after cross-check");
                    	DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected number of neighbor processes to a face");
        			}

        			gridstorage_.PB2PBNeighborRank_[FACE_CODIM][thisFaceLocalPBIndex].push_back(iProc);
        		}
        	}
        }


        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Finished computing face process boundary neighbors");
    }


    /** Communicates all process boundary face global Id's to the neighbors if owned
     *
     * Algorithm:
     *
     * 1) Loop over all process boundary faces, split faces on the ones to be sent and to be received
     * 1.1) Assemble a global index array to send to each process
     * 1.2) Note how many faces will be received from each process
     * 2) Assemble one big send array from small arrays (FaceKey + globalIndex)
     * 3) MPI_Alltoallv - communicate this array
     * 4) Save global indices for non-owned faces. Find the exact face by using the communicated FaceKey
     *
     * Optimization Proposal:
     * In principle communication of the FaceKey is not necessary. Instead, the natural FaceKey "<" operator
     * can be used to sort all communicated faces, thus allowing the receiving process to "figure out" what are
     * the faces sent to it by sorting its own faces.
     * 1) Sort all process boundary faces wrt FaceKey
     * 2) Fill arrays to send according to this sorted order
     * 3) Make map for each process from rank & received face to local face index accoridng to the sorted FaceKey order
     * 3) Communicate only globalIndices
     * 4) Use constructed map to fill in received global Indices
     *
     * Will decrease the global communication at the expense of increasing local computation time
     *
     * */
    void globalDistributeMissingFaceGlobalIndex()
    {
    	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Started distributing missing face GlobalIndices");

        typedef std::pair<FaceKey, GlobalIndexType>  FaceInfo;
        std::vector< std::vector< FaceInfo > > facesToSend (size_);

        int totalRecvSize = 0;
        int N_INTEGER_FACEINFO = 4;
        std::vector<int> sendbuf, sendcounts(size_), sdispls;
        std::vector<int> recvbuf, recvcounts(size_), rdispls;

        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: -- Checking which faces are missing");


        // 1) Loop over all process boundary faces, split faces on the ones to be sent and to be received
        // ********************************************************************************************
        for (FaceMapIterator iter = processBoundaryFaceKey2LocalIndexMap_.begin(); iter != processBoundaryFaceKey2LocalIndexMap_.end(); iter++)
        {
            LocalIndexType localFaceIndex = (*iter).second;
            LocalIndexType localFacePBIndex = gridstorage_.processBoundaryIndexMap_[FACE_CODIM][localFaceIndex];
            int neighborRank = gridstorage_.PB2PBNeighborRank_[FACE_CODIM][localFacePBIndex][0];

            // If the neighbor of this face has lower rank, then add it to recv, else to send
            if (neighborRank < rank_)  { recvcounts[neighborRank] += N_INTEGER_FACEINFO; }
            else
            {
            	GlobalIndexType thisGlobalIndex = gridstorage_.face_[localFaceIndex].globalIndex;
                facesToSend[neighborRank].push_back(FaceInfo((*iter).first, thisGlobalIndex ));
            }
        }

        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: -- Assembling arrays to send");


        // 2) Fill in communication arrays
        // ********************************************************************************************
        for (int i = 0; i < size_; i++)
        {
            sendcounts[i] = facesToSend[i].size() * N_INTEGER_FACEINFO;
            totalRecvSize += recvcounts[i];
            sdispls.push_back((i == 0) ? 0 : sdispls[i-1] + sendcounts[i-1] );
            rdispls.push_back((i == 0) ? 0 : rdispls[i-1] + recvcounts[i-1] );

            for (int j = 0; j < facesToSend[i].size(); j++)
            {
                sendbuf.push_back(facesToSend[i][j].second);
                sendbuf.push_back(facesToSend[i][j].first.node0);
                sendbuf.push_back(facesToSend[i][j].first.node1);
                sendbuf.push_back(facesToSend[i][j].first.node2);
            }

        }

        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: -- Sending  sendcounts=(" + Dune::VectorHelper::vector2string(sendcounts) + ") recvcounts=(" + Dune::VectorHelper::vector2string(recvcounts) + ")");

        // 3) MPI_alltoall key + globalId
        // ********************************************************************************************
        recvbuf.resize(totalRecvSize, 0);
        MPI_Comm comm = Dune::MPIHelper::getCommunicator();
        MPI_Alltoallv (sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(recvbuf.data()), recvcounts.data(), rdispls.data(), MPI_INT, comm );


        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: -- Extracting missing indices");

        // 4) Mark all missing faces
        // ********************************************************************************************8
        // Note that recvcounts should be 0 for this rank, as it was assembled by only considering the neighbor ranks
        int iData = 0;
        for (int iProc = 0; iProc < size_; iProc++)
        {
            int nThisFaceInfo = recvcounts[iProc] / N_INTEGER_FACEINFO;

            for (int iFace = 0; iFace < nThisFaceInfo; iFace++)
            {
                FaceKey thisKey;
                GlobalIndexType thisGlobalId = recvbuf[iData++];
                thisKey.node0 = recvbuf[iData++];
                thisKey.node1 = recvbuf[iData++];
                thisKey.node2 = recvbuf[iData++];

                FaceMapIterator faceIter = processBoundaryFaceKey2LocalIndexMap_.find(thisKey);

                if (faceIter == processBoundaryFaceKey2LocalIndexMap_.end()) {
                	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_ERROR>( __FILE__, __LINE__, "CurvilinearGridConstructor: Communicated FaceKey does not correspond to any face on this process");
                	DUNE_THROW(Dune::IOError, "CurvilinearGrid: Communicated FaceKey does not correspond to any face on this process ");
                }
                else
                {
                    LocalIndexType localFaceIndex = (*faceIter).second;
                    gridstorage_.face_[localFaceIndex].globalIndex = thisGlobalId;
                }
            }
        }

        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Finished distributing missing face GlobalIndices");
    }


    /** Communicates all process boundary face global indices to the neighbors if owned
     *
     * Algorithm:
     *
     * 1) Loop over all process boundary edges, split edges into ones to be sent and to be received
     * 1.1) If this edge rank lower than all other neighbor ranks, note to send it to all neighbors,
     *      Otherwise note which neighbor to receive it from
     * 1.1) Assemble a global index array to send to each process
     * 1.2) Note how many edges will be received from each process
     * 2) Assemble one big send array from small arrays (EdgeKey + globalIndex)
     * 3) MPI_Alltoallv - communicate this array
     * 4) Save global indices for non-owned edges. Find the exact edge by using the communicated EdgeKey
     *
     * Optimization Proposal:
     * Same as for FaceGlobalIndex algorithm
     *
     * */
    void globalDistributeMissingEdgeGlobalIndex()
    {
    	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Started distributing missing edge GlobalIndices");

        typedef std::pair<EdgeKey, GlobalIndexType>  EdgeInfo;
        std::vector< std::vector< EdgeInfo > > edgesToSend (size_);

        int totalRecvSize = 0;
        int N_INTEGER_EDGEINFO = 3;
        std::vector<int> sendbuf, sendcounts(size_), sdispls;
        std::vector<int> recvbuf, recvcounts(size_), rdispls;

        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: -- Checking which edges are missing");

        // 1) Loop over all process boundary faces, split faces on the ones to be sent and to be received
        // ********************************************************************************************
        for (Local2LocalIterator iter  = gridstorage_.processBoundaryIndexMap_[EDGE_CODIM].begin();
        		                 iter != gridstorage_.processBoundaryIndexMap_[EDGE_CODIM].end(); iter++)
        {
            LocalIndexType localEdgeIndex = (*iter).first;
            LocalIndexType localEdgePBIndex = (*iter).second;

            // Construct EdgeKey
            std::vector<LocalIndexType> thisCornerLocalIndices =  gridbase_.entityCornerLocalIndex(EDGE_CODIM, localEdgeIndex);
            EdgeKey thisEdgeKey;
            thisEdgeKey.node0 = gridstorage_.point_[thisCornerLocalIndices[0]].globalIndex;
            thisEdgeKey.node1 = gridstorage_.point_[thisCornerLocalIndices[1]].globalIndex;
            thisEdgeKey.sort();


            int candidateOwnerRank = gridstorage_.PB2PBNeighborRank_[EDGE_CODIM][localEdgePBIndex][0];

            //std::cout << "process_" << rank_ <<  " EdgeKey=(" << thisEdgeKey.node0 << "," << thisEdgeKey.node1 << ") localIndex=" << localEdgeIndex <<  std::endl;

            // If the one of the neighbors of this edge has lower rank, then note one more received edge from that process
            // else note to send it to all other neighbors
            if (candidateOwnerRank < rank_)  { recvcounts[candidateOwnerRank] += N_INTEGER_EDGEINFO; }
            else
            {
            	GlobalIndexType thisGlobalIndex = gridstorage_.edge_[localEdgeIndex].globalIndex;

                EdgeInfo thisEdgeInfo(thisEdgeKey, thisGlobalIndex);

                for (int iNeighbor = 0; iNeighbor < gridstorage_.PB2PBNeighborRank_[EDGE_CODIM][localEdgePBIndex].size(); iNeighbor++)
                {
                    int thisNeighborRank = gridstorage_.PB2PBNeighborRank_[EDGE_CODIM][localEdgePBIndex][iNeighbor];
                    edgesToSend[thisNeighborRank].push_back(thisEdgeInfo);
                };
            }
        }

        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: -- Assembling arrays to send");


        // 2) Fill in communication arrays
        // ********************************************************************************************
        for (int i = 0; i < size_; i++)
        {
            sendcounts[i] = edgesToSend[i].size() * N_INTEGER_EDGEINFO;
            totalRecvSize += recvcounts[i];
            sdispls.push_back((i == 0) ? 0 : sdispls[i-1] + sendcounts[i-1] );
            rdispls.push_back((i == 0) ? 0 : rdispls[i-1] + recvcounts[i-1] );

            for (int j = 0; j < edgesToSend[i].size(); j++)
            {
                sendbuf.push_back(edgesToSend[i][j].second);
                sendbuf.push_back(edgesToSend[i][j].first.node0);
                sendbuf.push_back(edgesToSend[i][j].first.node1);
            }

        }

        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: -- Sending  sendcounts=(" + Dune::VectorHelper::vector2string(sendcounts) + ") recvcounts=(" + Dune::VectorHelper::vector2string(recvcounts) + ")");



        // 3) MPI_alltoall key + globalId
        // ********************************************************************************************
        recvbuf.resize(totalRecvSize, 0);   // There are 3 integers per sent FaceInfo
        MPI_Comm comm = Dune::MPIHelper::getCommunicator();
        MPI_Alltoallv (sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(recvbuf.data()), recvcounts.data(), rdispls.data(), MPI_INT, comm );


        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: -- Extracting missing indices");

        // 4) Mark all missing faces
        // ********************************************************************************************8
        // Note that recvcounts should be 0 for this rank, as it was assembled by only considering the neighbor ranks
        int iData = 0;
        for (int iProc = 0; iProc < size_; iProc++)
        {
            int nThisEdgeInfo = recvcounts[iProc] / N_INTEGER_EDGEINFO;

            for (int iEdge = 0; iEdge < nThisEdgeInfo; iEdge++)
            {
                EdgeKey thisKey;
                GlobalIndexType thisGlobalId = recvbuf[iData++];
                thisKey.node0 = recvbuf[iData++];
                thisKey.node1 = recvbuf[iData++];

                EdgeMapIterator edgeIter = edgeKey2LocalIndexMap_.find(thisKey);

                if (edgeIter == edgeKey2LocalIndexMap_.end()) {
                	std::stringstream log_str;
                	log_str << "CurvilinearGridConstructor: Communicated EdgeKey (" << thisKey.node0 << ", " << thisKey.node1 << ") does not correspond to any edge on this process";
                	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_ERROR>( __FILE__, __LINE__, log_str.str());
                	DUNE_THROW(Dune::IOError, "CurvilinearGrid: Communicated EdgeKey does not correspond to any edge on this process "); }
                else
                {
                    LocalIndexType localEdgeIndex = (*edgeIter).second;
                    gridstorage_.edge_[localEdgeIndex].globalIndex = thisGlobalId;
                }
            }
        }

        loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearGridConstructor: Finished distributing missing edge GlobalIndices");
    }



private: // Private members

    LoggingMessage & loggingmessage_;

    // Temporary maps necessary to locate and communicate entities during grid base construction
    EdgeKey2EdgeIndexMap edgeKey2LocalIndexMap_;                    // (global edgeKey -> edge_ index)
    FaceKey2FaceIndexMap internalFaceKey2LocalIndexMap_;            // (global faceKey -> gridstorage_.face_ index)
    FaceKey2FaceIndexMap domainBoundaryFaceKey2LocalIndexMap_;      // (global faceKey -> gridstorage_.face_ index)
    FaceKey2FaceIndexMap processBoundaryFaceKey2LocalIndexMap_;     // (global faceKey -> gridstorage_.face_ index)

    // Curvilinear Grid Storage Class
    GridStorageType & gridstorage_;

    // Reference to Curvilinear Grid Base - necessary for OCTree construction
    GridBaseType & gridbase_;

    MPIHelper &mpihelper_;
    int rank_;
    int size_;
};

} // namespace Dune

#endif  // DUNE_CURVILINEARGRIDCONSTRUCTOR_HH
