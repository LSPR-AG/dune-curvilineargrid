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

#include <dune/curvilineargrid/curvilineargridbase/curvilineargridstorage.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilineargridbase.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilinearoctreenode.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilinearlooseoctree.hh>


namespace Dune {




// Forwards-declatation of the base class
// **********************************************
template <class ct>
class CurvilinearGridBase;


// Constructor class
// **********************************************
template <class ct>
class CurvilinearGridConstructor {
public:

    /* public types */
    typedef Dune::CurvilinearGridStorage<ct>        GridStorageType;
    typedef Dune::CurvilinearGridBase<ct>           GridBaseType;

    typedef typename GridStorageType::Vertex                 Vertex;
    typedef typename GridStorageType::VertexStorage          VertexStorage;
    typedef typename GridStorageType::EdgeStorage            EdgeStorage;
    typedef typename GridStorageType::FaceStorage            FaceStorage;
    typedef typename GridStorageType::EntityStorage          EntityStorage;

    typedef typename GridStorageType::EdgeKey                EdgeKey;
    typedef typename GridStorageType::FaceKey                FaceKey;

    typedef typename GridStorageType::EdgeKey2EdgeIdMap      EdgeKey2EdgeIdMap;
    typedef typename GridStorageType::FaceKey2FaceIdMap      FaceKey2FaceIdMap;
    typedef typename GridStorageType::Index2IndexMap         Index2IndexMap;

    typedef typename EdgeKey2EdgeIdMap::iterator             EdgeMapIterator;
    typedef typename FaceKey2FaceIdMap::iterator             FaceMapIterator;
    typedef typename Index2IndexMap::iterator                IndexMapIterator;


	typedef Dune::CurvilinearOctreeNode<ct>                       NodeType;
    typedef Dune::CurvilinearLooseOctree<ct, 3, NodeType>         CurvilinearLooseOctree;


    // Logging Message Typedefs
    static const unsigned int LOG_PHASE_DEV = Dune::LoggingMessage::Phase::DEVELOPMENT_PHASE;
    static const unsigned int LOG_CATEGORY_DEBUG = Dune::LoggingMessage::Category::DEBUG;



public: /* public methods */

    /** Parallel constructor - USE THIS CONSTRUCTOR*/
    CurvilinearGridConstructor(
    		bool withGhostElements,
    		bool verbose,
    		bool processVerbose,
    		GridStorageType & gridstorage,
    		GridBaseType & gridbase,
    		MPIHelper &mpihelper ) :
        withGhostElements_(withGhostElements),
        verbose_(verbose),
        processVerbose_(processVerbose),
        gridstorage_(gridstorage),
        gridbase_(gridbase),
        mpihelper_(mpihelper)
    {
        rank_ = mpihelper_.rank();
        size_ = mpihelper_.size();

        std::string log_string = "Initialized CurvilinearGridConstructor withGhostElements=" + std::to_string(withGhostElements);
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);
    }


public:

    /* ***************************************************************************
     * Section: Loading the mesh
     * ***************************************************************************/


    /** \brief Add a new vertex to the mesh
     * \param[in] globalIndex      global index of this vertex
     * \param[in] p                coordinate of this vertex
     * */
    void insertVertex(Vertex p, int globalIndex)
    {
        VertexStorage point;
        point.coord = p;
        point.globalIndex = globalIndex;

        gridstorage_.vertexGlobal2LocalMap_[globalIndex] = gridstorage_.point_.size();
        gridstorage_.point_.push_back(point);

        std::stringstream log_stream;
        log_stream << "CurvilinearGridConstructor: Inserted vertex LocalIndex=" << gridstorage_.point_.size()-1 << " GlobalIndex=" << globalIndex;
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_stream.str());
    }

    /** \brief Insert an element into the mesh
     * \param[in] gt               geometry type of the element (hopefully tetrahedron)
     * \param[in] globalId         the global index of this element
     * \param[in] vertexIndexSet   local indices of interpolatory vertices. Local with respect to the order the vertices were inserted
     * \param[in] order            interpolatory order of the element
     * \param[in] physicalTag      physical tag of the element
     *
     * Note: Even though we pass the globalId as a parameter from GMSH, it is a globalIndex for the set Elements+BoundaryFaces,
     * therefore obtaining globalIndex for elements from it is not possible, it will have to be communicated
     *
     * */
    void insertElement(Dune::GeometryType gt, int globalId, const std::vector<int> & vertexIndexSet, int order, int physicalTag)
    {
        if (!gt.isTetrahedron() || (vertexIndexSet.size() != Dune::CurvilinearGeometryHelper::dofPerOrder(gt, order)))  {
            DUNE_THROW(Dune::IOError, "CurvilinearGrid: insertElement() unexpected element type or number of interpolatory points");
        }

        EntityStorage thisElement;

        thisElement.geometryType = gt;
        thisElement.globalIndex = 0;        // At this stage globalIndex is not known yet
        thisElement.interpOrder = order;
        thisElement.physicalTag = physicalTag;
        thisElement.vertexIndexSet = vertexIndexSet;

        int thisLocalIndex = gridstorage_.element_.size();
        gridstorage_.elementGlobal2LocalMap_[globalId] = thisLocalIndex;
        gridstorage_.element_.push_back(thisElement);

        std::stringstream log_stream;
        log_stream << "CurvilinearGridConstructor: Inserted Element Type=" << Dune::CurvilinearGeometryHelper::geometryName(gt);
        log_stream << " LocalIndex=" << thisLocalIndex;
        log_stream << " Order=" << order;
        log_stream << " PhysicalTag=" << physicalTag;
        log_stream << " VertexIndices=(" << vector2string(vertexIndexSet) << ")";
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_stream.str());
    }

    /** Insert a boundary segment into the mesh
     *
     *     Note: It is expected that all faces - domain an process boundaries - are inserted by the factory before finalising
     *     Note: Only domain boundary faces have initial globalId given by GMSH. Therefore, we ignore it, and generate our own
     *     globalId for all faces at a later stage.
     *
     *  \param[in] gt                       geometry type of the face (should be a triangle)
     *  \param[in] globalId                 id of the element this face is associated to
     *  \param[in] associatedElementIndex   local index of the element this face is associated to
     *  \param[in] vertexIndexSet           local indices of the interpolatory vertices of this face
     *  \param[in] order                    interpolatory order of the face
     *  \param[in] physicalTag              physical tag of the element (material property)
     *
     * */

    void insertBoundarySegment(Dune::GeometryType gt, int globalId, int associatedElementIndex, const std::vector<int> & vertexIndexSet, int order, int physicalTag)
    {
        if (!gt.isTriangle() || (vertexIndexSet.size() != Dune::CurvilinearGeometryHelper::dofPerOrder(gt, order)))  {
            DUNE_THROW(Dune::IOError, "CurvilinearGridConstructor: insertBoundarySegment() unexpected number of interpolatory points");
        }


        // Get corners of this face
        // **********************************************************************************
        std::vector<int> faceCorners = entityVertexCornerSubset(gt, vertexIndexSet, order);
        FaceKey thisFaceKey;
        thisFaceKey.node0 = gridstorage_.point_[faceCorners[0]].globalIndex;
        thisFaceKey.node1 = gridstorage_.point_[faceCorners[1]].globalIndex;
        thisFaceKey.node2 = gridstorage_.point_[faceCorners[2]].globalIndex;

        // Sort in ascending order
        thisFaceKey.sort();


        // Take associated element, get all its corners, get all keys, compare to face key
        // **********************************************************************************
        std::vector<int> elementCorners = entityVertexCornerSubset(
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
                DUNE_THROW(Dune::IOError, "CurvilinearGrid: insertBoundarySegment() did not find the face in the associated element");
            }

            std::vector<int> internalLinearSubentityIndices =
            		Dune::CurvilinearGeometryHelper::linearElementSubentityCornerInternalIndexSet(gridstorage_.element_[associatedElementIndex].geometryType, 1, j);

            // Define (key = sorted localIndices of corners)
            FaceKey thisKey;
            thisKey.node0 = gridstorage_.point_[elementCorners[internalLinearSubentityIndices[0]]].globalIndex;
            thisKey.node1 = gridstorage_.point_[elementCorners[internalLinearSubentityIndices[1]]].globalIndex;
            thisKey.node2 = gridstorage_.point_[elementCorners[internalLinearSubentityIndices[2]]].globalIndex;


            // Sort in ascending order
            thisKey.sort();

            // By comparison find internalIndex of this face
            if (thisKey == thisFaceKey)
            {
                found_face = true;

                // Store domain and process boundaries separately for faster iterators
                // Store Map (key -> domainBoundaryIndex), Vector (domainBoundaryIndex -> faceIndex)
                int localDomainBoundaryIndex = gridstorage_.domainBoundaryFaceIndex_.size();
                int localFaceIndex = gridstorage_.face_.size();
                boundaryInternalMap_[thisKey] = localDomainBoundaryIndex;
                gridstorage_.domainBoundaryFaceIndex_.push_back(localFaceIndex);

                // Store Vector (faceId -> associated element)
                FaceStorage thisFaceAsSubentity;
                thisFaceAsSubentity.globalIndex  = 0;                  // At this stage the globalId is not known yet
                thisFaceAsSubentity.structuralType = Dune::CurvilinearGridFaceType::DomainBoundary;
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
                Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_stream.str());
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

    void generateMesh(int nVertexTotalMesh, int nElementTotalMesh) {

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Initializing mesh");


        // Temporary maps for Global Index construction
        // ************************************************************
        Index2IndexMap processBoundaryCornerMap;                            // (vertex global index -> processBoundaryCornerNeighborRank_ index)
        EdgeKey2EdgeIdMap processBoundaryEdgeMap;                           // (EdgeKey             -> processBoundaryEdgeNeighborRank_ index)


        // Construct missing parts of the mesh
        // ************************************************************
        gridstorage_.nVertexTotal_ = nVertexTotalMesh;
        gridstorage_.nElementTotal_ = nElementTotalMesh;

        computeProcessBoundingBox();
        generateEdges();
        generateFaces();

        if (size_ > 1)
        {
#if HAVE_MPI
        	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Assembling Parallel Grid");

        	// Parallel case
            generateProcessBoundaryCorners(processBoundaryCornerMap);
            generateProcessBoundaryEdges(processBoundaryEdgeMap);

            generateGlobalIndices(processBoundaryCornerMap, processBoundaryEdgeMap);
            if (withGhostElements_) { generateGhostElements(); }
#endif
        }
        else
        {
        	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Assembling Serial Grid");
            // Serial case:
            // * Boundary Neighbors not necessary, since all boundaries are domain boundaries
            // * No ghost elements, even if requested by user
            // * Fake globalIndex by making it equal to localIndex

            withGhostElements_ = false;
            for (int i = 0; i < gridstorage_.edge_.size();    i++)  { gridstorage_.edge_[i].globalIndex = i;     gridstorage_.edgeGlobal2LocalMap_[i] = i; }
            for (int i = 0; i < gridstorage_.face_.size();    i++)  { gridstorage_.face_[i].globalIndex = i;     gridstorage_.faceGlobal2LocalMap_[i] = i; }
            for (int i = 0; i < gridstorage_.element_.size(); i++)  { gridstorage_.element_[i].globalIndex = i;  gridstorage_.elementGlobal2LocalMap_[i] = i; }

            gridstorage_.nEdgeTotal_ = gridstorage_.edge_.size();
            gridstorage_.nFaceTotal_ = gridstorage_.face_.size();
        }

        constructOctree();


        // Deletes all temporary memory
        // ************************************************************
        edgemap_.clear();
        internalInternalMap_.clear();
        boundaryInternalMap_.clear();
        processInternalMap_.clear();
    }




protected:


    /* ***************************************************************************
     * Section: Auxiliary Methods
     * ***************************************************************************/

    // Converts an arbitrary vector into string by sticking all of the entries together
    // Whitespace separated
    template <class T>
    std::string vector2string(const T & V)
    {
        std::stringstream tmp_stream;

        int nEntry = V.size();
        if (nEntry == 0)  { tmp_stream << "Null"; }
        for (int i = 0; i < nEntry; i++) {
        	tmp_stream << V[i];
        	if (i != nEntry - 1) { tmp_stream << " "; }
        }
        return tmp_stream.str();
    }

    // Returns corner id's of this entity
    // TODO: Move functionality to the curvilinear geometry
    std::vector<int> entityVertexCornerSubset(
            Dune::GeometryType gt,
            const std::vector<int> & vertexIndexSet,
            int order) const
    {
        std::vector<int> corner;

        // Get corner number
        int cornerNo = gt.dim() + 1;
        //cornerNumber = ReferenceElements::general(geomType).size( geomType.dim() );

        // Get corners
        for (int j = 0; j < cornerNo; j++) {
            int internalId = Dune::CurvilinearGeometryHelper::cornerID(gt, order, j );
            corner.push_back(vertexIndexSet[internalId]);
        }

        return corner;
    }

    // Takes two sorted arrays with non-repeating entries
    // Returns an array which only has entries found in both input arrays
    std::vector<int> sortedSetIntersection(std::vector<int> A, std::vector<int> B)
    {
        std::vector<int>  rez;

        int indA = 0;
        int indB = 0;

        while ((indA < A.size()) && (indB < B.size()))
        {
                  if (A[indA] < B[indB]) { indA++; }
             else if (A[indA] > B[indB]) { indB++; }
             else {
                 rez.push_back(A[indA]);
                 indA++;
                 indB++;
             }
        }

        return rez;
    }


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

    // Generates all edges
    // FIXME: Original code seems to give edges some orientation
    // FIXME: Currently subentity orientation does not match the one of Dune
    // TODO: Use more generic functions when extending to any mesh other than tetrahedral
    void generateEdges()
    {
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Started constructing edges");

    	// Init the subentity index vector
    	int nEdgePerTetrahedron = 6;
    	int nElem = gridstorage_.element_.size();
    	gridstorage_.elementSubentityCodim1_ = std::vector<std::vector<int> > (nElem, std::vector<int>(nEdgePerTetrahedron));

        // Loop over all elements and their edges
        for (int iElem = 0; iElem < nElem; iElem++)
        {
        	EntityStorage & thisElem = gridstorage_.element_[iElem];
            std::vector<int> elementCornerLocalIndexSet = entityVertexCornerSubset(thisElem.geometryType, thisElem.vertexIndexSet, thisElem.interpOrder);

            for (int iEdge = 0; iEdge < nEdgePerTetrahedron; iEdge++)
            {
                std::vector<int> internalLinearSubentityIndices = Dune::CurvilinearGeometryHelper::linearElementSubentityCornerInternalIndexSet(thisElem.geometryType, 2, iEdge);

                // Define (key = sorted globalIndices of corners)
                EdgeKey thisKey;
                thisKey.node0 = gridstorage_.point_[elementCornerLocalIndexSet[internalLinearSubentityIndices[0]]].globalIndex;
                thisKey.node1 = gridstorage_.point_[elementCornerLocalIndexSet[internalLinearSubentityIndices[1]]].globalIndex;

                // Sort in ascending order
                thisKey.sort();
                //std::cout << "process_" << rank_ << "GenerateEdgeKey=(" << thisKey.node0 << "," << thisKey.node1 << ")" << std::endl;


                EdgeMapIterator edgeIter = edgemap_.find(thisKey);

                // If this edge has not been added already, add it to the map
                // Find local index of this edge and note it as this element subentity
                if (edgeIter == edgemap_.end())
                {
                    // Store map (key -> edgeIndex)
                	int localEdgeIndex = gridstorage_.edge_.size();
                    edgemap_[thisKey] = localEdgeIndex;

                    // Store vector (edgeId -> elemId + edgeElemIndex)
                    // Note: Edges do not have physical tag at all so we do not even store it
                    EdgeStorage thisEdge;
                    thisEdge.globalIndex = 0;        // GlobalId for edge determined later using global communication
                    thisEdge.elementIndex = iElem;
                    thisEdge.subentityIndex = iEdge;

                    // Log output
                    std::stringstream log_stream;
                    log_stream << "CurvilinearGridConstructor: Added Edge";
                    log_stream << " LocalIndex=" << localEdgeIndex;
                    log_stream << " AssociatedElementIndex=" << iElem;
                    log_stream << " InternalSubentityIndex=" << iEdge;
                    Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_stream.str());

                    gridstorage_.edge_.push_back(thisEdge);

                    gridstorage_.elementSubentityCodim1_[iElem][iEdge] = localEdgeIndex;
                } else {
                	gridstorage_.elementSubentityCodim1_[iElem][iEdge] = (*edgeIter).second;
                }
            }
        }
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Finished constructing edges");
    }

    // Generates Internal and ProcessBoundary Faces. (!!!) Assumes that all Domain Boundary Faces have been added.
    // FIXME: Currently subentity orientation does not match the one of Dune
    // TODO: If assume mesh with non-uniform p-refinement, it may make sense to point the face to the element which has the higher refinement
    void generateFaces()
    {
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Started generating faces");

    	// Init the subentity index vector
        int nFacePerTetrahedron = 4;
    	int nElem = gridstorage_.element_.size();
    	gridstorage_.elementSubentityCodim2_ = std::vector<std::vector<int> > (nElem, std::vector<int>(nFacePerTetrahedron));

        typedef std::map<FaceKey, std::vector<int>> tmpFace2InfoMap;
        typedef typename tmpFace2InfoMap::iterator  tmpMapIterator;
        tmpFace2InfoMap tmpFaceMap;

        // Loop over all elements and their faces
        for (int iElem = 0; iElem < nElem; iElem++)
        {
        	EntityStorage & thisElem = gridstorage_.element_[iElem];
            std::vector<int> elementCornerLocalIndexSet = entityVertexCornerSubset(thisElem.geometryType, thisElem.vertexIndexSet, thisElem.interpOrder);

            // Store info for all faces except of domain boundaries
            // Store it in a map, not to store internal faces twice
            for (int iFace = 0; iFace < nFacePerTetrahedron; iFace++)
            {
                std::vector<int> internalLinearSubentityIndices = Dune::CurvilinearGeometryHelper::linearElementSubentityCornerInternalIndexSet(thisElem.geometryType, 1, iFace);

                // Define (key = sorted globalIndices of corners)
                FaceKey thisKey;
                thisKey.node0 = gridstorage_.point_[elementCornerLocalIndexSet[internalLinearSubentityIndices[0]]].globalIndex;
                thisKey.node1 = gridstorage_.point_[elementCornerLocalIndexSet[internalLinearSubentityIndices[1]]].globalIndex;
                thisKey.node2 = gridstorage_.point_[elementCornerLocalIndexSet[internalLinearSubentityIndices[2]]].globalIndex;

                // Sort in ascending order
                thisKey.sort();


                FaceMapIterator faceIter = boundaryInternalMap_.find(thisKey);

                // Mark this face for creation if it is not an already existing Domain Boundary
                // Otherwise note its local index
                if (faceIter == boundaryInternalMap_.end())
                {
                    std::vector<int> connectedFaceInfo;
                    tmpMapIterator iter = tmpFaceMap.find(thisKey);

                    if (iter != tmpFaceMap.end()) { connectedFaceInfo = std::vector<int> ( (*iter).second ); }

                    connectedFaceInfo.push_back(iElem);
                    connectedFaceInfo.push_back(iFace);

                    tmpFaceMap[thisKey] = connectedFaceInfo;


                    std::stringstream log_stream;
                    log_stream << "CurvilinearGridConstructor: Adding FaceKey=(" << thisKey.node0 << ", " << thisKey.node1 << ", " << thisKey.node2 << ") attached to total of " << connectedFaceInfo.size() / 2 << " elements";
                    Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_stream.str());
                } else {
                	int localFaceIndex = gridstorage_.domainBoundaryFaceIndex_[(*faceIter).second];
                	gridstorage_.elementSubentityCodim2_[iElem][iFace] = localFaceIndex;
                }


            }
        }

        // Add internal and process boundary faces to the mesh
        for (tmpMapIterator iter = tmpFaceMap.begin();  iter != tmpFaceMap.end();  iter++)
        {
            FaceStorage thisFace;
            int localFaceIndex = gridstorage_.face_.size();
            std::vector<int> connectedFaceInfo = (*iter).second;

            // Store the face local index as element subentity
            gridstorage_.elementSubentityCodim2_[connectedFaceInfo[0]][connectedFaceInfo[1]] = localFaceIndex;

            // Add face to the mesh
            thisFace.globalIndex = 0;       // GlobalId is defined at a later stage
            thisFace.element1Index = connectedFaceInfo[0];
            thisFace.element1SubentityIndex = connectedFaceInfo[1];
            thisFace.physicalTag = -1;    // At the moment physicalTag of an internal face is not defined as it could be inbetween two different elements


            std::stringstream log_stream;
            log_stream << "CurvilinearGridConstructor: Added Face";
            log_stream << " LocalIndex=" << localFaceIndex;
            log_stream << " AssociatedElementIndex=" << thisFace.element1Index;
            log_stream << " InternalSubentityIndex=" << thisFace.element1SubentityIndex;

            // Store internal, domain and process boundaries separately for faster iterators
            if (connectedFaceInfo.size() == 2)
            {
                thisFace.structuralType = Dune::CurvilinearGridFaceType::ProcessBoundary;

                // Store Map (key -> processBoundaryFaceIndex), Vector (processBoundaryFaceIndex -> faceIndex)
                int localProcessBoundaryIndex = gridstorage_.processBoundaryFaceIndex_.size();
                processInternalMap_[(*iter).first] = localProcessBoundaryIndex;
                gridstorage_.processBoundaryFaceIndex_.push_back(localFaceIndex);

                thisFace.element2Index = 0;    // Eventually this will be the Ghost Element Index

                log_stream << " StructuralType=processBoundary";
            }
            else
            {
                thisFace.structuralType = Dune::CurvilinearGridFaceType::Internal;

                // Store Map (key -> internalFaceIndex), Vector (internalFaceIndex -> faceIndex)
                int localInternalFaceIndex = gridstorage_.internalFaceIndex_.size();
                internalInternalMap_[(*iter).first] = localInternalFaceIndex;
                gridstorage_.internalFaceIndex_.push_back(localFaceIndex);

                thisFace.element2Index = connectedFaceInfo[2];  // This is the 2nd neighbor of this internal face

                log_stream << " StructuralType=internal";
            }


            Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_stream.str());
            gridstorage_.face_.push_back(thisFace);
        }

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Finished generating faces");
    }


    // TODO: To make work with arbitrary geometries, must replace numbers, as well as make FaceKeys more abstract
    void generateProcessBoundaryCorners(Index2IndexMap & processBoundaryCornerMap)
    {
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Started generating BoundaryCorneers");

        // Construct the set of EdgeKeys corresponding to edges of processBoundaries
        // ********************************************************
        for (FaceMapIterator faceIter = processInternalMap_.begin(); faceIter != processInternalMap_.end(); faceIter++)
        {
            // Get global indices of the associated vertices from the map
            FaceKey thisFaceKey = (*faceIter).first;

            int thisVertexKey[3] = {thisFaceKey.node0, thisFaceKey.node1, thisFaceKey.node2};

            for (int i = 0; i < 3; i++)
            {
                // For each vertex, if it has not yet been added to the map, add it, mapping to a new entry
                // in the processBoundaryCornerNeighborRank_
                if (processBoundaryCornerMap.find(thisVertexKey[i]) == processBoundaryCornerMap.end())
                {
                	// DO NOT combine these two lines into one - then the map returns the size after insertion which is +1 to expected
                    int processBoundaryCornerIndex = processBoundaryCornerMap.size();
                	processBoundaryCornerMap[thisVertexKey[i]] = processBoundaryCornerIndex;

                    Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: -- Adding boundary corner GlobalIndex=" + std::to_string(thisVertexKey[i]));
                }
            }
        }

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Finished generating BoundaryCorneers");
    }


    void generateProcessBoundaryEdges(EdgeKey2EdgeIdMap & processBoundaryEdgeMap)
    {
        // Construct the set of process boundary corners - corners necessary to make process boundary faces on this process
        // ********************************************************
        for (FaceMapIterator faceIter = processInternalMap_.begin(); faceIter != processInternalMap_.end(); faceIter++)
        {
            // Get global indices of the associated vertices from the map
            FaceKey thisFaceKey = (*faceIter).first;
            EdgeKey thisEdgeKey[3];

            thisEdgeKey[0].node0 = thisFaceKey.node0;  thisEdgeKey[0].node1 = thisFaceKey.node1;
            thisEdgeKey[1].node0 = thisFaceKey.node0;  thisEdgeKey[1].node1 = thisFaceKey.node2;
            thisEdgeKey[2].node0 = thisFaceKey.node1;  thisEdgeKey[2].node1 = thisFaceKey.node2;

            for (int i = 0; i < 3; i++)
            {
                // For each vertex, if it has not yet been added to the map, add it, mapping to a new entry
                // in the processBoundaryCornerNeighborRank_
                if (processBoundaryEdgeMap.find(thisEdgeKey[i]) == processBoundaryEdgeMap.end())
                {
                	// DO NOT combine these two lines into one - then the map returns the size after insertion which is +1 to expected
                    int processBoundaryEdgeIndex = processBoundaryEdgeMap.size();
                    processBoundaryEdgeMap[thisEdgeKey[i]]  = processBoundaryEdgeIndex;

                    std::stringstream log_stream;
                    log_stream << "CurvilinearGridConstructor: -- Adding boundary EdgeKey= (" << thisEdgeKey[i].node0 << ", " << thisEdgeKey[i].node1 << ")";
                    Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_stream.str());
                }
            }
        }
    }



    // TODO: Inefficient Algorithm
    // The convention of owning an entity based on rank priority implies that processes with lower rank
    // have to do most of the work, and then perform a lot of communication. To balance out the workload
    // one would derive a more balanced owning paradigm
    //
    // Propose balanced paradigm:
    // Sum(Key[i]) mod nNeighbors
    void generateGlobalIndices(
        Index2IndexMap & processBoundaryCornerMap,
        EdgeKey2EdgeIdMap & processBoundaryEdgeMap
    )
    {
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Constructing Global Indices");
        Dune::CollectiveCommunication<MPI_Comm> collective_comm = mpihelper_.getCollectiveCommunication();


        // 1) Communicate process ranks associated with each process boundary corner
        // Then compute that for edges and faces using set intersection
        // *************************************************************************
        std::vector<std::vector<int> > processBoundaryCornerNeighborRank( processBoundaryCornerMap.size(), std::vector<int>() );   // List of ranks of all other processes sharing this vertex
        std::vector<std::vector<int> > processBoundaryEdgeNeighborRank( processBoundaryEdgeMap.size(), std::vector<int>() );       // List of ranks of all other processes sharing this edge
        gridstorage_.processBoundaryNeighborProcess_.resize(processInternalMap_.size());

        globalCommunicateVertexNeighborRanks(processBoundaryCornerMap, processBoundaryCornerNeighborRank);
        globalComputeEdgeNeighborRanks(processBoundaryCornerMap, processBoundaryEdgeMap, processBoundaryCornerNeighborRank, processBoundaryEdgeNeighborRank);
        globalComputeFaceNeighborRank(processBoundaryCornerMap, processBoundaryCornerNeighborRank);

        processBoundaryCornerMap.clear();
        processBoundaryCornerNeighborRank.clear();


        // 2) Get edges and faces on this process that are not owned by this process
        // *************************************************************************
        EdgeKey2EdgeIdMap edgeNonOwned;
        FaceKey2FaceIdMap faceNonOwned;

        // Edges
        for (EdgeMapIterator edgeIter = processBoundaryEdgeMap.begin(); edgeIter != processBoundaryEdgeMap.end(); edgeIter++ )
        {
            int edgeOwnerCandidateRank = processBoundaryEdgeNeighborRank[(*edgeIter).second][0];
            if (edgeOwnerCandidateRank < rank_) { edgeNonOwned[(*edgeIter).first] = edgeOwnerCandidateRank; }
        }

        // Faces
        for (FaceMapIterator faceIter = processInternalMap_.begin(); faceIter != processInternalMap_.end(); faceIter++ )
        {
            int faceOwnerCandidateRank = gridstorage_.processBoundaryNeighborProcess_[(*faceIter).second];
            if (faceOwnerCandidateRank < rank_) { faceNonOwned[(*faceIter).first] = faceOwnerCandidateRank; }
        }

        int nEdgeOwned = gridstorage_.edge_.size() - edgeNonOwned.size();
        int nFaceOwned = gridstorage_.face_.size() - faceNonOwned.size();
        int elementsOwned = gridstorage_.element_.size();


        // 3) Communicate number of edges and faces owned by each process to all
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
        	gridstorage_.nEdgeTotal_ += edgesOnProcess[iProc];
        	gridstorage_.nFaceTotal_ += facesOnProcess[iProc];

            if (iProc < rank_)
            {
                edgesBeforeMe += edgesOnProcess[iProc];
                facesBeforeMe += facesOnProcess[iProc];
                elementsBeforeMe += elementsOnProcess[iProc];
            }
        }


        // 3) Enumerate all edges, faces and elements that you own
        // *************************************************************************

        int iEdgeGlobalId = edgesBeforeMe;
        int iFaceGlobalId = facesBeforeMe;

        // Enumerating elements is simply shifting the local index, since all elements on this process are owned by it
        for (int i = 0; i < gridstorage_.element_.size(); i++)  { gridstorage_.element_[i].globalIndex = elementsBeforeMe + i; }

        // Faces that are not shared with other processes are automatically owned by this process
        for (int i = 0; i < gridstorage_.internalFaceIndex_.size(); i++)        { gridstorage_.face_[gridstorage_.internalFaceIndex_[i]].globalIndex = iFaceGlobalId++; }
        for (int i = 0; i < gridstorage_.domainBoundaryFaceIndex_.size(); i++)  { gridstorage_.face_[gridstorage_.domainBoundaryFaceIndex_[i]].globalIndex = iFaceGlobalId++; }

        for (FaceMapIterator iter = processInternalMap_.begin(); iter != processInternalMap_.end(); iter++)
        {
            // This process owns this face if it is in the processInternalMap and if it is not in the non-owned face map
            if (faceNonOwned.find((*iter).first) == faceNonOwned.end())  { gridstorage_.face_[(*iter).second].globalIndex = iFaceGlobalId++; }
        }

        for (EdgeMapIterator iter = edgemap_.begin(); iter != edgemap_.end(); iter++)
        {
            // This process owns this edge if it is in the edgemap and if it is not in the non-owned edge map
            if (edgeNonOwned.find((*iter).first) == edgeNonOwned.end())  { gridstorage_.edge_[(*iter).second].globalIndex = iEdgeGlobalId++; }
        }


        // 4) Communicate missing edge and face global indices to their corresponding neighbors. Receive them and assign.
        // *************************************************************************

        globalDistributeMissingEdgeGlobalIndex(processBoundaryEdgeMap, processBoundaryEdgeNeighborRank);
        globalDistributeMissingFaceGlobalIndex();

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Finished Constructing Global Indices");
    }


    /** Communicates the Ghost Elements
     *
     * Prerequisites:
     * * Requires all processBoundaries to know neighboring process rank
     * * Requires existence of face global index
     *
     * Aspects:
     * * Ghost elements can have different interpolation order
     *
     * Algorithm:
     * 1) Compute number of process boundaries shared with each process, as well as and the exact interpolation orders of Ghost Elements
     * 2) MPI_alltoallv - communicate to each process a list of interpolation orders of the Ghost Elements it is going to receive
     * 3) MPI_alltoallv - Package and send element globalIndex + elementPhysicalTag + all interpolatory Vertex global indices
     * 4) Add all ghost elements to the mesh. Calculate which vertices are missing from which processes
     * 5) MPI_alltoall - Communicates to each process the number of missing vertices out of the ones it had communicated.
     * 5.1) Then communicate the globalIndices's of all missing vertices
     * 5.2) Then communicate the vertex coordinates corresponding to received global indices
     * 6) Distrubute vertex coordinates and add received coordinates to the mesh
     *
     *
     * TODO: Only supports tetrahedral ghost elements at the moment
     * [FIXME] Check if everything makes sense in terms of self-to-self communication
     * [FIXME] Check if there is balance between send and receive
     *
     * */
    void generateGhostElements()
    {
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Generating Ghost Elements");

        // 1) Compute number of process boundaries shared with each process, as well as and the exact interpolation orders of Ghost Elements
        // *************************************************************************************

        // For each other process stores the set of element indices local to this process. These elements will be Ghost Elements on the other processes
        std::vector< std::vector<int> > thisProcessGhostElementIndexSet (size_, std::vector<int>() );
        // For each process store the set of face global indices, to which the corresponding Ghost Elements are associated
        std::vector< std::vector<int> > thisProcessGhostFaceGlobalIndexSet (size_, std::vector<int>() );

        for (int iPBFace = 0; iPBFace < gridstorage_.processBoundaryFaceIndex_.size(); iPBFace++ )  {
            int thisFaceIndex = gridstorage_.processBoundaryFaceIndex_[iPBFace];
            int thisNeighborRank = gridstorage_.processBoundaryNeighborProcess_[iPBFace];

            int thisGhostIndex = gridstorage_.face_[thisFaceIndex].element1Index;
            int thisFaceGlobalIndex = gridstorage_.face_[thisFaceIndex].globalIndex;
            thisProcessGhostElementIndexSet[thisNeighborRank].push_back(thisGhostIndex);
            thisProcessGhostFaceGlobalIndexSet[thisNeighborRank].push_back(thisFaceGlobalIndex);
        }

        std::vector<int> processNeighborSize;
        for (int i = 0; i < size_; i++) { processNeighborSize.push_back(thisProcessGhostElementIndexSet[i].size()); }
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Ghost elements shared with neighbor processes (" + vector2string(processNeighborSize) + ")" );


        // 2) MPI_alltoallv - communicate to each process a list of interpolation orders of the ghost elements it is going to receive
        // *************************************************************************************

        // For each other process stores the set of interpolation orders of Ghost Elements that process wishes to communicate to this process
        std::vector< std::vector<int> > neighborProcessGhostOrder (size_, std::vector<int>() );
        ghostDistributeInterpolationOrders(thisProcessGhostElementIndexSet, neighborProcessGhostOrder);


        // 3) MPI_alltoallv - Package element globalIndex + elementPhysicalTag + all interpVertex globalIds
        // *************************************************************************************
        std::vector<int> packageGhostElementData;
        ghostDistributeGhostElements(thisProcessGhostElementIndexSet, thisProcessGhostFaceGlobalIndexSet, neighborProcessGhostOrder, packageGhostElementData );

        thisProcessGhostElementIndexSet.clear();
        thisProcessGhostFaceGlobalIndexSet.clear();


        // 4) Add all ghost elements to the mesh. Calculate which vertices are missing from which processes
        // *************************************************************************************
        std::vector<std::set<int> > missingVertices (size_, std::set<int>());
        ghostInsertGhostElements(neighborProcessGhostOrder, packageGhostElementData, missingVertices);

        neighborProcessGhostOrder.clear();
        packageGhostElementData.clear();


        // 5) Communicates to each process the number of missing vertices out of the ones it had communicated
        // Then communicate the globalId's of all missing vertices
        // *************************************************************************************
        std::vector<int> packageMissingVertexGlobalIndices;
        std::vector<int> verticesRequested;
        std::vector<int> verticesToSend;

        ghostCommunicateMissingVertexGlobalIndices(missingVertices, packageMissingVertexGlobalIndices, verticesRequested, verticesToSend);


        // 6) Distrubute vertex coordinates and add received coordinates to the mesh
        // *************************************************************************************
        ghostCommunicateMissingVertexCoordinates (missingVertices, packageMissingVertexGlobalIndices, verticesRequested, verticesToSend);

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Finished Generating Ghost Elements");
    }


    /** Construct OCTree for locating tetrahedrons in mesh */

    // TODO: Use standard logging message
    // TODO: Original octree has diagnostics output under #if 0, can append at later stage

    // FIXME: Replace OCTree pointer by just an instance
    void constructOctree() {
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Started OCTree construction");

        // bounding box of whole mesh
        Vertex center = gridstorage_.boundingBoxCenter_;
        Vertex extent = gridstorage_.boundingBoxExtent_;

        // octree length is the largest component of extent
        double length = extent[0];
        if (extent[1] > length)  { length = extent[1]; }
        if (extent[2] > length)  { length = extent[2]; }

        // construct LooseOctree with large max depth
        gridstorage_.octree_ = new CurvilinearLooseOctree(center, length, 100, verbose_, processVerbose_, mpihelper_);

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
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, outputString.str());
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
     * */

    void globalCommunicateVertexNeighborRanks (
            Index2IndexMap & processBoundaryCornerMap,
            std::vector<std::vector<int> > & processBoundaryCornerNeighborRank
    )
    {
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Started communicating process boundary neighbors");

        // 1) collective_comm.max() - find the maximal number of process boundary corners per process
        // ********************************************************

        Dune::CollectiveCommunication<MPI_Comm> collective_comm = mpihelper_.getCollectiveCommunication();

        // Reserve memory for saving ranks associated to process boundary corners
        int thisProcessBoundarySize = processBoundaryCornerMap.size();

        int maxProcessBoundarySize = collective_comm.max(thisProcessBoundarySize);

        // 2) collective_comm.allgather() - communicate global index of your process boundary corner to all other processes
        // Repeat this process until every process has communicated all its corners
        // If you run out of corners, communicate fake corners
        // ********************************************************
        IndexMapIterator procCornerIter = processBoundaryCornerMap.begin();

        for (int iCorner = 0; iCorner < maxProcessBoundarySize; iCorner++)
        {
            // If all process boundary corners have been sent start sending fake corners
            int thisCornerInd = (iCorner < thisProcessBoundarySize) ? (*(procCornerIter++)).first : -1;

            // Communicate
            std::vector<int> procCornerIndexSet (size_);
            collective_comm.allgather(&thisCornerInd, 1, reinterpret_cast<int*> (procCornerIndexSet.data()) );

            // Loop over corners sent by other processes. If this corner present, note its sender rank
            for (int iProc = 0; iProc < size_; iProc++)
            {
                // Only consider non-fake corners sent by other processes
                if ((iProc != rank_) && (procCornerIndexSet[iProc] >= 0))
                {
                    // Attempt to find this corner global id among process boundary corners of this process
                	IndexMapIterator tmpIter = processBoundaryCornerMap.find(procCornerIndexSet[iProc]);

                    // If this corner is present, note its sender process
                    if (tmpIter != processBoundaryCornerMap.end()) {
                    	processBoundaryCornerNeighborRank[(*tmpIter).second].push_back(iProc);
                    }
                }
            }
        }

        // 3) Sort all neighbor rank sets, to accelerate set intersection algorithm in future
        // ********************************************************
        for (int i = 0; i < processBoundaryCornerNeighborRank.size(); i++)
        {
            std::sort(processBoundaryCornerNeighborRank[i].begin(), processBoundaryCornerNeighborRank[i].end());
        }


        // Testing output
        std::stringstream log_stream;
        log_stream << "CurvilinearGridConstructor: -- Process boundary corner";
        for (IndexMapIterator cornerIter = processBoundaryCornerMap.begin(); cornerIter != processBoundaryCornerMap.end(); cornerIter++)
        {
        	log_stream << " GlobalIndex=" << (*cornerIter).first;
        	log_stream << " has Neighbors=(" << vector2string(processBoundaryCornerNeighborRank[(*cornerIter).second]) << ")";
        }
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_stream.str());
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Finished process boundary neighbors");
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
    void globalComputeEdgeNeighborRanks(
            Index2IndexMap & processBoundaryCornerMap,
            EdgeKey2EdgeIdMap & processBoundaryEdgeMap,
            std::vector<std::vector<int> > & processBoundaryCornerNeighborRank,
            std::vector<std::vector<int> > & processBoundaryEdgeNeighborRank)
    {
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Started computing edge process boundary neighbors");

        for (EdgeMapIterator edgeIter = processBoundaryEdgeMap.begin(); edgeIter != processBoundaryEdgeMap.end(); edgeIter++ )
        {
            // Get corners of the edge
            EdgeKey thisEdgeKey = (*edgeIter).first;

            // Get neighbor processes associated with each corner
            std::vector<int> corner0neighborset = processBoundaryCornerNeighborRank[processBoundaryCornerMap[thisEdgeKey.node0]];
            std::vector<int> corner1neighborset = processBoundaryCornerNeighborRank[processBoundaryCornerMap[thisEdgeKey.node1]];

            // Find neighbors common to both edge corners
            std::vector<int> edgeneighborset = sortedSetIntersection(corner0neighborset, corner1neighborset);

            std::stringstream log_stream;
            log_stream << "CurvilinearGridConstructor: --";
            log_stream << "Neighbors[0]=(" << vector2string(corner0neighborset) << ")";
            log_stream << " Neighbors[1]=(" << vector2string(corner1neighborset) << ")";
            log_stream << " Intersection=" << vector2string(edgeneighborset);
            Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_stream.str());

            if (edgeneighborset.size() < 1) { DUNE_THROW(Dune::IOError, "CurvilinearGrid: Found no neighbor processes to an edge "); }

            // Store the edge neighbor rank set
            edgeneighborset.swap(processBoundaryEdgeNeighborRank[(*edgeIter).second]);
        }
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Finished computing edge process boundary neighbors");
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
    void globalComputeFaceNeighborRank(
            Index2IndexMap & processBoundaryCornerMap,
            std::vector<std::vector<int> > & processBoundaryCornerNeighborRank)
    {
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Started computing face process boundary neighbors");

        for (FaceMapIterator faceIter = processInternalMap_.begin(); faceIter != processInternalMap_.end(); faceIter++ )
        {
            // Get corners of the edge
            FaceKey thisFaceKey = (*faceIter).first;

            // Get neighbor processes associated with each corner
            std::vector<int> corner0neighborset = processBoundaryCornerNeighborRank[processBoundaryCornerMap[thisFaceKey.node0]];
            std::vector<int> corner1neighborset = processBoundaryCornerNeighborRank[processBoundaryCornerMap[thisFaceKey.node1]];
            std::vector<int> corner2neighborset = processBoundaryCornerNeighborRank[processBoundaryCornerMap[thisFaceKey.node2]];

            // Find neighbors common to all 3 face corners. Need to intersect sets twice
            std::vector<int> faceneighborset;
            faceneighborset = sortedSetIntersection(corner0neighborset, corner1neighborset);
            faceneighborset = sortedSetIntersection(faceneighborset,    corner2neighborset);

            std::stringstream log_stream;
            log_stream << "CurvilinearGridConstructor: --";
            log_stream << "Neighbors[0]=(" << vector2string(corner0neighborset) << ")";
            log_stream << " Neighbors[1]=(" << vector2string(corner1neighborset) << ")";
            log_stream << " Neighbors[2]=(" << vector2string(corner2neighborset) << ")";
            log_stream << " Intersection=" << vector2string(faceneighborset);
            Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_stream.str());

            if (faceneighborset.size() != 1) { DUNE_THROW(Dune::IOError, "CurvilinearGrid: Found wrong number of neighbor processes to a face"); }

            // Store the face neighbor rank. Face is only allowed to have exactly one neighbor
            gridstorage_.processBoundaryNeighborProcess_[(*faceIter).second] = faceneighborset[0];
        }

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Finished computing face process boundary neighbors");
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
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Started distributing missing face GlobalIndices");

        typedef std::pair<FaceKey, int>  FaceInfo;
        std::vector< std::vector< FaceInfo > > facesToSend (size_);

        int totalRecvSize = 0;
        int N_INTEGER_FACEINFO = 4;
        std::vector<int> sendbuf, sendcounts(size_), sdispls;
        std::vector<int> recvbuf, recvcounts(size_), rdispls;

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: -- Checking which faces are missing");


        // 1) Loop over all process boundary faces, split faces on the ones to be sent and to be received
        // ********************************************************************************************
        for (FaceMapIterator iter = processInternalMap_.begin(); iter != processInternalMap_.end(); iter++)
        {
            int localProcessBoundaryIndex = (*iter).second;
            int localFaceIndex = gridstorage_.processBoundaryFaceIndex_[localProcessBoundaryIndex];
            int neighborRank = gridstorage_.processBoundaryNeighborProcess_[localProcessBoundaryIndex];

            // If the neighbor of this face has lower rank, then add it to recv, else to send
            if (neighborRank < rank_)  { recvcounts[neighborRank] += N_INTEGER_FACEINFO; }
            else
            {
                int thisGlobalIndex = gridstorage_.face_[localFaceIndex].globalIndex;
                facesToSend[neighborRank].push_back(FaceInfo((*iter).first, thisGlobalIndex ));
            }
        }

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: -- Assembling arrays to send");


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

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: -- Sending  sendcounts=(" + vector2string(sendcounts) + ") recvcounts=(" + vector2string(recvcounts) + ")");

        // 3) MPI_alltoall key + globalId
        // ********************************************************************************************
        recvbuf.resize(totalRecvSize, 0);
        MPI_Comm comm = Dune::MPIHelper::getCommunicator();
        MPI_Alltoallv (sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(recvbuf.data()), recvcounts.data(), rdispls.data(), MPI_INT, comm );


        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: -- Extracting missing indices");

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
                int thisGlobalId = recvbuf[iData++];
                thisKey.node0 = recvbuf[iData++];
                thisKey.node1 = recvbuf[iData++];
                thisKey.node2 = recvbuf[iData++];

                FaceMapIterator faceIter = processInternalMap_.find(thisKey);

                if (faceIter == processInternalMap_.end()) { DUNE_THROW(Dune::IOError, "CurvilinearGrid: Communicated FaceKey does not correspond to any face on this process "); }
                else
                {
                    int localFaceIndex = gridstorage_.processBoundaryFaceIndex_[(*faceIter).second];
                    gridstorage_.face_[localFaceIndex].globalIndex = thisGlobalId;
                }
            }
        }

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Finished distributing missing face GlobalIndices");
    }


    /** Communicates all process boundary face global Id's to the neighbors if owned
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
    void globalDistributeMissingEdgeGlobalIndex(
            EdgeKey2EdgeIdMap & processBoundaryEdgeMap,
            std::vector<std::vector<int> > & processBoundaryEdgeNeighborRank)
    {
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Started distributing missing edge GlobalIndices");

        typedef std::pair<EdgeKey, int>  EdgeInfo;
        std::vector< std::vector< EdgeInfo > > edgesToSend (size_);

        int totalRecvSize = 0;
        int N_INTEGER_EDGEINFO = 3;
        std::vector<int> sendbuf, sendcounts(size_), sdispls;
        std::vector<int> recvbuf, recvcounts(size_), rdispls;

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: -- Checking which edges are missing");

        // 1) Loop over all process boundary faces, split faces on the ones to be sent and to be received
        // ********************************************************************************************
        for (EdgeMapIterator iter = processBoundaryEdgeMap.begin(); iter != processBoundaryEdgeMap.end(); iter++)
        {
            EdgeKey thisEdgeKey = (*iter).first;
            int localProcessBoundaryIndex = (*iter).second;

            EdgeMapIterator thisEdgeMapIter = edgemap_.find(thisEdgeKey);
            if (thisEdgeMapIter == edgemap_.end()) {
            	std::cout << ":((" << std::endl;
            	//std::cout << "process_" << rank_ << "ThisEdgeKey=(" << thisEdgeKey.node0 << "," << thisEdgeKey.node1 << ")" << std::endl;
            	DUNE_THROW(Dune::IOError, "CurvilinearGrid: Process boundary edge not found among all edges");
            }

            int localEdgeIndex = (*thisEdgeMapIter).second;
            int candidateOwnerRank = processBoundaryEdgeNeighborRank[localProcessBoundaryIndex][0];

            //std::cout << "process_" << rank_ <<  " EdgeKey=(" << thisEdgeKey.node0 << "," << thisEdgeKey.node1 << ") localIndex=" << localEdgeIndex <<  std::endl;

            // If the one of the neighbors of this edge has lower rank, then note one more received edge from that process
            // else note to send it to all other neighbors
            if (candidateOwnerRank < rank_)  { recvcounts[candidateOwnerRank] += N_INTEGER_EDGEINFO; }
            else
            {
                int thisGlobalIndex = gridstorage_.edge_[localEdgeIndex].globalIndex;

                EdgeInfo thisEdgeInfo(thisEdgeKey, thisGlobalIndex);

                for (int iNeighbor = 0; iNeighbor < processBoundaryEdgeNeighborRank[localProcessBoundaryIndex].size(); iNeighbor++)
                {
                    int thisNeighborRank = processBoundaryEdgeNeighborRank[localProcessBoundaryIndex][iNeighbor];
                    edgesToSend[thisNeighborRank].push_back(thisEdgeInfo);
                };
            }
        }

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: -- Assembling arrays to send");


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

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: -- Sending  sendcounts=(" + vector2string(sendcounts) + ") recvcounts=(" + vector2string(recvcounts) + ")");



        // 3) MPI_alltoall key + globalId
        // ********************************************************************************************
        recvbuf.resize(totalRecvSize, 0);   // There are 3 integers per sent FaceInfo
        MPI_Comm comm = Dune::MPIHelper::getCommunicator();
        MPI_Alltoallv (sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(recvbuf.data()), recvcounts.data(), rdispls.data(), MPI_INT, comm );


        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: -- Extracting missing indices");

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
                int thisGlobalId = recvbuf[iData++];
                thisKey.node0 = recvbuf[iData++];
                thisKey.node1 = recvbuf[iData++];

                EdgeMapIterator edgeIter = edgemap_.find(thisKey);

                if (edgeIter == edgemap_.end()) {
                	std::cout << " 0_o" << std::endl;
                	DUNE_THROW(Dune::IOError, "CurvilinearGrid: Communicated EdgeKey does not correspond to any edge on this process "); }
                else
                {
                    int localEdgeIndex = (*edgeIter).second;
                    gridstorage_.edge_[localEdgeIndex].globalIndex = thisGlobalId;
                }
            }
        }

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Finished distributing missing edge GlobalIndices");
    }


    /* ***************************************************************************
     * SubSection: Auxiliary methods of generateGhostElements()
     * ***************************************************************************/


    /** Communicate to each process a list of interpolation orders of the ghost elements it is going to receive
     *
     * Algorithm:
     *
     * 1) Put interpolation orders of elements of this process that will become Ghost Elements on other processes
     *    into send array in the correct order
     * 2) Communicate
     * 3) Store interpolation orders of Ghost Elements that will be sent to this process from each other process
     *
     * Optimization Proposal:
     * Same as for FaceGlobalIndex algorithm
     *
     * */
    void ghostDistributeInterpolationOrders(
            const std::vector< std::vector<int> > & thisProcessGhostElementIndexSet,
            std::vector< std::vector<int> > & neighborProcessGhostOrder
    )
    {
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Communicate ghost element interpolatory orders");
        std::vector<int> sendbuf, sendcounts, sdispls;
        std::vector<int> recvbuf, recvcounts, rdispls;

        // We should receive the number interpolation orders equal to the number of process boundaries
        recvbuf.resize(gridstorage_.processBoundaryFaceIndex_.size(), 0);

        for (int iProc = 0; iProc < size_; iProc++)
        {
            for (int iElem = 0; iElem < thisProcessGhostElementIndexSet[iProc].size(); iElem++) {
                sendbuf.push_back(gridstorage_.element_[thisProcessGhostElementIndexSet[iProc][iElem]].interpOrder);
            }
            sendcounts.push_back(thisProcessGhostElementIndexSet[iProc].size());
            sdispls.push_back((iProc == 0) ? 0 : sdispls[iProc-1] + sendcounts[iProc-1] );

            // For ghost elements we send/receive same amount to/from each processor
            recvcounts.push_back(sendcounts[iProc]);
            rdispls.push_back(sdispls[iProc]);
        }

        MPI_Comm comm = Dune::MPIHelper::getCommunicator();
        MPI_Alltoallv (sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(recvbuf.data()), recvcounts.data(), rdispls.data(), MPI_INT, comm );

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Interpolatory orders of neighbor ghosts = (" + vector2string(recvbuf) + ")");


        // Parse the received data
        // Note that with Ghost Elements, we receive from each neighbor process as many entities as we send to it
        int iData = 0;
        for (int iProc = 0; iProc < size_; iProc++)
        {
            for (int iElem = 0; iElem < thisProcessGhostElementIndexSet[iProc].size(); iElem++)
            {
                neighborProcessGhostOrder[iProc].push_back(recvbuf[iData++]);
            }
        }
    }


    /** Communicate to all ghost element information except of the explicit vertex coordinates
     * MPI_alltoallv - Package element globalIndex + elementPhysicalTag + assocFaceGlobalIndex + all interpVertex globalIds
     *
     * Optimization Proposal:
     * It is possible not to send the associated face global index, and instead figure out the order of elements
     * according to the order of associated FaceKeys. The gain however is very small
     *
     * TODO: If planning to use with non-tetrahedral meshes, need to pass element type as well
     *
     * */
    void ghostDistributeGhostElements(
            std::vector< std::vector<int> > & thisProcessGhostElementIndexSet,
            std::vector< std::vector<int> > & thisProcessGhostFaceGlobalIndexSet,
            std::vector< std::vector<int> > & neighborProcessGhostOrder,
            std::vector<int> & recvbuf )
    {
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Communicate ghost element data");

        std::vector<int> sendbuf, sendcounts, sdispls;
        std::vector<int> recvcounts, rdispls;

        // Calculates total amount of integers to receive during DoF communication stage
        int totalRecvSize = 0;

        Dune::GeometryType meshGeometryType;
        meshGeometryType.makeTetrahedron();

        for (int i = 0; i < thisProcessGhostElementIndexSet.size(); i++)
        {
            int thisSendCounts = 0;
            int thisRecvCounts = 0;

            for (int j = 0; j < thisProcessGhostElementIndexSet[i].size(); j++)
            {
                int ghostElementIndex = thisProcessGhostElementIndexSet[i][j];
                int ghostElementFaceGlobalIndex = thisProcessGhostFaceGlobalIndexSet[i][j];
                int thisDofNum = gridstorage_.element_[ghostElementIndex].vertexIndexSet.size();
                thisSendCounts += 3 + thisDofNum;
                thisRecvCounts += 3 + Dune::CurvilinearGeometryHelper::dofPerOrder(meshGeometryType, neighborProcessGhostOrder[i][j]);
                totalRecvSize += thisRecvCounts;

                sendbuf.push_back(gridstorage_.element_[ghostElementIndex].globalIndex);
                sendbuf.push_back(gridstorage_.element_[ghostElementIndex].physicalTag);
                sendbuf.push_back(ghostElementFaceGlobalIndex);
                for (int iDof = 0; iDof < thisDofNum; iDof++)
                {
                    int localVertexIndex = gridstorage_.element_[thisProcessGhostElementIndexSet[i][j]].vertexIndexSet[iDof];
                	sendbuf.push_back(gridstorage_.point_[localVertexIndex].globalIndex);
                }
            }

            sendcounts.push_back(thisSendCounts);
            recvcounts.push_back(thisRecvCounts);
            sdispls.push_back((i == 0) ? 0 : sdispls[i-1] + sendcounts[i-1] );
            rdispls.push_back((i == 0) ? 0 : rdispls[i-1] + recvcounts[i-1] );
        }

        std::cout << "process_" << rank_ << " sendbuf=(" << vector2string(sendbuf) << ")" << std::endl;

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Communicating ghost elements sendcounts=(" + vector2string(sendcounts) + ") recvcounts=(" + vector2string(recvcounts) + ")" );

        recvbuf.resize(totalRecvSize, 0);

        MPI_Comm comm = Dune::MPIHelper::getCommunicator();
        MPI_Alltoallv (sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(recvbuf.data()), recvcounts.data(), rdispls.data(), MPI_INT, comm );

        std::cout << "process_" << rank_ << " recvbuf=(" << vector2string(recvbuf) << ")" << std::endl;
    }


    /** Add received elements to the mesh. For each vertex global index, find if coordinate is already present on this process
     *  If not, mark this vertex as a missing vertex for further communication.
     *
     *
     * */
    void ghostInsertGhostElements (
            std::vector< std::vector<int> > & neighborProcessGhostOrder,
            std::vector< int > & packageGhostElementData,
            std::vector<std::set<int> > & missingVertices
    )
    {
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Inserting communicated ghost elements");

        int iData = 0;
        Dune::GeometryType meshGeometryType;
        meshGeometryType.makeTetrahedron();

        for (int iProc = 0; iProc < neighborProcessGhostOrder.size(); iProc++) {

            for (int iGhost = 0; iGhost < neighborProcessGhostOrder[iProc].size(); iGhost++)
            {
                EntityStorage thisElement;
                thisElement.geometryType = meshGeometryType;
                thisElement.globalIndex = packageGhostElementData[iData++];
                thisElement.interpOrder = neighborProcessGhostOrder[iProc][iGhost];
                thisElement.physicalTag = packageGhostElementData[iData++];

                int associatedFaceGlobalIndex = packageGhostElementData[iData++];

                int thisElementDof = Dune::CurvilinearGeometryHelper::dofPerOrder(meshGeometryType, thisElement.interpOrder);
                for (int iDof = 0; iDof < thisElementDof; iDof++)
                {
                    int thisVertexGlobalIndex = packageGhostElementData[iData++];
                    IndexMapIterator vertexIter = gridstorage_.vertexGlobal2LocalMap_.find(thisVertexGlobalIndex);

                    // If this vertex already exists, just reuse it. Otherwise, create new vertex, and later request its coordinate
                    if (vertexIter != gridstorage_.vertexGlobal2LocalMap_.end()) {
                        int thisVertexLocalIndex = (*vertexIter).second;
                        thisElement.vertexIndexSet.push_back(thisVertexLocalIndex);

                        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Ghost Vertex already on this process GlobalIndex=" + std::to_string(thisVertexGlobalIndex));
                    }
                    else
                    {
                    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Ghost Vertex missing on this process GlobalIndex=" + std::to_string(thisVertexGlobalIndex));

                        // Create a new vertex with local index pointing to the end of current vertex array
                        int localVertexIndex = gridstorage_.point_.size();
                    	thisElement.vertexIndexSet.push_back(localVertexIndex);

                        // Insert the fake vertex into the mesh
                        Vertex fakeCoord;
                        insertVertex(fakeCoord, thisVertexGlobalIndex);

                        // Note that this vertex needs communicating
                        missingVertices[iProc].insert(thisVertexGlobalIndex);
                    }
                }

                // Associate a (hopefully processBoundary) face with this ghost element
                IndexMapIterator faceIndexIter = gridstorage_.faceGlobal2LocalMap_.find(associatedFaceGlobalIndex);
                if (faceIndexIter == gridstorage_.faceGlobal2LocalMap_.end())  {
                	std::cout << ":<  " << associatedFaceGlobalIndex << std::endl;

                	DUNE_THROW(Dune::IOError, "CurvilinearGridConstructor: Received Ghost process boundary face not found among faces of this process");
                }

                int localFaceIndex = (*faceIndexIter).second;
                int localGhostIndex = gridstorage_.ghostElement_.size();
                gridstorage_.face_[localFaceIndex].element2Index = localGhostIndex;
                gridstorage_.ghostElement_.push_back(thisElement);
            }
        }
    }


    //
    //

    /** Communicate to each process the number of missing vertices out of the ones it had provided with Ghost Elements
     *  Then communicate the globalIndices of all missing vertices
     *
     * */
    void ghostCommunicateMissingVertexGlobalIndices(
            std::vector<std::set<int> > & missingVertices,
            std::vector<int> & recvbuf,
            std::vector<int> & verticesRequested,
            std::vector<int> & verticesToSend
    )
    {
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Communicating ghost element missing vertex indices");

        std::vector<int> sendbuf;
        recvbuf.resize(size_, 0);

        // 4.1) MPI_alltoallv - tell each process the number of coordinates you want from it
        for (int iProc = 0; iProc < size_; iProc++)  { sendbuf.push_back(missingVertices[iProc].size()); }

        MPI_Comm comm = Dune::MPIHelper::getCommunicator();
        MPI_Alltoall(sendbuf.data(), 1, MPI_INT, reinterpret_cast<int*>(recvbuf.data()), 1, MPI_INT, comm);


        // 4.2) MPI_alltoallv - tell each process the list of global Indices of coordinates you want from it
        // Cleanup
        std::vector<int> sendcounts, sdispls, recvcounts, rdispls;
        sendbuf.clear();
        int totalRecvSize = 0;
        int iData = 0;

        for (int iProc = 0; iProc < size_; iProc++)
        {
            int thisSendSize = missingVertices[iProc].size();
            int thisRecvSize = recvbuf[iData++];

            recvcounts.push_back(thisRecvSize);
            sendcounts.push_back(thisSendSize);
            totalRecvSize += thisRecvSize;

            for (std::set<int>::iterator mVertIter = missingVertices[iProc].begin(); mVertIter != missingVertices[iProc].end(); mVertIter++ )  { sendbuf.push_back(*mVertIter); }

            sdispls.push_back((iProc == 0) ? 0 : sdispls[iProc-1] + sendcounts[iProc-1] );
            rdispls.push_back((iProc == 0) ? 0 : rdispls[iProc-1] + recvcounts[iProc-1] );
        }

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Missing vertex numbers sendcounts=(" + vector2string(sendcounts) + ") recvcounts=(" + vector2string(recvcounts) + ")" );

        recvbuf.clear(); recvbuf.resize(totalRecvSize, 0);
        MPI_Alltoallv (sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(recvbuf.data()), recvcounts.data(), rdispls.data(), MPI_INT, comm );

        // We will require the information about requested and sent vertices when we communicate the coordinates
        verticesRequested.swap(sendcounts);
        verticesToSend.swap(recvcounts);
    }


    // Distrubute vertex coordinates and add received coordinates to the mesh
    void ghostCommunicateMissingVertexCoordinates (
            std::vector<std::set<int> > & missingVertices,
            std::vector<int> & packageMissingVertexGlobalIndices,
            std::vector<int> & verticesRequested,
            std::vector<int> & verticesToReceive
    )
    {
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: Communicating ghost element missing vertex coordinates");

        // 4.3) MPI_alltoallv - package list of globalId+coordinate for each process and send it
        std::vector<int> sendcounts(size_), sdispls;
        std::vector<int> recvcounts(size_), rdispls;
        std::vector<double> recvbuf, sendbuf;

        int iData = 0;
        int totalRecvSize = 0;

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridConstructor: verticesToReceive=(" + vector2string(verticesToReceive) + ") verticesRequestedPerProcess=(" + vector2string(verticesRequested) + ")" + " missingGlobalIndicesPackage=(" + vector2string(packageMissingVertexGlobalIndices) + ")" );


        for (int i = 0; i < size_; i++)
        {
            // Go through all vertices requested from this process. Package coordinates
            for (int j = 0; j < verticesRequested[i]; j++)
            {
                int thisVertexGlobalIndex = packageMissingVertexGlobalIndices[iData++];
                int thisVertexLocalIndex = gridstorage_.vertexGlobal2LocalMap_[thisVertexGlobalIndex];

                Vertex p = gridstorage_.point_[thisVertexLocalIndex].coord;

                for (int iDim = 0; iDim < 3; iDim++)  { sendbuf.push_back(p[iDim]); }
            }

            // We communicate (coord = 3 doubles) for each sent/received vertex
            // We now receive the amount we sent before, and send the amount we received before
            int thisSendSize = 3 * verticesRequested[i];
            int thisRecvSize = 3 * verticesToReceive[i];

            sendcounts[i] = thisSendSize;
            recvcounts[i] = thisRecvSize;
            totalRecvSize += thisRecvSize;

            sdispls.push_back((i == 0) ? 0 : sdispls[i-1] + sendcounts[i-1] );
            rdispls.push_back((i == 0) ? 0 : rdispls[i-1] + recvcounts[i-1] );
        }

        recvbuf.resize(totalRecvSize, 0);

        MPI_Comm comm = Dune::MPIHelper::getCommunicator();
        MPI_Alltoallv (sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_DOUBLE, reinterpret_cast<double*>(recvbuf.data()), recvcounts.data(), rdispls.data(), MPI_DOUBLE, comm );

        // Assign coordinates to all missing vertices
        iData = 0;

        for (int i = 0; i < missingVertices.size(); i++)
        {
            for (std::set<int>::iterator mVertIter = missingVertices[i].begin(); mVertIter != missingVertices[i].end(); mVertIter++ )
            {
                Vertex thisCoord;
                thisCoord[0] = recvbuf[iData++];
                thisCoord[1] = recvbuf[iData++];
                thisCoord[2] = recvbuf[iData++];

                int thisVertexGlobalIndex = *mVertIter;
                int thisVertexLocalIndex = gridstorage_.vertexGlobal2LocalMap_[thisVertexGlobalIndex];
                gridstorage_.point_[thisVertexLocalIndex].coord = thisCoord;
            }
        }
    }



private: // Private members

    bool verbose_;
    bool processVerbose_;
    bool withGhostElements_;

    // Temporary maps necessary to locate and communicate entities during grid base construction
    EdgeKey2EdgeIdMap edgemap_;                // (global edgeKey -> edge_ index)
    FaceKey2FaceIdMap internalInternalMap_;    // (global faceKey -> internalFaceIndex_)
    FaceKey2FaceIdMap boundaryInternalMap_;    // (global faceKey -> domainBoundaryFaceIndex_)
    FaceKey2FaceIdMap processInternalMap_;     // (global faceKey -> processBoundaryFaceIndex_)

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
