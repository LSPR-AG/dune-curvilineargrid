#ifndef DUNE_GLOBALBOUNDARYCONTAINER_HH
#define DUNE_GLOBALBOUNDARYCONTAINER_HH

#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <iostream>
#include <assert.h>
#include <numeric>

#include <fstream>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>


namespace Dune {

namespace CurvGrid {




/** \brief Storage class for boundary segments. Minimal info needed for communication */
template<class Grid, int NVERTEX, int NEDGE, int NCORNER>
struct BoundarySegmentContainer {
	static const int dimension = Grid::dimension;
	static const int ELEMENT_CODIM = 0;
	static const int FACE_CODIM = 1;

	typedef typename Grid::LeafGridView   LeafGridView;
	typedef typename LeafGridView::template Codim<FACE_CODIM>::Geometry	GeometryFace;
	typedef typename GeometryFace::LocalCoordinate					LocalCoordinate;
	typedef typename GeometryFace::GlobalCoordinate				GlobalCoordinate;

    typedef typename Grid::GlobalIdSet             GlobalIdSet;
    typedef typename GlobalIdSet::IdType           GlobalId;

    typedef unsigned int UInt;

	/************************
	 * Storage
	 ************************/

	GlobalCoordinate p_[NVERTEX];	// Vertices of the boundary segment in the original order
	GlobalCoordinate n_;						// Unit outer normal of this boundary segment, pointing to the outside of the domain
	UInt order_;										// Curvilinear order of the entity
	UInt gindface_;								// Global index of the face associated with the boundary segment
	UInt gindedge_[NEDGE];				// Global index of all edges of this face in the correct subentity order
	UInt gindcorner_[NCORNER]	;		// Global index of all corners of this face in the correct subentity order
};



/** \brief This class uses all-to-all communication to gather geometric information
 * of all boundary segments (faces) on all processes, except those normally contained on each process (to save time)
 *
 * This method is not expected to be used on large meshes, as all-to-all communication is too expensive for that.
 * It is designed to enable the implementation of all-to-all coupling boundary methods such as Boundary Integral in electromagnetics.
 *
 * Algorithm:
 * 1) Loop over all boundary segments on all processes
 * 1.1) Assemble all data necessary for communication into a vector
 * 2) Mpi-all-to-all communicate the array to all processes except self
 * 3) Store this array for further access by an iterator
 *
 * */
template <class Grid>
class GlobalBoundaryContainer
{
	static const int dimension = Grid::dimension;
	static const int ELEMENT_CODIM = 0;
	static const int FACE_CODIM = 1;
	static const int EDGE_CODIM = 2;
	static const int VERTEX_CODIM = 3;

	// Geometry Typedefs



	// Grid Typedefs
	typedef typename Grid::ctype CoordinateType;
	typedef typename Grid::LeafGridView   LeafGridView;
	typedef typename Grid::template Codim<FACE_CODIM>::EntityGeometryMappingImpl  BaseGeometryFace;
	typedef typename BaseGeometryFace::LocalCoordinate					LocalCoordinate;
	typedef typename BaseGeometryFace::GlobalCoordinate				GlobalCoordinate;
	typedef unsigned int UInt;
	typedef typename LeafGridView::template Codim<FACE_CODIM>::Entity			EntityFace;
	typedef typename LeafGridView::template Codim<EDGE_CODIM>::Entity			EntityEdge;
	typedef typename LeafGridView::template Codim<VERTEX_CODIM>::Entity		EntityVertex;

	typedef Dune::ReferenceElements<CoordinateType, dimension> ReferenceElements3d;
	typedef Dune::ReferenceElements<CoordinateType, dimension-1> ReferenceElements2d;

	// BoundaryContainer Typedefs

	typedef BoundarySegmentContainer <Grid, 3, 3, 3> ContainerTriangleLinear;  // [TODO] Extend this to curvilinear case some day
	typedef std::vector<ContainerTriangleLinear>	ContainerVector;
	typedef typename ContainerVector::const_iterator  ContainerVectorConstIter;

public:

	GlobalBoundaryContainer(const Grid & grid) :
		grid_(grid),
		rank_(grid.mpihelper().rank()),
		size_(grid.mpihelper().size())
	{


		/****************************************************
		 *   Stage 1: Construct Local Domain Boundary Data
		 ****************************************************/

  		/** \brief Iterate ove all elements of Interior Border partition */
		UInt elemCount = 0;
  		UInt nElementInterior = grid_.numInternal(ELEMENT_CODIM);
  		UInt nBoundarySegment = grid_.numBoundarySegments();
  		LocalCoordinate centerFaceLocal = ReferenceElements2d::simplex().position( 0, 0 );

  		std::vector<ContainerTriangleLinear> boundaryContainerThis(nBoundarySegment);

		LeafGridView leafView = grid_.leafGridView();
		for (auto&& elemThis : elements(leafView, Dune::Partitions::interiorBorder))
		{
			LoggingMessage::writePatience("[[ Preparing domain boundary for global communication...", elemCount++, nElementInterior);

			for (auto&& interThis : intersections(leafView, elemThis))
			{
				/*************************************************************************/
				/** Domain Boundary                                                                                       */
				/*************************************************************************/
				if ( (interThis.boundary() == true) && (interThis.neighbor() == false) )
				{
					UInt boundaryIndex = interThis.boundarySegmentIndex();
					UInt intrIndexInInside = interThis.indexInInside();
					EntityFace faceThis = elemThis.template subEntity<FACE_CODIM>(intrIndexInInside);
					BaseGeometryFace faceGeomBase =  grid_.template entityBaseGeometry<FACE_CODIM>(faceThis);

					UInt nEdge = 3;  // [TODO] Replace me with correct expression from a reference element
					UInt nVertex = faceGeomBase.nVertex();
					UInt nCorner = faceGeomBase.nCorner();

					ContainerTriangleLinear & thisCont = boundaryContainerThis[boundaryIndex];

					// [TODO] This paradigm does not scale up in Curvilinear, as normal changes value over the face.
					// Need to design a method that computes the normal just from the face coordinates, without the
					// need for containing element. Then, only need communicate the sign of the outer normal
					thisCont.n_ = interThis.unitOuterNormal(centerFaceLocal);
					thisCont.gindface_ = grid_.template entityGlobalIndex<FACE_CODIM>(faceThis);
					thisCont.order_ = grid_.template entityInterpolationOrder<ELEMENT_CODIM>(elemThis);

					for (UInt iEdge = 0; iEdge < nEdge; iEdge++)  {
						UInt indexEdgeInElem = ReferenceElements3d::general(elemThis.type()).subEntity(intrIndexInInside, FACE_CODIM, 0, EDGE_CODIM);
						EntityEdge edgeThis = elemThis.template subEntity<EDGE_CODIM>(indexEdgeInElem);
						thisCont.gindedge_[iEdge] = grid_.template entityGlobalIndex<EDGE_CODIM>(edgeThis);
					}
					for (UInt iCorner = 0; iCorner < nCorner; iCorner++)  {
						UInt indexCornerInElem = ReferenceElements3d::general(elemThis.type()).subEntity(intrIndexInInside, FACE_CODIM, 0, VERTEX_CODIM);
						EntityVertex cornerThis = elemThis.template subEntity<VERTEX_CODIM>(indexCornerInElem);
						thisCont.gindcorner_[iCorner] = grid_.template entityGlobalIndex<VERTEX_CODIM>(cornerThis);
					}
					for (UInt iVertex = 0; iVertex < nVertex; iVertex++)  { thisCont.p_[iVertex] = faceGeomBase.vertex(iVertex);}
				}
			}
		}


		/********************************************************************
		 *   Stage 2: MPI communicate boundary data to all other processes
		 ********************************************************************/
		MPI_Comm comm = Dune::MPIHelper::getCommunicator();
		int mpi_err;

		// Algorithm

		/*********************************************************
		 *  1) Preliminary steps - compute sizes for communication
		 *********************************************************/
		int structSize = sizeof(ContainerTriangleLinear);	// the size of the communiated structure
		int nCommThis = boundaryContainerThis.size();		// the number of structures to be communicated
		int sizeCommThis = nCommThis * structSize;			// the total communication size in bytes

		/*********************************************************************************
		 *  2) Send number of structures this proc will communicate to all other processes
		 *********************************************************************************/
		std::vector<int> nCommRecv(size_, 0);
		mpi_err = MPI_Allgather(&nCommThis, 1, MPI_INT, reinterpret_cast<int *>(nCommRecv.data()), size_, MPI_INT, comm);
		int nCommRecvTot = std::accumulate(nCommRecv.begin(), nCommRecv.end(), 0);		// The total number of structrures that will be received by any one process
		int sizeCommRecvTot = nCommRecvTot * structSize;														// The total number of bytes that will be received by any one process

		/*********************************************************************************
		 *  3) MPI_COMM local struct array to all faces except this one
		 *********************************************************************************/
		boundaryContainerGlobal_.resize(nCommRecvTot);   // Reserve memory in boundary container for structures that will be received
		mpi_err = MPI_Allgather(
				reinterpret_cast<ContainerTriangleLinear *>(boundaryContainerThis.data()),
				sizeCommThis,
				MPI_BYTE,
				reinterpret_cast<ContainerTriangleLinear *>(boundaryContainerGlobal_.data()),
				sizeCommRecvTot,
				MPI_BYTE,
				comm
		);

		/*********************************************************************************
		 *  4) Explicitly remove all boundary segments already stored on this process
		 *********************************************************************************/
		// Find the index of data associated with this same rank within the boundaryContainerGlobal
		int startThisData = 0;
		for (int i = 0; i < rank_ - 1; i++) { startThisData += nCommRecv[i]; }

		// Swap all data associated with this process with the data at the end of the array
		if (rank_ != size_ - 1) {  // No point in doing this operation if this process data is already at the end
			for (int i = 0; i < nCommThis; i++) {
				std::swap(boundaryContainerGlobal_[startThisData + i], boundaryContainerGlobal_[nCommRecvTot - 1 - i]);
			}
		}

		// Truncate the array to get rid of the unnecessary entities
		boundaryContainerGlobal_.resize(nCommRecvTot - nCommThis);
	}


	ContainerVectorConstIter begin()  { return boundaryContainerGlobal_.begin(); }
	ContainerVectorConstIter end()  { return boundaryContainerGlobal_.end(); }


private:
	const Grid & grid_;
	int rank_;
	int size_;

	std::vector<ContainerTriangleLinear> boundaryContainerGlobal_;


};



/** \brief This class is an iterator over GlobalBoundaryContainer.
 * It is designed to extract useful information about the boundary faces, such as geometry and global indices
 * */
template <class Grid>
class GlobalBoundaryIterator {

	static const int ELEMENT_CODIM = 0;
	static const int FACE_CODIM = 1;
	static const int EDGE_CODIM = 2;
	static const int VERTEX_CODIM = 3;

	static const int dimension = Grid::dimension;
	typedef typename Grid::ctype  CoordinateType;
	typedef unsigned int UInt;

	typedef Dune::ReferenceElements<CoordinateType, dimension> ReferenceElements2D;

	typedef GlobalBoundaryIterator<Grid> This;
	typedef GlobalBoundaryContainer<Grid> BaseContainer;
	typedef typename BaseContainer::ContainerVectorConstIter  BaseConstIter;

	typedef typename Grid::template Codim<FACE_CODIM>::EntityGeometryMappingImpl  BaseGeometryFace;
	typedef typename Grid::template Codim<EDGE_CODIM>::EntityGeometryMappingImpl  BaseGeometryEdge;

	typedef typename BaseGeometryFace::GlobalCoordinate  GlobalCoordinate;

public:

	GlobalBoundaryIterator(const BaseContainer & baseContainer) :
		baseContainer_(baseContainer),
		baseiter_(baseContainer.begin())
	{

	}


	/** \brief Extract GlobalIndex of a subentity of the boundary segment
	 *  Allowed codimensions:  (1 : Face itself), (2 : SubEdge), (3 : SubCorner) */
	template <int codim>
	UInt globalIndex(UInt subIndex) {
		Dune::GeometryType gt;   gt.makeTriangle();
		assert(subIndex < ReferenceElements2D::general(gt).size(codim));  // Do not allow subentity index out of bounds

		switch(codim) {
		case FACE_CODIM : return baseiter_ -> gindface_;
		case EDGE_CODIM : return baseiter_ -> gindedge_[subIndex];
		case VERTEX_CODIM : return baseiter_ -> gindcorner_[subIndex];
		default : DUNE_THROW(Dune::IOError, "BoundaryContainer: Unexpected Subentity Index");
		}
	}

	/** \brief Extract curvilinear order of the boundary segment */
	UInt order()  { return baseiter_ -> order_; }

	/** \brief Extract unit outer normal of the boundary segment, pointing to the outside of the domain */
	GlobalCoordinate unitOuterNormal() { return baseiter_ -> n_; }

	/** \brief Extract geometry of the associated face */
	BaseGeometryFace geometry() {
		Dune::GeometryType gt;   gt.makeTriangle();

		GlobalCoordinate (&pFaceArr)[] = baseiter_ -> p_;
		std::vector<GlobalCoordinate> pFace(pFaceArr, pFaceArr + sizeof(pFaceArr) / sizeof(pFaceArr[0]));
		return BaseGeometryFace(gt, pFace, order());
	}

	/** \brief Extract geometry of a subentity edge of the associated face, given by the subentity index */
	BaseGeometryEdge geometryEdge(UInt subIndex) {
		Dune::GeometryType gtEdge;   gtEdge.makeLine();
		std::vector<int> edgeVertexInteriorInd = CurvilinearGeometryHelper::subentityInternalCoordinateSet(gtEdge, order(), EDGE_CODIM, subIndex);
		std::vector<GlobalCoordinate> pEdge;
		for (UInt i = 0; i < edgeVertexInteriorInd.size(); i++)  { pEdge.push_back((baseiter_->p_)[edgeVertexInteriorInd[i]]); }
		return BaseGeometryEdge(gtEdge, pEdge, order());
	}

	/** \brief Increase this iterator */
	This & operator ++() {
		baseiter_++;
		return *this;
	}

	/** \brief Check if this iterator has reached its end */
	bool end() { return (baseiter_ == baseContainer_.end()); }


private:

	const BaseContainer & baseContainer_;
	BaseConstIter baseiter_;

};











} // namespace CurvGrid

} // namespace Dune

#endif  // DUNE_GLOBALBOUNDARYCONTAINER_HH
