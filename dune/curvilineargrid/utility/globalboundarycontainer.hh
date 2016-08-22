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




/** \brief Storage class for arbitrary global closed boundary within the domain. Minimal info needed for communication */
template<class Grid>
struct BoundarySegmentContainer {
	static const int dimension = Grid::dimension;
	typedef typename Grid::ctype  CoordinateType;
	typedef Dune::FieldVector<CoordinateType, dimension>  GlobalCoordinate;

    typedef unsigned int UInt;

	/************************
	 * Storage
	 ************************/
    UInt order_;															// Curvilinear order of the entity
    std::vector<UInt> subBoundarySegInd_;			// Subentity boundary segment index set
    std::vector<GlobalCoordinate> p_;					// Vertices of the boundary face in the original order
	std::vector<std::vector<UInt> > subGInd_;		// Global index of all subentities of this element by codim, in correct subentity order
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

public:

	// Grid Typedefs
	typedef typename Grid::ctype CoordinateType;
	typedef typename Grid::LeafGridView   LeafGridView;
	typedef typename Grid::PhysicalTagType  PhysicalTagType;
	typedef typename Grid::template Codim<ELEMENT_CODIM>::EntityGeometryMappingImpl  BaseGeometryElement;
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

	typedef BoundarySegmentContainer<Grid> BoundaryContainer;
	typedef std::vector<BoundaryContainer>	ContainerVector;
	typedef typename ContainerVector::const_iterator  ContainerVectorConstIter;

public:

	GlobalBoundaryContainer(const Grid & grid, bool findDomainBoundary, PhysicalTagType volumeTag = 0, PhysicalTagType surfaceTag = 0) :
		grid_(grid),
		rank_(grid.mpihelper().rank()),
		size_(grid.mpihelper().size()),
		findDB_(findDomainBoundary),
		volumeTag_(volumeTag),
		surfaceTag_(surfaceTag)
	{

		// Unless the user wishes to obtain a domain boundary container, user must specify the tags associated with the boundary
		assert(findDB_ || ((volumeTag_ > 0) && (surfaceTag_ > 0)));

		{
			std::stringstream logstr;
			logstr << "Initializing GlobalBoundaryContainer with ";
			logstr << "volumeTag=" << volumeTag_ << ", ";
			logstr << "surfaceTag=" << surfaceTag_ << ", ";
			logstr << "findDomainBoundary=" << findDomainBoundary;
			LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, logstr.str());
		}

		// MPI_COMM Data
		std::vector<UInt> elemOrder;
		std::vector<UInt> nBS;
		std::vector<UInt> bsSubInd;
		std::vector<UInt> nPoint;
		std::vector<GlobalCoordinate> point;
		std::vector<UInt> gInd[4];

		// [TODO] Hard-coded for Tetrahedron - can extend to other geometry types if necessary
		int NELEM = 1;
		int NFACE = 4;
		int NEDGE = 6;
		int NCORNER = 4;



		/****************************************************
		 *   Stage 1: Construct Local Domain Boundary Data
		 ****************************************************/

  		/** \brief Iterate ove all elements of Interior Border partition */
		UInt elemCount = 0;
  		UInt nElementInterior = grid_.numInternal(ELEMENT_CODIM);
  		LocalCoordinate centerFaceLocal = ReferenceElements2d::simplex().position( 0, 0 );

		LeafGridView leafView = grid_.leafGridView();
		for (auto&& elemThis : elements(leafView, Dune::Partitions::interiorBorder))
		{
			LoggingMessage::writePatience("Preparing domain boundary for global communication...", elemCount++, nElementInterior);
			PhysicalTagType elementThisTag = grid_.template entityPhysicalTag<ELEMENT_CODIM>(elemThis);

			int nBSThis = 0;

			for (auto&& interThis : intersections(leafView, elemThis))
			{
				// Obtain the entity associated with this intersection
				unsigned int intrBraIndexInInside = interThis.indexInInside();
				EntityFace faceBra = elemThis.template subEntity<FACE_CODIM>(intrBraIndexInInside);
				PhysicalTagType boundaryTagBra = grid_.template entityPhysicalTag<FACE_CODIM>(faceBra);

				// Determine if the intersection corresponds to the boundary in question
				bool isBoundary = findDB_
						? (interThis.boundary() == true) && (interThis.neighbor() == false)
						:  (elementThisTag == volumeTag_) && (boundaryTagBra == surfaceTag_);

				if (isBoundary)  {
					nBSThis++;
					bsSubInd.push_back(interThis.indexInInside());
				}
			}

			// If this element contains any boundary segments at all, fill in other information and add it to the set
			if (nBSThis > 0) {
				elemOrder.push_back(grid_.template entityInterpolationOrder<ELEMENT_CODIM>(elemThis));

				nBS.push_back(nBSThis);

				BaseGeometryElement elemGeomBase = grid_.template entityBaseGeometry<ELEMENT_CODIM>(elemThis);
				std::vector<GlobalCoordinate> pThis = elemGeomBase.vertexSet();
				nPoint.push_back(pThis.size());
				for (int ip = 0; ip < pThis.size(); ip++) { point.push_back(pThis[ip]); }

				gInd[ELEMENT_CODIM].push_back(grid_.template entityGlobalIndex<ELEMENT_CODIM>(elemThis));
				for (int ip = 0; ip < NFACE; ip++)  { gInd[FACE_CODIM].push_back(grid_.template subentityGlobalIndex<ELEMENT_CODIM, FACE_CODIM>(elemThis, ip)); }
				for (int ip = 0; ip < NEDGE; ip++)  { gInd[EDGE_CODIM].push_back(grid_.template subentityGlobalIndex<ELEMENT_CODIM, EDGE_CODIM>(elemThis, ip)); }
				for (int ip = 0; ip < NCORNER; ip++)  { gInd[VERTEX_CODIM].push_back(grid_.template subentityGlobalIndex<ELEMENT_CODIM, VERTEX_CODIM>(elemThis, ip)); }



				/* Test for Grid self-consistency between globalIndex of subentity and subentityGlobalIndex
				bool BBS = false;  // Bloody BullShit :D
				for (UInt iCorner = 0; iCorner < nCorner; iCorner++)  {
					UInt indexCornerInElem = ReferenceElements3d::general(elemThis.type()).subEntity(intrIndexInInside, FACE_CODIM, iCorner, VERTEX_CODIM);
					EntityVertex cornerThis = elemThis.template subEntity<VERTEX_CODIM>(indexCornerInElem);

					int gindcornerv1 = grid_.template entityGlobalIndex<VERTEX_CODIM>(cornerThis);
					int gindcornerv2 = grid_.template subentityGlobalIndex<FACE_CODIM, VERTEX_CODIM>(faceThis, iCorner);
					int gindcornerv3 = grid_.template subentityGlobalIndex<ELEMENT_CODIM, VERTEX_CODIM>(elemThis, indexCornerInElem);

					// All different ways of obtaining the global index of a corner should give an equivalent result
					//std::cout << "Test: " << boundaryIndex << " " << iCorner << " " << indexCornerInElem << " " << gindcornerv1 << " " << gindcornerv2 << " " << gindcornerv3 << std::endl;
					if(gindcornerv1 != gindcornerv2) { BBS = true; }
					if(gindcornerv2 != gindcornerv3) { BBS = true; }
				}
				assert(!BBS);*/
			}
		}

		std::vector<UInt> elemOrderTotal;
		std::vector<UInt> nBSTotal;
		std::vector<UInt> bsSubIndTotal;
		std::vector<UInt> nPointTotal;
		std::vector<GlobalCoordinate> pointTotal;
		std::vector<UInt> gIndTotal[4];

		std::vector<int> commSize;

		/*
		std::cout << rank_ << " BeforeComm: "
				<< elemOrder.size() << " "
				<< nBS.size() << " "
				<< bsSubInd.size() << " "
				<< nPoint.size() << " "
				<< point.size() << " "
				<< gInd[0].size() << " "
				<< gInd[1].size() << " "
				<< gInd[2].size() << " "
				<< gInd[3].size() << std::endl;
		*/

		LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "Started Global Boundary Communication");

		mpiAllGatherVwrapper(nBS, nBSTotal, commSize);
		mpiAllGatherVwrapper(bsSubInd, bsSubIndTotal, commSize);
		mpiAllGatherVwrapper(nPoint, nPointTotal, commSize);
		mpiAllGatherVwrapper(point, pointTotal, commSize);
		mpiAllGatherVwrapper(gInd[0], gIndTotal[0], commSize);
		mpiAllGatherVwrapper(gInd[1], gIndTotal[1], commSize);
		mpiAllGatherVwrapper(gInd[2], gIndTotal[2], commSize);
		mpiAllGatherVwrapper(gInd[3], gIndTotal[3], commSize);

		// Have order last, so that we know from commSize, how many boundary elements there are on each process
		mpiAllGatherVwrapper(elemOrder, elemOrderTotal, commSize);

		LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "Finished Global Boundary Communication");

		/*
		std::cout << rank_  << " AfterComm: "
				<< elemOrderTotal.size() << " "
				<< nBSTotal.size() << " "
				<< bsSubIndTotal.size() << " "
				<< nPointTotal.size() << " "
				<< pointTotal.size() << " "
				<< gIndTotal[0].size() << " "
				<< gIndTotal[1].size() << " "
				<< gIndTotal[2].size() << " "
				<< gIndTotal[3].size() << std::endl;
		*/

		// Test if communication produced arrays of expected size
		assert(commSize.size() == size_);
		UInt nFaceCommTotal = nBSTotal.size();
		assert(nPointTotal.size() == nFaceCommTotal);
		assert(gIndTotal[0].size() == NELEM * nFaceCommTotal);
		assert(gIndTotal[1].size() == NFACE * nFaceCommTotal);
		assert(gIndTotal[2].size() == NEDGE * nFaceCommTotal);
		assert(gIndTotal[3].size() == NCORNER * nFaceCommTotal);

		/*********************************************************************************
		 *  4) Transfer the communicated arrays to a useable storage
		 *  Note: Explicitly remove all boundary segments already stored on this process
		 *********************************************************************************/

		int iElem = 0;
		int iRank = 0;
		int iSubInt = 0;
		int iPoint = 0;
		std::vector<int> igInd {0, 0, 0, 0};

	    UInt order_;															// Curvilinear order of the entity
	    std::vector<UInt> subBoundarySegInd_;			// Subentity boundary segment index set
	    std::vector<GlobalCoordinate> p_;					// Vertices of the boundary face in the original order
		std::vector<std::vector<UInt> > subGInd_;		// Global index of all subentities of this element by codim, in correct subentity order

		while (iRank < size_) {

			//std::cout << rank_ << " IterTest: " << iRank << " " << iElem << " " << iSubInt << " " << iPoint << " " << igInd[0] << " " << igInd[1] << " " << igInd[2] << " " << igInd[3] << std::endl;

			BoundaryContainer thisCont;							// Make boundary container
			thisCont.order_ = elemOrderTotal[iElem];		// Insert order

			//std::cout << "aaa" << std::endl;

			// Insert subInt
			assert((nBSTotal[iElem] > 0) && (nBSTotal[iElem] < 4));
			for (int i = iSubInt; i < iSubInt + nBSTotal[iElem]; i++)  { thisCont.subBoundarySegInd_.push_back(bsSubIndTotal[i]); }
			iSubInt += nBSTotal[iElem];

			//std::cout << "bbb" << std::endl;

			// Insert points
			assert(iPoint + nPointTotal[iElem] <= pointTotal.size());
			for (int i = iPoint; i < iPoint + nPointTotal[iElem]; i++)  { thisCont.p_.push_back(pointTotal[i]); }
			iPoint += nPointTotal[iElem];

			//std::cout << "ccc" << std::endl;

			// insert global indices
			thisCont.subGInd_ = std::vector<std::vector<UInt> > (4);
			for (int i = igInd[0]; i < igInd[0] + NELEM; i++)  { thisCont.subGInd_[0].push_back(gIndTotal[0][i]); }
			for (int i = igInd[1]; i < igInd[1] + NFACE; i++)  { thisCont.subGInd_[1].push_back(gIndTotal[1][i]); }
			for (int i = igInd[2]; i < igInd[2] + NEDGE; i++)  { thisCont.subGInd_[2].push_back(gIndTotal[2][i]); }
			for (int i = igInd[3]; i < igInd[3] + NCORNER; i++)  { thisCont.subGInd_[3].push_back(gIndTotal[3][i]); }
			igInd[0] += NELEM;
			igInd[1] += NFACE;
			igInd[2] += NEDGE;
			igInd[3] += NCORNER;

			//std::cout << "ddd" << std::endl;

			// Ignore entities on this process, they are already located in the local grid, save space
			iElem++;
			if (iRank != rank_) { boundaryContainerGlobal_.push_back(thisCont); }

			//std::cout << "eee" << std::endl;

			// If traversed all elements sent by this rank, jump to next one
			// Skip ranks that have communicated nothing
			while ((iRank < size_) && (iElem == commSize[iRank])) {
				iElem = 0;
				iRank++;
			}
		}

		//std::cout << "rank=" << rank_ << " Finished Loop" << std::endl;
	}


	ContainerVectorConstIter begin()  const { return boundaryContainerGlobal_.begin(); }
	ContainerVectorConstIter end()  const { return boundaryContainerGlobal_.end(); }


	template<class DataType>
	void mpiAllGatherVwrapper(std::vector <DataType> & src, std::vector<DataType> & rez, std::vector<int> & nCommRecv) {

		/********************************************************************
		 *   Stage 2: MPI communicate boundary data to all other processes
		 ********************************************************************/
		MPI_Comm comm = Dune::MPIHelper::getCommunicator();
		int mpi_err;

		// Algorithm

		/*********************************************************
		 *  1) Preliminary steps - compute sizes for communication
		 *********************************************************/
		int structSize = sizeof(DataType);	// the size of the communiated structure
		int nCommThis = src.size();		// the number of structures to be communicated
		int sizeCommThis = nCommThis * structSize;			// the total communication size in bytes

		/*********************************************************************************
		 *  2) Send number of structures this proc will communicate to all other processes
		 *********************************************************************************/

		LoggingMessage::template write<CurvGrid::LOG_MSG_PRODUCTION>( __FILE__, __LINE__, "--GlobalBoundaryContainer: Communicating Domain Size. ThisSize=" + std::to_string(nCommThis));

		nCommRecv = std::vector<int>(size_, 0);		// Number of elements that are sent by each of the processes
		std::vector<int> sizeCommRecv(size_, 0);		// Total size in bytes sent by each process
		std::vector<int> displCommRecv(size_, 0);		// Displacement in bytes
		mpi_err = MPI_Allgather(&nCommThis, 1, MPI_INT, reinterpret_cast<int *>(nCommRecv.data()), 1, MPI_INT, comm);
		int nCommRecvTot = std::accumulate(nCommRecv.begin(), nCommRecv.end(), 0);		// The total number of structrures that will be received by any one process
		for (int i = 0; i < size_; i++) { sizeCommRecv[i] = nCommRecv[i] * structSize; }
		for (int i = 1; i < size_; i++) { displCommRecv[i] = displCommRecv[i-1] + sizeCommRecv[i-1]; }

		//std::cout << "nCommThis=" << nCommThis << std::endl;
		//std::cout << "nCommRecv=" << VectorHelper::vector2string(nCommRecv) << std::endl;

		/*********************************************************************************
		 *  3) MPI_COMM local struct array to all faces except this one
		 *********************************************************************************/

		LoggingMessage::template write<CurvGrid::LOG_MSG_PRODUCTION>( __FILE__, __LINE__, "--GlobalBoundaryContainer: Communicating Domain. nComm=" + std::to_string(nCommRecvTot) + ", single structure size=" + std::to_string(structSize));

		// Reserve memory in boundary container for structures that will be received
		rez = std::vector<DataType>(nCommRecvTot);

		mpi_err = MPI_Allgatherv(
			reinterpret_cast<DataType *>(src.data()),
			sizeCommThis,
			MPI_BYTE,
			reinterpret_cast<BoundaryContainer *>(rez.data()),
			reinterpret_cast<int *>(sizeCommRecv.data()),
			reinterpret_cast<int *>(displCommRecv.data()),
			MPI_BYTE,
			comm
		);
	}



private:
	const Grid & grid_;
	int rank_;
	int size_;

	bool findDB_;
	PhysicalTagType volumeTag_;
	PhysicalTagType surfaceTag_;

	std::vector<BoundaryContainer> boundaryContainerGlobal_;


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

	typedef Dune::ReferenceElement<CoordinateType, dimension>  ReferenceElement3D;
	typedef Dune::ReferenceElement<CoordinateType, dimension-1>  ReferenceElement2D;
	typedef Dune::ReferenceElements<CoordinateType, dimension> ReferenceElements3D;
	typedef Dune::ReferenceElements<CoordinateType, dimension-1> ReferenceElements2D;

	typedef GlobalBoundaryIterator<Grid> This;
	typedef GlobalBoundaryContainer<Grid> BaseContainer;
	typedef typename BaseContainer::ContainerVectorConstIter  BaseConstIter;
	typedef typename BaseContainer::BoundaryContainer SegmentContainer;

	typedef typename Grid::template Codim<ELEMENT_CODIM>::EntityGeometryMappingImpl  BaseGeometryElement;
	typedef typename Grid::template Codim<FACE_CODIM>::EntityGeometryMappingImpl  BaseGeometryFace;
	typedef typename Grid::template Codim<EDGE_CODIM>::EntityGeometryMappingImpl  BaseGeometryEdge;

	typedef typename BaseGeometryFace::LocalCoordinate   LocalCoordinate2d;
	typedef typename BaseGeometryElement::LocalCoordinate   LocalCoordinate3d;
	typedef typename BaseGeometryElement::GlobalCoordinate  GlobalCoordinate;

	static const int NUMBER_VERTEX = SegmentContainer::NUMBER_VERTEX;

public:

	GlobalBoundaryIterator(const BaseContainer & baseContainer) :
		baseContainer_(baseContainer),
		baseiter_(baseContainer.begin()),
		subIndexIter_(0),
		thisElemGeom_(nullptr)
	{
		init();
	}

	~GlobalBoundaryIterator()
	{
		if (thisElemGeom_ != nullptr) { delete thisElemGeom_; }
	}

	/** \brief Extract GlobalIndex of a subentity of the boundary segment
	 *  Allowed codimensions:  (1 : Face itself), (2 : SubEdge), (3 : SubCorner) */

	// [TODO] Currently the storage allows for all global indices of element and its subentities
	// However we only ever require the global indices of the boundary faces of the element and their subentities
	// Can optimize by storing and communicating less stuff
	template <int codim>
	UInt globalIndex(UInt subIndexInFace) const {
		Dune::GeometryType gt2d;   gt2d.makeTriangle();
		Dune::GeometryType gt3d;   gt3d.makeTetrahedron();
		assert((codim > 0)&&(codim <= dimension));
		assert(subIndexInFace < ReferenceElements2D::general(gt2d).size(codim - FACE_CODIM));  // Do not allow subentity index out of bounds
		UInt subIndexInElem = ReferenceElements3D::general(gt3d).subEntity(indexInInside(), FACE_CODIM, subIndexInFace, codim);

		return globalIndexInParent<codim>(subIndexInElem);
	}

	template <int codim>
	UInt globalIndexInParent(UInt subIndexInElem) const {
		assert(subIndexInElem < baseiter_->subGInd_[codim].size() );
		return baseiter_->subGInd_[codim][subIndexInElem];
	}

	/** \brief Extract curvilinear order of the boundary segment */
	UInt order() const { return baseiter_ -> order_; }


	GlobalCoordinate unitOuterNormal(const LocalCoordinate2d & x2d) const {
		LocalCoordinate3d x3d = thisElemGeom_->coordinateInParent(x2d, indexInInside());
		return thisElemGeom_->subentityUnitNormal(indexInInside(), x3d, nullptr);
	}


	UInt indexInInside() const { return baseiter_->subBoundarySegInd_[subIndexIter_]; }


	BaseGeometryElement geometryInInside() const { return *thisElemGeom_; }

	/** \brief Extract geometry of the associated face */
	BaseGeometryFace geometry() const {
		return thisElemGeom_->template subentityGeometry<dimension - 1>(baseiter_->subBoundarySegInd_[subIndexIter_]);
	}

	/** \brief Extract geometry of a subentity edge of the associated face, given by the subentity index */
	BaseGeometryEdge geometryEdge(UInt subIndex) const {
		Dune::GeometryType gt;   gt.makeTetrahedron();
		UInt edgeSubIndex =  ReferenceElements3D::general(gt).subEntity(indexInInside(), FACE_CODIM, subIndex, EDGE_CODIM);
		return thisElemGeom_->template subentityGeometry<dimension - 2>(edgeSubIndex);
	}

	/** \brief Increase this iterator */
	This & operator ++() {
		if (subIndexIter_ + 1 < baseiter_->subBoundarySegInd_.size())  { subIndexIter_++; }
		else {
			subIndexIter_ = 0;
			baseiter_++;
			init();
		}
		return *this;
	}

	/** \brief Dereference the iterator */
	const SegmentContainer & operator *() const { return *baseiter_; }

	/** \brief Check if this iterator has reached its end */
	bool end() const { return (baseiter_ == baseContainer_.end()); }

private:

	void init() {
		if (!end()) {
			if (thisElemGeom_ != nullptr) { delete thisElemGeom_; }
			thisElemGeom_ = new BaseGeometryElement(Dune::GeometryType(Dune::GeometryType::simplex, dimension), baseiter_->p_, baseiter_->order_);
		}
	}


private:

	const BaseContainer & baseContainer_;
	BaseConstIter baseiter_;
	UInt subIndexIter_;

	BaseGeometryElement * thisElemGeom_;

};











} // namespace CurvGrid

} // namespace Dune

#endif  // DUNE_GLOBALBOUNDARYCONTAINER_HH
