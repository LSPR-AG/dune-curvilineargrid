/********************************************
 * Dune-CurvilinearGrid
 * Tutorial 5a: Communication: Constant Communication
 *
 * Author: Aleksejs Fomins
 *
 * Description: This tutorial integrates the vector unit outer normal over the domain boundary.
 * According to the Divergence theorem, this integral should amount to zero vector.
 ********************************************/

#include <config.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/curvilineargrid/curvilineargrid/grid.hh>
#include <dune/curvilineargrid/curvilineargridhowto/creategrid.hh>


using namespace Dune;
using namespace Dune::CurvGrid;

const bool isCached = true;
const int DIM0D = 0;   const int VERTEX_CODIM = 3;
const int DIM1D = 1;   const int EDGE_CODIM = 2;
const int DIM2D = 2;   const int FACE_CODIM = 1;
const int DIM3D = 3;   const int ELEMENT_CODIM = 0;



/* A DataHandle class that communicates a fixed constant for all entities participating in the communication */
template<class GridType, class IndexSet, class ValueType> // mapper type and vector type
class DataHandleConstant
  : public Dune::CommDataHandleIF<DataHandleConstant<GridType, IndexSet,ValueType>, ValueType>
{
	typedef typename  IndexSet::IndexType             IndexType;
	typedef std::map<IndexType, ValueType>  DataMap;
	typedef typename  DataMap::iterator               DataIter;
	typedef std::vector<Dune::GeometryType> TypeVector;

public:

	//! constructor
	DataHandleConstant (const GridType& grid, const IndexSet& indexset, DataMap & in, DataMap & out, int codim)
	  : grid_(grid), indexset_(indexset), in_(in), out_(out), codim_(codim)
	{
		TypeVector gtVec = indexset_.types(codim_);
		assert(gtVec.size() == 1);

		std::cout << "expecting to communicate " << in_.size() << " of " << indexset_.size(gtVec[0]) << " entities" << std::endl;
		//std::cout << "the expected send-indices are " << VectorHelper::map2string(in_) << std::endl;
	}

	//! returns true if data for this codim should be communicated
	bool contains (int dim, int codim) const  { return haveDim(dim)&&(codim==codim_); }

	//! returns true if size per entity of given dim and codim is a constant
	bool fixedsize (int dim, int codim) const { return haveDim(dim) ; }

	/*! how many objects of type DataType have to be sent for a given entity
	*  Note: Only the sender side needs to know this size. */
	template<class EntityType>
	size_t size (EntityType& e) const  { return 1 ; }

	/*! pack data from user to message buffer */
	template<class MessageBuffer, class EntityType>
	void gather (MessageBuffer& buff, const EntityType& e) const
	{
		IndexType ind = indexset_.index(e);
		DataIter iter = in_.find(ind);


		if (iter == in_.end())
		{
			std::cout << " Error: DataHandleConstant: sending data from unexpected entity index=" << ind
					<< " of ptype=" << Dune::PartitionName(e.partitionType()) << std::endl;

			DUNE_THROW( Dune::GridError, " DataHandleConst: sending data from entity that is not supposed to communicate" );
		}


		if (codim_ == FACE_CODIM) {
			bool isPB = e.partitionType() == Dune::PartitionType::BorderEntity;
			bool isPeriodic = (e.partitionType() == Dune::PartitionType::InteriorEntity) &&
					grid_.gridbase().faceBoundaryType(iter->second) != GridType::GridStorageType::FaceBoundaryType::PeriodicBoundary;

			if (!(isPB || isPeriodic))
			{
				std::cout << " Error: DataHandleConstant: sending data from unexpected entity index=" << ind
						<< " of ptype=" << Dune::PartitionName(e.partitionType()) << std::endl;

				DUNE_THROW( Dune::GridError, " DataHandleConst: sending data from entity that is not supposed to communicate" );
			}

			//if (isPeriodic) { std::cout << " sending from periodic face " << std::endl; }
		}

		buff.write(iter->second);
	}

	/** \brief Unpack data from message buffer to user
	*
	* \param n The number of objects sent by the sender
	*/
	template<class MessageBuffer, class EntityType>
	void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
	{
		IndexType ind = indexset_.index(e);
		DataIter iter = out_.find(ind);

		if (iter == out_.end())
		{
			//if ((e.partitionType() == Dune::PartitionType::GhostEntity) && gridbase_.)


			std::cout << " Error: DataHandleConstant: received data from entity index=" << ind
					<< " of ptype=" << e.partitionType()
					<< " which was not among the expected receiving entities" << std::endl;
			DUNE_THROW( Dune::GridError, " DataHandleConst: received data from entity not present on this process" );
		}

		ValueType x;
		buff.read(x);
		out_[ind] = x;
	}

protected:
	bool haveDim(int dim) const  { return dim == 3; }

private:
	const GridType& grid_;
	const IndexSet& indexset_;
	DataMap& in_;
	DataMap& out_;
	int codim_;
};







template <class GridType>
class ConstCommTest {

private:
	static const int cdim = GridType::dimension;

	typedef typename GridType::LeafGridView   LeafGridView;
	typedef typename GridType::LeafIndexSet   LeafIndexSet;
	typedef typename LeafIndexSet::IndexType  IndexType;

	typedef typename LeafGridView::template Codim< ELEMENT_CODIM >::Entity Element;
	typedef typename LeafGridView::template Codim< FACE_CODIM >::Entity Face;
	typedef typename GridType::Traits::LeafIntersection               Intersection;

	typedef int                               ValueType;
	typedef typename std::map<IndexType, ValueType>    DataMap;
	typedef typename DataMap::iterator                 DataIter;

	typedef DataHandleConstant<GridType, LeafIndexSet, ValueType>        DHConstImpl;
	typedef typename Dune::CommDataHandleIF< DHConstImpl, ValueType >    DHConst;

	typedef Dune::ReferenceElements<double, cdim>  ReferenceElements3D;


public:

	/* This class runs simple grid communication of selected codim and interface
	 * It then checks that exactly the entities that were supposed to participate in the communication participated in it
	 *
	 * NOTE: This procedure is much more complicated than necessary for regular communication and is excessive for user.
	 *   For a simpler example look at CommunicationGlobalIndex tutorial
	 *   */
	template<int codim>
	static void communicateConst(Dune::InterfaceType iftype, Dune::CommunicationDirection dir, Dune::MPIHelper & mpihelper, GridType & grid)
	{
		if (mpihelper.size() <= 1)  { std::cout << "--skipping codim " << codim << " communication test because there is only 1 process" << std::endl; }
		else
		{
			std::cout << "Started const-communication example for codim=" << codim << std::endl;

			const LeafIndexSet & indset = grid.leafIndexSet();
			LeafGridView leafView = grid.leafGridView();

			DataMap in;
			DataMap out;
			int TMP = 5; // Just some number to be communicated

			// Iterate over entities of this codimension
			bool withGhosts = grid.gridbase().withGhostElements();
			bool withPeriodic = grid.withPeriodic();

			//std::cout << " -- creating send-arrays" << std::endl;
			for (auto&& elementThis : elements(leafView, Dune::Partitions::interiorBorder)) {
				bool marked_internal = false;

				for (auto&& intersection : intersections(leafView, elementThis)) {

					bool isDB = intersection.boundary() && !intersection.neighbor();
					bool isPB = !intersection.boundary() && intersection.neighbor() && (intersection.outside().partitionType() == Dune::PartitionType::GhostEntity);
					bool isPeriodic = intersection.boundary() && intersection.neighbor();

					assert(!isPeriodic || withPeriodic);  // Can only find periodic if grid has periodic, otherwise something weird happens

					// We are interested in marking communicating boundaries, and their neighbor entities
					if (isPB || isPeriodic) {
						//std::cout << "process=" << mpihelper.rank() << " visiting face " << indset.subIndex(entity, intersection.indexInInside(), 1) << std::endl;

						const Face thisFace = elementThis.template subEntity<FACE_CODIM>(intersection.indexInInside());
						const Element & elementNeighbor = intersection.outside();

						// Mark subentities of the shared face for communication
						// This includes process and periodic boundary communication
						if (codim > ELEMENT_CODIM) {  // Element can't be a subentity of a face
							mark_subentity<Face, codim>(thisFace, in, out, iftype, dir, grid, TMP);
						}

						// Mark subentities of Ghosts and corresponding Interior elements, that would communicate if ghosts are present in the mesh
						if (withGhosts)        {
							if (!marked_internal)  { mark_subentity<Element, codim>(elementThis, in, out, iftype, dir, grid, TMP); }
							mark_subentity<Element, codim>(elementNeighbor, in, out, iftype, dir, grid, TMP);
							marked_internal = true;
						}
					}
				}
			}
			int recvSize = out.size();

			// Create datahandle
			std::cout << " -- Initialising datahandle" << std::endl;
			DHConstImpl dhimpl(grid, indset, in, out, codim);

			// Perform communication
			std::cout << " -- started communicate" << std::endl;
			grid.communicate (dhimpl, iftype, dir);

			// assert all values of out array to be what they need to be
			assert(out.size() == recvSize);

			/*
			std::cout << "after communication: " << std::endl;
			for (DataIter datait = out.begin(); datait != out.end(); datait++)
			{
			  std::cout << "     " << mpihelper.rank() << " " << (*datait).first << " " << (*datait).second << " type " << grid.gridbase().entityStructuralType(codim, (*datait).first) << std::endl;
			}
			*/


			for (const auto & datait : out)
			{
				if (datait.second != TMP) { DUNE_THROW(Dune::IOError, "received unexpected data = " + std::to_string(datait.second)); }
			}
		}
	}




private:

	static bool usePartitionType(Dune::PartitionType ptype, Dune::InterfaceType iftype, Dune::CommunicationDirection dir, bool send)
	{
		switch (iftype)
		{
			case Dune::InterfaceType::InteriorBorder_InteriorBorder_Interface : return false;
			case Dune::InterfaceType::All_All_Interface :                       return true;
			case Dune::InterfaceType::InteriorBorder_All_Interface  :
			{
				if ((dir == Dune::CommunicationDirection::ForwardCommunication) == send)  { return ptype == Dune::PartitionType::InteriorEntity; }
				else                                                                      { return ptype == Dune::PartitionType::GhostEntity; }
			} break;
		}
	}


	// [TODO] Construct subentity number without explicitly specifying the GeometryType
	template<class Entity, int codim>
	static void mark_subentity(
		const Entity & entity,
		DataMap & in,
		DataMap & out,
		Dune::InterfaceType iftype,
		Dune::CommunicationDirection dir,
		GridType & grid,
		int TMP)
	{
		const int entityCodim = Entity::codimension;

		assert(entityCodim <= codim);

		// Add subentities from all PB faces
		bool commSharedFace = entityCodim == FACE_CODIM;

		// Add subentities of entities if they are to be communicated
		bool commElementSend = ((entityCodim == ELEMENT_CODIM)&&(usePartitionType(entity.partitionType(), iftype, dir, true) ));
		bool commElementRecv = ((entityCodim == ELEMENT_CODIM)&&(usePartitionType(entity.partitionType(), iftype, dir, false) ));

		if (commSharedFace || commElementSend || commElementRecv)
		{
			const LeafIndexSet & indset = grid.leafIndexSet();

			// Find out how many subentities of the given codim are in this entity
			// [TODO] generalize for generic geometry types
			Dune::GeometryType gt;  gt.makeSimplex(cdim);
			int nSub = ReferenceElements3D::general(gt).size(0, entityCodim, codim);


			for (int i = 0; i < nSub; i++) {
				 int ind = indset.subIndex(entity, i, codim);

				 if (commSharedFace || commElementSend)  { in[ind] = TMP;  }	// Note that this index will send
				 if (commSharedFace || commElementRecv)  { out[ind] = 0;   }	// Note that this index will receive, and that it has not received yet
			}
		}
	}
};















int main (int argc , char **argv) {
	static Dune::MPIHelper & mpihelper = Dune::MPIHelper::instance(argc, argv);

	// Define curvilinear grid
	const int dim = 3;
	typedef  double    ctype;
	typedef Dune::CurvilinearGrid<ctype, dim, isCached> GridType;

	// Create Grid
	GridType * grid = createGrid<GridType>(mpihelper, argc, argv);

	// Construct all interfaces we are interested to test
	// That is, IB->IB, IB->ALL, ALL->IB and ALL->ALL
	Dune::InterfaceType interfType[4] {
		Dune::InterfaceType::InteriorBorder_InteriorBorder_Interface,
		Dune::InterfaceType::InteriorBorder_All_Interface,
		Dune::InterfaceType::InteriorBorder_All_Interface,
		Dune::InterfaceType::All_All_Interface
	};

	Dune::CommunicationDirection interfDir[4] {
		Dune::CommunicationDirection::ForwardCommunication,
		Dune::CommunicationDirection::ForwardCommunication,
		Dune::CommunicationDirection::BackwardCommunication,
		Dune::CommunicationDirection::ForwardCommunication
	};


	/*
	std::cout << " Brief self-test";
	std::cout << " " << usePartitionType(Dune::PartitionType::InteriorEntity, Dune::InterfaceType::InteriorBorder_All_Interface, Dune::CommunicationDirection::ForwardCommunication, true);
	std::cout << " " << usePartitionType(Dune::PartitionType::GhostEntity, Dune::InterfaceType::InteriorBorder_All_Interface, Dune::CommunicationDirection::ForwardCommunication, true);
	std::cout << " " << usePartitionType(Dune::PartitionType::InteriorEntity, Dune::InterfaceType::InteriorBorder_All_Interface, Dune::CommunicationDirection::ForwardCommunication, false);
	std::cout << " " << usePartitionType(Dune::PartitionType::GhostEntity, Dune::InterfaceType::InteriorBorder_All_Interface, Dune::CommunicationDirection::ForwardCommunication, false);
	std::cout << std::endl;
	*/


	// For each interface and codimension, perform simple communication test given by communicateConst
	typedef ConstCommTest<GridType> TestClass;

	std::vector<bool> doInterface = grid->withPeriodic()
			? std::vector<bool> {true, false, false, false}
			: std::vector<bool> {true, false, false, false};

	std::vector<bool> doCodim = grid->withPeriodic()
			? std::vector<bool> {true, true, false, false}
			: std::vector<bool> {true, true, true, true};

	// For each interface and codimension, perform simple communication test given by communicateConst
	for (int iInterface = 0; iInterface <= 3; iInterface++)
	{
		if (doInterface[iInterface]) {
			for (int iCodim = 0; iCodim <= 3; iCodim++)
			{
				if (doCodim[iCodim]) {
					if (iCodim == ELEMENT_CODIM) { TestClass::communicateConst<ELEMENT_CODIM>(interfType[iInterface], interfDir[iInterface], mpihelper, *grid); }
					if (iCodim == FACE_CODIM) { TestClass::communicateConst<FACE_CODIM>(interfType[iInterface], interfDir[iInterface], mpihelper, *grid); }
					if (iCodim == EDGE_CODIM) { TestClass::communicateConst<EDGE_CODIM>(interfType[iInterface], interfDir[iInterface], mpihelper, *grid); }
					if (iCodim == VERTEX_CODIM) { TestClass::communicateConst<VERTEX_CODIM>(interfType[iInterface], interfDir[iInterface], mpihelper, *grid); }
				} else {
					if (grid->comm().rank() == MPI_MASTER_RANK) {
						std::cout << "Skipping codim " << iCodim << " because not defined for periodic case" << std::endl;
					}
				}
			}
		} else {
			if (grid->comm().rank() == MPI_MASTER_RANK) {
				std::cout << "Skipping interface " << iInterface << " because not defined for periodic case" << std::endl;
			}
		}


	}



	typedef LoggingTimer<LoggingMessage>                 LoggingTimerDev;
	LoggingTimerDev::reportParallel();


    // Delete the grid
    delete grid;


    return 0;
}
