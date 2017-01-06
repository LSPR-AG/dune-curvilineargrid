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
		std::cout << "the expected send-indices are " << VectorHelper::map2string(in_) << std::endl;
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
			std::cout << " Error: DataHandleConstant: sending data from unexpected entity index=" << ind;

			/*
			std::cout << " of structural type " << Dune::PartitionName(e.partitionType());
			int parentIndex = grid_.gridbase().edgeNeighbor(ind);
			std::cout << " of parent index " << parentIndex << " of structural type " << grid_.gridbase().entityStructuralType(0, parentIndex);

			const auto & storage = grid_.gridbase().gridstorage();
			for (int i = 0; i < storage.elementSubentityCodim1_[parentIndex].size(); i++)
			{
				std::cout << " has neighbor face " << storage.elementSubentityCodim1_[parentIndex][i] << " of type " << grid_.gridbase().entityStructuralType(0, parentIndex);
			}
			std::cout << std::endl;
			*/


			DUNE_THROW( Dune::GridError, " DataHandleConst: sending data from entity that is not supposed to communicate" );
		}

		buff.write((*iter).second);
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
			std::cout << " Error: DataHandleConstant: received data from entity index=" << ind << " which was not found" << std::endl;
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
			int TMP = 5;

			// Iterate over entities of this codimension
			bool withGhosts = grid.gridbase().withGhostElements();

			//std::cout << " -- creating send-arrays" << std::endl;
			for (auto&& elementThis : elements(leafView, Dune::Partitions::interiorBorder)) {
				bool marked_internal = false;

				for (auto&& intersection : intersections(leafView, elementThis)) {
					// We are interested in marking process boundaries, and their neighbor entities
					if (isProcessBoundary(intersection, withGhosts)) {
						//std::cout << "process=" << mpihelper.rank() << " visiting face " << indset.subIndex(entity, intersection.indexInInside(), 1) << std::endl;

						const Face thisFace = elementThis.template subEntity<FACE_CODIM>(intersection.indexInInside());
						const Element & elementNeighbor = intersection.outside();

						mark_subentity<Face, codim>(thisFace, in, out, iftype, dir, grid, TMP);
						if (!marked_internal)  { mark_subentity<Element, codim>(elementThis, in, out, iftype, dir, grid, TMP); }
						if (withGhosts)        { mark_subentity<Element, codim>(elementNeighbor, in, out, iftype, dir, grid, TMP); }
						marked_internal = true;
					}
				}
			}
			int recvSize = out.size();

			// Create datahandle
			std::cout << " -- Initialising datahandle" << std::endl;
			DHConstImpl dhimpl(grid, indset, in, out, codim);
			//DHConst dh(DHConstImpl(indset, in, out, codim));

			// Perform communication
			std::cout << " -- started communicate" << std::endl;
			//grid.communicate (dh, iftype, dir);
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


			for (DataIter datait = out.begin(); datait != out.end(); datait++)
			{
				if ((*datait).second != TMP) { std::cout << "received unexpected data = " << (*datait).second << std::endl; }
				assert((*datait).second == TMP);
			}
		}
	}




private:

	// Warning, this procedure can not tell apart Ghosts and ProcessBoundaries
	static bool isProcessBoundary(const Intersection & intr, bool withGhosts)
	{
	    if (intr.boundary()) { return false; }
	    if (intr.neighbor()) {
	    	return intr.outside().partitionType() == Dune::PartitionType::GhostEntity;
	    } else {
	    	if (!withGhosts)  { return true; }
	    	else
	    	{
	        	std::cout << "Found non-boundary entity with no neighbour, given that a mesh has ghosts" << std::endl;
	        	DUNE_THROW( Dune::IOError, "Found non-boundary entity with no neighbour" );
	    	}
	    }
	}

	static bool usePartitionType(Dune::PartitionType ptype, Dune::InterfaceType iftype, Dune::CommunicationDirection dir, bool send)
	{
		switch (iftype)
		{
			case Dune::InterfaceType::InteriorBorder_InteriorBorder_Interface : return false;
			case Dune::InterfaceType::All_All_Interface :                       return true;
			case Dune::InterfaceType::InteriorBorder_All_Interface  :
			{
				if ((dir == Dune::CommunicationDirection::ForwardCommunication) == send)  { return ptype == Dune::PartitionType::InteriorEntity; }
				else                                                                      { return ptype != Dune::PartitionType::InteriorEntity; }
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

		// Add subentities from all PB faces
		bool useFace = ((entityCodim == 1)&&(codim >= 1));

		// Add subentities of entities if they are to be communicated
		bool useEntitySend = ((entityCodim == 0)&&(usePartitionType(entity.partitionType (), iftype, dir, true) ));
		bool useEntityRecv = ((entityCodim == 0)&&(usePartitionType(entity.partitionType (), iftype, dir, false) ));

		if (useFace || useEntitySend || useEntityRecv)
		{
			const LeafIndexSet & indset = grid.leafIndexSet();
			// Get index set

			Dune::GeometryType gt;  gt.makeSimplex(cdim);
			int nSub = ReferenceElements3D::general(gt).size(0, entityCodim, codim);

			for (int i = 0; i < nSub; i++) {
				 int ind = indset.subIndex(entity, i, codim);

				 if (useFace || useEntitySend)  { in[ind] = TMP;  }
				 if (useFace || useEntityRecv)  { out[ind] = 0;   }
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

	for (int i = 0; i < 4; i++)
	{
		TestClass::communicateConst<ELEMENT_CODIM>(interfType[i], interfDir[i], mpihelper, *grid);  // Elements
		TestClass::communicateConst<FACE_CODIM>(interfType[i], interfDir[i], mpihelper, *grid);  // Faces
		TestClass::communicateConst<EDGE_CODIM>(interfType[i], interfDir[i], mpihelper, *grid);  // Edges
		TestClass::communicateConst<VERTEX_CODIM>(interfType[i], interfDir[i], mpihelper, *grid);  // Vertices
	}


	typedef LoggingTimer<LoggingMessage>                 LoggingTimerDev;
	LoggingTimerDev::reportParallel();


    // Delete the grid
    delete grid;


    return 0;
}
