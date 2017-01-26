/********************************************
 * Dune-CurvilinearGrid
 * Tutorial 5b: Communication: Global Index Communication
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
template<class Grid> // mapper type and vector type
class DataHandleGlobalIndex
		: public Dune::CommDataHandleIF<DataHandleGlobalIndex<Grid>, typename Grid::GlobalIdSet::IdType>
{
	typedef typename Grid::LeafIndexSet   LeafIndexSet;
	typedef typename Grid::GlobalIdSet    GlobalIdSet;
	typedef typename Grid::LocalIndexType  LocalIndexType;
	typedef typename  GlobalIdSet::IdType          IdType;
//	typedef typename  std::map<IdType, IdType>  GlobalIdMap;
//	typedef typename  GlobalIdMap::iterator        DataIter;

public:

	//! constructor
	DataHandleGlobalIndex (Grid & grid, int codim)
		: grid_(grid), codim_(codim),
		  indexset_(grid.leafIndexSet()),
		  idset_(grid.globalIdSet())
	{  }


	//! returns true if data for this codim should be communicated
	bool contains (int dim, int codim) const  { return haveDim(dim)&&(codim==codim_); }

	//! returns true if size per entity of given dim and codim is a constant
	bool fixedsize (int dim, int codim) const  { return true; }

	/*! how many objects of type DataType have to be sent for a given entity
	*  Note: Only the sender side needs to know this size. */
	template<class EntityType>
	size_t size (EntityType& e) const  { return 1; }

	/*! Get globalId of this entity and write it to buffer */
	template<class MessageBuffer, class EntityType>
	void gather (MessageBuffer& buff, const EntityType& e) const { buff.write(idset_.id(e)); }

	/** \brief  Get globalId of this entity, Read globalId from buffer, Assert that the Id's are equal */
	template<class MessageBuffer, class EntityType>
	void scatter (MessageBuffer& buff, const EntityType& e, size_t n) {
		const int codim = EntityType::codimension;

		IdType thisProcId = idset_.id(e);
		IdType thisRecvId;
		buff.read(thisRecvId);

		bool directMatch = thisProcId == thisRecvId;
		bool periodicMatch = false;

		if (codim == FACE_CODIM) {
			auto boundaryType = grid_.template entityBoundaryType<codim>(e);

			if (boundaryType == Grid::PERIODIC_BOUNDARY_TYPE) {
				LocalIndexType periodicFaceLocalIndex = indexset_.index(e);
				LocalIndexType neighborFaceLocalIndex = grid_.gridbase().intersection().periodicNeighborFace(periodicFaceLocalIndex);
				auto neighborBoundaryType = grid_.gridbase().intersection().boundaryType(neighborFaceLocalIndex);

				if			(neighborBoundaryType == Grid::PERIODIC_BOUNDARY_TYPE)  { assert(directMatch); }
				else {
					IdType expectedId = grid_.gridbase().entity().globalId(codim, neighborFaceLocalIndex);

					 // Periodic boundary face indices should be different
					if ((directMatch)||(expectedId != thisRecvId)) {
						std::cout << " Error: DataHandleGlobalId: received periodic id="
								<< thisRecvId << " does not match the expected id="
								<< expectedId << " which was sent to the true id=" <<
								thisProcId << std::endl;
						DUNE_THROW( Dune::GridError, " DataHandleConst: received data from entity not present on this process" );
					}
				}
				periodicMatch = true;
			}
		}

		if (!(directMatch || periodicMatch))
		{
			std::cout << " Error: DataHandleGlobalId: received id=" << thisRecvId << " does not match true id=" << thisProcId << std::endl;
			DUNE_THROW( Dune::GridError, " DataHandleConst: received data from entity not present on this process" );
		}
	}

protected:
	bool haveDim(int dim) const  { return dim == 3; }

private:
	Grid & grid_;
	const LeafIndexSet & indexset_;
	const GlobalIdSet& idset_;
	int codim_;
};


/* This class runs simple grid communication of selected codim and interface
 * To each entity it sends that entity's globalId, which should be the same on the receiving and sending side
 * On the receiving size, we check if indeed the received and owned GlobalId match
 *   */
template<class GridType, int codim>
void communicateGlobalIndex(Dune::InterfaceType iftype, Dune::CommunicationDirection dir, Dune::MPIHelper & mpihelper, GridType & grid)
{
	if (mpihelper.size() <= 1)  { std::cout << "--skipping codim " << codim << " communication test because there is only 1 process" << std::endl; }
	else {
		std::cout << "Started const-communication example for codim=" << codim << std::endl;

		typedef typename GridType::LeafGridView   LeafGridView;
		typedef typename LeafGridView::IntersectionIterator IntersectionIterator;
		typedef typename IntersectionIterator :: Intersection Intersection;

		typedef typename GridType::GridBaseType    GridBaseType;
		typedef typename GridBaseType::LocalIndexType  LocalIndexType;
;
		LeafGridView leafView = grid.leafGridView();

		// Create datahandle
		std::cout << " -- Initialising datahandle" << std::endl;
		typedef DataHandleGlobalIndex<GridType> DHGlobalIndexImpl;
		DHGlobalIndexImpl dhimpl(grid, codim);

		// Perform communication
		std::cout << " -- started communicate" << std::endl;
		grid.communicate (dhimpl, iftype, dir);
	}
}



int main (int argc , char **argv) {
	static Dune::MPIHelper & mpihelper = Dune::MPIHelper::instance(argc, argv);

	// Define curvilinear grid
	const int dim = 3;
	typedef  double    ctype;
	typedef Dune::CurvilinearGrid<ctype, dim, isCached> GridType;
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

	std::vector<bool> doInterface = grid->withPeriodic()
			? std::vector<bool> {true, true, false, false}
			: std::vector<bool> {true, true, true, true};

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
					if (iCodim == ELEMENT_CODIM) { communicateGlobalIndex<GridType, ELEMENT_CODIM>(interfType[iInterface], interfDir[iInterface], mpihelper, *grid); }
					if (iCodim == FACE_CODIM) { communicateGlobalIndex<GridType, FACE_CODIM>(interfType[iInterface], interfDir[iInterface], mpihelper, *grid); }
					if (iCodim == EDGE_CODIM) { communicateGlobalIndex<GridType, EDGE_CODIM>(interfType[iInterface], interfDir[iInterface], mpihelper, *grid); }
					if (iCodim == VERTEX_CODIM) { communicateGlobalIndex<GridType, VERTEX_CODIM>(interfType[iInterface], interfDir[iInterface], mpihelper, *grid); }
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
