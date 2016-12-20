// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

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


/* A DataHandle class that communicates a fixed constant for all entities participating in the communication */
template<class GlobalIdSet, class IndexType> // mapper type and vector type
class DataHandleGlobalIndex
  : public Dune::CommDataHandleIF<DataHandleGlobalIndex<GlobalIdSet, IndexType>, typename GlobalIdSet::IdType>
{
  typedef typename  GlobalIdSet::IdType          IdType;
  typedef typename  std::map<IndexType, IdType>  GlobalIdMap;
  typedef typename  GlobalIdMap::iterator        DataIter;

public:

  //! constructor
  DataHandleGlobalIndex (
    const GlobalIdSet& idset,
    int codim)
      : idset_(idset),
        codim_(codim)
  {

  }

  //! returns true if data for this codim should be communicated
  bool contains (int dim, int codim) const
  {
    return haveDim(dim)&&(codim==codim_);
  }

  //! returns true if size per entity of given dim and codim is a constant
  bool fixedsize (int dim, int codim) const
  {
    return true;
  }

  /*! how many objects of type DataType have to be sent for a given entity
   *
   *  Note: Only the sender side needs to know this size. */
  template<class EntityType>
  size_t size (EntityType& e) const
  {
    return 1;
  }

  /*! Get globalId of this entity and write it to buffer */
  template<class MessageBuffer, class EntityType>
  void gather (MessageBuffer& buff, const EntityType& e) const
  {
  	buff.write(idset_.id(e));
  }

  /** \brief  Get globalId of this entity, Read globalId from buffer, Assert that the Id's are equal
   */
  template<class MessageBuffer, class EntityType>
  void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
  {
	IdType thisProcId = idset_.id(e);
	IdType thisRecvId;
	buff.read(thisRecvId);

    if (thisProcId != thisRecvId)
    {
    	std::cout << " Error: DataHandleGlobalId: received id=" << thisRecvId << " does not match true id=" << thisProcId << std::endl;
    	DUNE_THROW( Dune::GridError, " DataHandleConst: received data from entity not present on this process" );
    }
  }

protected:
  bool haveDim(int dim) const  { return dim == 3; }

private:
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
  else
  {
	  std::cout << "Started const-communication example for codim=" << codim << std::endl;

	  typedef typename GridType::LeafGridView   LeafGridView;
	  typedef typename GridType::GlobalIdSet    GlobalIdSet;
	  typedef typename LeafGridView::IntersectionIterator IntersectionIterator;
	  typedef typename IntersectionIterator :: Intersection Intersection;

	  typedef typename GlobalIdSet::IdType                     IdType;
	  typedef typename GridType::GridBaseType::LocalIndexType  LocalIndexType;

	  const GlobalIdSet & idset = grid.globalIdSet();
	  LeafGridView leafView = grid.leafGridView();

	  // Create datahandle
	  std::cout << " -- Initialising datahandle" << std::endl;
	  typedef DataHandleGlobalIndex<GlobalIdSet, LocalIndexType> DHGlobalIndexImpl;
	  DHGlobalIndexImpl dhimpl(idset, codim);

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

	// For each interface and codimension, perform simple communication test given by communicateConst
	for (int i = 0; i < 4; i++)
	{
		communicateGlobalIndex<GridType, 0>(interfType[i], interfDir[i], mpihelper, *grid);  // Elements
		communicateGlobalIndex<GridType, 1>(interfType[i], interfDir[i], mpihelper, *grid);  // Faces
		communicateGlobalIndex<GridType, 2>(interfType[i], interfDir[i], mpihelper, *grid);  // Edges
		communicateGlobalIndex<GridType, 3>(interfType[i], interfDir[i], mpihelper, *grid);  // Vertices
	}

	typedef LoggingTimer<LoggingMessage>                 LoggingTimerDev;
	LoggingTimerDev::reportParallel();

    // Delete the grid
    delete grid;


    return 0;
}
