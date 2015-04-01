// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

#include <config.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/curvilineargrid/curvilineargrid/grid.hh>
#include <dune/curvilineargrid/curvilineargrid/factory.hh>

#include <dune/curvilineargrid/io/file/curvilineargmshreader.hh>


const bool isCached = true;


template <class GridType>
GridType * createGrid(Dune::MPIHelper & mpihelper)
{
	// Obtain path to the mesh folder from a CMake constant
    const std::string CURVILINEARGRID_TEST_GRID_PATH = std::string(DUNE_CURVILINEARGRID_EXAMPLE_GRIDS_PATH) + "curvilinear/";

    // A choice of 32 tetrahedra spheres with interpolation orders 1 to 5.
    const std::string GMSH_FILE_NAME[5] {"sphere32.msh", "sphere32ord2.msh", "sphere32ord3.msh", "sphere32ord4.msh", "sphere32ord5.msh"};

    // Choice of file name
    int interpOrder = 1;
    std::string filename = CURVILINEARGRID_TEST_GRID_PATH + "sphere2000ord3.msh";//GMSH_FILE_NAME[interpOrder - 1];

    // Additional constants
    bool insertBoundarySegment = true;  // If boundary segments will be inserted from GMSH. At the moment MUST BE true
    bool withGhostElements = true;      // to create Ghost elements
    bool writeReaderVTKFile = false;    // to write mesh to VTK during reading stage

    // Construct the grid factory
    Dune::CurvilinearGridFactory<GridType> factory(withGhostElements, mpihelper);

    // Factory requires total vertex and element number in the mesh for faster performance
    int nVertexTotal;
    int nElementTotal;

    // Read the mesh into the factory using Curvilinear GMSH Reader
    Dune::CurvilinearGmshReader< GridType >::read(factory,
                                                  filename,
                                                  mpihelper,
                                                  writeReaderVTKFile,
                                                  insertBoundarySegment);

    // Create the grid
    return factory.createGrid();
}


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

    // Instantiation of the logging message and loggingtimer
    typedef Dune::LoggingMessage<Dune::LoggingMessageHelper::Phase::DEVELOPMENT_PHASE>   LoggingMessageDev;
    typedef Dune::LoggingTimer<LoggingMessageDev>                                        LoggingTimerDev;
    LoggingMessageDev::getInstance().init(mpihelper, true, true);
    LoggingTimerDev::getInstance().init(false);

	typedef Dune::CurvilinearGrid<ctype, dim, isCached, LoggingMessageDev> GridType;


	// Create Grid
	GridType * grid = createGrid<GridType>(mpihelper);

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


    // Delete the grid
    delete grid;


    return 0;
}
