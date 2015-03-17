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
    std::string filename = CURVILINEARGRID_TEST_GRID_PATH + GMSH_FILE_NAME[interpOrder - 1];

    // Additional constants
    bool insertBoundarySegment = true;  // If boundary segments will be inserted from GMSH. At the moment MUST BE true
    bool withGhostElements = true;      // to create Ghost elements
    bool verbose = true;                // to write logging output on master process
    bool processVerbose = true;         // to write logging output on all processes
    bool writeReaderVTKFile = false;    // to write mesh to VTK during reading stage

    // Construct the grid factory
    Dune::CurvilinearGridFactory<GridType> factory(withGhostElements, verbose, processVerbose, mpihelper);

    // Factory requires total vertex and element number in the mesh for faster performance
    int nVertexTotal;
    int nElementTotal;

    // Read the mesh into the factory using Curvilinear GMSH Reader
    Dune::CurvilinearGmshReader< GridType >::read(factory,
                                                  filename,
                                                  mpihelper,
                                                  verbose,
                                                  processVerbose,
                                                  writeReaderVTKFile,
                                                  insertBoundarySegment);

    // Create the grid
    return factory.createGrid();
}


/* A DataHandle class that communicates a fixed constant for all entities participating in the communication */
template<class IndexSet, class ValueType> // mapper type and vector type
class DataHandleConstant
  : public Dune::CommDataHandleIF<DataHandleConstant<IndexSet,ValueType>, ValueType>
{
  typedef typename  IndexSet::IndexType             IndexType;
  typedef typename  std::map<IndexType, ValueType>  DataMap;
  typedef typename  DataMap::iterator               DataIter;

public:

  //! constructor
  DataHandleConstant (
    const IndexSet& indexset,
    DataMap & in,
    DataMap & out,
    int codim)
      : indexset_(indexset),
        in_(in),
        out_(out),
        codim_(codim)
  {
	std::vector<Dune::GeometryType>  gtVec = indexset_.types(codim_);
	assert(gtVec.size() == 1);

	std::cout << "expecting to communicate " << in_.size() << " of " << indexset_.size(gtVec[0]) << " entities" << std::endl;
  }

  //! returns true if data for this codim should be communicated
  bool contains (int dim, int codim) const
  {
    return haveDim(dim)&&(codim==codim_);
  }

  //! returns true if size per entity of given dim and codim is a constant
  bool fixedsize (int dim, int codim) const
  {
    return haveDim(dim) ;
  }

  /*! how many objects of type DataType have to be sent for a given entity
   *
   *  Note: Only the sender side needs to know this size. */
  template<class EntityType>
  size_t size (EntityType& e) const
  {
    return 1 ;
  }

  /*! pack data from user to message buffer */
  template<class MessageBuffer, class EntityType>
  void gather (MessageBuffer& buff, const EntityType& e) const
  {
	IndexType ind = indexset_.index(e);
	DataIter iter = in_.find(ind);

    if (iter == in_.end())
    {
    	std::cout << " Error: DataHandleConstant: sending data from unexpected entity index=" << ind << std::endl;
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
  const IndexSet& indexset_;
  DataMap& in_;
  DataMap& out_;
  int codim_;
};


// [FIXME] Check if the entity has neighbour of that comm type before adding it
template<class IndexType, class GridType>
bool checkSend(int codim, IndexType ind, Dune::InterfaceType iftype, Dune::CommunicationDirection dir, GridType & grid)
{
  typedef typename GridType::GridStorageType  GridStorageType;
  typedef typename GridType::GridBaseType     GridBaseType;
  const GridBaseType & gridbase = grid.gridbase();

  typedef typename GridStorageType::StructuralType  StructuralType;
  const int InternalType        = GridStorageType::PartitionType::Internal;
  const int ProcessBoundaryType = GridStorageType::PartitionType::ProcessBoundary;
  const int DomainBoundaryType  = GridStorageType::PartitionType::DomainBoundary;
  const int GhostType           = GridStorageType::PartitionType::Ghost;

  StructuralType entityType = gridbase.entityStructuralType(codim, ind);

  switch (iftype)
  {
  case Dune::InterfaceType::All_All_Interface  :
	return (entityType == ProcessBoundaryType) || (entityType == InternalType) || (entityType == GhostType);  break;
  case Dune::InterfaceType::InteriorBorder_InteriorBorder_Interface  :  return entityType == ProcessBoundaryType;  break;
  case Dune::InterfaceType::InteriorBorder_All_Interface  :
  {
	if (dir == Dune::CommunicationDirection::ForwardCommunication)  { return (entityType == InternalType) || (entityType == ProcessBoundaryType); }
	else                                                            { return (entityType == ProcessBoundaryType) || (entityType == GhostType); }
  } break;
  default :  DUNE_THROW( Dune::IOError, "Attempt to use illegal comm protocol" ); break;
  }
}


template<class IndexType, class GridType>
bool checkRecv(int codim, IndexType ind, Dune::InterfaceType iftype, Dune::CommunicationDirection dir, GridType & grid)
{
	if (dir == Dune::CommunicationDirection::ForwardCommunication)   { return checkSend(codim, ind, iftype, Dune::CommunicationDirection::BackwardCommunication, grid); }
	if (dir == Dune::CommunicationDirection::BackwardCommunication)  { return checkSend(codim, ind, iftype, Dune::CommunicationDirection::ForwardCommunication, grid); }
}


/* This class runs simple grid communication of selected codim and interface
 * It then checks that exactly the entities that were supposed to participate in the communication participated in it
 *   */
template<class GridType, int codim>
void communicateConst(Dune::InterfaceType iftype, Dune::CommunicationDirection dir, Dune::MPIHelper & mpihelper, GridType & grid)
{
  if (mpihelper.size() <= 1)  { std::cout << "--skipping codim " << codim << " communication test because there is only 1 process" << std::endl; }
  else
  {
	  std::cout << "Started const-communication example for codim=" << codim << std::endl;

	  typedef typename GridType::LeafGridView   LeafGridView;
	  typedef typename GridType::LeafIndexSet   LeafIndexSet;

	  typedef typename LeafIndexSet::IndexType  IndexType;
	  typedef int                               ValueType;
	  typedef typename std::map<IndexType, ValueType>    DataMap;
	  typedef typename DataMap::iterator                 DataIter;

	  const LeafIndexSet & indset = grid.leafIndexSet();
	  LeafGridView leafView = grid.leafGridView();

	  DataMap in;
	  DataMap out;
	  int TMP = 5;

	  // Iterate over entities of this codimension
	  typedef typename LeafGridView::template Codim<codim>::Iterator EntityLeafIterator;
	  EntityLeafIterator ibegin = leafView.template begin<codim>();
	  EntityLeafIterator iend   = leafView.template end<codim>();


	  std::cout << " -- creating send-arrays" << std::endl;
	  for (EntityLeafIterator it = ibegin; it != iend; ++it)
	  {
		// Get entities that participate in comm - thus assemble in and out maps
		IndexType ind = indset.index(*it);
		//std::cout << "--- checking entity index=" << ind << " of type " << grid.gridbase().entityStructuralType(codim, ind) << std::endl;
		if (checkSend(codim, ind, iftype, dir, grid))  { in[ind] = TMP; }// std::cout << "-----mark send" << std::endl; }
		if (checkRecv(codim, ind, iftype, dir, grid))  { out[ind] = 0;  }// std::cout << "-----mark recv" << std::endl; }
	  }
	  int recvSize = out.size();

	  // Create datahandle
	  std::cout << " -- Initialising datahandle" << std::endl;
	  typedef DataHandleConstant<LeafIndexSet, ValueType>        DHConstImpl;
	  typedef typename Dune::CommDataHandleIF< DHConstImpl, ValueType >   DHConst;
	  DHConstImpl dhimpl(indset, in, out, codim);
	  //DHConst dh(DHConstImpl(indset, in, out, codim));

	  // Perform communication
	  std::cout << " -- started communicate" << std::endl;
	  //grid.communicate (dh, iftype, dir);
	  grid.communicate (dhimpl, iftype, dir);

	  // assert all values of out array to be what they need to be
	  assert(out.size() == recvSize);

	  std::cout << "after communication: " << std::endl;
	  for (DataIter datait = out.begin(); datait != out.end(); datait++)
	  {
		  std::cout << "     " << (*datait).first << " " << (*datait).second << std::endl;
	  }


	  for (DataIter datait = out.begin(); datait != out.end(); datait++)
	  {
		  assert((*datait).second == TMP);
	  }
  }
}






int main (int argc , char **argv) {
	static Dune::MPIHelper & mpihelper = Dune::MPIHelper::instance(argc, argv);

	// Define curvilinear grid
	const int dim = 3;
	const int dimworld = 3;
	typedef  double    ctype;
	typedef Dune::CurvilinearGrid<dim, dimworld, ctype, isCached> GridType;


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
		communicateConst<GridType, 0>(interfType[i], interfDir[i], mpihelper, *grid);  // Elements
		communicateConst<GridType, 1>(interfType[i], interfDir[i], mpihelper, *grid);  // Faces
		communicateConst<GridType, 2>(interfType[i], interfDir[i], mpihelper, *grid);  // Edges
		communicateConst<GridType, 3>(interfType[i], interfDir[i], mpihelper, *grid);  // Vertices
	}


    // Delete the grid
    delete grid;


    return 0;
}
