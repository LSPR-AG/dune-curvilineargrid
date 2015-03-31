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
    std::string filename = CURVILINEARGRID_TEST_GRID_PATH + "sphere2000ord3.msh"; // GMSH_FILE_NAME[interpOrder - 1]; //

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
template<class GridType, class IndexSet, class ValueType> // mapper type and vector type
class DataHandleConstant
  : public Dune::CommDataHandleIF<DataHandleConstant<GridType, IndexSet,ValueType>, ValueType>
{
  typedef typename  IndexSet::IndexType             IndexType;
  typedef typename  std::map<IndexType, ValueType>  DataMap;
  typedef typename  DataMap::iterator               DataIter;

public:

  //! constructor
  DataHandleConstant (
	const GridType& grid,
    const IndexSet& indexset,
    DataMap & in,
    DataMap & out,
    int codim)
      : grid_(grid),
    	indexset_(indexset),
        in_(in),
        out_(out),
        codim_(codim)
  {
	std::vector<Dune::GeometryType>  gtVec = indexset_.types(codim_);
	assert(gtVec.size() == 1);

	std::cout << "expecting to communicate " << in_.size() << " of " << indexset_.size(gtVec[0]) << " entities" << std::endl;
	std::cout << "the expected send-indices are " << Dune::VectorHelper::map2string(in_) << std::endl;
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
    	std::cout << " Error: DataHandleConstant: sending data from unexpected entity index=" << ind;

    	std::cout << " of structural type " << grid_.gridbase().entityStructuralType(EntityType::codimension, ind);
    	int parentIndex = grid_.gridbase().edgeNeighbor(ind);
    	std::cout << " of parent index " << parentIndex << " of structural type " << grid_.gridbase().entityStructuralType(0, parentIndex);

    	const auto & storage = grid_.gridbase().gridstorage();
    	for (int i = 0; i < storage.elementSubentityCodim1_[parentIndex].size(); i++)
    	{
    		std::cout << " has neighbor face " << storage.elementSubentityCodim1_[parentIndex][i] << " of type " << grid_.gridbase().entityStructuralType(0, parentIndex);
    	}
    	std::cout << std::endl;


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


// Warning, this procedure can not tell apart Ghosts and ProcessBoundaries
template <class Intersection>
bool isProcessBoundary(Intersection intr, bool withGhosts)
{
    if (intr.boundary()) { return false; }
    if (intr.neighbor())
    {
    	if (intr.outside().partitionType() == Dune::PartitionType::GhostEntity)  { return true; }
    	else                                                                     { return false; }
    } else {
    	if (!withGhosts)  { return true; }
    	else
    	{
        	std::cout << "Found non-boundary entity with no neighbour, given that a mesh has ghosts" << std::endl;
        	DUNE_THROW( Dune::IOError, "Found non-boundary entity with no neighbour" );
    	}
    }
}


bool usePartitionType(
  Dune::PartitionType ptype,
  Dune::InterfaceType iftype,
  Dune::CommunicationDirection dir,
  bool send)
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
template<class GridType, class Entity, class DataMap, int codim>
void mark_subentity(
  const Entity & entity,
  DataMap & in,
  DataMap & out,
  Dune::InterfaceType iftype,
  Dune::CommunicationDirection dir,
  GridType & grid,
  int TMP)
{
	const int cdim = GridType::dimension;
	const int entityCodim = Entity::codimension;

	// Add subentities from all PB faces
	bool useFace = ((entityCodim == 1)&&(codim >= 1));

	// Add subentities of entities if they are to be communicated
	bool useEntitySend = ((entityCodim == 0)&&(usePartitionType(entity.partitionType (), iftype, dir, true) ));
	bool useEntityRecv = ((entityCodim == 0)&&(usePartitionType(entity.partitionType (), iftype, dir, false) ));

	if (useFace || useEntitySend || useEntityRecv)
	{
      typedef typename GridType::LeafIndexSet   LeafIndexSet;
      const LeafIndexSet & indset = grid.leafIndexSet();
      // Get index set

      Dune::GeometryType gt;  gt.makeSimplex(cdim);
      int nSub = Dune::ReferenceElements<double, cdim>::general(gt).size(0, entityCodim, codim);

      for (int i = 0; i < nSub; i++)
      {
         int ind = indset.subIndex(entity, i, codim);

         if (useFace || useEntitySend)  { in[ind] = TMP;  }
         if (useFace || useEntityRecv)  { out[ind] = 0;   }
      }
	}
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
	  typedef typename LeafGridView::IntersectionIterator IntersectionIterator;
	  typedef typename IntersectionIterator :: Intersection Intersection;

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
	  typedef typename LeafGridView::template Codim<0>::Iterator EntityLeafIterator;
	  EntityLeafIterator ibegin = leafView.template begin<0>();
	  EntityLeafIterator iend   = leafView.template end<0>();

	  bool withGhosts = grid.gridbase().withGhostElements();

	  //std::cout << " -- creating send-arrays" << std::endl;
	  for (EntityLeafIterator it = ibegin; it != iend; ++it)
	  {
		typedef typename LeafGridView::template Codim< 0 >::Entity Element;
		typedef typename LeafGridView::template Codim< 1 >::Entity Face;

		const Element &entity = *it;
		bool marked_internal = false;

		// Only iterate over internal elements for simplicity
		if (entity.partitionType() != Dune::PartitionType::GhostEntity)
		{
	        const IntersectionIterator nend = leafView.iend(entity);
			for( IntersectionIterator nit = leafView.ibegin(entity); nit != nend; ++nit )
			{
			  const Intersection &intersection = *nit;

			  // We are interested in marking process boundaries, and their neighbor entities
			  if (isProcessBoundary(intersection, withGhosts))
			  {
				  //std::cout << "process=" << mpihelper.rank() << " visiting face " << indset.subIndex(entity, intersection.indexInInside(), 1) << std::endl;

				  const Face thisFace = entity.template subEntity<1>(intersection.indexInInside());
				  const Element & entityOut = intersection.outside();

				  mark_subentity<GridType, Face, DataMap, codim>(thisFace, in, out, iftype, dir, grid, TMP);
				  if (!marked_internal)  { mark_subentity<GridType, Element, DataMap, codim>(entity, in, out, iftype, dir, grid, TMP); }
				  if (withGhosts)        { mark_subentity<GridType, Element, DataMap, codim>(entityOut, in, out, iftype, dir, grid, TMP); }
				  marked_internal = true;
			  }
			}
		}


	  }
	  int recvSize = out.size();

	  // Create datahandle
	  std::cout << " -- Initialising datahandle" << std::endl;
	  typedef DataHandleConstant<GridType, LeafIndexSet, ValueType>        DHConstImpl;
	  typedef typename Dune::CommDataHandleIF< DHConstImpl, ValueType >    DHConst;
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






int main (int argc , char **argv) {
	static Dune::MPIHelper & mpihelper = Dune::MPIHelper::instance(argc, argv);

	// Define curvilinear grid
	const int dim = 3;
	typedef  double    ctype;

    // Instantiation of the logging message
    typedef Dune::LoggingMessage<Dune::LoggingMessageHelper::Phase::DEVELOPMENT_PHASE>   LoggingMessageDev;
    LoggingMessageDev::getInstance().verbose(true);
    LoggingMessageDev::getInstance().processVerbose(true);

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


	/*
	std::cout << " Brief self-test";
	std::cout << " " << usePartitionType(Dune::PartitionType::InteriorEntity, Dune::InterfaceType::InteriorBorder_All_Interface, Dune::CommunicationDirection::ForwardCommunication, true);
	std::cout << " " << usePartitionType(Dune::PartitionType::GhostEntity, Dune::InterfaceType::InteriorBorder_All_Interface, Dune::CommunicationDirection::ForwardCommunication, true);
	std::cout << " " << usePartitionType(Dune::PartitionType::InteriorEntity, Dune::InterfaceType::InteriorBorder_All_Interface, Dune::CommunicationDirection::ForwardCommunication, false);
	std::cout << " " << usePartitionType(Dune::PartitionType::GhostEntity, Dune::InterfaceType::InteriorBorder_All_Interface, Dune::CommunicationDirection::ForwardCommunication, false);
	std::cout << std::endl;
	*/


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
