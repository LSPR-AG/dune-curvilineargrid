// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CURVGRID_ENTITYSEED_HH
#define DUNE_CURVGRID_ENTITYSEED_HH

#include <dune/common/typetraits.hh>

#include <dune/grid/common/entityseed.hh>
#include <dune/curvilineargrid/curvilineargrid/capabilities.hh>

namespace Dune
{

  namespace CurvGrid
  {


    template< int codim, class Grid >
    class CurvEntitySeed
    {
    	typedef typename remove_const< Grid >::type::Traits  Traits;
    	typedef typename remove_const< Grid >::type::ctype   ctype;

    public:
      static const int codimension    = codim;
      static const int dimension      = remove_const< Grid >::type::dimension;
      static const int mydimension    = dimension - codimension;
      static const int dimensionworld = remove_const< Grid >::type::dimensionworld;

      typedef Dune::CurvilinearGridBase<ctype,dimension>    GridBaseType;
      typedef typename GridBaseType::LocalIndexType         LocalIndexType;

      //! default construct an invalid entity seed
      CurvEntitySeed (LocalIndexType index, Dune::PartitionIteratorType pitype, const GridBaseType & gridbase)
  	  	  : index_(index), pitype_(pitype), gridbase_(&gridbase)
      {}


	  //! Copy assignment operator from an existing seed.
      CurvEntitySeed& operator=(const CurvEntitySeed& other)
	  {
    	  gridbase_ = other.gridbase_;
    	  index_ = other.index_;
    	  pitype_ = other.pitype_;
	      return *this;
	  }


      //! check whether the EntitySeed refers to a valid Entity
      bool isValid() const
      {
    	  return !(gridbase_->entityIndexIterator(codim, localIndex) == gridbase_->entityIndexEnd(codim) );
      }

      // Returns local index of the associated entity
      LocalIndexType localIndex() const  { return index_; }

      // Returns partition type of the associated entity.
      // Not needed for entity itself, but needed for the range the entity would iterate over if necessary
      Dune::PartitionIteratorType partitionIteratorType() const  { return pitype_; }

      // Returns a reference to the grid base class
      GridBaseType & gridBase() const  { return *gridbase_; }

    private:
      GridBaseType * gridbase_;
      LocalIndexType index_;
      Dune::PartitionIteratorType pitype_;
    };

  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_ENTITYSEED_HH
