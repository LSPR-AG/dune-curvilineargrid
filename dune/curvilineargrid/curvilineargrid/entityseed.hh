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
    	typedef typename remove_const< Grid >::type::Traits Traits;

    public:
      static const int codimension    = codim;
      static const int dimension      = Traits::dimension;
      static const int mydimension    = dimension - codimension;
      static const int dimensionworld = Traits::dimensionworld;

	  typedef typename Traits::GridStorageType      GridStorageType;
	  typedef typename Traits::GridBaseType         GridBaseType;


      typedef typename Traits::GlobalIndexType           GlobalIndexType;
      typedef typename Traits::LocalIndexType            LocalIndexType;

      //! default construct an invalid entity seed
      CurvEntitySeed (LocalIndexType index, Dune::PartitionType pitype, GridBaseType & gridbase)
  	  	  : index_(index), pitype_(pitype), gridbase_(gridbase)
      {}

      //! check whether the EntitySeed refers to a valid Entity
      bool isValid() const
      {
    	  return !(gridbase_.entityIndexIterator(codim, localIndex) == gridbase_.entityIndexEnd(codim) );
      }

      // Returns local index of the associated entity
      LocalIndexType localIndex() const  { return index_; }

      // Returns partition type of the associated entity.
      // Not needed for entity itself, but needed for the range the entity would iterate over if necessary
      Dune::PartitionType partitionType() const  { return pitype_; }

      // Returns a reference to the grid base class
      GridBaseType & gridBase() const  { return gridbase_; }

    private:
      GridBaseType & gridbase_;
      LocalIndexType index_;
      Dune::PartitionType pitype_;
    };

  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_ENTITYSEED_HH
