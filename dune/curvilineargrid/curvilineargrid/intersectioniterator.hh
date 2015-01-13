// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CURVGRID_INTERSECTIONITERATOR_HH
#define DUNE_CURVGRID_INTERSECTIONITERATOR_HH

#include <dune/curvilineargrid/curvilineargrid/entitypointer.hh>
#include <dune/curvilineargrid/curvilineargrid/intersection.hh>

namespace Dune
{

  namespace CurvGrid
  {

    // IntersectionIterator
    // --------------------

    template< class Grid >
    class IntersectionIterator
    {
      typedef typename remove_const< Grid >::type::Traits Traits;

      typedef typename Traits::IntersectionImpl IntersectionImpl;

      typedef typename Traits::template Codim< 0 >::EntityPointerImpl EntityPointerImpl;
      typedef typename Traits::template Codim< 0 >::Geometry ElementGeometry;

    public:
      typedef Dune::Intersection< Grid, IntersectionImpl > Intersection;

      typedef typename Traits::template Codim< 0 >::EntityPointer EntityPointer;

      template< class Entity >
      IntersectionIterator (LocalIndexType localIndexInside,
    		  InternalIndexType subIndexInside,
    		  GridBaseType & gridbase)
        : intersection_(localIndexInside, subIndexInside, gridbase)
      {}

      IntersectionIterator ( const IntersectionIterator &other )
        : intersection_( IntersectionImpl( Grid::getRealImplementation( other.intersection_ ) ) )
      {}

      IntersectionIterator &operator= ( const IntersectionIterator &other )
      {
        hostIterator_ = other.hostIterator_;
        Grid::getRealImplementation( intersection_ ) = Grid::getRealImplementation( other.intersection_ );
        return *this;
      }

      bool equals ( const IntersectionIterator &other ) const
      {
        return (hostIterator_ == other.hostIterator_);
      }

      void increment ()
      {
    	  // If (SubentityIndex >= SubentitySize)  { THROW ERROR }
    	  // Increase iterator, while
    	  // 1) SubentityIndex < SubentitySize
    	  // 2) FacePartitionType == Ghost


        ++hostIterator_;
        intersectionImpl().invalidate();
      }

      const Intersection &dereference () const
      {
        if( !intersectionImpl() )
          intersectionImpl().initialize( *hostIterator_ );
        return intersection_;
      }

    private:
      IntersectionImpl &intersectionImpl () const
      {
        return Grid::getRealImplementation( intersection_ );
      }

      mutable Intersection intersection_;
    };

  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_INTERSECTIONITERATOR_HH
