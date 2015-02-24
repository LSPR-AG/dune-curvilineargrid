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

    // CurvIntersectionIterator
    // --------------------

    template< class Grid >
    class CurvIntersectionIterator
    {
      typedef typename remove_const< Grid >::type::Traits  Traits;
      typedef typename remove_const< Grid >::type::ctype   ctype;

      static const int dimension   = remove_const< Grid >::type::dimension;

      typedef Dune::CurvilinearGridBase<ctype,dimension>    GridBaseType;

      typedef typename GridBaseType::LocalIndexType             LocalIndexType;
      typedef typename GridBaseType::InternalIndexType          InternalIndexType;
      typedef typename GridBaseType::StructuralType             StructuralType;
      typedef typename GridBaseType::InterpolatoryOrderType     InterpolatoryOrderType;

      typedef Dune::CurvGrid::CurvIntersection<Grid>  IntersectionImpl;


    public:
      typedef Dune::Intersection< Grid, IntersectionImpl > Intersection;

      typedef typename Traits::template Codim< 0 >::EntityPointer EntityPointer;

      CurvIntersectionIterator (LocalIndexType localIndexInside,
    		  InternalIndexType subIndexInside,
    		  GridBaseType & gridbase)
        : intersection_( IntersectionImpl(localIndexInside, subIndexInside, gridbase))
      {}

      CurvIntersectionIterator ( const CurvIntersectionIterator &other )
        : intersection_( IntersectionImpl( Grid::getRealImplementation( other.intersection_ ) ) )
      {}

      CurvIntersectionIterator &operator= ( const CurvIntersectionIterator &other )
      {
        Grid::getRealImplementation( intersection_ ) = Grid::getRealImplementation( other.intersection_ );
        return *this;
      }

      bool equals ( const CurvIntersectionIterator &other ) const  { return (intersection_ == other.intersection_); }

      void increment ()  { intersectionImpl().next(); }

      const Intersection &dereference () const  { return intersection_; }

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
