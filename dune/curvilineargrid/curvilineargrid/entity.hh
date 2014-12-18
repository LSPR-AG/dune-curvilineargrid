// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CURVGRID_ENTITY_HH
#define DUNE_CURVGRID_ENTITY_HH

#include <dune/common/nullptr.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/grid.hh>
#include <dune/curvilineargrid/curvilineargrid/capabilities.hh>

namespace Dune
{

  namespace CurvGrid
  {


    template<int cd, int dim, class GridImp, template<int,int,class> class EntityImp>
    class Entity
    {
    protected:
      // type of underlying implementation, for internal use only
      typedef EntityImp< cd, dim, GridImp > Implementation;

      //! Return reference to the real implementation
      Implementation &impl () { return realEntity; }
      //! Return const reference to the real implementation
      const Implementation &impl () const { return realEntity; }

    protected:
      Implementation realEntity;

    public:

      //! \brief The corresponding geometry type
      typedef typename GridImp::template Codim<cd>::Geometry Geometry;

      //! \brief The corresponding entity seed (for storage of entities)
      typedef typename GridImp::template Codim<cd>::EntitySeed EntitySeed;

      enum  { codimension=cd };
      enum  { dimension=dim  };
      enum  { mydimension=dim-cd };


      //! The level of this entity
      int level () const { return realEntity.level(); }

      //! Partition type of this entity
      PartitionType partitionType () const { return realEntity.partitionType(); }

      Geometry geometry () const { return realEntity.geometry(); }

      GeometryType type () const { return realEntity.type(); }

      EntitySeed seed () const { return realEntity.seed(); }

      bool operator==(const Entity& other) const  { return realEntity.equals(other.realEntity); }
      bool operator!=(const Entity& other) const  { return !realEntity.equals(other.realEntity); }


      Entity(const Entity& other) : realEntity(other.realEntity)  { }

      Entity(Entity&& other) : realEntity(std::move(other.realEntity))  { }

      Entity& operator=(const Entity& other)
      {
        realEntity = other.realEntity;
        return *this;
      }

      //! Move assignment operator from an existing entity.
      Entity& operator=(Entity&& other)
      {
        realEntity = std::move(other.realEntity);
        return *this;
      }


      //! Copy constructor from EntityImp
      explicit Entity(const EntityImp<cd,dim,GridImp> & e) : realEntity(e) {}
    };




    template<int dim, class GridImp, template<int,int,class> class EntityImp>
    class Entity <0,dim,GridImp,EntityImp>
    {

      public:
        /** \name Attributes
         *  \{ */

        //! codimensioon of the entity
        static const int codimension;
        //! dimension of the grid
        static const int dimension;
        //! dimension of the entity
        static const int mydimension;
        //! dimension of the world
        static const int dimensionworld;

        //! \b true, if the entity is faked, i.e., if there is no corresponding host entity
        static const bool fake;
        /** \} */

        /** \name Types Required by DUNE
         *  \{ */

        //! type of corresponding local geometry
        typedef typename Traits::template Codim< codimension >::LocalGeometry LocalGeometry;
        //! type of corresponding entity pointer
        typedef typename Traits::template Codim< codimension >::EntityPointer EntityPointer;

        //! type of hierarchic iterator
        typedef typename Traits::HierarchicIterator HierarchicIterator;
        //! type of leaf intersection iterator
        typedef typename Traits::LeafIntersectionIterator LeafIntersectionIterator;
        //! type of level intersection iterator
        typedef typename Traits::LevelIntersectionIterator LevelIntersectionIterator;

        /** \} */

        explicit Entity ( const Grid &grid )  {}

        explicit Entity ( const GeometryImpl &geo )  {}

        template< int codim >
        int count () const  { }

        unsigned int subEntities (unsigned int codim) const  {  }

        template< int codim >
        EntityPointer<codim> subEntity ( int i ) const  { }

        LevelIntersectionIterator ilevelbegin () const  { }

        LevelIntersectionIterator ilevelend () const  {  }

        LeafIntersectionIterator ileafbegin () const  { }

        LeafIntersectionIterator ileafend () const  {  }

        bool hasBoundaryIntersections () const  {  }

        bool isLeaf () const  {  }

        EntityPointer father () const  {  }

        bool hasFather () const  {  }

        LocalGeometry geometryInFather () const  {  }

        HierarchicIterator hbegin ( int maxLevel ) const  {  }

        HierarchicIterator hend ( int maxLevel ) const  {  }

        bool isRegular () const  {  }

        bool isNew () const  {  }

        bool mightVanish () const  {}
    };

  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_ENTITY_HH
