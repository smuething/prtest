// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_POWERCOMPOSITEGRIDFUNCTIONSPACEBASE_HH
#define DUNE_PDELAB_POWERCOMPOSITEGRIDFUNCTIONSPACEBASE_HH

#include <cstddef>
#include <functional>
#include <vector>

#include <dune/pdelab/common/typetree/fixedcapacitystack.hh>
#include <dune/pdelab/common/typetree/traversal.hh>
#include <dune/pdelab/common/typetree/traversalutilities.hh>
#include <dune/pdelab/common/typetree/utility.hh>
#include <dune/pdelab/common/typetree/transformation.hh>
#include <dune/pdelab/common/typetree/visitor.hh>
#include <dune/pdelab/constraints/constraintstransformation.hh>

#include <dune/pdelab/gridfunctionspace/tags.hh>
#include <dune/pdelab/gridfunctionspace/utility.hh>
#include <dune/pdelab/ordering/lexicographicordering.hh>
#include <dune/pdelab/ordering/entityblockedlocalordering.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace grid function space
    //! \ingroup PDELab
    //! \{

    //! Trait class for the multi component grid function spaces
    template<typename G, typename B, typename M, std::size_t k>
    struct PowerCompositeGridFunctionSpaceTraits
    {
      enum{
        //! \brief True if this grid function space is composed of others.
        isComposite = 1,
        //! \brief number of child spaces
        noChilds = k
      };

      const static std::size_t CHILDREN = k;

      //! \brief the grid view where grid function is defined upon
      typedef G GridViewType;

      typedef G GridView;

      //! \brief vector backend
      typedef B BackendType;

      typedef B Backend;

      //! \brief mapper
      typedef M MapperType;

      //! \brief short cut for size type exported by Backend
      typedef typename B::size_type SizeType;
    };

    //! Mixin class providing common functionality of PowerGridFunctionSpace and CompositeGridFunctionSpace
    template<typename GridFunctionSpace, typename GV, typename B, typename Mapper, std::size_t k>
    class PowerCompositeGridFunctionSpaceBase
    {

#ifndef DOXYGEN

      const GridFunctionSpace& gfs() const
      {
        return static_cast<const GridFunctionSpace&>(*this);
      }

      GridFunctionSpace& gfs()
      {
        return static_cast<GridFunctionSpace&>(*this);
      }

#endif // DOXYGEN

    public:

      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<GV,B,Mapper,k> Traits;

      typedef Mapper OrderingTag;

      // TODO: Do not just use constraints from child 0!
      //! extract type for storing constraints
      template<typename E>
      struct ConstraintsContainer
      {
        typedef typename SelectType<
          is_same<
            typename GridFunctionSpace::template Child<0>::type::template ConstraintsContainer<E>::Type,
            EmptyTransformation
            >::value,
          EmptyTransformation,
          ConstraintsTransformation<
            typename GridFunctionSpace::Ordering::Traits::DOFIndex,
            typename GridFunctionSpace::Ordering::Traits::ContainerIndex,
            E
            >
          >::Type Type;
      };

      //! recalculate sizes
      void update ()
      {
        gfs().ordering().update();
      }

      //! get dimension of root finite element space
      typename Traits::SizeType globalSize () const
      {
        return gfs().ordering().size();
      }

      //! get dimension of this finite element space
      typename Traits::SizeType size () const
      {
        return gfs().ordering().size();
      }

      //! get max dimension of shape function space
      typename Traits::SizeType maxLocalSize () const
      {
        // this is bullshit !
        return gfs().ordering().maxLocalSize();
      }

      //! Returns whether this GridFunctionSpace contains entities with PartitionType partition.
      bool containsPartition(PartitionType partition) const
      {
        return gfs().ordering().containsPartition(partition);
      }

      //! get grid view
      const typename Traits::GridViewType& gridview () const DUNE_DEPRECATED_MSG("Use gridView() instead of gridview()")
      {
        return gfs().template child<0>().gridView();
      }

      //! get grid view
      const typename Traits::GridViewType& gridView () const
      {
        return gfs().template child<0>().gridView();
      }

      B& backend()
      {
        return _backend;
      }

      const B& backend() const
      {
        return _backend;
      }

      PowerCompositeGridFunctionSpaceBase(const B& backend)
        : _backend(backend)
      {}

      const std::string& name() const
      {
        return _name;
      }

      void name(const std::string& name)
      {
        _name = name;
      }

    private:

      B _backend;
      std::string _name;

    };

  }

}
//! \}
#endif // DUNE_PDELAB_POWERCOMPOSITEGRIDFUNCTIONSPACEBASE_HH
