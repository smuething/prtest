// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_BACKEND_ISTLMATRIXBACKEND_HH
#define DUNE_PDELAB_BACKEND_ISTLMATRIXBACKEND_HH

#include<utility>
#include<vector>
#include<set>

#include <dune/pdelab/common/unordered_map.hh>
#include <dune/pdelab/common/unordered_set.hh>
#include <dune/pdelab/ordering/orderingbase.hh>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <dune/pdelab/backend/istl/tags.hh>
#include <dune/pdelab/backend/backendselector.hh>

namespace Dune {
  namespace PDELab {

    template<typename RV, typename CV, typename block_type>
    struct matrix_for_vectors;

    template<typename B1, typename A1, typename B2, typename A2, typename block_type>
    struct matrix_for_vectors<Dune::BlockVector<B1,A1>,Dune::BlockVector<B2,A2>,block_type>
    {
      typedef Dune::BCRSMatrix<block_type> type;
    };

    template<typename B1, int n1, typename B2, int n2, typename block_type>
    struct matrix_for_vectors<Dune::FieldVector<B1,n1>,Dune::FieldVector<B2,n2>,block_type>
    {
      typedef Dune::FieldMatrix<block_type,n1,n2> type;
    };

    template<typename E, typename RV, typename CV, std::size_t blocklevel>
    struct recursive_build_matrix_type
    {
      typedef typename matrix_for_vectors<RV,CV,typename recursive_build_matrix_type<E,typename RV::block_type,typename CV::block_type,blocklevel-1>::type>::type type;
    };

    template<typename E, typename RV, typename CV>
    struct recursive_build_matrix_type<E,RV,CV,1>
    {
      typedef Dune::FieldMatrix<E,RV::dimension,CV::dimension> type;
    };


    template<typename E, typename RV, typename CV>
    struct build_matrix_type
    {

      dune_static_assert(RV::blocklevel == CV::blocklevel,"Both vectors must have identical blocking depth");

      typedef typename recursive_build_matrix_type<E,RV,CV,RV::blocklevel>::type type;

    };

    template<typename RowOrdering, typename ColOrdering, typename SubPattern_ = void>
    class Pattern
      : public std::vector<unordered_map<std::size_t,SubPattern_> >
    {

    public:

      typedef SubPattern_ SubPattern;

      template<typename RI, typename CI>
      void add_link(const RI& ri, const CI& ci)
      {
        recursive_add_entry(ri.view(),ci.view());
      }

      template<typename RI, typename CI>
      void recursive_add_entry(const RI& ri, const CI& ci)
      {
        this->resize(_row_ordering.blockCount());
        std::pair<typename unordered_map<std::size_t,SubPattern>::iterator,bool> r = (*this)[ri.back()].insert(make_pair(ci.back(),SubPattern(_row_ordering.childOrdering(ri.back()),_col_ordering.childOrdering(ci.back()))));
        r.first->second.recursive_add_entry(ri.back_popped(),ci.back_popped());
      }

      Pattern(const RowOrdering& row_ordering, const ColOrdering& col_ordering)
        : _row_ordering(row_ordering)
        , _col_ordering(col_ordering)
      {}

    private:

      const RowOrdering& _row_ordering;
      const ColOrdering& _col_ordering;

    };

    template<typename RowOrdering, typename ColOrdering>
    class Pattern<RowOrdering,ColOrdering,void>
      : public std::vector<unordered_set<std::size_t> >
    {

    public:

      typedef void SubPattern;

      template<typename RI, typename CI>
      void add_link(const RI& ri, const CI& ci)
      {
        recursive_add_entry(ri,ci);
      }

      template<typename RI, typename CI>
      void recursive_add_entry(const RI& ri, const CI& ci)
      {
        this->resize(_row_ordering.blockCount());
        (*this)[ri.back()].insert(ci.back());
      }

      Pattern(const RowOrdering& row_ordering, const ColOrdering& col_ordering)
        : _row_ordering(row_ordering)
        , _col_ordering(col_ordering)
      {}

    private:

      const RowOrdering& _row_ordering;
      const ColOrdering& _col_ordering;

    };

    template<typename M, int blocklevel = M::blocklevel>
    struct requires_pattern
    {
      static const bool value = requires_pattern<typename M::block_type,blocklevel-1>::value;
    };

    template<typename M>
    struct requires_pattern<M,0>
    {
      static const bool value = false;
    };

    template<typename B, typename A, int blocklevel>
    struct requires_pattern<BCRSMatrix<B,A>,blocklevel>
    {
      static const bool value = true;
    };

    template<typename M, typename RowOrdering, typename ColOrdering, bool pattern>
    struct _build_pattern_type
    {
      typedef void type;
    };

    template<typename M, typename RowOrdering, typename ColOrdering>
    struct _build_pattern_type<M,RowOrdering,ColOrdering,true>
    {
      typedef Pattern<RowOrdering,ColOrdering,typename _build_pattern_type<typename M::block_type,RowOrdering,ColOrdering,requires_pattern<typename M::block_type>::value>::type> type;
    };

    template<typename M, typename GFSV, typename GFSU, typename Tag>
    struct build_pattern_type
    {

      typedef OrderingBase<
        typename GFSV::Ordering::Traits::DOFIndex,
        typename GFSV::Ordering::Traits::GlobalDOFIndex,
        typename GFSV::Ordering::Traits::ContainerIndex
        > RowOrdering;

      typedef OrderingBase<
        typename GFSU::Ordering::Traits::DOFIndex,
        typename GFSU::Ordering::Traits::GlobalDOFIndex,
        typename GFSU::Ordering::Traits::ContainerIndex
        > ColOrdering;

      typedef typename _build_pattern_type<M,RowOrdering,ColOrdering,requires_pattern<M>::value>::type type;
    };

    template<typename M, typename GFSV, typename GFSU>
    struct build_pattern_type<M,GFSV,GFSU,FlatContainerAllocationTag>
    {
      typedef Pattern<typename GFSV::Ordering, typename GFSU::Ordering> type;
    };


    template<typename RI, typename CI, typename Block>
    typename Block::field_type&
    access_istl_matrix_element(istl::tags::bcrs_matrix, Block& b, const RI& ri, const CI& ci, int i, int j)
    {
      return access_istl_matrix_element(istl::container_tag(b[ri[i]][ci[j]]),b[ri[i]][ci[j]],ri,ci,i-1,j-1);
    }

    template<typename RI, typename CI, typename Block>
    typename Block::field_type&
    access_istl_matrix_element(istl::tags::field_matrix_1_1, Block& b, const RI& ri, const CI& ci, int i, int j)
    {
      assert(i == -1);
      assert(j == -1);
      return b[0][0];
    }

    template<typename RI, typename CI, typename Block>
    typename Block::field_type&
    access_istl_matrix_element(istl::tags::field_matrix_n_m, Block& b, const RI& ri, const CI& ci, int i, int j)
    {
      assert(i == 0);
      assert(j == 0);
      return b[ri[0]][ci[0]];
    }

    template<typename RI, typename CI, typename Block>
    typename Block::field_type&
    access_istl_matrix_element(istl::tags::field_matrix_1_m, Block& b, const RI& ri, const CI& ci, int i, int j)
    {
      assert(i == -1);
      assert(j == 0);
      return b[0][ci[0]];
    }

    template<typename RI, typename CI, typename Block>
    typename Block::field_type&
    access_istl_matrix_element(istl::tags::field_matrix_n_1, Block& b, const RI& ri, const CI& ci, int i, int j)
    {
      assert(i == 0);
      assert(j == -1);
      return b[ri[0]][0];
    }


    template<typename RI, typename CI, typename Block>
    const typename Block::field_type&
    access_istl_matrix_element(istl::tags::bcrs_matrix, const Block& b, const RI& ri, const CI& ci, int i, int j)
    {
      return access_istl_matrix_element(istl::container_tag(b[ri[i]][ci[j]]),b[ri[i]][ci[j]],ri,ci,i-1,j-1);
    }

    template<typename RI, typename CI, typename Block>
    const typename Block::field_type&
    access_istl_matrix_element(istl::tags::field_matrix_1_1, const Block& b, const RI& ri, const CI& ci, int i, int j)
    {
      assert(i == -1);
      assert(j == -1);
      return b[0][0];
    }

    template<typename RI, typename CI, typename Block>
    const typename Block::field_type&
    access_istl_matrix_element(istl::tags::field_matrix_n_m, const Block& b, const RI& ri, const CI& ci, int i, int j)
    {
      assert(i == 0);
      assert(j == 0);
      return b[ri[0]][ci[0]];
    }

    template<typename RI, typename CI, typename Block>
    const typename Block::field_type&
    access_istl_matrix_element(istl::tags::field_matrix_1_m, const Block& b, const RI& ri, const CI& ci, int i, int j)
    {
      assert(i == -1);
      assert(j == 0);
      return b[0][ci[0]];
    }

    template<typename RI, typename CI, typename Block>
    const typename Block::field_type&
    access_istl_matrix_element(istl::tags::field_matrix_n_1, const Block& b, const RI& ri, const CI& ci, int i, int j)
    {
      assert(i == 0);
      assert(j == -1);
      return b[ri[0]][0];
    }


    template<typename RI, typename Block>
    void clear_istl_matrix_row(istl::tags::bcrs_matrix, Block& b, const RI& ri, int i)
    {
      typedef typename Block::ColIterator col_iterator_type;
      const col_iterator_type end = b[ri[i]].end();
      for(col_iterator_type cit = b[ri[i]].begin(); cit != end; ++cit)
        clear_istl_matrix_row(istl::container_tag(*cit),*cit,ri,i-1);
    }

    template<typename RI, typename Block>
    void clear_istl_matrix_row(istl::tags::field_matrix_1_any, Block& b, const RI& ri, int i)
    {
      assert(i == -1);
      b[0] = 0;
    }

    template<typename RI, typename Block>
    void clear_istl_matrix_row(istl::tags::field_matrix_n_any, Block& b, const RI& ri, int i)
    {
      assert(i == 0);
      b[ri[0]] = 0;
    }


    template<typename OrderingV, typename OrderingU, typename Pattern, typename Container>
    typename enable_if<
      !is_same<typename Pattern::SubPattern,void>::value &&
      requires_pattern<Container>::value
      >::type
    allocate_istl_matrix(const OrderingV& ordering_v,
                         const OrderingU& ordering_u,
                         const Pattern& p,
                         Container& c)
    {
      c.setSize(ordering_v.blockCount(),ordering_u.blockCount(),false);
      c.setBuildMode(Container::random);

      for (std::size_t i = 0; i < c.N(); ++i)
        c.setrowsize(i,p[i].size());
      c.endrowsizes();

      for (std::size_t i = 0; i < c.N(); ++i)
        for (typename Pattern::value_type::const_iterator cit = p[i].begin(); cit != p[i].end(); ++cit)
          c.addindex(i,cit->first);
      c.endindices();

      for (std::size_t i = 0; i < c.N(); ++i)
        for (typename Pattern::value_type::const_iterator cit = p[i].begin(); cit != p[i].end(); ++cit)
          {
            allocate_istl_matrix(ordering_v.childOrdering(i),
                                 ordering_u.childOrdering(cit->first),
                                 cit->second,
                                 c[i][cit->first]);
          }
    }

    template<typename OrderingV, typename OrderingU, typename Pattern, typename Container>
    typename enable_if<
      !is_same<typename Pattern::SubPattern,void>::value &&
      !requires_pattern<Container>::value
      >::type
    allocate_istl_matrix(const OrderingV& ordering_v,
                         const OrderingU& ordering_u,
                         const Pattern& p,
                         Container& c)
    {
      for (std::size_t i = 0; i < c.N(); ++i)
        for (typename Pattern::value_type::iterator cit = p[i].begin(); cit != p[i].end(); ++cit)
          {
            allocate_istl_matrix(ordering_v.childOrdering(i),
                                 ordering_u.childOrdering(cit->first),
                                 cit->second,
                                 c[i][cit->first]);
          }
    }

    template<typename OrderingV, typename OrderingU, typename Pattern, typename Container>
    typename enable_if<
      is_same<typename Pattern::SubPattern,void>::value
      >::type
    allocate_istl_matrix(const OrderingV& ordering_v,
                         const OrderingU& ordering_u,
                         const Pattern& p,
                         Container& c)
    {
      c.setSize(ordering_v.blockCount(),ordering_u.blockCount(),false);
      c.setBuildMode(Container::random);

      for (std::size_t i = 0; i < c.N(); ++i)
        c.setrowsize(i,p[i].size());
      c.endrowsizes();

      for (std::size_t i = 0; i < c.N(); ++i)
        for (typename Pattern::value_type::const_iterator cit = p[i].begin(); cit != p[i].end(); ++cit)
          c.addindex(i,*cit);
      c.endindices();
    }


    template<typename GFSV, typename GFSU, typename C>
    class ISTLMatrixContainer
    {

    public:

      typedef typename C::field_type ElementType;
      typedef ElementType E;
      typedef C Container;
      typedef C BaseT;
      typedef typename C::field_type field_type;
      typedef typename C::block_type block_type;
      typedef typename C::size_type size_type;

      typedef typename GFSV::Ordering::Traits::ContainerIndex RowIndex;
      typedef typename GFSU::Ordering::Traits::ContainerIndex ColIndex;

      typedef typename build_pattern_type<C,GFSV,GFSU,typename GFSV::Ordering::ContainerAllocationTag>::type Pattern;

      template<typename RowCache, typename ColCache>
      class LocalView
      {

        dune_static_assert((is_same<typename RowCache::LocalFunctionSpace::Traits::GridFunctionSpace,GFSV>::value),
                           "The RowCache passed to LocalView must belong to the underlying GFSV");

        dune_static_assert((is_same<typename ColCache::LocalFunctionSpace::Traits::GridFunctionSpace,GFSU>::value),
                           "The ColCache passed to LocalView must belong to the underlying GFSU");

      public:

        typedef typename C::field_type E;
        typedef E ElementType;

        typedef RowCache RowIndexCache;
        typedef ColCache ColIndexCache;

        typedef typename RowCache::LocalFunctionSpace LFSV;
        typedef typename ColCache::LocalFunctionSpace LFSU;

        typedef typename LFSV::Traits::DOFIndex RowDOFIndex;
        typedef typename LFSV::Traits::GridFunctionSpace::Ordering::Traits::ContainerIndex RowContainerIndex;

        typedef typename LFSU::Traits::DOFIndex ColDOFIndex;
        typedef typename LFSU::Traits::GridFunctionSpace::Ordering::Traits::ContainerIndex ColContainerIndex;

        LocalView()
          : _container(nullptr)
          , _row_cache(nullptr)
          , _col_cache(nullptr)
        {}

        LocalView(ISTLMatrixContainer& container)
          : _container(&container)
          , _row_cache(nullptr)
          , _col_cache(nullptr)
        {}

        const RowIndexCache& rowIndexCache() const
        {
          assert(_row_cache);
          return *_row_cache;
        }

        const ColIndexCache& colIndexCache() const
        {
          assert(_col_cache);
          return *_col_cache;
        }

        void attach(ISTLMatrixContainer& container)
        {
          _container = &container;
        }

        void detach()
        {
          _container = nullptr;
        }

        void bind(const RowCache& row_cache, const ColCache& col_cache)
        {
          _row_cache = &row_cache;
          _col_cache = &col_cache;
        }

        void unbind()
        {
        }

        size_type N() const
        {
          return _row_cache->size();
        }

        size_type M() const
        {
          return _col_cache->size();
        }

        template<typename LC>
        void read(LC& local_container) const
        {
          for (size_type i = 0; i < N(); ++i)
            for (size_type j = 0; j < M(); ++j)
              local_container.getEntry(i,j) = (*_container)(_row_cache->containerIndex(i),_col_cache->containerIndex(j));
        }

        template<typename LC>
        void write(const LC& local_container) const
        {
          for (size_type i = 0; i < N(); ++i)
            for (size_type j = 0; j < M(); ++j)
              (*_container)(_row_cache->containerIndex(i),_col_cache->containerIndex(j)) = local_container.getEntry(i,j);
        }

        template<typename LC>
        void add(const LC& local_container) const
        {
          for (size_type i = 0; i < N(); ++i)
            for (size_type j = 0; j < M(); ++j)
              (*_container)(_row_cache->containerIndex(i),_col_cache->containerIndex(j)) += local_container.getEntry(i,j);
        }

        void commit()
        {
        }


        ElementType& operator()(size_type i, size_type j)
        {
          return (*_container)(_row_cache->containerIndex(i),_col_cache->containerIndex(j));
        }

        const ElementType& operator()(size_type i, size_type j) const
        {
          return (*_container)(_row_cache->containerIndex(i),_col_cache->containerIndex(j));
        }

        ElementType& operator()(const RowDOFIndex& i, const ColDOFIndex& j)
        {
          return (*_container)(_row_cache->containerIndex(i),_col_cache->containerIndex(j));
        }

        const ElementType& operator()(const RowDOFIndex& i, const ColDOFIndex& j) const
        {
          return (*_container)(_row_cache->containerIndex(i),_col_cache->containerIndex(j));
        }

        ElementType& operator()(const RowContainerIndex& i, const ColContainerIndex& j)
        {
          return (*_container)(i,j);
        }

        const ElementType& operator()(const RowContainerIndex& i, const ColContainerIndex& j) const
        {
          return (*_container)(i,j);
        }

        ElementType& operator()(const RowContainerIndex& i, size_type j)
        {
          return (*_container)(i,_col_cache->containerIndex(j));
        }

        const ElementType& operator()(const RowContainerIndex& i, size_type j) const
        {
          return (*_container)(i,_col_cache->containerIndex(j));
        }

        ElementType& operator()(size_type i, const ColContainerIndex& j)
        {
          return (*_container)(_row_cache->containerIndex(i),j);
        }

        const ElementType& operator()(size_type i, const ColContainerIndex& j) const
        {
          return (*_container)(_row_cache->containerIndex(i),j);
        }


        void add(size_type i, size_type j, const ElementType& v)
        {
          (*_container)(_row_cache->containerIndex(i),_col_cache->containerIndex(j)) += v;
        }

        void add(const RowDOFIndex& i, const ColDOFIndex& j, const ElementType& v)
        {
          (*_container)(_row_cache->containerIndex(i),_col_cache->containerIndex(j)) += v;
        }

        void add(const RowContainerIndex& i, const ColContainerIndex& j, const ElementType& v)
        {
          (*_container)(i,j) += v;
        }

        void add(const RowContainerIndex& i, size_type j, const ElementType& v)
        {
          (*_container)(i,_col_cache->containerIndex(j)) += v;
        }

        void add(size_type i, const ColContainerIndex& j, const ElementType& v)
        {
          (*_container)(_row_cache->containerIndex(i),j) += v;
        }

        ISTLMatrixContainer& global_container()
        {
          return *_container;
        }

        const ISTLMatrixContainer& global_container() const
        {
          return *_container;
        }

      private:

        ISTLMatrixContainer* _container;
        const RowCache* _row_cache;
        const ColCache* _col_cache;

      };

      template<typename GO>
      ISTLMatrixContainer (const GO& go)
      {
        Pattern pattern(go.testGridFunctionSpace().ordering(),go.trialGridFunctionSpace().ordering());
        go.fill_pattern(pattern);
        allocate_istl_matrix(go.testGridFunctionSpace().ordering(),
                             go.trialGridFunctionSpace().ordering(),
                             pattern,
                             _container);
      }

      template<typename GO>
      ISTLMatrixContainer (const GO& go, const E& e)
      {
        Pattern pattern(go.testGridFunctionSpace().ordering(),go.trialGridFunctionSpace().ordering());
        go.fill_pattern(pattern);
        allocate_istl_matrix(go.testGridFunctionSpace().ordering(),
                             go.trialGridFunctionSpace().ordering(),
                             pattern,
                             _container);
        _container = e;
      }

      ISTLMatrixContainer () : _container() {}

      size_type N() const
      {
        return _container.N();
      }

      size_type M() const
      {
        return _container.M();
      }

      ISTLMatrixContainer& operator= (const E& e)
      {
        _container = e;
        return *this;
      }

      ISTLMatrixContainer& operator*= (const E& e)
      {
        _container *= e;
        return *this;
      }

      E& operator()(const RowIndex& ri, const ColIndex& ci)
      {
        return access_istl_matrix_element(istl::container_tag(_container),_container,ri,ci,ri.size()-1,ci.size()-1);
      }

      const E& operator()(const RowIndex& ri, const ColIndex& ci) const
      {
        return access_istl_matrix_element(istl::container_tag(_container),_container,ri,ci,ri.size()-1,ci.size()-1);
      }

      const Container& base() const
      {
        return _container;
      }

      Container& base()
      {
        return _container;
      }

      void flush()
      {}

      void finalize()
      {}

      void clear_row(const RowIndex& ri, const E& diagonal_entry)
      {
        clear_istl_matrix_row(istl::container_tag(_container),_container,ri,ri.size()-1);
        (*this)(ri,ri) = diagonal_entry;
      }

    private:

      Container _container;

    };

    struct ISTLMatrixBackend
    {

      typedef std::size_t size_type;

      template<typename VV, typename VU, typename E>
      struct MatrixHelper
      {
        typedef ISTLMatrixContainer<typename VV::GridFunctionSpace,typename VU::GridFunctionSpace,typename build_matrix_type<E,typename VV::Container,typename VU::Container>::type > type;
      };
    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTLMATRIXBACKEND_HH
