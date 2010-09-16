// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_ISTLSOLVERBACKEND_HH
#define DUNE_ISTLSOLVERBACKEND_HH

#include <dune/common/deprecated.hh>
#include <dune/common/mpihelper.hh>

#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/solvercategory.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/io.hh>
#include <dune/istl/superlu.hh>

#include "../gridfunctionspace/constraints.hh"
#include "../gridfunctionspace/genericdatahandle.hh"
#include "../newton/newton.hh"
#include "istlvectorbackend.hh"

namespace Dune {
  namespace PDELab {

    template<typename X, typename Y, typename GOS>
    class OnTheFlyOperator : public Dune::LinearOperator<X,Y>
    {
    public:
      typedef X domain_type;
      typedef Y range_type;
      typedef typename X::field_type field_type;

      enum {category=Dune::SolverCategory::sequential};

      OnTheFlyOperator (GOS& gos_)
        : gos(gos_)
      {}

      virtual void apply (const X& x, Y& y) const
      {
        y = 0.0;
        gos.jacobian_apply(x,y);
      }

      virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
      {
        Y temp(y);
        temp = 0.0;
        gos.jacobian_apply(x,temp);
        y.axpy(alpha,temp);
      }

    private:
      GOS& gos;
    };


    //========================================================
    // A parallel helper class providing a nonoverlapping
    // decomposition of all degrees of freedom
    //========================================================

    // operator that resets result to zero at constrained DOFS
    template<typename GFS>
    class ParallelISTLHelper
    {
      /**
       * @brief Writes 1<<24 to each data item (of the container) that is gathered or scattered
       * and is neither interior nor border.
       *
       * Can be used to mark ghost cells.
       */
      class GhostGatherScatter
      {
      public:
        template<class MessageBuffer, class EntityType, class DataType>
        void gather (MessageBuffer& buff, const EntityType& e, DataType& data)
        {
          if (e.partitionType()!=Dune::InteriorEntity && e.partitionType()!=Dune::BorderEntity)
            data = (1<<24);
          buff.write(data);
        }

        template<class MessageBuffer, class EntityType, class DataType>
        void scatter (MessageBuffer& buff, const EntityType& e, DataType& data)
        {
          DataType x;
          buff.read(x);
          if (e.partitionType()!=Dune::InteriorEntity && e.partitionType()!=Dune::BorderEntity)
            data = (1<<24);
        }
      };

      /**
       * @brief GatherScatter handle that sets 1<<24 for data items neither associated to
       * the interior or border and take the minimum when scattering.
       *
       * Used to compute an owner rank for each unknown.
       */
      class InteriorBorderGatherScatter
      {
      public:
        template<class MessageBuffer, class EntityType, class DataType>
        void gather (MessageBuffer& buff, const EntityType& e, DataType& data)
        {
          if (e.partitionType()!=Dune::InteriorEntity && e.partitionType()!=Dune::BorderEntity)
            data = (1<<24);
          buff.write(data);
        }

        template<class MessageBuffer, class EntityType, class DataType>
        void scatter (MessageBuffer& buff, const EntityType& e, DataType& data)
        {
          DataType x;
          buff.read(x);
          if (e.partitionType()!=Dune::InteriorEntity && e.partitionType()!=Dune::BorderEntity)
            data = x;
          else
            data = std::min(data,x);
        }
      };

      /**
       * @brief GatherScatter handle for finding out about neighbouring processor ranks.
       *
       */
      template<typename T>
      struct NeighbourGatherScatter
      {
        NeighbourGatherScatter(int rank_, std::set<int>& neighbours_)
          : myrank(rank_), neighbours(neighbours_)
        {}

        template<class MessageBuffer, class DataType>
        void gather (MessageBuffer& buff, DataType& data)
        {
          buff.write(myrank);
        }

        template<class MessageBuffer, class DataType>
        void scatter (MessageBuffer& buff, DataType& data)
        {
          DataType x;
          buff.read(x);
          neighbours.insert((int)x);
        }

        T myrank;
        std::set<int>& neighbours;
      };


      /**
       * @brief GatherScatter handle for finding out about neighbouring processor ranks.
       *
       */
      struct SharedGatherScatter
      {
        template<class MessageBuffer, class DataType>
        void gather (MessageBuffer& buff, DataType& data)
        {
          data=true;
          buff.write(data);
        }

        template<class MessageBuffer, class DataType>
        void scatter (MessageBuffer& buff, DataType& data)
        {
          bool x;
          buff.read(x);
          data = data || x;
        }
      };

      /**
       * @brief GatherScatter handle for finding out about neighbouring processor ranks.
       *
       */
      template<typename B, typename V1>
      struct GlobalIndexGatherScatter
      {
        GlobalIndexGatherScatter(const V1& mask_)
          : mask(mask_)
        {}

        template<class MessageBuffer, class DataType>
        void gather (MessageBuffer& buff, typename B::size_type i, DataType& data)
        {
          //if(B::access(mask, i)>0)
          if(data < std::numeric_limits<DataType>::max())
            // We now the global index and therefore write it
            buff.write(data);
          else
            buff.write(std::numeric_limits<DataType>::max());
        }

        template<class MessageBuffer, class DataType>
        void scatter (MessageBuffer& buff, typename B::size_type i, DataType& data)
        {
          DataType x;
          buff.read(x);
          data = std::min(data, x);
        }
        V1 mask;
      };

      typedef typename GFS::template VectorContainer<double>::Type V;

    public:

      ParallelISTLHelper (const GFS& gfs_)
        : gfs(gfs_), v(gfs,(double)gfs.gridview().comm().rank())
      {
        // find out about ghosts
        Dune::PDELab::GenericDataHandle2<GFS,V,GhostGatherScatter> gdh(gfs,v,GhostGatherScatter());
        if (gfs.gridview().comm().size()>1)
          gfs.gridview().communicate(gdh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);

        // partition interior/border
        Dune::PDELab::GenericDataHandle2<GFS,V,InteriorBorderGatherScatter> dh(gfs,v,InteriorBorderGatherScatter());
        if (gfs.gridview().comm().size()>1)
          gfs.gridview().communicate(dh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);

        // convert vector into mask vector
        for (typename V::size_type i=0; i<v.N(); ++i)
          for (typename V::size_type j=0; j<v[i].N(); ++j)
            if (v[i][j]==gfs.gridview().comm().rank())
              v[i][j] = 1.0;
            else
              v[i][j] = 0.0;
      }

      // keep only DOFs assigned to this processor
      template<typename W>
      void mask (W& w) const
      {
        for (typename V::size_type i=0; i<v.N(); ++i)
          for (typename V::size_type j=0; j<v[i].N(); ++j)
            w[i][j] *= v[i][j];
      }

      // access to mask vector
      double mask (typename V::size_type i, typename V::size_type j) const
      {
        return v[i][j];
      }

#if HAVE_MPI

      /**
       * @brief Creates a matrix suitable for parallel AMG and the parallel information
       *
       * It is silently assumed that the unknows are associated with vertices.
       *
       * @tparam MatrixType The type of the ISTL matrix used.
       * @tparam Comm The type of the OwnerOverlapCopyCommunication
       * @param m The local matrix.
       * @param c The parallel information object providing index set, interfaces and
       * communicators.
       */
      template<typename MatrixType, typename Comm>
      void createIndexSetAndProjectForAMG(MatrixType& m, Comm& c);
#endif
    private:
      const GFS& gfs;
      V v;
    };


    namespace
    {
      template<typename GFS, int k>
      struct BlockwiseIndicesHelper
      {
        enum{ value = false };
      };

      template<typename M, typename B, int k>
      struct BlockSizeIsEqual
      {
        enum{ value = false };
      };

      template<int k>
      struct BlockSizeIsEqual<GridFunctionSpaceBlockwiseMapper,ISTLVectorBackend<k>,k>
      {
        enum{ value = true };
      };

      template<typename GFS>
      struct BlockwiseIndicesHelper<GFS,1>
      {
        enum{ value =  BlockSizeIsEqual<typename GFS::Traits::MapperType,
              typename GFS::Traits::BackendType, GFS::Traits::noChilds>::value};
      };

      template<typename GFS>
      struct BlockwiseIndices
      {
        enum{
          value = BlockwiseIndicesHelper<GFS,GFS::Traits::isComposite>::value
        };
      };

      template<typename GFS, bool b, int k>
      struct BlockProcessorHelper
      {};

      template<typename GFS>
      struct BlockProcessorHelper<GFS,false,1>
      {

        template<typename T>
        struct AMGVectorTypeSelector
        {
          typedef typename T::BaseT  Type;
        };

        template<typename T>
        static typename AMGVectorTypeSelector<T>::Type& getVector(T& t)
        {
          return t.base();
        }

        template<typename G>
        static void postProcessCount(G& g)
        {}

        template<typename G>
        static void increment(G& g, std::size_t i)
        {
          ++g;
        }
        template<typename M, typename TI>
        static void addIndex(const typename TI::GlobalIndex& gi, std::size_t i,
                             typename TI::LocalIndex::Attribute attr, M& m, TI& idxset)
        {
          // Add index
          idxset.add(gi, typename TI::LocalIndex(i, attr));
        }
      };

      template<typename GFS>
      struct BlockProcessorHelper<GFS,true,1>
        : public BlockProcessorHelper<GFS,false,1>
      {};

      template<typename GFS, int k>
      struct BlockProcessorHelper<GFS, true, k>
      {
        template<typename T>
        struct AMGVectorTypeSelector
        {
          typedef typename T::BaseT Type;
        };

        template<typename T>
        static typename AMGVectorTypeSelector<T>::Type& getVector(T& t)
        {
          return t.base();
        }
        template<typename G>
        static void postProcessCount(G& g)
        {
          g=g;
        }

        template<typename G>
        static void increment(G& g, std::size_t i)
        {
          if((i+1)%GFS::Traits::noChilds==0)
            ++g;
        }
        template<typename M, typename TI>
        static void addIndex(const typename TI::GlobalIndex& gi, std::size_t i,
                             typename TI::LocalIndex::Attribute attr, M& m, TI& idxset)
        {
          if(i%GFS::Traits::noChilds==0)
            BlockProcessorHelper<GFS, false, 1>::addIndex(gi, i/GFS::Traits::noChilds,
                                                          attr, m, idxset);
        }

      };

      template<typename GFS>
      struct BlockProcessor
        : public BlockProcessorHelper<GFS, BlockwiseIndices<GFS>::value,
                                      GFS::Traits::BackendType::BlockSize>
      {};

    } // end anonymous namspace


#if HAVE_MPI
    template<typename GFS>
    template<typename M, typename C>
    void ParallelISTLHelper<GFS>::createIndexSetAndProjectForAMG(M& m, C& c)
    {
      typedef typename GFS::Traits::GridViewType GV;
      const GV& gv = gfs.gridview();
      static const std::size_t dim = GV::Grid::dimension;
      if(gv.comm().size()>1 && gv.grid().overlapSize(dim)<1)
        DUNE_THROW(Dune::InvalidStateException, "ParallelISTLHelper::createIndexSetAndProjectForAMG: "
                   <<"Only grids with at least one layer of overlap cells are supported");

      // First find out which dofs we share with other processors
      typedef typename GFS::template VectorContainer<bool>::Type BoolVector;
      BoolVector sharedDOF(gfs, false);
      Dune::PDELab::GenericDataHandle<GFS,BoolVector,SharedGatherScatter> gdh(gfs,sharedDOF,SharedGatherScatter());

      if (gfs.gridview().comm().size()>1)
        gfs.gridview().communicate(gdh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);

      // Count shared dofs that we own
      typedef typename C::ParallelIndexSet::GlobalIndex GlobalIndex;
      GlobalIndex count=0;
      std::size_t noScalars=0;

      for (typename V::size_type i=0; i<v.N(); ++i)
        for (typename V::size_type j=0; j<v[i].N(); ++j, ++noScalars)
          if(v[i][j]==1.0 && sharedDOF[i][j])
            ++count;

      std::cout<<gv.comm().rank()<<": shared count is"<< count.touint()<<std::endl;

      // Maybe divide by block size?
      BlockProcessor<GFS>::postProcessCount(count);
      std::cout<<gv.comm().rank()<<": shared block count is"<< count.touint()<<std::endl;

      std::vector<GlobalIndex> counts(gfs.gridview().comm().size());
      MPI_Allgather(&count, 1, MPITraits<GlobalIndex>::getType(), &(counts[0]),
                    1, MPITraits<GlobalIndex>::getType(),
                    gfs.gridview().comm());

      // Compute start index start_p = \sum_{i=0}^{i<p} counts_i
      GlobalIndex start=0;
      for(int i=0; i<gfs.gridview().comm().rank(); ++i)
        start=start+counts[i];
      //std::cout<<gv.comm().rank()<<": start index = "<<start.touint()<<std::endl;

      typedef typename GFS::template VectorContainer<GlobalIndex>::Type GIVector;
      GIVector scalarIndices(gfs, std::numeric_limits<GlobalIndex>::max());


      for (typename V::size_type i=0, ii=0; i<v.N(); ++i)
        for (typename V::size_type j=0; j<v[i].N(); ++j)
          if(v[i][j]==1.0 && sharedDOF[i][j]){
            scalarIndices[i][j]=start;
            BlockProcessor<GFS>::increment(start, ii++);
          }

      // publish global indices for the shared DOFS to other processors.
      typedef GlobalIndexGatherScatter<typename GFS::Traits::BackendType,V> GIGS;
      Dune::PDELab::GenericDataHandle3<GFS,GIVector,GIGS> gdhgi(gfs, scalarIndices, GIGS(v));
      if (gfs.gridview().comm().size()>1)
        gfs.gridview().communicate(gdhgi,Dune::All_All_Interface,Dune::ForwardCommunication);

      // Setup the index set
      c.indexSet().beginResize();
      for (typename V::size_type i=0, ii=0; i<v.N(); ++i)
        for (typename V::size_type j=0; j<v[i].N(); ++j, ++ii){
          Dune::OwnerOverlapCopyAttributeSet::AttributeSet attr;
          if(scalarIndices[i][j]!=std::numeric_limits<GlobalIndex>::max()){
            // global index exist in index set
            if(v[i][j]>0){
              // This dof is managed by us.
              attr = Dune::OwnerOverlapCopyAttributeSet::owner;
            }else{
              attr = Dune::OwnerOverlapCopyAttributeSet::copy;
            }
            BlockProcessor<GFS>::
              addIndex(scalarIndices[i][j], ii, attr, m, c.indexSet());
          }
        }
      c.indexSet().endResize();
      //std::cout<<gv.comm().rank()<<": index set size = "<<c.indexSet().size()<<std::endl;
      //std::cout<<gv.comm().rank()<<": "<<c.indexSet()<<std::endl;

      // Compute neighbours using communication
      typedef NeighbourGatherScatter<typename V::ElementType> NeighbourGS;
      std::set<int> neighbours;
      Dune::PDELab::GenericDataHandle<GFS,V,NeighbourGS> gdhn(gfs, v, NeighbourGS(gfs.gridview().comm().rank(),
                                                                                  neighbours));
      if (gfs.gridview().comm().size()>1)
        gfs.gridview().communicate(gdhn,Dune::All_All_Interface,Dune::ForwardCommunication);
      c.remoteIndices().setNeighbours(neighbours);
      //std::cout<<gv.comm().rank()<<": no neighbours="<<neighbours.size()<<std::endl;

      c.remoteIndices().template rebuild<false>();
      //std::cout<<c.remoteIndices()<<std::endl;
    }
#endif

    //========================================================
    // Generic support for nonoverlapping grids
    //========================================================

    //! Operator for the non-overlapping parallel case
    /**
     * Calculate \f$y:=Ax\f$.
     *
     * \tparam GFS The GridFunctionSpace the vectors apply to.
     * \tparam M   Type of the matrix.  Should be one of the ISTL matrix types.
     * \tparam X   Type of the vectors the matrix is applied to.
     * \tparam Y   Type of the result vectors.
     */
    template<class GFS, class M, class X, class Y>
    class NonoverlappingOperator : public Dune::AssembledLinearOperator<M,X,Y>
    {
    public:
      //! export type of matrix
      typedef M matrix_type;
      //! export type of vectors the matrix is applied to
      typedef X domain_type;
      //! export type of result vectors
      typedef Y range_type;
      //! export type of the entries for x
      typedef typename X::field_type field_type;

      //redefine the category, that is the only difference
      enum {category=Dune::SolverCategory::nonoverlapping};

      //! Construct a non-overlapping operator
      /**
       * \param gfs_    GridFunctionsSpace for the vectors.
       * \param A       Matrix for this operator.  This should be the locally
       *                assembled matrix.
       * \param helper_ Helper for parallel communication (not used).
       *
       * \note The constructed object stores references to all the objects
       *       given as parameters here.  They should be valid for as long as
       *       the constructed object is used.  They are not needed to
       *       destruct the constructed object.
       *
       * \deprecated The helper_ parameter is unused.  Use the constructor
       *             without the helper_ parameter instead.
       */
      NonoverlappingOperator (const GFS& gfs_, const M& A,
                              const ParallelISTLHelper<GFS>& helper_)
        DUNE_DEPRECATED
        : gfs(gfs_), _A_(A)
      {
      }

      //! Construct a non-overlapping operator
      /**
       * \param gfs_ GridFunctionsSpace for the vectors.
       * \param A    Matrix for this operator.  This should be the locally
       *             assembled matrix.
       *
       * \note The constructed object stores references to all the objects
       *       given as parameters here.  They should be valid for as long as
       *       the constructed object is used.  They are not needed to
       *       destruct the constructed object.
       */
      NonoverlappingOperator (const GFS& gfs_, const M& A)
        : gfs(gfs_), _A_(A)
      { }

      //! apply operator
      /**
       * Compute \f$y:=A(x)\f$ on this process, then make y consistent (sum up
       * corresponding entries of y on the different processes and store the
       * result back in y on each process).
       */
      virtual void apply (const X& x, Y& y) const
      {
        // apply local operator; now we have sum y_p = sequential y
        _A_.mv(x,y);

        // accumulate y on border
        Dune::PDELab::AddDataHandle<GFS,Y> adddh(gfs,y);
        if (gfs.gridview().comm().size()>1)
          gfs.gridview().communicate(adddh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
      }

      //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
      /**
       * Compute \f$y:=\alpha A(x)\f$ on this process, then make y consistent
       * (sum up corresponding entries of y on the different processes and
       * store the result back in y on each process).
       */
      virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
      {
        // apply local operator; now we have sum y_p = sequential y
        _A_.usmv(alpha,x,y);

        // accumulate y on border
        Dune::PDELab::AddDataHandle<GFS,Y> adddh(gfs,y);
        if (gfs.gridview().comm().size()>1)
          gfs.gridview().communicate(adddh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
      }

      //! extract the matrix
      virtual const M& getmat () const
      {
        return _A_;
      }

    private:
      const GFS& gfs;
      const M& _A_;
    };


    // parallel scalar product assuming no overlap
    template<class GFS, class X>
    class NonoverlappingScalarProduct : public Dune::ScalarProduct<X>
    {
    public:
      //! export types
      typedef X domain_type;
      typedef typename X::ElementType field_type;

      //! define the category
      enum {category=Dune::SolverCategory::nonoverlapping};

      /*! \brief Constructor needs to know the grid function space
       */
      NonoverlappingScalarProduct (const GFS& gfs_, const ParallelISTLHelper<GFS>& helper_)
        : gfs(gfs_), helper(helper_)
      {}

      /*! \brief Dot product of two vectors.
        It is assumed that the vectors are consistent on the interior+border
        partition.
      */
      virtual field_type dot (const X& x, const X& y)
      {
        // do local scalar product on unique partition
        field_type sum = 0;
        for (typename X::size_type i=0; i<x.N(); ++i)
          for (typename X::size_type j=0; j<x[i].N(); ++j)
            sum += (x[i][j]*y[i][j])*helper.mask(i,j);

        // do global communication
        return gfs.gridview().comm().sum(sum);
      }

      /*! \brief Norm of a right-hand side vector.
        The vector must be consistent on the interior+border partition
      */
      virtual double norm (const X& x)
      {
        return sqrt(static_cast<double>(this->dot(x,x)));
      }

      /*! \brief make additive vector consistent
       */
      void make_consistent (X& x) const
      {
        Dune::PDELab::AddDataHandle<GFS,X> adddh(gfs,x);
        if (gfs.gridview().comm().size()>1)
          gfs.gridview().communicate(adddh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
      }

    private:
      const GFS& gfs;
      const ParallelISTLHelper<GFS>& helper;
    };

    // parallel Richardson preconditioner
    template<class GFS, class X, class Y>
    class NonoverlappingRichardson : public Dune::Preconditioner<X,Y>
    {
    public:
      //! \brief The domain type of the preconditioner.
      typedef X domain_type;
      //! \brief The range type of the preconditioner.
      typedef Y range_type;
      //! \brief The field type of the preconditioner.
      typedef typename X::ElementType field_type;

      // define the category
      enum {
        //! \brief The category the preconditioner is part of.
        category=Dune::SolverCategory::nonoverlapping
      };

      //! \brief Constructor.
      NonoverlappingRichardson (const GFS& gfs_, const ParallelISTLHelper<GFS>& helper_)
        : gfs(gfs_), helper(helper_)
      {
      }

      /*!
        \brief Prepare the preconditioner.

        \copydoc Preconditioner::pre(X&,Y&)
      */
      virtual void pre (X& x, Y& b) {}

      /*!
        \brief Apply the precondioner.

        \copydoc Preconditioner::apply(X&,const Y&)
      */
      virtual void apply (X& v, const Y& d)
      {
        v = d;
        // no communication is necessary here because defect is already consistent!
//         helper.mask(v);
//         Dune::PDELab::AddDataHandle<GFS,X> adddh(gfs,v);
//         if (gfs.gridview().comm().size()>1)
//           gfs.gridview().communicate(adddh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
      }

      /*!
        \brief Clean up.

        \copydoc Preconditioner::post(X&)
      */
      virtual void post (X& x) {}

    private:
      const GFS& gfs;
      const ParallelISTLHelper<GFS>& helper;
    };


    //========================================================
    // Generic support for overlapping grids
    // (need to be used with appropriate constraints)
    //========================================================

    // operator that resets result to zero at constrained DOFS
    template<class CC, class M, class X, class Y>
    class OverlappingOperator : public Dune::AssembledLinearOperator<M,X,Y>
    {
    public:
      //! export types
      typedef M matrix_type;
      typedef X domain_type;
      typedef Y range_type;
      typedef typename X::field_type field_type;

      //redefine the category, that is the only difference
      enum {category=Dune::SolverCategory::overlapping};

      OverlappingOperator (const CC& cc_, const M& A)
        : cc(cc_), _A_(A)
      {}

      //! apply operator to x:  \f$ y = A(x) \f$
      virtual void apply (const X& x, Y& y) const
      {
        _A_.mv(x,y);
        Dune::PDELab::set_constrained_dofs(cc,0.0,y);
      }

      //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
      virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
      {
        _A_.usmv(alpha,x,y);
        Dune::PDELab::set_constrained_dofs(cc,0.0,y);
      }

      //! get matrix via *
      virtual const M& getmat () const
      {
        return _A_;
      }

    private:
      const CC& cc;
      const M& _A_;
    };

    // new scalar product assuming at least overlap 1
    // uses unique partitioning of nodes for parallelization
    template<class GFS, class X>
    class OverlappingScalarProduct : public Dune::ScalarProduct<X>
    {
    public:
      //! export types
      typedef X domain_type;
      typedef typename X::ElementType field_type;

      //! define the category
      enum {category=Dune::SolverCategory::overlapping};

      /*! \brief Constructor needs to know the grid function space
       */
      OverlappingScalarProduct (const GFS& gfs_, const ParallelISTLHelper<GFS>& helper_)
        : gfs(gfs_), helper(helper_)
      {}


      /*! \brief Dot product of two vectors.
        It is assumed that the vectors are consistent on the interior+border
        partition.
      */
      virtual field_type dot (const X& x, const X& y)
      {
        // do local scalar product on unique partition
        field_type sum = 0;
        for (typename X::size_type i=0; i<x.N(); ++i)
          for (typename X::size_type j=0; j<x[i].N(); ++j)
            sum += (x[i][j]*y[i][j])*helper.mask(i,j);

        // do global communication
        return gfs.gridview().comm().sum(sum);
      }

      /*! \brief Norm of a right-hand side vector.
        The vector must be consistent on the interior+border partition
      */
      virtual double norm (const X& x)
      {
        return sqrt(static_cast<double>(this->dot(x,x)));
      }

    private:
      const GFS& gfs;
      const ParallelISTLHelper<GFS>& helper;
    };

    // wrapped sequential preconditioner
    template<class CC, class GFS, class P>
    class OverlappingWrappedPreconditioner
      : public Dune::Preconditioner<typename P::domain_type,typename P::range_type>
    {
    public:
      //! \brief The domain type of the preconditioner.
      typedef typename P::domain_type domain_type;
      //! \brief The range type of the preconditioner.
      typedef typename P::range_type range_type;

      // define the category
      enum {
        //! \brief The category the preconditioner is part of.
        category=Dune::SolverCategory::overlapping
      };

      //! Constructor.
      OverlappingWrappedPreconditioner (const GFS& gfs_, P& prec_, const CC& cc_,
                                        const ParallelISTLHelper<GFS>& helper_)
        : gfs(gfs_), prec(prec_), cc(cc_), helper(helper_)
      {}

      /*!
        \brief Prepare the preconditioner.

        \copydoc Preconditioner::pre(domain_type&,range_type&)
      */
      virtual void pre (domain_type& x, range_type& b)
      {
        prec.pre(x,b);
      }

      /*!
        \brief Apply the precondioner.

        \copydoc Preconditioner::apply(domain_type&,const range_type&)
      */
      virtual void apply (domain_type& v, const range_type& d)
      {
        range_type dd(d);
        set_constrained_dofs(cc,0.0,dd);
        prec.apply(v,dd);
        Dune::PDELab::AddDataHandle<GFS,domain_type> adddh(gfs,v);
        if (gfs.gridview().comm().size()>1)
          gfs.gridview().communicate(adddh,Dune::All_All_Interface,Dune::ForwardCommunication);
      }

      /*!
        \brief Clean up.

        \copydoc Preconditioner::post(domain_type&)
      */
      virtual void post (domain_type& x)
      {
        prec.post(x);
      }

    private:
      const GFS& gfs;
      P& prec;
      const CC& cc;
      const ParallelISTLHelper<GFS>& helper;
    };


#if HAVE_SUPERLU
    // exact subdomain solves with SuperLU as preconditioner
    template<class GFS, class M, class X, class Y>
    class SuperLUSubdomainSolver : public Dune::Preconditioner<X,Y>
    {
      typedef typename M::BaseT ISTLM;

    public:
      //! \brief The domain type of the preconditioner.
      typedef X domain_type;
      //! \brief The range type of the preconditioner.
      typedef Y range_type;
      //! \brief The field type of the preconditioner.
      typedef typename X::ElementType field_type;


      // define the category
      enum {
        //! \brief The category the preconditioner is part of.
        category=Dune::SolverCategory::overlapping
      };

      /*! \brief Constructor.

        Constructor gets all parameters to operate the prec.
        \param A The matrix to operate on.
        \param n The number of iterations to perform.
        \param w The relaxation factor.
      */
      SuperLUSubdomainSolver (const GFS& gfs_, const M& A_)
        : gfs(gfs_), A(A_), solver(A_,false) // this does the decomposition
      {}

      /*!
        \brief Prepare the preconditioner.

        \copydoc Preconditioner::pre(X&,Y&)
      */
      virtual void pre (X& x, Y& b) {}

      /*!
        \brief Apply the precondioner.

        \copydoc Preconditioner::apply(X&,const Y&)
      */
      virtual void apply (X& v, const Y& d)
      {
        Dune::InverseOperatorResult stat;
        Y b(d); // need copy, since solver overwrites right hand side
        solver.apply(v,b,stat);
        Dune::PDELab::AddDataHandle<GFS,X> adddh(gfs,v);
        if (gfs.gridview().comm().size()>1)
          gfs.gridview().communicate(adddh,Dune::All_All_Interface,Dune::ForwardCommunication);
      }

      /*!
        \brief Clean up.

        \copydoc Preconditioner::post(X&)
      */
      virtual void post (X& x) {}

    private:
      const GFS& gfs;
      const M& A;
      Dune::SuperLU<ISTLM> solver;
    };

    // exact subdomain solves with SuperLU as preconditioner
    template<class GFS, class M, class X, class Y>
    class RestrictedSuperLUSubdomainSolver : public Dune::Preconditioner<X,Y>
    {
      typedef typename M::BaseT ISTLM;

    public:
      //! \brief The domain type of the preconditioner.
      typedef X domain_type;
      //! \brief The range type of the preconditioner.
      typedef Y range_type;
      //! \brief The field type of the preconditioner.
      typedef typename X::ElementType field_type;


      // define the category
      enum {
        //! \brief The category the preconditioner is part of.
        category=Dune::SolverCategory::overlapping
      };

      /*! \brief Constructor.

        Constructor gets all parameters to operate the prec.
        \param A The matrix to operate on.
        \param n The number of iterations to perform.
        \param w The relaxation factor.
      */
      RestrictedSuperLUSubdomainSolver (const GFS& gfs_, const M& A_,
                                        const ParallelISTLHelper<GFS>& helper_)
        : gfs(gfs_), A(A_), solver(A_,false), helper(helper_) // this does the decomposition
      {}

      /*!
        \brief Prepare the preconditioner.

        \copydoc Preconditioner::pre(X&,Y&)
      */
      virtual void pre (X& x, Y& b) {}

      /*!
        \brief Apply the precondioner.

        \copydoc Preconditioner::apply(X&,const Y&)
      */
      virtual void apply (X& v, const Y& d)
      {
        Dune::InverseOperatorResult stat;
        Y b(d); // need copy, since solver overwrites right hand side
        solver.apply(v,b,stat);
        helper.mask(v);
        Dune::PDELab::AddDataHandle<GFS,X> adddh(gfs,v);
        if (gfs.gridview().comm().size()>1)
          gfs.gridview().communicate(adddh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
      }

      /*!
        \brief Clean up.

        \copydoc Preconditioner::post(X&)
      */
      virtual void post (X& x) {}

    private:
      const GFS& gfs;
      const M& A;
      Dune::SuperLU<ISTLM> solver;
      const ParallelISTLHelper<GFS>& helper;
    };
#endif



    //==============================================================================
    // Here we add some standard linear solvers conforming to the linear solver
    // interface required to solve linear and nonlinear problems.
    //==============================================================================

    struct SequentialNorm
    {/*! \brief compute global norm of a vector

        \param[in] v the given vector
      */
      template<class V>
      typename V::ElementType norm(const V& v) const
      {
        return v.two_norm();
      }
    };

    class LinearResultStorage
    {
    public:
      /*! \brief Return access to result data */
      const Dune::PDELab::LinearSolverResult<double>& result() const
      {
        return res;
      }

    protected:
      Dune::PDELab::LinearSolverResult<double> res;
    };
    

    template<template<class,class,class,int> class Preconditioner,
             template<class> class Solver>
    class ISTLBackend_SEQ_Base
      : public SequentialNorm, public LinearResultStorage
    {
    public:
      /*! \brief make a linear solver object

        \param[in] maxiter maximum number of iterations to do
        \param[in] verbose print messages if true
      */
      explicit ISTLBackend_SEQ_Base(unsigned maxiter_=5000, bool verbose_=true)
        : maxiter(maxiter_), verbose(verbose_)
      {}
      
      

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename W::ElementType reduction)
      {
        Dune::MatrixAdapter<M,V,W> opa(A);
        Preconditioner<M,V,W,1> ssor(A, 3, 1.0);
        Solver<V> solver(opa, ssor, reduction, maxiter, verbose);
        Dune::InverseOperatorResult stat;
        solver.apply(z, r, stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
      }

    private:
      unsigned maxiter;
      bool verbose;
    };
    
    template<template<typename> class Solver>
    class ISTLBackend_SEQ_ILU0 
      :  public SequentialNorm, public LinearResultStorage
    {
    public:
      explicit ISTLBackend_SEQ_ILU0 (unsigned maxiter_=5000, bool verbose_=true)
        : maxiter(maxiter_), verbose(verbose_)
       {}
      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename W::ElementType reduction)
      {
        Dune::MatrixAdapter<M,V,W> opa(A);
        Dune::SeqILU0<M,V,W> ilu0(A, 1.0);
        Solver<V> solver(opa, ilu0, reduction, maxiter, verbose);
        Dune::InverseOperatorResult stat;
        solver.apply(z, r, stat);
        res.converged  = stat.converged;
        res.reduction  = stat.reduction;
       }
    private:
      unsigned maxiter;
      bool verbose;
    };
    

    class ISTLBackend_SEQ_BCGS_SSOR
      : public ISTLBackend_SEQ_Base<Dune::SeqSSOR, Dune::BiCGSTABSolver>
    {
    public:
      /*! \brief make a linear solver object

        \param[in] maxiter maximum number of iterations to do
        \param[in] verbose print messages if true
      */
      explicit ISTLBackend_SEQ_BCGS_SSOR (unsigned maxiter_=5000, bool verbose_=true)
        : ISTLBackend_SEQ_Base<Dune::SeqSSOR, Dune::BiCGSTABSolver>(maxiter_, verbose_)
      {}
    };

    class ISTLBackend_SEQ_BCGS_ILU0
      : public ISTLBackend_SEQ_ILU0<Dune::BiCGSTABSolver>
    {
    public:
      /*! \brief make a linear solver object

        \param[in] maxiter maximum number of iterations to do
        \param[in] verbose print messages if true
      */
      explicit ISTLBackend_SEQ_BCGS_ILU0 (unsigned maxiter_=5000, bool verbose_=true)
        : ISTLBackend_SEQ_ILU0<Dune::BiCGSTABSolver>(maxiter_, verbose_)
      {}
    };
    

    class ISTLBackend_SEQ_CG_ILU0
      : public ISTLBackend_SEQ_ILU0<Dune::CGSolver>
    {
    public:
      /*! \brief make a linear solver object

        \param[in] maxiter maximum number of iterations to do
        \param[in] verbose print messages if true
      */
      explicit ISTLBackend_SEQ_CG_ILU0 (unsigned maxiter_=5000, bool verbose_=true)
        : ISTLBackend_SEQ_ILU0<Dune::CGSolver>(maxiter_, verbose_)
      {}
    };

    class ISTLBackend_SEQ_CG_SSOR
      : public ISTLBackend_SEQ_Base<Dune::SeqSSOR, Dune::CGSolver>
    {
    public:
      /*! \brief make a linear solver object

        \param[in] maxiter maximum number of iterations to do
        \param[in] verbose print messages if true
      */
      explicit ISTLBackend_SEQ_CG_SSOR (unsigned maxiter_=5000, bool verbose_=true)
        : ISTLBackend_SEQ_Base<Dune::SeqSSOR, Dune::CGSolver>(maxiter_, verbose_)
      {}
    };
    

    class ISTLBackend_SEQ_SuperLU
      : public SequentialNorm, public LinearResultStorage
    {
    public:
      /*! \brief make a linear solver object

        \param[in] maxiter maximum number of iterations to do
        \param[in] verbose print messages if true
      */
      explicit ISTLBackend_SEQ_SuperLU (bool verbose_=true)
        : verbose(verbose_)
      {}

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename W::ElementType reduction)
      {
#if HAVE_SUPERLU
        typedef typename M::BaseT ISTLM;
        Dune::SuperLU<ISTLM> solver(A, verbose);
        Dune::InverseOperatorResult stat;
        solver.apply(z, r, stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
#else
        std::cout << "No superLU support, please install and configure it." << std::endl;
#endif
      }

    private:
      bool verbose;
    };

    //! Solver to be used for explicit time-steppers with (block-)diagonal mass matrix
    class ISTLBackend_SEQ_ExplicitDiagonal
      : public SequentialNorm, public LinearResultStorage
    {
    public:
      /*! \brief make a linear solver object

        \param[in] maxiter maximum number of iterations to do
        \param[in] verbose print messages if true
      */
      explicit ISTLBackend_SEQ_ExplicitDiagonal ()
      {}

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename W::ElementType reduction)
      {
        Dune::SeqJac<M,V,W> jac(A,1,1.0);
        jac.pre(z,r);
        jac.apply(z,r);
        jac.post(z);
        res.converged  = true;
        res.iterations = 1;
        res.elapsed    = 0.0;
        res.reduction  = reduction;
      }
    };

    template<class GFS>
    class ISTLBackend_NOVLP_BCGS_NOPREC
    {
      typedef Dune::PDELab::ParallelISTLHelper<GFS> PHELPER;

    public:
      /*! \brief make a linear solver object

        \param[in] gfs a grid function space
        \param[in] maxiter maximum number of iterations to do
        \param[in] verbose print messages if true
      */
      explicit ISTLBackend_NOVLP_BCGS_NOPREC (const GFS& gfs_, unsigned maxiter_=5000, int verbose_=1)
        : gfs(gfs_), phelper(gfs), maxiter(maxiter_), verbose(verbose_)
      {}

      /*! \brief compute global norm of a vector

        \param[in] v the given vector
      */
      template<class V>
      typename V::ElementType norm (const V& v) const
      {
        V x(v); // make a copy because it has to be made consistent
        typedef Dune::PDELab::NonoverlappingScalarProduct<GFS,V> PSP;
        PSP psp(gfs,phelper);
        psp.make_consistent(x);
        return psp.norm(x);
      }

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename V::ElementType reduction)
      {
        typedef Dune::PDELab::NonoverlappingOperator<GFS,M,V,W> POP;
        POP pop(gfs,A);
        typedef Dune::PDELab::NonoverlappingScalarProduct<GFS,V> PSP;
        PSP psp(gfs,phelper);
        typedef Dune::PDELab::NonoverlappingRichardson<GFS,V,W> PRICH;
        PRICH prich(gfs,phelper);
        int verb=0;
        if (gfs.gridview().comm().rank()==0) verb=verbose;
        Dune::BiCGSTABSolver<V> solver(pop,psp,prich,reduction,maxiter,verb);
        Dune::InverseOperatorResult stat;
        solver.apply(z,r,stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
      }

      /*! \brief Return access to result data */
      const Dune::PDELab::LinearSolverResult<double>& result() const
      {
        return res;
      }

    private:
      const GFS& gfs;
      PHELPER phelper;
      Dune::PDELab::LinearSolverResult<double> res;
      unsigned maxiter;
      int verbose;
    };

    //! Solver to be used for explicit time-steppers with (block-)diagonal mass matrix
    template<typename GFS>
    class ISTLBackend_NOVLP_ExplicitDiagonal
    {
      typedef Dune::PDELab::ParallelISTLHelper<GFS> PHELPER;

      const GFS& gfs;
      PHELPER phelper;
      Dune::PDELab::LinearSolverResult<double> res;

    public:
      /*! \brief make a linear solver object

        \param[in] gfs GridFunctionSpace, used to identify DoFs for parallel
        communication
      */
      explicit ISTLBackend_NOVLP_ExplicitDiagonal(const GFS& gfs_)
        : gfs(gfs_), phelper(gfs)
      {}

      /*! \brief compute global norm of a vector

        \param[in] v the given vector
      */
      template<class V>
      typename V::ElementType norm (const V& v) const
      {
        dune_static_assert
          (AlwaysFalse<V>::value,
           "ISTLBackend_NOVLP_ExplicitDiagonal::norm() should not be "
           "neccessary, so we skipped the testing.  If you have a scenario "
           "where you need it, please verify the implementation and remove "
           "this assert or report back to us so we can remove it.");

        typedef Dune::PDELab::NonoverlappingScalarProduct<GFS,V> PSP;
        V x(v); // make a copy because it has to be made consistent
        PSP psp(gfs,phelper);
        psp.make_consistent(x);
        return psp.norm(x);
      }

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename W::ElementType reduction)
      {
        Dune::SeqJac<M,V,W> jac(A,1,1.0);
        jac.pre(z,r);
        jac.apply(z,r);
        jac.post(z);
        if (gfs.gridview().comm().size()>1)
        {
          Dune::PDELab::AddDataHandle<GFS,V> adddh(gfs,z);
          gfs.gridview().communicate(adddh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
        }
        res.converged  = true;
        res.iterations = 1;
        res.elapsed    = 0.0;
        res.reduction  = reduction;
      }

      /*! \brief Return access to result data */
      const Dune::PDELab::LinearSolverResult<double>& result() const
      {
        return res;
      }
    };

    template<typename GFS>
    class OVLPScalarProductImplementation
    {
    public:
      OVLPScalarProductImplementation(const GFS& gfs_)
        : gfs(gfs_), helper(gfs_)
      {}

      /*! \brief Dot product of two vectors.
        It is assumed that the vectors are consistent on the interior+border
        partition.
      */
      template<typename X>
      typename X::ElementType dot (const X& x, const X& y) const
      {
        // do local scalar product on unique partition
        typename X::ElementType sum = 0;
        for (typename X::size_type i=0; i<x.N(); ++i)
          for (typename X::size_type j=0; j<x[i].N(); ++j)
            sum += (x[i][j]*y[i][j])*helper.mask(i,j);

        // do global communication
        return gfs.gridview().comm().sum(sum);
      }

      /*! \brief Norm of a right-hand side vector.
        The vector must be consistent on the interior+border partition
      */
       template<typename X>
      typename X::ElementType norm (const X& x) const
      {
        return sqrt(static_cast<double>(this->dot(x,x)));
      }

      const  ParallelISTLHelper<GFS>& parallelHelper()
      {
        return helper;
      }
      
    private:
      const GFS& gfs;
      ParallelISTLHelper<GFS> helper;
    };
    

    template<typename GFS, typename X>
    class OVLPScalarProduct
      : public ScalarProduct<X>
    {
    public:
      enum {category=Dune::SolverCategory::overlapping};
      OVLPScalarProduct(const OVLPScalarProductImplementation<GFS>& implementation_)
        : implementation(implementation_)
      {}
      virtual typename X::ElementType dot(const X& x, const X& y)
      {
        return implementation.dot(x,y);
      }
      
       virtual typename X::ElementType norm (const X& x)
      {
        return sqrt(static_cast<double>(this->dot(x,x)));
      }

    private:
      const OVLPScalarProductImplementation<GFS>& implementation;
    };
    
    template<class GFS, class C,
             template<class,class,class,int> class Preconditioner,
             template<class> class Solver>
    class ISTLBackend_OVLP_Base
      : public OVLPScalarProductImplementation<GFS>, public LinearResultStorage
    {
    public:
      /*! \brief make a linear solver object

        \param[in] gfs a grid function space
        \param[in] c a constraints object
        \param[in] maxiter maximum number of iterations to do
        \param[in] kssor number of SSOR steps to apply as inner iteration
        \param[in] verbose print messages if true
      */
      explicit ISTLBackend_OVLP_Base (const GFS& gfs_, const C& c_, unsigned maxiter_=5000,
                                            int steps_=5, int verbose_=1)
        : OVLPScalarProductImplementation<GFS>(gfs_), gfs(gfs_), c(c_), maxiter(maxiter_), steps(steps_), verbose(verbose_)
      {}
      
      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename V::ElementType reduction)
      {
        typedef Dune::PDELab::OverlappingOperator<C,M,V,W> POP;
        POP pop(c,A);
        typedef OVLPScalarProduct<GFS,V> PSP;
        PSP psp(*this);
        typedef Preconditioner<M,V,W,1> SeqPrec;
        SeqPrec seqprec(A,steps,1.0);
        typedef Dune::PDELab::OverlappingWrappedPreconditioner<C,GFS,SeqPrec> WPREC;
        WPREC wprec(gfs,seqprec,c,this->parallelHelper());
        int verb=0;
        if (gfs.gridview().comm().rank()==0) verb=verbose;
        Solver<V> solver(pop,psp,wprec,reduction,maxiter,verb);
        Dune::InverseOperatorResult stat;
        solver.apply(z,r,stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
      }
    private:
      const GFS& gfs;
      const C& c;
      unsigned maxiter;
      int steps;
      bool verbose;
    };
    

    template<class GFS, class C>
    class ISTLBackend_OVLP_BCGS_SSORk
      : public ISTLBackend_OVLP_Base<GFS,C,Dune::SeqSSOR, Dune::BiCGSTABSolver>
    {
      typedef Dune::PDELab::ParallelISTLHelper<GFS> PHELPER;

    public:
      /*! \brief make a linear solver object

        \param[in] gfs a grid function space
        \param[in] c a constraints object
        \param[in] maxiter maximum number of iterations to do
        \param[in] steps number of SSOR steps to apply as inner iteration
        \param[in] verbose print messages if true
      */
      explicit ISTLBackend_OVLP_BCGS_SSORk (const GFS& gfs, const C& c, unsigned maxiter=5000,
                                            int steps=5, int verbose=1)
        : ISTLBackend_OVLP_Base<GFS,C,Dune::SeqSSOR, Dune::BiCGSTABSolver>(gfs, c, maxiter, steps, verbose)
      {}
    };

    template<class GFS, class C>
    class ISTLBackend_OVLP_CG_SSORk
      : public ISTLBackend_OVLP_Base<GFS,C,Dune::SeqSSOR, Dune::CGSolver>
    {
      typedef Dune::PDELab::ParallelISTLHelper<GFS> PHELPER;

    public:
      /*! \brief make a linear solver object

        \param[in] gfs a grid function space
        \param[in] c a constraints object
        \param[in] maxiter maximum number of iterations to do
        \param[in] steps number of SSOR steps to apply as inner iteration
        \param[in] verbose print messages if true
      */
      explicit ISTLBackend_OVLP_CG_SSORk (const GFS& gfs, const C& c, unsigned maxiter=5000,
                                            int steps=5, int verbose=1)
        : ISTLBackend_OVLP_Base<GFS,C,Dune::SeqSSOR, Dune::CGSolver>(gfs, c, maxiter, steps, verbose)
      {}
    };

    template<class GFS, class C, template<typename> class Solver>
    class ISTLBackend_OVLP_SuperLU_Base
      : public OVLPScalarProductImplementation<GFS>, public LinearResultStorage
    {
      typedef Dune::PDELab::ParallelISTLHelper<GFS> PHELPER;

    public:
      /*! \brief make a linear solver object

        \param[in] gfs a grid function space
        \param[in] c a constraints object
        \param[in] maxiter maximum number of iterations to do
        \param[in] kssor number of SSOR steps to apply as inner iteration
        \param[in] verbose print messages if true
      */
      explicit ISTLBackend_OVLP_SuperLU_Base (const GFS& gfs_, const C& c_, unsigned maxiter_=5000,
                                              int verbose_=1)
        : OVLPScalarProductImplementation<GFS>(gfs_), gfs(gfs_), c(c_), maxiter(maxiter_), verbose(verbose_)
      {}

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename V::ElementType reduction)
      {
        typedef Dune::PDELab::OverlappingOperator<C,M,V,W> POP;
        POP pop(c,A);
        typedef OVLPScalarProduct<GFS,V> PSP;
        PSP psp(*this);
#if HAVE_SUPERLU
        typedef Dune::PDELab::SuperLUSubdomainSolver<GFS,M,V,W> PREC;
        PREC prec(gfs,A);
        int verb=0;
        if (gfs.gridview().comm().rank()==0) verb=verbose;
        Solver<V> solver(pop,psp,prec,reduction,maxiter,verbose);
        Dune::InverseOperatorResult stat;
        solver.apply(z,r,stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
#else
        std::cout << "No superLU support, please install and configure it." << std::endl;
#endif
      }

    private:
      const GFS& gfs;
      const C& c;
      unsigned maxiter;
      int verbose;
    };

    template<class GFS, class C>
    class ISTLBackend_OVLP_BCGS_SuperLU
      : public ISTLBackend_OVLP_SuperLU_Base<GFS,C,Dune::BiCGSTABSolver>
    {
    public:
      
      /*! \brief make a linear solver object

        \param[in] gfs a grid function space
        \param[in] c a constraints object
        \param[in] maxiter maximum number of iterations to do
        \param[in] kssor number of SSOR steps to apply as inner iteration
        \param[in] verbose print messages if true
      */
      explicit ISTLBackend_OVLP_BCGS_SuperLU (const GFS& gfs_, const C& c_, unsigned maxiter_=5000,
                                              int verbose_=1)
        : ISTLBackend_OVLP_SuperLU_Base<GFS,C,Dune::BiCGSTABSolver>(gfs_,c_,maxiter_,verbose_)
      {}
    };
    
    template<class GFS, class C>
    class ISTLBackend_OVLP_CG_SuperLU
      : public ISTLBackend_OVLP_SuperLU_Base<GFS,C,Dune::CGSolver>
    {
    public:
      
      /*! \brief make a linear solver object

        \param[in] gfs a grid function space
        \param[in] c a constraints object
        \param[in] maxiter maximum number of iterations to do
        \param[in] kssor number of SSOR steps to apply as inner iteration
        \param[in] verbose print messages if true
      */
      explicit ISTLBackend_OVLP_CG_SuperLU (const GFS& gfs_, const C& c_, 
                                              unsigned maxiter_=5000,
                                              int verbose_=1)
        : ISTLBackend_OVLP_SuperLU_Base<GFS,C,Dune::CGSolver>(gfs_,c_,maxiter_,verbose_)
      {}
    };


    //! Solver to be used for explicit time-steppers with (block-)diagonal mass matrix
    template<class GFS>
    class ISTLBackend_OVLP_ExplicitDiagonal
      : public LinearResultStorage
    {
    public:
      /*! \brief make a linear solver object

        \param[in] maxiter maximum number of iterations to do
        \param[in] verbose print messages if true
      */
      explicit ISTLBackend_OVLP_ExplicitDiagonal (const GFS& gfs_)
        : gfs(gfs_)
      {}

      /*! \brief compute global norm of a vector

        \param[in] v the given vector
      */
      template<class V>
      typename V::ElementType norm(const V& v) const
      {
        dune_static_assert
          (AlwaysFalse<V>::value,
           "ISTLBackend_OVLP_ExplicitDiagonal::norm() should not be "
           "neccessary, so we skipped the implementation.  If you have a "
           "scenario where you need it, please implement it or report back to "
           "us.");
      }

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename W::ElementType reduction)
      {
        Dune::SeqJac<M,V,W> jac(A,1,1.0);
        jac.pre(z,r);
        jac.apply(z,r);
        jac.post(z);
        if (gfs.gridview().comm().size()>1)
        {
          Dune::PDELab::CopyDataHandle<GFS,V> copydh(gfs,z);
          gfs.gridview().communicate(copydh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
        }
        res.converged  = true;
        res.iterations = 1;
        res.elapsed    = 0.0;
        res.reduction  = reduction;
      }

    private:
      const GFS& gfs;
    };


    template<int s, bool isFakeMPIHelper>
    struct CommSelector
    {
      typedef Dune::Amg::SequentialInformation type;
    };

    // Need MPI for OwnerOverlapCopyCommunication
#if HAVE_MPI
    template<int s>
    struct CommSelector<s,false>
    {
      typedef OwnerOverlapCopyCommunication<bigunsignedint<s>,int> type;
    };
#endif

    template<class GFS, int s, template<class,class,class,int> class SMI, template<class> class SOI>
    class ISTLBackend_AMG
    {
      typedef Dune::PDELab::ParallelISTLHelper<GFS> PHELPER;

    public:
      ISTLBackend_AMG(const GFS& gfs_, int smoothsteps=2,
                      unsigned maxiter_=5000, int verbose_=1)
        : gfs(gfs_), phelper(gfs), maxiter(maxiter_), steps(smoothsteps), verbose(verbose_)
      {}


      /*! \brief compute global norm of a vector

        \param[in] v the given vector
      */
      template<class V>
      typename V::ElementType norm (const V& v) const
      {
        typedef Dune::PDELab::OverlappingScalarProduct<GFS,V> PSP;
        PSP psp(gfs,phelper);
        return psp.norm(v);
      }

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V>
      void apply(M& A, V& z, V& r, typename V::ElementType reduction)
      {
        typedef typename M::BaseT MatrixType;
        typedef typename BlockProcessor<GFS>::template AMGVectorTypeSelector<V>::Type
          VectorType;
        typedef typename CommSelector<s,Dune::MPIHelper::isFake>::type Comm;

        Comm oocc(gfs.gridview().comm());
        MatrixType& mat=A.base();
        phelper.createIndexSetAndProjectForAMG(mat, oocc);
        typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<MatrixType,
          Dune::Amg::FirstDiagonal> >
          Criterion;
        typedef SMI<MatrixType,VectorType,VectorType,1> Smoother;
        typedef Dune::BlockPreconditioner<VectorType,VectorType,Comm,Smoother> ParSmoother;
        typedef typename Dune::Amg::SmootherTraits<ParSmoother>::Arguments SmootherArgs;
        typedef Dune::OverlappingSchwarzOperator<MatrixType,VectorType,VectorType,Comm> Operator;
        typedef Dune::Amg::AMG<Operator,VectorType,ParSmoother,Comm> AMG;
        SmootherArgs smootherArgs;
        smootherArgs.iterations = 1;
        smootherArgs.relaxationFactor = 1;

        Criterion criterion(15,2000);
        criterion.setDefaultValuesIsotropic(GFS::Traits::GridViewType::Traits::Grid::dimension);
        criterion.setDebugLevel(verbose);
        Dune::OverlappingSchwarzScalarProduct<VectorType,Comm> sp(oocc);
        Operator oop(mat, oocc);
        //oocc.copyOwnerToAll(BlockProcessor<GFS>::getVector(r), BlockProcessor<GFS>::getVector(r));
        AMG amg=AMG(oop, criterion, smootherArgs, 1, steps, steps, false, oocc);

        Dune::InverseOperatorResult stat;
        int verb=0;
        if (gfs.gridview().comm().rank()==0) verb=verbose;

        SOI<VectorType> solver(oop,sp,amg,reduction,maxiter,verb);
        solver.apply(BlockProcessor<GFS>::getVector(z),BlockProcessor<GFS>::getVector(r),stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        //oocc.copyOwnerToAll(BlockProcessor<GFS>::getVector(z), BlockProcessor<GFS>::getVector(z));
      }

      /*! \brief Return access to result data */
      const Dune::PDELab::LinearSolverResult<double>& result() const
      {
        return res;
      }

    private:
      const GFS& gfs;
      PHELPER phelper;
      Dune::PDELab::LinearSolverResult<double> res;
      unsigned maxiter;
      int steps;
      int verbose;
    };

    /**
     * @brief Parallel cojugate gradient solver preconditioned with AMG smoothed by SSOR
     * @tparam GFS The type of the grid functions space.
     * @tparam s The bits to use for the globale index.
     */
    template<class GFS, int s=96>
    class ISTLBackend_CG_AMG_SSOR
      : public ISTLBackend_AMG<GFS, s, Dune::SeqSSOR, Dune::CGSolver>
    {

    public:
      /**
       * @brief Constructor
       * @param gfs_ The grid function space used.
       * @param smoothsteps The number of steps to use for both pre and post smoothing.
       * @param maxiter_ The maximum number of iterations allowed.
       * @param verbose_ The verbosity level to use.
       */
      ISTLBackend_CG_AMG_SSOR(const GFS& gfs_,int smoothsteps=2,
                              unsigned maxiter_=5000, int verbose_=1)
        : ISTLBackend_AMG<GFS, s, Dune::SeqSSOR, Dune::CGSolver>(gfs_,smoothsteps, maxiter_,verbose_)
      {}
    };

    /**
     * @brief Parallel BiCGStab solver preconditioned with AMG smoothed by SSOR
     * @tparam GFS The type of the grid functions space.
     * @tparam s The bits to use for the globale index.
     */
    template<class GFS, int s=96>
    class ISTLBackend_BCGS_AMG_SSOR
      : public ISTLBackend_AMG<GFS, s, Dune::SeqSSOR, Dune::BiCGSTABSolver>
    {

    public:
      /**
       * @brief Constructor
       * @param gfs_ The grid function space used.
       * @param smoothsteps The number of steps to use for both pre and post smoothing.
       * @param maxiter_ The maximum number of iterations allowed.
       * @param verbose_ The verbosity level to use.
       */
      ISTLBackend_BCGS_AMG_SSOR(const GFS& gfs_, int smoothsteps=2,
                                unsigned maxiter_=5000, int verbose_=1)
        : ISTLBackend_AMG<GFS, s, Dune::SeqSSOR, Dune::BiCGSTABSolver>(gfs_,smoothsteps, maxiter_,verbose_)
      {}
    };

    template<class GFS, template<class,class,class,int> class SMI, template<class> class SOI>
    class ISTLBackend_SEQ_AMG
    {

    public:
      ISTLBackend_SEQ_AMG(int smoothsteps=2,
                          unsigned maxiter_=5000, int verbose_=1)
        : maxiter(maxiter_), steps(smoothsteps), verbose(verbose_)
      {}


      /*! \brief compute global norm of a vector

        \param[in] v the given vector
      */
      template<class V>
      typename V::ElementType norm (const V& v) const
      {
        return v.two_norm();
      }

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V>
      void apply(M& A, V& z, V& r, typename V::ElementType reduction)
      {
        typedef typename M::BaseT MatrixType;
        typedef typename BlockProcessor<GFS>::template AMGVectorTypeSelector<V>::Type
          VectorType;
        MatrixType& mat=A.base();
        typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<MatrixType,
          Dune::Amg::FirstDiagonal> >
          Criterion;
        typedef SMI<MatrixType,VectorType,VectorType,1> Smoother;
        typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments
          SmootherArgs;
        typedef Dune::MatrixAdapter<MatrixType,VectorType,VectorType> Operator;
        typedef Dune::Amg::AMG<Operator,VectorType,Smoother> AMG;
        SmootherArgs smootherArgs;
        smootherArgs.iterations = 1;
        smootherArgs.relaxationFactor = 1;

        Criterion criterion(15,2000);
        criterion.setDebugLevel(verbose?3:0);
        Operator oop(mat);
        AMG amg=AMG(oop, criterion, smootherArgs, 1, steps, steps);

        Dune::InverseOperatorResult stat;

        SOI<VectorType> solver(oop,amg,reduction,maxiter,verbose);
        solver.apply(BlockProcessor<GFS>::getVector(z),BlockProcessor<GFS>::getVector(r),stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
      }

      /*! \brief Return access to result data */
      const Dune::PDELab::LinearSolverResult<double>& result() const
      {
        return res;
      }

    private:
      Dune::PDELab::LinearSolverResult<double> res;
      unsigned maxiter;
      int steps;
      int verbose;
    };

    /**
     * @brief Sequential conjugate gradient solver preconditioned with AMG smoothed by SSOR
     * @tparam GFS The type of the grid functions space.
     */
    template<class GFS>
    class ISTLBackend_SEQ_CG_AMG_SSOR
      : public ISTLBackend_SEQ_AMG<GFS, Dune::SeqSSOR, Dune::CGSolver>
    {

    public:
      /**
       * @brief Constructor
       * @param smoothsteps The number of steps to use for both pre and post smoothing.
       * @param maxiter_ The maximum number of iterations allowed.
       * @param verbose_ The verbosity level to use.
       */
      ISTLBackend_SEQ_CG_AMG_SSOR(int smoothsteps=2,
                                  unsigned maxiter_=5000, int verbose_=1)
        : ISTLBackend_SEQ_AMG<GFS, Dune::SeqSSOR, Dune::CGSolver>(smoothsteps, maxiter_,verbose_)
      {}
    };

    /**
     * @brief Sequential BiCGStab solver preconditioned with AMG smoothed by SSOR
     * @tparam GFS The type of the grid functions space.
     */
    template<class GFS>
    class ISTLBackend_SEQ_BCGS_AMG_SSOR
      : public ISTLBackend_SEQ_AMG<GFS, Dune::SeqSSOR, Dune::BiCGSTABSolver>
    {

    public:
      /**
       * @brief Constructor
       * @param smoothsteps The number of steps to use for both pre and post smoothing.
       * @param maxiter_ The maximum number of iterations allowed.
       * @param verbose_ The verbosity level to use.
       */
      ISTLBackend_SEQ_BCGS_AMG_SSOR(int smoothsteps=2,
                                    unsigned maxiter_=5000, int verbose_=1)
        : ISTLBackend_SEQ_AMG<GFS, Dune::SeqSSOR, Dune::BiCGSTABSolver>(smoothsteps, maxiter_, verbose_)
      {}
    };
    
    template<class GFS>
    class ISTLBackend_SEQ_BCGS_AMG_SOR
      : public ISTLBackend_SEQ_AMG<GFS, Dune::SeqSOR, Dune::BiCGSTABSolver>
    {

    public:
      /**
       * @brief Constructor
       * @param smoothsteps The number of steps to use for both pre and post smoothing.
       * @param maxiter_ The maximum number of iterations allowed.
       * @param verbose_ The verbosity level to use.
       */
      ISTLBackend_SEQ_BCGS_AMG_SOR(int smoothsteps=2,
                                    unsigned maxiter_=5000, int verbose_=1)
        : ISTLBackend_SEQ_AMG<GFS, Dune::SeqSOR, Dune::BiCGSTABSolver>(smoothsteps, maxiter_, verbose_)
      {}
    };

    template<class GFS>
    class ISTLBackend_SEQ_LS_AMG_SSOR
      : public ISTLBackend_SEQ_AMG<GFS, Dune::SeqSSOR, Dune::LoopSolver>
    {

    public:
      /**
       * @brief Constructor
       * @param smoothsteps The number of steps to use for both pre and post smoothing.
       * @param maxiter_ The maximum number of iterations allowed.
       * @param verbose_ The verbosity level to use.
       */
      ISTLBackend_SEQ_LS_AMG_SSOR(int smoothsteps=2,
                                    unsigned maxiter_=5000, int verbose_=1)
        : ISTLBackend_SEQ_AMG<GFS, Dune::SeqSSOR, Dune::LoopSolver>(smoothsteps, maxiter_, verbose_)
      {}
    };

    template<class GFS>
    class ISTLBackend_SEQ_LS_AMG_SOR
      : public ISTLBackend_SEQ_AMG<GFS, Dune::SeqSOR, Dune::LoopSolver>
    {

    public:
      /**
       * @brief Constructor
       * @param smoothsteps The number of steps to use for both pre and post smoothing.
       * @param maxiter_ The maximum number of iterations allowed.
       * @param verbose_ The verbosity level to use.
       */
      ISTLBackend_SEQ_LS_AMG_SOR(int smoothsteps=2,
                                    unsigned maxiter_=5000, int verbose_=1)
        : ISTLBackend_SEQ_AMG<GFS, Dune::SeqSOR, Dune::LoopSolver>(smoothsteps, maxiter_, verbose_)
      {}
    };
  } // namespace PDELab
} // namespace Dune

#endif
