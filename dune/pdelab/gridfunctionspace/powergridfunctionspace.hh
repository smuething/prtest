// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_POWERGRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_POWERGRIDFUNCTIONSPACE_HH

#include "gridfunctionspace.hh"

namespace Dune {
  namespace PDELab {

    //=======================================
    // power grid function space
    //=======================================

    template<typename T, int k, typename P>
    class PowerGridFunctionSpace;

    //! product of identical grid function spaces
    //! base class that holds implementation of the methods
    //! this is the default version with lexicographic ordering
    //!
    //! \tparam T PLEASE DOCUMENT
    //! \tparam k PLEASE DOCUMENT
    //! \tparam P PLEASE DOCUMENT
    template<typename T, int k, typename P>
    class PowerGridFunctionSpaceBase 
      : public PowerNode<T,k,CountingPointerStoragePolicy>,
        public Countable
    {
    public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType,
                                                    P, k>
      Traits;

      //! extract type of container storing Es
      template<typename E>
      struct VectorContainer
      {
        //! \brief define Type as the Type of a container of E's
        typedef typename Traits::BackendType::template VectorContainer<PowerGridFunctionSpaceBase,E> Type;	
      private:
        VectorContainer () {}
      };

      //! extract type for storing constraints
      template<typename E>
      struct ConstraintsContainer
      {
        //! \brief define Type as the Type of a container of E's
        typedef ConstraintsTransformation<typename Traits::SizeType,E> Type;	
      private:
        ConstraintsContainer () {}
      };

      // define local function space parametrized by self 
      typedef Dune::PDELab::PowerLocalFunctionSpaceNode<PowerGridFunctionSpaceBase> LocalFunctionSpace;


      // it is part of the trick to have a constructor without arguments
      // setting of the children is then done by the constructors
      // of the specialized derived classes
      PowerGridFunctionSpaceBase ()
      {
      }

      // get grid view
      const typename Traits::GridViewType& gridview () const
      {
        return this->template getChild<0>().gridview();
      }

      //! get dimension of root finite element space
      typename Traits::SizeType globalSize () const
      {
        // this is bullshit all children may have different
        // size although they have the same type ...
        // [Jö] well, it does happen for the elements I use for the Yee FDTD
        //      scheme, at least
        return offset[k];
      }

      //! get dimension of this finite element space
      typename Traits::SizeType size () const
      {
        return offset[k];
      }

      // get max dimension of shape function space
      typename Traits::SizeType maxLocalSize () const
      {
        // this is bullshit !
        return maxlocalsize;
      }

      //! map index from our index set [0,size()-1] to root index set
      typename Traits::SizeType upMap (typename Traits::SizeType i) const
      {
        return i;
      }

      //! map index from child i's index set into our index set
      template<int i>
      typename Traits::SizeType subMap (typename Traits::SizeType j) const
      {
        return offset[i]+j;
      }

      //------------------------------
      // generic data handle interface
      //------------------------------

      //! returns true if data for this codim should be communicated
      bool dataHandleContains (int dim, int codim) const
      {
        for (int i=0; i<k; i++)
          if (this->getChild(i).dataHandleContains(dim,codim))
            return true;
        return false;
      }
      
      //! returns true if size per entity of given dim and codim is a constant
      bool dataHandleFixedSize (int dim, int codim) const
      {
        for (int i=0; i<k; i++)
          if (!this->getChild(i).dataHandleFixedSize(dim,codim))
            return false;
        return true;
      }
      
      /*! how many objects of type DataType have to be sent for a given entity
        
        Note: Only the sender side needs to know this size. 
      */
      template<class EntityType>
      size_t dataHandleSize (const EntityType& e) const
      {
        size_t n=0;
        for (int i=0; i<k; i++)
          n += this->getChild(i).dataHandleSize(e);
        return n;
      }

      //! return vector of global indices associated with the given entity
      template<class EntityType>
      void dataHandleGlobalIndices (const EntityType& e, 
                                    std::vector<typename Traits::SizeType>& global) const
      {
        size_t n=0;
        for (int i=0; i<k; i++)
          n += this->getChild(i).dataHandleSize(e);
        global.resize(n);
        n = 0;
        for (int i=0; i<k; i++)
          {
            this->getChild(i).dataHandleGlobalIndices(e,childglobal);
            for (size_t j=0; j<childglobal.size(); j++)
              global[n+j] = offset[i]+childglobal[j];
            n += childglobal.size();
          }          
      }

      //------------------------------

      // recalculate sizes
      void update ()
      {
        for (int i=0; i<k; i++)
          (*this)[i].update();
        setup();
      }

    protected:
      void setup ()
      {
        Dune::dinfo << "PowerGridFunctionSpace(lexicographic version):"
                    << std::endl;
        Dune::dinfo << "( ";
        offset[0] = 0;
        maxlocalsize = 0;
        for (int i=0; i<k; i++)
          {
            childSize[i] = this->getChild(i).globalSize();
            Dune::dinfo << childSize[i] << " ";
            offset[i+1] = offset[i]+childSize[i];
            maxlocalsize += this->getChild(i).maxLocalSize();
          }
        Dune::dinfo << ") total size = " << offset[k]
                    << " max local size = " << maxlocalsize
                    << std::endl;
        childglobal.resize(maxlocalsize);
      }

    private:
      typename Traits::SizeType childSize[k];
      typename Traits::SizeType offset[k+1];
      typename Traits::SizeType maxlocalsize;
      mutable std::vector<typename Traits::SizeType> childglobal;
    };


    // product of identical grid function spaces
    // base class that holds implementation of the methods
    // specialization for blockwise ordering
    template<typename T, int k, int s>
    class PowerGridFunctionSpaceBase<T,k,GridFunctionSpaceComponentBlockwiseMapper<s> >
      : public PowerNode<T,k,CountingPointerStoragePolicy>,
        public Countable
    {
    public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType,
                                                    GridFunctionSpaceComponentBlockwiseMapper<s>, k>
      Traits;

      //! extract type of container storing Es
      template<typename E>
      struct VectorContainer
      {
        //! \brief define Type as the Type of a container of E's
        typedef typename Traits::BackendType::template VectorContainer<PowerGridFunctionSpaceBase,E> Type;	
      private:
        VectorContainer () {}
      };

      //! extract type for storing constraints
      template<typename E>
      struct ConstraintsContainer
      {
        //! \brief define Type as the Type of a container of E's
        typedef ConstraintsTransformation<typename Traits::SizeType,E> Type;	
      private:
        ConstraintsContainer () {}
      };

      // define local function space parametrized by self 
      typedef Dune::PDELab::PowerLocalFunctionSpaceNode<PowerGridFunctionSpaceBase> LocalFunctionSpace;
      
      // it is part of the trick to have a constructor without arguments
      // setting of the children is then done by the constructors
      // of the specialized derived classes
      PowerGridFunctionSpaceBase ()
      {
      }

      // get grid view
      const typename Traits::GridViewType& gridview () const
      {
        return this->template getChild<0>().gridview();
      }

      //! get dimension of root finite element space
      typename Traits::SizeType globalSize () const
      {
        // this is bullshit all children may have different
        // size although they have the same type ...
        // [Jö] well, it does happen for the elements I use for the Yee FDTD
        //      scheme, at least
        return offset[k];
      }

      //! get dimension of this finite element space
      typename Traits::SizeType size () const
      {
        return offset[k];
      }

      // get max dimension of shape function space
      typename Traits::SizeType maxLocalSize () const
      {
        // this is bullshit !
        return maxlocalsize;
      }

      //! map index from our index set [0,size()-1] to root index set
      typename Traits::SizeType upMap (typename Traits::SizeType i) const
      {
        return i;
      }

      //! map index from child i's index set into our index set
      template<int i>
      typename Traits::SizeType subMap (typename Traits::SizeType j) const
      {
        return (j%s)+(j/s)*k*s+i*s;
      }

      //------------------------------
      // generic data handle interface
      //------------------------------

      //! returns true if data for this codim should be communicated
      bool dataHandleContains (int dim, int codim) const
      {
        for (int i=0; i<k; i++)
          if (this->getChild(i).dataHandleContains(dim,codim))
            return true;
        return false;
      }
      
      //! returns true if size per entity of given dim and codim is a constant
      bool dataHandleFixedSize (int dim, int codim) const
      {
        for (int i=0; i<k; i++)
          if (!this->getChild(i).dataHandleFixedSize(dim,codim))
            return false;
        return true;
      }
      
      /*! how many objects of type DataType have to be sent for a given entity
        
        Note: Only the sender side needs to know this size. 
      */
      template<class EntityType>
      size_t dataHandleSize (const EntityType& e) const
      {
        size_t n=0;
        for (int i=0; i<k; i++)
          n += this->getChild(i).dataHandleSize(e);
        return n;
      }

      //! return vector of global indices associated with the given entity
      template<class EntityType>
      void dataHandleGlobalIndices (const EntityType& e, 
                                    std::vector<typename Traits::SizeType>& global) const
      {
        size_t n=0;
        for (int i=0; i<k; i++)
          n += this->getChild(i).dataHandleSize(e);
        global.resize(n);
        n = 0;
        for (int i=0; i<k; i++)
          {
            this->getChild(i).dataHandleGlobalIndices(e,childglobal);
            for (size_t j=0; j<childglobal.size(); j++)
              global[n+j] = childglobal[j]*k+i;
            n += childglobal.size();
          }          
      }

      //------------------------------

      // recalculate sizes
      void update ()
      {
        for (int i=0; i<k; i++)
          (*this)[i].update();
        setup();
      }

    protected:
      void setup ()
      {
        Dune::dinfo << "PowerGridFunctionSpace(blockwise version):"
                    << std::endl;
        Dune::dinfo << "( ";
        offset[0] = 0;
        maxlocalsize = 0;
        for (int i=0; i<k; i++)
          {
            childSize[i] = this->getChild(i).globalSize();
            offset[i+1] = offset[i]+childSize[i];
            Dune::dinfo << childSize[i] << "[" << offset[i] << "] ";
            maxlocalsize += this->getChild(i).maxLocalSize();
          }
        Dune::dinfo << ") total size = " << offset[k]
                    << " max local size = " << maxlocalsize
                    << std::endl;
        /* check the local block size */
        if (this->getChild(0).maxLocalSize()%s != 0)
          DUNE_THROW(Exception,
                     "number of DOFs (" << this->getChild(0).maxLocalSize() << ") per component "
                     "must be a multiple of the BlockSize (" << s << ")");
        for (int i=1; i<k; i++)
          if (childSize[i]!=childSize[0])
            DUNE_THROW(Exception, "components must be of equal size");
        childglobal.resize(maxlocalsize);
      }

    private:
      typename Traits::SizeType childSize[k];
      typename Traits::SizeType offset[k+1];
      typename Traits::SizeType maxlocalsize;
      mutable std::vector<typename Traits::SizeType> childglobal;
    };

    template<typename T, int k>
    class PowerGridFunctionSpaceBase<T,k,GridFunctionSpaceBlockwiseMapper >
      : public PowerGridFunctionSpaceBase<T,k,GridFunctionSpaceComponentBlockwiseMapper<1> >
    {
    protected:
      using PowerGridFunctionSpaceBase<T,k,GridFunctionSpaceComponentBlockwiseMapper<1> >::setup;
    };


    /** 
        \brief Tupel of grid function spaces base class that holds
        implementation of the methods specialization for dynamic
        blockwise ordering
    */
    template<typename T, int k>
    class PowerGridFunctionSpaceBase<T,k,GridFunctionSpaceDynamicBlockwiseMapper >
      : public PowerNode<T,k,CountingPointerStoragePolicy>,
        public Countable
    {
    public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType,
                                                    GridFunctionSpaceDynamicBlockwiseMapper, k>
      Traits;

      typedef PowerNode<T,k,CountingPointerStoragePolicy> BaseT;

      //! extract type of container storing Es
      template<typename E>
      struct VectorContainer
      {
        //! \brief define Type as the Type of a container of E's
        typedef typename Traits::BackendType::template VectorContainer<PowerGridFunctionSpaceBase,E> Type;	
      private:
        VectorContainer () {}
      };

      //! extract type for storing constraints
      template<typename E>
      struct ConstraintsContainer
      {
        //! \brief define Type as the Type of a container of E's
        typedef ConstraintsTransformation<typename Traits::SizeType,E> Type;	
      private:
        ConstraintsContainer () {}
      };

      // define local function space parametrized by self 
      typedef Dune::PDELab::PowerLocalFunctionSpaceNode<PowerGridFunctionSpaceBase> LocalFunctionSpace;
      
      // it is part of the trick to have a constructor without arguments
      // setting of the children is then done by the constructors
      // of the specialized derived classes
      PowerGridFunctionSpaceBase ()
      {
      }

      // get grid view
      const typename Traits::GridViewType& gridview () const
      {
        return this->template getChild<0>().gridview();
      }

      //! get dimension of root finite element space
      typename Traits::SizeType globalSize () const
      {
        // this is bullshit all children may have different
        // size although they have the same type ...
        // [Jö] well, it does happen for the elements I use for the Yee FDTD
        //      scheme, at least
        return offset[k];
      }

      //! get dimension of this finite element space
      typename Traits::SizeType size () const
      {
        return offset[k];
      }

      // get max dimension of shape function space
      typename Traits::SizeType maxLocalSize () const
      {
        // this is bullshit !
        return maxlocalsize;
      }

      //! map index from our index set [0,size()-1] to root index set
      typename Traits::SizeType upMap (typename Traits::SizeType i) const
      {
        return i;
      }

      //! map index from child i's index set into our index set
      template<int i>
      typename Traits::SizeType subMap (typename Traits::SizeType j) const
      {
        BlockIndexRangeIterator & it = blockIndexIterators[i];
        while(it->first < j)++it;
        while(it->first > j)--it;
        const typename Traits::SizeType global_index = it->second + j - it->first;
        return global_index;
      }

      //------------------------------
      // generic data handle interface
      //------------------------------

      //! returns true if data for this codim should be communicated
      bool dataHandleContains (int dim, int codim) const
      {
        for (int i=0; i<k; i++)
          if (this->getChild(i).dataHandleContains(dim,codim))
            return true;
        return false;
      }
      
      //! returns true if size per entity of given dim and codim is a constant
      bool dataHandleFixedSize (int dim, int codim) const
      {
        for (int i=0; i<k; i++)
          if (!this->getChild(i).dataHandleFixedSize(dim,codim))
            return false;
        return true;
      }
      
      /*! how many objects of type DataType have to be sent for a given entity
        
        Note: Only the sender side needs to know this size. 
      */
      template<class EntityType>
      size_t dataHandleSize (const EntityType& e) const
      {
        size_t n=0;
        for (int i=0; i<k; i++)
          n += this->getChild(i).dataHandleSize(e);
        return n;
      }

      //! return vector of global indices associated with the given entity
      template<class EntityType>
      void dataHandleGlobalIndices (const EntityType& e, 
                                    std::vector<typename Traits::SizeType>& global) const
      {
        size_t n=0;
        for (int i=0; i<k; i++)
          n += this->getChild(i).dataHandleSize(e);
        global.resize(n);
        n = 0;
        for (int i=0; i<k; i++)
          {
            this->getChild(i).dataHandleGlobalIndices(e,childglobal);
            for (size_t j=0; j<childglobal.size(); j++)
              global[n+j] = childglobal[j]*k+i;
            n += childglobal.size();
          }          
      }

      //------------------------------

      // recalculate sizes
      void update ()
      {
        for (int i=0; i<k; i++)
          (*this)[i].update();
        setup();
      }

    protected:
      void setup ()
      {
        
        typedef typename Traits::GridViewType GridView;
        const GridView & gv = gridview();

        // Initialize offset array for each child
        blockIndices.clear();
        blockIndices.resize(BaseT::CHILDREN);

        std::vector<std::vector<typename Traits::SizeType> > childOffsets;
        childOffsets.resize(BaseT::CHILDREN);
        for(int i=0; i<BaseT::CHILDREN; ++i)
          childOffsets[i].resize(gv.size(0)+1);


        // Iterate grid and determine the offsets for each child
        typedef typename GridView::template Codim<0>::Iterator Iterator;
        Iterator it = gv.template begin<0>();
        const Iterator eit = gv.template end<0>();
        typename Traits::SizeType running_index(0);
        for(; it!=eit; ++it){
          typename Traits::SizeType e_index = gv.indexSet().index(*it);

          // Loop over children (realized by meta-program)
          DynamicBlockwiseMapperImp::GetChildOffsetsMetaProgram<PowerGridFunctionSpaceBase,BaseT::CHILDREN,0>::
            getChildOffsets(*this,*it,childOffsets);
          
          for(int i=0; i<BaseT::CHILDREN; ++i){
            if(childOffsets[i][e_index+1] != childOffsets[i][e_index]){
              // Add new block index range element to list
              blockIndices[i].push_back(SizeTypePair(childOffsets[i][e_index],running_index));

              // Update running index
              running_index += childOffsets[i][e_index+1] - childOffsets[i][e_index];
            }
          }
        }
        
        // Insert a "stop entry" at end of every list
        for(int i=0; i<BaseT::CHILDREN; ++i)
          blockIndices[i].push_back(SizeTypePair(childOffsets[i][gv.size(0)],running_index));


        blockIndexIterators.clear();
        for(int i=0; i<BaseT::CHILDREN; ++i)
          blockIndexIterators.push_back(blockIndices[i].begin());


        Dune::dinfo << "PowerGridFunctionSpace(blockwise version):"
                    << std::endl;
        Dune::dinfo << "( ";
        offset[0] = 0;
        maxlocalsize = 0;
        for (int i=0; i<k; i++)
          {
            childSize[i] = this->getChild(i).globalSize();
            offset[i+1] = offset[i]+childSize[i];
            Dune::dinfo << childSize[i] << "[" << offset[i] << "] ";
            maxlocalsize += this->getChild(i).maxLocalSize();
          }
        Dune::dinfo << ") total size = " << offset[k]
                    << " max local size = " << maxlocalsize
                    << std::endl;

        childglobal.resize(maxlocalsize);
      }

    private:
      typename Traits::SizeType childSize[k];
      typename Traits::SizeType offset[k+1];
      typename Traits::SizeType maxlocalsize;
      mutable std::vector<typename Traits::SizeType> childglobal;
      typedef std::pair<typename Traits::SizeType,typename Traits::SizeType> SizeTypePair;
      typedef std::list<SizeTypePair> BlockIndexRangeList;
      std::vector<BlockIndexRangeList> blockIndices;
      typedef typename BlockIndexRangeList::const_iterator BlockIndexRangeIterator;
      mutable std::vector<BlockIndexRangeIterator> blockIndexIterators;
    };


    //! \addtogroup GridFunctionSpace
    //! \{

    /// \brief product of identical grid function spaces
    ///
    /// the specializations of this class just set the members
    /// all the methods are generic in the implementation
    /// \tparam T The type of the underlying grid function space
    /// \tparam k how many identical function space to use in the composition
    /// \tparam P The type of the mapper used. The mapper maps each degree of freedom
    /// of each function space to a unique index. Use e.g. 
    /// \link GridFunctionSpaceLexicographicMapper GridFunctionSpaceLexicographicMapper \endlink
    /// or \link  GridFunctionSpaceComponentBlockwiseMapper  GridFunctionSpaceComponentBlockwiseMapper \endlink
    /// or \link  GridFunctionSpaceBlockwiseMapper  GridFunctionSpaceBlockwiseMapper \endlink
    template<typename T, int k, typename P=GridFunctionSpaceLexicographicMapper>
    class PowerGridFunctionSpace 
      : public PowerGridFunctionSpaceBase<T,k,P>
    {
    public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType,
                                                    P, k>
      Traits;

      //! \brief Construct a PowerGridFunction with k clones of the function t
      //! \param t grid function space to clone.
      PowerGridFunctionSpace (T& t)
      {
        for (int i=0; i<k; i++)
          setChild(i,t);
        PowerGridFunctionSpaceBase<T,k,P>::setup();
      }

      //! \brief Construct a PowerGridFunction with k clones of the function t
      //! \param t grid function space to clone.
      PowerGridFunctionSpace (T** t)  
      {
        for (int i=0; i<k; i++)
          setChild(i,*(t[i]));
        PowerGridFunctionSpaceBase<T,k,P>::setup();
      }
    };

    //! \}

#ifndef DOXYGEN
    template<typename T, typename P>
    class PowerGridFunctionSpace<T,2,P> 
      : public PowerGridFunctionSpaceBase<T,2,P>
    {
    public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType,
                                                    P, 2>
      Traits;

      //! Construct a PowerGridFunction with k clones of the function t
      PowerGridFunctionSpace (T& t)
      {
        this->template setChild<0>(t);
        this->template setChild<1>(t);
        PowerGridFunctionSpaceBase<T,2,P>::setup();
      }

      //! Construct a PowerGridFunction with k clones of the function t
      PowerGridFunctionSpace (T& t0, T& t1)
      {
        this->template setChild<0>(t0);
        this->template setChild<1>(t1);
        PowerGridFunctionSpaceBase<T,2,P>::setup();
      }
    };

    template<typename T, typename P>
    class PowerGridFunctionSpace<T,3,P> 
      : public PowerGridFunctionSpaceBase<T,3,P>
    {
    public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType,
                                                    P, 3>
      Traits;

      //! Construct a PowerGridFunction with k clones of the function t
      PowerGridFunctionSpace (T& t)
      {
        this->template setChild<0>(t);
        this->template setChild<1>(t);
        this->template setChild<2>(t);
        PowerGridFunctionSpaceBase<T,3,P>::setup();
      }

      //! Construct a PowerGridFunction with k clones of the function t
      PowerGridFunctionSpace (T& t0, T& t1, T& t2)
      {
        this->template setChild<0>(t0);
        this->template setChild<1>(t1);
        this->template setChild<2>(t2);
        PowerGridFunctionSpaceBase<T,3,P>::setup();
      }
    };

    template<typename T, typename P>
    class PowerGridFunctionSpace<T,4,P> 
      : public PowerGridFunctionSpaceBase<T,4,P>
    {
    public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType,
                                                    P, 4>
      Traits;

      //! Construct a PowerGridFunction with k clones of the function t
      PowerGridFunctionSpace (T& t)
      {
        this->template setChild<0>(t);
        this->template setChild<1>(t);
        this->template setChild<2>(t);
        this->template setChild<3>(t);
        PowerGridFunctionSpaceBase<T,4,P>::setup();
      }

      //! Construct a PowerGridFunction with k clones of the function t
      PowerGridFunctionSpace (T& t0, T& t1, T& t2, T& t3)
      {
        this->template setChild<0>(t0);
        this->template setChild<1>(t1);
        this->template setChild<2>(t2);
        this->template setChild<3>(t3);
        PowerGridFunctionSpaceBase<T,4,P>::setup();
      }
    };

    template<typename T, typename P>
    class PowerGridFunctionSpace<T,5,P> 
      : public PowerGridFunctionSpaceBase<T,5,P>
    {
    public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType,
                                                    P, 5>
      Traits;

      //! Construct a PowerGridFunction with k clones of the function t
      PowerGridFunctionSpace (T& t)
      {
        this->template setChild<0>(t);
        this->template setChild<1>(t);
        this->template setChild<2>(t);
        this->template setChild<3>(t);
        this->template setChild<4>(t);
        PowerGridFunctionSpaceBase<T,5,P>::setup();
      }

      //! Construct a PowerGridFunction with k clones of the function t
      PowerGridFunctionSpace (T& t0, T& t1, T& t2, T& t3, T& t4)
      {
        this->template setChild<0>(t0);
        this->template setChild<1>(t1);
        this->template setChild<2>(t2);
        this->template setChild<3>(t3);
        this->template setChild<4>(t4);
        PowerGridFunctionSpaceBase<T,5,P>::setup();
      }
    };

    template<typename T, typename P>
    class PowerGridFunctionSpace<T,6,P> 
      : public PowerGridFunctionSpaceBase<T,6,P>
    {
    public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType,
                                                    P, 6>
      Traits;

      //! Construct a PowerGridFunction with k clones of the function t
      PowerGridFunctionSpace (T& t)
      {
        this->template setChild<0>(t);
        this->template setChild<1>(t);
        this->template setChild<2>(t);
        this->template setChild<3>(t);
        this->template setChild<4>(t);
        this->template setChild<5>(t);
        PowerGridFunctionSpaceBase<T,6,P>::setup();
      }

      //! Construct a PowerGridFunction with k clones of the function t
      PowerGridFunctionSpace (T& t0, T& t1, T& t2, T& t3, T& t4, T& t5)
      {
        this->template setChild<0>(t0);
        this->template setChild<1>(t1);
        this->template setChild<2>(t2);
        this->template setChild<3>(t3);
        this->template setChild<4>(t4);
        this->template setChild<5>(t5);
        PowerGridFunctionSpaceBase<T,6,P>::setup();
      }
    };

    template<typename T, typename P>
    class PowerGridFunctionSpace<T,7,P> 
      : public PowerGridFunctionSpaceBase<T,7,P>
    {
    public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType,
                                                    typename T::Traits::BackendType,
                                                    P, 7>
      Traits;

      //! Construct a PowerGridFunction with k clones of the function t
      PowerGridFunctionSpace (T& t)
      {
        this->template setChild<0>(t);
        this->template setChild<1>(t);
        this->template setChild<2>(t);
        this->template setChild<3>(t);
        this->template setChild<4>(t);
        this->template setChild<5>(t);
        this->template setChild<6>(t);
        PowerGridFunctionSpaceBase<T,7,P>::setup();
      }

      //! Construct a PowerGridFunction with k clones of the function t
      PowerGridFunctionSpace (T& t0, T& t1, T& t2, T& t3, T& t4, T& t5, T& t6)
      {
        this->template setChild<0>(t0);
        this->template setChild<1>(t1);
        this->template setChild<2>(t2);
        this->template setChild<3>(t3);
        this->template setChild<4>(t4);
        this->template setChild<5>(t5);
        this->template setChild<6>(t6);
        PowerGridFunctionSpaceBase<T,7,P>::setup();
      }
    };

    template<typename T, typename P>
    class PowerGridFunctionSpace<T,8,P> 
      : public PowerGridFunctionSpaceBase<T,8,P>
    {
    public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType,
                                                    P, 8>
      Traits;

      //! Construct a PowerGridFunction with k clones of the function t
      PowerGridFunctionSpace (T& t)
      {
        this->template setChild<0>(t);
        this->template setChild<1>(t);
        this->template setChild<2>(t);
        this->template setChild<3>(t);
        this->template setChild<4>(t);
        this->template setChild<5>(t);
        this->template setChild<6>(t);
        this->template setChild<7>(t);
        PowerGridFunctionSpaceBase<T,8,P>::setup();
      }

      //! Construct a PowerGridFunction with k clones of the function t
      PowerGridFunctionSpace (T& t0, T& t1, T& t2, T& t3, T& t4, T& t5, T& t6, T& t7)
      {
        this->template setChild<0>(t0);
        this->template setChild<1>(t1);
        this->template setChild<2>(t2);
        this->template setChild<3>(t3);
        this->template setChild<4>(t4);
        this->template setChild<5>(t5);
        this->template setChild<6>(t6);
        this->template setChild<7>(t7);
        PowerGridFunctionSpaceBase<T,8,P>::setup();
      }
    };

    template<typename T, typename P>
    class PowerGridFunctionSpace<T,9,P> 
      : public PowerGridFunctionSpaceBase<T,9,P>
    {
    public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType,
                                                    P, 9>
      Traits;

      //! Construct a PowerGridFunction with k clones of the function t
      PowerGridFunctionSpace (T& t)
      {
        this->template setChild<0>(t);
        this->template setChild<1>(t);
        this->template setChild<2>(t);
        this->template setChild<3>(t);
        this->template setChild<4>(t);
        this->template setChild<5>(t);
        this->template setChild<6>(t);
        this->template setChild<7>(t);
        this->template setChild<8>(t);
        PowerGridFunctionSpaceBase<T,9,P>::setup();
      }

      //! Construct a PowerGridFunction with k clones of the function t
      PowerGridFunctionSpace (T& t0, T& t1, T& t2, T& t3, T& t4, T& t5, T& t6, T& t7, T& t8)
      {
        this->template setChild<0>(t0);
        this->template setChild<1>(t1);
        this->template setChild<2>(t2);
        this->template setChild<3>(t3);
        this->template setChild<4>(t4);
        this->template setChild<5>(t5);
        this->template setChild<6>(t6);
        this->template setChild<7>(t7);
        this->template setChild<8>(t8);
        PowerGridFunctionSpaceBase<T,9,P>::setup();
      }
    };

    template<typename T, typename P>
    class PowerGridFunctionSpace<T,10,P> 
      : public PowerGridFunctionSpaceBase<T,10,P>
    {
    public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType,
                                                    P, 10>
      Traits;

      //! Construct a PowerGridFunction with k clones of the function t
      PowerGridFunctionSpace (T& t)
      {
        this->template setChild<0>(t);
        this->template setChild<1>(t);
        this->template setChild<2>(t);
        this->template setChild<3>(t);
        this->template setChild<4>(t);
        this->template setChild<5>(t);
        this->template setChild<6>(t);
        this->template setChild<7>(t);
        this->template setChild<8>(t);
        this->template setChild<9>(t);
        PowerGridFunctionSpaceBase<T,10,P>::setup();
      }

      //! Construct a PowerGridFunction with k clones of the function t
      PowerGridFunctionSpace (T& t0, T& t1, T& t2, T& t3, T& t4, T& t5, T& t6, T& t7, T& t8, T& t9)
      {
        this->template setChild<0>(t0);
        this->template setChild<1>(t1);
        this->template setChild<2>(t2);
        this->template setChild<3>(t3);
        this->template setChild<4>(t4);
        this->template setChild<5>(t5);
        this->template setChild<6>(t6);
        this->template setChild<7>(t7);
        this->template setChild<8>(t8);
        this->template setChild<9>(t9);
        PowerGridFunctionSpaceBase<T,10,P>::setup();
      }
    };
#endif // DOXYGEN

  }
}    
#endif
