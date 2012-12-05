#ifndef DUNE_PDELAB_COMMON_PROJECTION_HH
#define DUNE_PDELAB_COMMON_PROJECTION_HH

#include <dune/common/dynmatrix.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>

namespace Dune {
  namespace PDELab {


    template<typename M>
    void transposeMatrix(const M& A, M& B)
    {
      B.setSize(A.M(),A.N(),A.nonzeroes());
      B.setBuildMode(M::random);
      for (auto rowit = A.begin(); rowit != A.end(); ++rowit)
        for (auto colit = rowit->begin(); colit != rowit->end(); ++colit)
          B.incrementrowsize(colit.index());

      B.endrowsizes();

      for (auto rowit = A.begin(); rowit != A.end(); ++rowit)
        for (auto colit = rowit->begin(); colit != rowit->end(); ++colit)
          B.addindex(colit.index(),rowit.index());

      B.endindices();

      for (auto rowit = A.begin(); rowit != A.end(); ++rowit)
        for (auto colit = rowit->begin(); colit != rowit->end(); ++colit)
          B[colit.index()][rowit.index()] = *colit;
    }


    template<typename Pattern>
    struct ScalarVolumePattern
      : public TypeTree::TreePairVisitor
      , public TypeTree::DynamicTraversal
      , public FullVolumePattern
    {

      template<typename LFSU, typename LFSV, typename TreePath>
      void leaf(const LFSU& lfsu, const LFSV& lfsv, TreePath trePath) const
      {
        pattern_volume(lfsu,lfsv,_pattern);
      }

      ScalarVolumePattern(Pattern& pattern)
        : _pattern(pattern)
      {}

      Pattern& _pattern;

    };



    template<typename GV, typename EG, typename J>
    struct ScalarProjection
      : public TypeTree::TreePairVisitor
      , public TypeTree::DynamicTraversal
    {

      template<typename LFSU, typename LFSV, typename TreePath>
      void leaf(const LFSU& lfsu, const LFSV& lfsv, TreePath treePath) const
      {
        const bool forward = true;

        // Switches between local and global interface
        typedef FiniteElementInterfaceSwitch<
          typename LFSU::Traits::FiniteElementType
          > LFSU_FESwitch;
        typedef BasisInterfaceSwitch<
          typename LFSU_FESwitch::Basis
          > LFSU_BasisSwitch;

        typedef FiniteElementInterfaceSwitch<
          typename LFSV::Traits::FiniteElementType
          > LFSV_FESwitch;
        typedef BasisInterfaceSwitch<
          typename LFSV_FESwitch::Basis
          > LFSV_BasisSwitch;

        // domain and range field type
        typedef typename LFSU_BasisSwitch::DomainField DF;
        typedef typename LFSU_BasisSwitch::RangeField RF;
        typedef typename LFSU_BasisSwitch::Range RangeType;

        typedef typename LFSU::Traits::SizeType size_type;

        // dimensions
        const int dim = EG::Geometry::dimension;

        const int integration_order = 2 * std::max(
          LFSU_FESwitch::basis(lfsu.finiteElement()).order(),
          LFSV_FESwitch::basis(lfsv.finiteElement()).order()
          );

        const size_type u_size = lfsu.size();
        const size_type v_size = lfsv.size();

        DynamicMatrix<RF> mass_matrix_f(forward ? v_size : 0, forward ? v_size : 0, 0);
        DynamicMatrix<RF> mass_matrix_b(forward ? 0 : u_size, forward ? 0 : u_size, 0);
        DynamicMatrix<RF> proj_matrix(forward ? u_size : v_size, forward ? v_size : u_size, 0);

        // select quadrature rule
        GeometryType gt = eg.geometry().type();
        const QuadratureRule<DF,dim>& rule = QuadratureRules<DF,dim>::rule(gt,integration_order);

        std::vector<RangeType> phi(u_size);
        std::vector<RangeType> psi(v_size);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(),
               endit = rule.end();
             it != endit;
             ++it)
          {
            LFSU_FESwitch::basis(lfsu.finiteElement()).evaluateFunction(it->position(),phi);
            LFSV_FESwitch::basis(lfsv.finiteElement()).evaluateFunction(it->position(),psi);

            const RF factor = it->weight();
            for (size_type i = 0; i < v_size; ++i)
              {
                if (forward)
                  for (size_type j = 0; j < u_size; ++j)
                    proj_matrix[j][i] += psi[i]*phi[j]*factor;
                else
                  for (size_type j = 0; j < u_size; ++j)
                    proj_matrix[i][j] += psi[i]*phi[j]*factor;

                if (forward)
                  for (size_type j = 0; j < v_size; ++j)
                    mass_matrix_f[i][j] += psi[i]*psi[j]*factor;
              }

            if (!forward)
              for (size_type i = 0; i < u_size; ++i)
                for (size_type j = 0; j < u_size; ++j)
                  mass_matrix_b[i][j] += phi[i]*phi[j]*factor;
          }


        if (_adaptive_weights && !forward)
          for (size_type i = 0; i < v_size; ++i)
            {
              const auto& key = LFSV_FESwitch::coefficients(lfsv.finiteElement()).localKey(i);
              if (key.codim() != dim)
                continue;

              const double weight = _cell_counts[_gv.indexSet().subIndex(eg.entity(),key.subEntity(),dim)];

              if (weight != 1.0)
                //mass_matrix_f[i][i] *= weight;
                for (int j = 0; j < u_size; ++j)
                  proj_matrix[j][i] *= weight;
            }

        //for (size_type i = 0; i < u_size; ++i)
        //  mass_matrix_up[i][i] = 1.0;

        if (forward)
          mass_matrix_f.invert();
        else
          mass_matrix_b.invert();

        /*
          for (size_type i = 0; i < v_size; ++i)
            {
              const RF factor = mass_matrix_f[i][i];
              for (size_type j = 0; j < v_size; ++j)
              mass_matrix_f[i][j] /= factor;
            }
        */

        /*
          for (size_type i = 0; i < u_size; ++i)
            {
              const RF factor = mass_matrix_up[i][i];
              for (size_type j = 0; j < u_size; ++j)
              mass_matrix_up[i][j] /= factor;
            }
        */
        //typename DynamicMatrix<RF>::row_type row(u_size,0);

        for (size_type i = 0; i < v_size; ++i)
          {
            const auto& key = LFSV_FESwitch::coefficients(lfsv.finiteElement()).localKey(i);
            double weight = 1.0;

            if (_adaptive_weights)
              {
                weight = _cell_counts[_gv.indexSet().subIndex(eg.entity(),key.subEntity(),key.codim())];
                weight = weight > 0 ? 1.0/weight : 1.0;
              }

            //mass_matrix.mv(proj_matrix_T[i],row);
            for (size_type j = 0; j < u_size; ++j)
              {
                //jac.accumulate(lfsv,i,lfsu,j,row[j]);
                jac.accumulate(lfsv,i,lfsu,j,weight * (mass_matrix_f[i]*proj_matrix[j]));
                //jac.accumulate(lfsv,i,lfsu,j,proj_matrix[j][i]);
              }
          }

      }

      ScalarProjection(const EG& eg_, J& j, const GV& gv, const std::vector<unsigned char>& cell_counts, bool adaptive_weights)
        : eg(eg_)
        , jac(j)
        , _gv(gv)
        , _cell_counts(cell_counts)
        , _adaptive_weights(adaptive_weights)
      {}

      const EG& eg;
      J& jac;
      GV _gv;
      const std::vector<unsigned char>& _cell_counts;
      const bool _adaptive_weights;

    };


    template<typename GV>
    class LocalProjectionOperator
      : public LocalOperatorDefaultFlags
    {

    public:

      static const bool doPatternVolume = true;
      static const bool doAlphaVolume = true;

      template<typename LFSU, typename LFSV>
      void pattern_volume (const LFSU& lfsu, const LFSV& lfsv,
                           LocalSparsityPattern& pattern) const
      {
        ScalarVolumePattern<LocalSparsityPattern> svp(pattern);
        TypeTree::applyToTreePair(lfsu,lfsv,svp);
      }


      template<typename EG, typename LFSU, typename X, typename LFSV, typename J>
      void jacobian_volume(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, J& jac) const
      {
        ScalarProjection<GV,EG,J> sp(eg,jac,_gv,_cell_counts,_adaptive_weights);
        TypeTree::applyToTreePair(lfsu,lfsv,sp);
      }

      LocalProjectionOperator(const GV& gv, bool adaptive_weights = false)
        : _gv(gv)
        , _cell_counts(adaptive_weights ? gv.indexSet().size(GV::dimension) : 0)
        , _adaptive_weights(adaptive_weights)
      {
        if (!_adaptive_weights)
          return;
        const int dim = GV::dimension;
        for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it)
          {
            const auto& ref_el = GenericReferenceElements<typename GV::ctype,dim>::general(it->type());
            for (std::size_t i = 0; i < ref_el.size(dim); ++i)
              ++_cell_counts[gv.indexSet().subIndex(*it,i,dim)];
          }
      }

      bool weighted() const
      {
        return _adaptive_weights;
      }

      void setWeighted(bool weighted)
      {
        _adaptive_weights = weighted;
      }

      GV _gv;
      std::vector<unsigned char> _cell_counts;
      bool _adaptive_weights;

    };


  } // namespace PDELab
} // namespace Dune


#endif // DUNE_PDELAB_COMMON_PROJECTION_HH
