#include "config.h"

#include <dune/common/parametertreeparser.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/finiteelementmap/q1fem.hh>
#include <dune/pdelab/finiteelementmap/p1fem.hh>
#include <dune/pdelab/finiteelementmap/q22dfem.hh>
#include <dune/pdelab/finiteelementmap/q12dfem.hh>
#include <dune/pdelab/finiteelementmap/pk2dfem.hh>
#include <dune/pdelab/finiteelementmap/opbfem.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include <dune/pdelab/constraints/constraintsparameters.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/common/projection.hh>

#include <dune/pdelab/test/gridexamples.hh>


template<typename GV, typename RF>
class InitialFunction :
  public Dune::PDELab::AnalyticGridFunctionBase<
  Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
  InitialFunction<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, InitialFunction<GV,RF> > BaseT;

  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeType RangeType;

  InitialFunction(const GV & gv) : BaseT(gv) {}

  inline void evaluateGlobal(const DomainType & pos, RangeType & val) const
  {
    const RF x = pos[0];
    const RF y = pos[1];
    val = std::exp(-10.0*((0.5-x)*(0.5-x) + (0.5-y)*(0.5-y))) * std::sin(20.0*(x-0.5)*(y-0.5)) + 0.2;
  }
};


int main(int argc, char** argv) {

  try {

    Dune::MPIHelper::instance(argc,argv);

    Dune::ParameterTree parameters;
    Dune::ParameterTreeParser::readINITree(argc > 1 ? argv[1] : "testprojection.ini",parameters);

    const int dim = 2;

    std::vector<int> boundaryIndexToPhysicalGroup, elementIndexToPhysicalGroup;

    ALUUnitSquare unit_square;
    typedef Dune::ALUGrid<dim,dim,Dune::simplex,Dune::nonconforming> Grid;

    Grid& grid = unit_square;

    typedef Grid::ctype ctype;
    typedef Grid::LeafGridView GV;
    GV gv = grid.leafView();

    grid.globalRefine(parameters.get<int>("meshrefine"));

    const int darcy_k = 3;

    typedef ctype DF;
    typedef double RF;

    typedef Dune::PDELab::OPBLocalFiniteElementMap<DF,RF,darcy_k,dim,Dune::GeometryType::simplex,Dune::GMPField<512> > DGFEM;
    DGFEM dg_fem;

    typedef Dune::PDELab::NoConstraints NOCON;

    typedef Dune::PDELab::ISTLVectorBackend<1> VectorBackend;
    typedef VectorBackend::MatrixBackend MBE;

    // DG GridFunctionSpace
    typedef Dune::PDELab::GridFunctionSpace<
      GV,DGFEM,
      NOCON,
      VectorBackend
      > GFS;
    GFS gfs(gv,dg_fem);

    typedef Dune::PDELab::BackendVectorSelector<GFS,RF>::Type V;

    typedef InitialFunction<GV,RF> Func;
    Func func(gv);

    V v(gfs,0.0);
    Dune::PDELab::interpolate(func,gfs,v);


    const int reduced_order = 1;
    typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV,DF,RF,reduced_order> ReducedFEM;
    ReducedFEM reduced_fem(gv);

    // Reduced (P1) GridFunctionSpace
    typedef Dune::PDELab::GridFunctionSpace<
      GV,ReducedFEM,
      NOCON,
      VectorBackend
      > ReducedGFS;
    ReducedGFS reduced_gfs(gv,reduced_fem);

    // Reduction operator
    typedef Dune::PDELab::LocalProjectionOperator<GV> LPO;
    LPO lpo(gv,parameters.get<bool>("adaptiveweights",false));

    typedef Dune::PDELab::GridOperator<
      GFS,ReducedGFS,LPO,
      MBE,RF,RF,RF
      > ReductionGridOperator;
    ReductionGridOperator reduction_gridoperator(gfs,reduced_gfs,lpo);

    typedef ReductionGridOperator::Traits::Jacobian RM;
    typedef ReductionGridOperator::Traits::Range RV;

    RV rv(reduced_gfs,0.0);
    RM rm(reduction_gridoperator);
    rm = 0.0;

    reduction_gridoperator.jacobian(v,rm);

    rm.base().mv(v.base(),rv.base());


#if 0

    /*
    typedef Dune::PDELab::GridOperator<
      ReducedGFS,ReducedGFS,LL2,
      MBE,RF,RF,RF
      > L2GridOperator;
    L2GridOperator l2_gridoperator(reduced_gfs,reduced_gfs,ll2);

    typedef L2GridOperator::Traits::Jacobian L2M;
    typedef L2GridOperator::Traits::Domain L2V;

    L2V l2v(reduced_gfs,0.0);
    L2M l2m(l2_gridoperator,0.0);

    l2_gridoperator.jacobian(l2v,l2m);

    typedef Dune::PDELab::LocalProjectionOperator<GV,L2M> LPO;
    LPO lpo(gv,l2m);

    typedef Dune::PDELab::GridOperator<
      GFS,ReducedGFS,LPO,
      MBE,RF,RF,RF
      > ReductionGridOperator;
    ReductionGridOperator reduction_gridoperator(gfs,reduced_gfs,lpo);

    typedef ReductionGridOperator::Traits::Jacobian RM;
    typedef ReductionGridOperator::Traits::Range RV;

    RV rv(reduced_gfs,0.0);
    RM rm(reduction_gridoperator,0.0);

    reduction_gridoperator.jacobian(u,rm);

    rm.base().mv(u.base(),rv.base());

    typedef Dune::PDELab::StationaryLinearProblemSolver<
      L2GridOperator,
      LS,
      L2V
      > ReductionSolver;

    ReductionSolver reduction_solver(l2_gridoperator,ls,1e-10);

    reduction_solver.apply(rv);

    V ur(gfs,0);

    rm.base().mtv(rv.base(),ur.base());

    */

    typedef Dune::PDELab::LocalProjectionOperator2<GFS,V> LPO;
    LPO lpo(gfs,u);

    typedef Dune::PDELab::GridOperator<
      ReducedGFS,ReducedGFS,LPO,
      MBE,RF,RF,RF
      > ReductionGridOperator;
    ReductionGridOperator reduction_gridoperator(reduced_gfs,reduced_gfs,lpo);

    typedef ReductionGridOperator::Traits::Range RV;

    RV rv(reduced_gfs,0.0);

    typedef Dune::PDELab::StationaryLinearProblemSolver<
      ReductionGridOperator,
      LS,
      RV
      > ReductionSolver;

    ReductionSolver reduction_solver(reduction_gridoperator,ls,1e-10);

    reduction_solver.apply(rv);

    //V ur(gfs,0);

    //rm.base().mtv(rv.base(),ur.base());

    typedef Dune::PDELab::LocalProjectionOperator3<GV> PO;
    PO po(gv);

    typedef Dune::PDELab::GridOperator<
      ReducedGFS,GFS,PO,
      MBE,RF,RF,RF
      > ProlongationGridOperator;
    ProlongationGridOperator prolongation_gridoperator(reduced_gfs,gfs,po);

    typedef ProlongationGridOperator::Traits::Jacobian PM;

    PM pm(prolongation_gridoperator,0.0);

    prolongation_gridoperator.jacobian(rv,pm);

    V ru(gfs,0.0);
    pm.base().mv(rv.base(),ru.base());

    RV rv2(reduced_gfs,0.0);

    pm.base().mtv(u.base(),rv2.base());


#endif

    /*
     * Output
     */

    {
      Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,parameters.get("refineoutput",0));

      typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
      typedef Dune::PDELab::DiscreteGridFunction<ReducedGFS,RV> RDGF;

      DGF dgf(gfs,v);
      RDGF rdgf(reduced_gfs,rv);

      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"v"));
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<RDGF>(rdgf,"rv"));

      vtkwriter.write(parameters["vtkfile"],Dune::VTK::appendedraw);
    }


  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (std::exception &e){
    std::cerr << "STL reported error: " << e.what() << std::endl;
    return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }

}
