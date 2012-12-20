#include "config.h"

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/dgfparser/gridptr.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>

#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/finiteelementmap/hangingnodeconstraints.hh>
#include <dune/pdelab/finiteelementmap/q1fem.hh>
#include <dune/pdelab/constraints/constraints.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/localoperator/l2.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>

#include <dune/istl/io.hh>


class ConstraintsParameters
  : public Dune::PDELab::DirichletConstraintsParameters
{

public:

  template<typename I>
  bool isDirichlet(const I & ig, const Dune::FieldVector<typename I::ctype, I::dimension-1> & x) const
  {
    return false;
  }

};

int main(int argc, char** argv)
{
  try {
    Dune::MPIHelper& mpi_helper = Dune::MPIHelper::instance(argc,argv);

    typedef Dune::UGGrid<2> Grid;

    Dune::GridPtr<Grid> grid_ptr("unitsquare.dgf");
    Grid& grid = *grid_ptr;

    grid.globalRefine(1);

    grid.mark(1,*grid.leafView().begin<0>());

    grid.preAdapt();
    grid.adapt();
    grid.postAdapt();

    typedef Grid::LeafGridView GV;
    GV gv = grid.leafView();

    typedef GV::ctype DF;
    typedef double RF;

    typedef Dune::PDELab::Q1LocalFiniteElementMap<DF,RF,2> FEM;
    FEM fem;

    typedef Dune::PDELab::ISTLVectorBackend<> VBE;
    typedef Dune::PDELab::HangingNodesDirichletConstraints<
      Grid,
      Dune::PDELab::HangingNodesConstraintsAssemblers::CubeGridQ1Assembler,
      ConstraintsParameters
      > Constraints;
    ConstraintsParameters constraintsParameters;
    Constraints constraints(grid,false,constraintsParameters);

    typedef Dune::PDELab::GridFunctionSpace<GV,FEM,Constraints,VBE> GFS;
    GFS gfs(gv,fem,constraints);
    gfs.name("x");

    typedef Dune::PDELab::BackendVectorSelector<GFS,int>::Type IV;
    IV x(gfs);

    IV::iterator vit = x.begin();
    int i = 0;
    for (GV::template Codim<2>::Iterator it = gv.begin<2>(),
           end = gv.end<2>();
         it != end;
         ++it, ++vit, ++i)
      {
        *vit = i;
      }

    Dune::VTKWriter<GV> vtk_writer(gv);
    Dune::PDELab::addSolutionToVTKWriter(vtk_writer,gfs,x);
    vtk_writer.write("hangingnodes");

    typedef GFS::ConstraintsContainer<RF>::Type CC;
    CC cc;

    Dune::PDELab::constraints(constraintsParameters,gfs,cc);

    for (auto& c : cc)
      {
        std::cout << c.first << ": ";
        for (auto& contrib : c.second)
          std::cout << contrib.first << "=" << contrib.second << " ";
        std::cout << std::endl;
      }

    typedef Dune::PDELab::L2 LOP;
    LOP lop;

    typedef Dune::PDELab::GridOperator<
      GFS,GFS,
      LOP,
      Dune::PDELab::ISTLMatrixBackend,
      RF,RF,RF,CC,CC
      > GridOperator;

    GridOperator grid_operator(gfs,cc,gfs,cc,lop);

    typedef GridOperator::Traits::Domain V;
    typedef GridOperator::Traits::Jacobian M;

    V v(gfs,0.0);
    M m(grid_operator,0.0);

    grid_operator.jacobian(v,m);

    Dune::printmatrix(std::cout,Dune::PDELab::istl::raw(m),"","");

    return 0;
  }
  catch (Dune::Exception& e) {
    std::cerr << "Dune exception: " << e << std::endl;
  }
  catch (std::exception& e) {
    std::cerr << "std::exception: " << e.what() << std::endl;
  }
  catch (...) {
    std::cerr << "unknown exception" << std::endl;
  }
  return 1;
}
