#include "Matrix_solvers.hpp"

Matrix_solvers::Matrix_solvers(SparseMatrix<double> A,VectorXd b,VectorXd x0,matrix_solver params)
{
  this->A = A;
  this->b = b;
  this->x0 = x0;
  this->maxIter = params.max_iter;
  this->tol = params.tol;

  this->solver = params.solver;

}

solution_info Matrix_solvers::solve()
{
  solution_info info;
  if((solver.compare("SOR")==0)||(solver.compare("GS")==0))
  info = gauss_seidel();

  if(solver.compare("BICGSTAB")==0)
  info = bicgstab();


  return info;

}

solution_info Matrix_solvers::gauss_seidel()
{
  SparseMatrix<double> U = A.triangularView<StrictlyUpper>();
  VectorXd x_new,x_old;
  x_old=x0;
  double diff=1;int iter=0;

  while((diff>tol)&&(iter<maxIter))
  {
    x_new = A.triangularView<Lower>().solve(b-U*x_old);
    diff = (x_old-x_new).norm()/x_new.norm();

    x_old = x_new;
    iter++;

  }

  solution_info info;
  info.x = x_new;
  info.iterations = iter;
  info.residual = diff;

  return info;

}

solution_info Matrix_solvers::bicgstab()
{
  BiCGSTAB<SparseMatrix<double>,IncompleteLUT<double> > solver;
  solver.setMaxIterations(maxIter);
  solver.setTolerance(tol);

  solver.compute(A);
  VectorXd x_new = solver.solveWithGuess(b,x0);


  solution_info info;
  info.x = x_new;
  info.iterations = solver.iterations();
  info.residual = solver.error();

  return info;


}
