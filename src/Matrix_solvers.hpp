#ifndef Matrix_solvers_H
#define Matrix_solvers_H

#include "headers.hpp"
#include "Schemes.hpp"

struct solution_info
{
  int iterations;
  double residual;
  VectorXd x;
};

class Matrix_solvers
{
  SparseMatrix<double> A;
  VectorXd b;
  VectorXd x0;
  int maxIter;
  double tol;
  string solver;

public:

  Matrix_solvers(SparseMatrix<double> A,VectorXd b,VectorXd x0,matrix_solver params);

  solution_info solve(); // General function which calls other solvers

  solution_info gauss_seidel();

  solution_info bicgstab();

};























#endif
