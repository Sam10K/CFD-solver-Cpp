#ifndef Pressure_poisson_H
#define Pressure_poisson_H

#include "headers.hpp"
#include "Schemes.hpp"
#include "Solver_setup.hpp"
#include "SIMPLE_params.hpp"
#include "Mesh.hpp"
#include "Pressure_field.hpp"
#include "Velocity_field.hpp"
#include "Matrix_solvers.hpp"

struct Matrices_poisson
{
  SparseMatrix<double> P_mat;
  VectorXd Bc;
};

struct poisson_eq_coeff
{
  double ac,af,bc;
};

class Pressure_poisson_eq
{
  Pressure Pdash;
  Velocity U;
  Mesh mesh;
  Solver_params solver_setup;
  Schemes schemes;
  SIMPLE_params SIMPLE;

public:

  Pressure_poisson_eq(Pressure p,Velocity u);

  poisson_eq_coeff poisson_matrix(int face_id,int cell_id);

  poisson_eq_coeff rhs_matrix(int face_id,int cell_id);

  Matrices_poisson matrix_assemble();

  Pressure matrix_solve();

};

#endif
