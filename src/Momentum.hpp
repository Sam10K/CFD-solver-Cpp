#ifndef Momentum_H
#define Momentum_H


#include "headers.hpp"
#include "Schemes.hpp"
#include "Solver_setup.hpp"
#include "SIMPLE_params.hpp"
#include "Mesh.hpp"
#include "Pressure_field.hpp"
#include "Velocity_field.hpp"
#include "Matrix_solvers.hpp"


struct Matrices
{
  SparseMatrix<double> U_mat,V_mat;
  MatrixXd Bc;
};

struct eq_coeff
{
  Matrix<double,1,2> ac,af,bc;
};

class Momentum_eq
{
  Pressure P;
  Velocity U;
  Mesh mesh;
  Solver_params solver_setup;
  Schemes schemes;
  SIMPLE_params SIMPLE;

public:

  Momentum_eq(Pressure p,Velocity u);

  eq_coeff diffusion_matrix(int face_id,int cell_id);

  eq_coeff convection_matrix(int face_id,int cell_id);

  eq_coeff pressure_gradient(int cell_id);

  Matrices matrix_assemble();

  Velocity matrix_solve();

};

#endif
