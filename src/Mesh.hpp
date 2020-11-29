#ifndef Mesh_H
#define Mesh_H

#include "headers.hpp"

struct scalarBC
{
  string type;
  double value;
};

struct vectorBC
{
  string type;
  Matrix<double,1,2> value;
};

struct Boundary
{
  string name;
  int face_id_start,face_id_end;
  int nfaces;
  scalarBC pbc;
  vectorBC ubc;
};


class Mesh
{
public:
  int n_faces,n_nodes,nd,n_cells;
  MatrixXi face_node_id,face_cell_id,ubound_type,pbound_type;//size is n_faces      //bound_type:
  MatrixXi cell_face_id,cell_node_id;//size is n_cells                              // fixed=1,zeroGradient=2,wall=3
  vector<vector<int>> node_cell_id,node_face_id; //size is n_nodes
  MatrixXd Nodes,centroid,ubound_value;
  VectorXd volume,pbound_value;
  vector <Boundary> Boundaries;
  Boundary Interior;

  void Read_Mesh();

  void Read_P_BC();

  void Read_U_BC();

  Matrix<double,1,2> mid(int face_id);

  double area(int face_id);

  Matrix<double,1,2> normal(int face_id,int owner);

  int neighb(int face_id,int owner);

};

#endif
