#include "Post_processing.hpp"

Post_process::Post_process(Pressure p,Velocity u)
{
  mesh = p.get_mesh();
  solver_setup = p.get_solver_params();
  schemes = p.get_schemes();
  SIMPLE = p.get_simple_params();
  P = p;
  U = u;
}

void Post_process::writeVTK(int folder)
{
  ofstream fout;

  //char file[20];
  //sprintf(file,"%d",folder);
  int status = mkdir("Data",0777);


  fout.open("./Data/Data_"+to_string(folder)+".vtk");
  fout<<"# vtk DataFile Version 2.0"<<endl;
  fout<<"Data\nASCII"<<endl<<endl;
  fout<<"DATASET UNSTRUCTURED_GRID"<<endl;
  fout<<"POINTS "<<mesh.n_nodes<<" double"<<endl;
  MatrixXd Nodes_data(mesh.n_nodes,3);Nodes_data.setZero();
  Nodes_data.col(0)=mesh.Nodes.col(0);
  Nodes_data.col(1)=mesh.Nodes.col(1);
  fout<<Nodes_data<<endl;

  fout<<"CELLS "<<mesh.n_cells<<" "<<mesh.n_cells*(mesh.nd+1)<<endl;
  MatrixXi Cell_data(mesh.n_cells,mesh.nd+1);
  VectorXi dummy(mesh.n_cells);dummy.setOnes();dummy = dummy*mesh.nd;
  Cell_data << dummy,mesh.cell_node_id;
  fout<<Cell_data<<endl;

  fout<<"CELL_TYPES "<<mesh.n_cells<<endl;
  VectorXi Cell_type(mesh.n_cells);Cell_type.setOnes();
  Cell_type = (mesh.nd==4)?Cell_type*9:Cell_type*5; // 9 is quad, 5 is tri for VTK
  fout<<Cell_type<<endl;

  fout<<"CELL_DATA "<<mesh.n_cells<<endl;
  fout<<"SCALARS P double"<<endl;
  fout<<"LOOKUP_TABLE default"<<endl;
  fout<<P.get_pressure()<<endl;

  MatrixXd vel(mesh.n_cells,3);vel.setZero();
  vel.col(0) = U.get_vel().col(0);
  vel.col(1) = U.get_vel().col(1);
  fout<<"VECTORS Velocity double"<<endl;
  fout<<vel<<endl;

  fout.close();
}
