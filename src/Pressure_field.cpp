#include "Pressure_field.hpp"

Pressure::Pressure(Solver_params solv,Schemes sch,SIMPLE_params simp,Mesh msh)
{
  mesh = msh;
  solver_setup = solv;
  schemes = sch;
  SIMPLE = simp;
  initialize();
}


void Pressure::initialize()
{

  ifstream fin;
  fin.open(file_bc_p);

  string line;

  char dummy[20];
  double dummy_f;

  while(fin)
  {
    getline(fin,line);

    if(line.compare("internal field")==0)
    {
      while(fin)
      {
        getline(fin,line);
        char buffer[line.length()+1];
        strcpy(buffer,line.c_str());

        if(sscanf(buffer, "list_type %[^;]", dummy) && strlen(buffer)!=0)
        {
          if(strcmp(dummy,"uniform")==0)
          {
            getline(fin,line);
            char buffer[line.length()+1];
            strcpy(buffer,line.c_str());

            if(sscanf(buffer, "list_values %lf;", &dummy_f) && strlen(buffer)!=0)
            {
              VectorXd p;
              p.setOnes(mesh.n_cells);
              p = p*dummy_f;
              set_c(p);
              break;
            }
          }
        }
      }
      break;
    }
  }

}


void Pressure::set_c(VectorXd p)
{
  this->p_c = p;

  MatrixXd face_p(mesh.n_faces,2);
  for(int i=0;i<mesh.n_faces;i++)
  {
    int cell1 = mesh.face_cell_id(i,0);
    int cell2 = mesh.face_cell_id(i,1);
    if(cell2==-1)
    cell2=cell1;

    double dummy3 = calc_p_f(i,cell1);
    double dummy4 = calc_p_f(i,cell2);
    face_p.row(i)<<dummy3,dummy4;
  }

  this->p_f=face_p;

  MatrixXd pgrad(mesh.n_cells,2);
  for(int i=0;i<mesh.n_cells;i++)
  pgrad.row(i) = calc_grad_p_c(i);

  this->grad_p_c = pgrad;

  MatrixXd pgrad_f(mesh.n_faces,4);
  for(int i=0;i<mesh.n_faces;i++)
  {
    int cell1 = mesh.face_cell_id(i,0);
    int cell2 = mesh.face_cell_id(i,1);
    if(cell2==-1)
    cell2=cell1;

    Matrix<double,1,2> dummy1 = calc_grad_p_f(i,cell1);
    Matrix<double,1,2> dummy2 = calc_grad_p_f(i,cell2);
    pgrad_f.row(i)<<dummy1,dummy2;

  }

  this->grad_p_f = pgrad_f;

}


double Pressure::get_f(int face_id,int cell_id)
{
  //double p_f = calc_p_f(face_id,cell_id);
  //return p_f;
  RowVectorXi ids = mesh.face_cell_id.row(face_id);
  double dummy;

  if(cell_id==ids(0))
  dummy = p_f(face_id,0);
  else
  dummy = p_f(face_id,1);

  return dummy;

}


Matrix<double,1,2> Pressure::grad_f(int face_id,int cell_id)
{
  //RowVectorXd grad_p_f = calc_grad_p_f(face_id,cell_id);
  //return grad_p_f;
  RowVectorXi ids = mesh.face_cell_id.row(face_id);
  RowVectorXd dummy(2);

  if(cell_id==ids(0))
  dummy<<grad_p_f(face_id,0),grad_p_f(face_id,1);
  else
  dummy<<grad_p_f(face_id,2),grad_p_f(face_id,3);

  return dummy;
}


double Pressure::calc_p_f(int face_id,int cell_id)
{
  Matrix<double,1,2> mid = mesh.mid(face_id);
  Matrix<double,1,2> centroid = mesh.centroid.row(cell_id);
  int neighb_id = mesh.neighb(face_id,cell_id);
  Matrix<double,1,2> neighb_nodes = (neighb_id!=-1)?mesh.centroid.row(neighb_id):mid;
  double val_c = p_c(cell_id);
  double val_n;

  if(neighb_id!=-1)
  val_n = p_c(neighb_id);

  else
  {
    int pbound_type = mesh.pbound_type(face_id);

    if(pbound_type==1)
    val_n = mesh.pbound_value(face_id);

    else if(pbound_type==2)
    val_n = val_c;

    else if(pbound_type==3)
    val_n = val_c;

  }

  Matrix<double,1,2> d_f_n = mid-neighb_nodes;
  Matrix<double,1,2> d_c_n = centroid-neighb_nodes;

  double gc = d_f_n.norm()/d_c_n.norm();
  double gn = 1.0-gc;

  double val_f = gc*val_c + gn*val_n;

  return val_f;
}


Matrix<double,1,2> Pressure::calc_grad_p_c(int cell_id)
{
  string grad_scheme = schemes.grad_scheme;
  RowVectorXi face_id = mesh.cell_face_id.row(cell_id);
  Matrix<double,1,2> centroid = mesh.centroid.row(cell_id);
  double vol = mesh.volume(cell_id);
  int nd = mesh.nd;

  Matrix<double, 1, 2> grad_val_c;

  if(grad_scheme.compare("LSQ")==0)
  {
    double a11,a12,a21,a22,b1,b2;
    a11=0.0;a12=0.0;a21=0.0;a22=0.0;b1=0.0;b2=0.0;
    for(int i=0;i<nd;i++)
    {
      int id = face_id(i);
      Matrix<double,1,2> mid = mesh.mid(id);
      int neighb_id = mesh.neighb(id,cell_id);
      Matrix<double,1,2> neighb_nodes = (neighb_id!=-1)?mesh.centroid.row(neighb_id):mid;
      double val_c = p_c(cell_id);
      double val_n = (neighb_id!=-1)?p_c(neighb_id):get_f(id,cell_id);

      double delta_x = neighb_nodes(0)-centroid(0);
      double delta_y = neighb_nodes(1)-centroid(1);
      double w = 1/(sqrt(pow(delta_x,2)+pow(delta_y,2)));
      double delta_val = val_n-val_c;

      a11 = a11 + w*delta_x*delta_x;
      a12 = a12 + w*delta_x*delta_y;
      a21 = a21 + w*delta_y*delta_x;
      a22 = a22 + w*delta_y*delta_y;

      b1 = b1 + w*delta_x*delta_val;
      b2 = b2 + w*delta_y*delta_val;

    }

    grad_val_c(1) = (b2*a11 - a21*b1) / (a22*a11 - a21*a12);
    grad_val_c(0) = (b1 - a12*grad_val_c(1))/a11;

    return grad_val_c;

  }
  if(grad_scheme.compare("GN")==0)
  {
    grad_val_c.setZero();
    for(int i=0;i<nd;i++)
    {
      int fid = face_id(i);
      Matrix<int,1,2> node_id = mesh.face_node_id.row(fid);
      double val_node_1 = get_n(node_id(0));
      double val_node_2 = get_n(node_id(1));
      double val_f = 0.5*(val_node_1 + val_node_2);
      double ds = mesh.area(fid);
      Matrix<double,1,2> n = mesh.normal(fid,cell_id);
      grad_val_c += val_f*n*ds;

    }

    grad_val_c /= vol;

    return grad_val_c;

  }
}

Matrix<double,1,2> Pressure::calc_grad_p_f(int face_id,int cell_id)
{
  Matrix<double,1,2> mid = mesh.mid(face_id);
  Matrix<double,1,2> centroid = mesh.centroid.row(cell_id);
  int neighb_id = mesh.neighb(face_id,cell_id);
  Matrix<double,1,2> neighb_nodes = (neighb_id!=-1)?mesh.centroid.row(neighb_id):mid;
  Matrix<double,1,2> grad_val_c = grad_c(cell_id);
  Matrix<double,1,2> grad_val_n = (neighb_id!=-1)?grad_c(neighb_id):grad_val_c;

  Matrix<double,1,2> d_f_n = mid-neighb_nodes;
  Matrix<double,1,2> d_c_n = neighb_nodes-centroid;
  double gc = d_f_n.norm()/d_c_n.norm();
  double gn = 1.0-gc;

  Matrix<double,1,2> dummy_zero; dummy_zero<<0.0,0.0;

  Matrix<double,1,2> grad_p_f_avg = (neighb_id!=-1)?(gc*grad_val_c + gn*grad_val_n):dummy_zero;

  Matrix<double,1,2> n = mesh.normal(face_id,cell_id);
  double ds = mesh.area(face_id);
  Matrix<double,1,2> S = n*ds;
  Matrix<double,1,2> d = (neighb_id!=-1)?d_c_n:(S*S.dot(mid-centroid)/S.squaredNorm());
  Matrix<double,1,2> e = d/d.norm();
  double val_c = get_c(cell_id);
  double val_n = (neighb_id!=-1)?p_c(neighb_id):get_f(face_id,cell_id);

  Matrix<double,1,2> grad_p_f_corr;

  grad_p_f_corr(0) = grad_p_f_avg(0) + e(0)*((val_n-val_c)/d.norm() - e.dot(grad_p_f_avg));
  grad_p_f_corr(1) = grad_p_f_avg(1) + e(1)*((val_n-val_c)/d.norm() - e.dot(grad_p_f_avg));

  return grad_p_f_corr;

}


double Pressure::calc_nodal_pressure(int node_id)
{
  vector<int> cell_ids = mesh.node_cell_id[node_id];

  double numerator=0.0;
  double denominator = 0.0;

  for(int i=0;i<cell_ids.size();i++)
  {
    Matrix<double,1,2> centroid = mesh.centroid.row(cell_ids[i]);
    Matrix<double,1,2> node = mesh.Nodes.row(node_id);
    double inv_dist = 1/(centroid-node).norm();
    double pressure_c = get_c(cell_ids[i]);

    numerator += pressure_c*inv_dist;
    denominator += inv_dist;

  }

  double nodal_pressure = numerator/denominator;

  vector<int> face_ids = mesh.node_face_id[node_id];

  for(int i=0;i<face_ids.size();i++)
  {
    if(mesh.pbound_type(face_ids[i])!=0)
    {
      if(mesh.pbound_type(face_ids[i])==1)
      {nodal_pressure = mesh.pbound_value(face_ids[i]);break;}

      else if(mesh.pbound_type(face_ids[i])==2)
      {nodal_pressure = get_c(mesh.neighb(face_ids[i],-1));break;}

      else if(mesh.pbound_type(face_ids[i])==3)
      {nodal_pressure = get_c(mesh.neighb(face_ids[i],-1));break;}

    }
  }

  return nodal_pressure;

}
