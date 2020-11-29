#include "Velocity_field.hpp"

Velocity::Velocity(Pressure p)
{
  mesh = p.get_mesh();
  solver_setup = p.get_solver_params();
  schemes = p.get_schemes();
  SIMPLE = p.get_simple_params();
  P = p;
  MatrixXd dummy(p.get_mesh().n_cells,2);
  dummy<<p.get_mesh().volume,p.get_mesh().volume;
  Ap = dummy;
  initialize();
}

void Velocity::initialize()
{

  ifstream fin;
  fin.open(file_bc_u);

  string line;

  char dummy[20];
  double dummy_f1,dummy_f2;

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

            if(sscanf(buffer, "list_values (%lf %lf);", &dummy_f1, &dummy_f2) && strlen(buffer)!=0)
            {
              VectorXd u,v;
              u.setOnes(mesh.n_cells);
              v.setOnes(mesh.n_cells);
              u = u*dummy_f1;
              v = v*dummy_f2;
              MatrixXd vel(mesh.n_cells,2);
              vel<<u,v;
              this->velocity = vel;
              break;
            }
          }
        }
      }
      break;
    }
  }

  // Routine to initially calculate the facial velocities
  MatrixXd vel_f(mesh.n_faces,4);
  for(int i=0;i<mesh.n_faces;i++)
  {
    int cell1 = mesh.face_cell_id(i,0);
    int cell2 = mesh.face_cell_id(i,1);
    if(cell2==-1)
    cell2=cell1;

    RowVectorXd val_f1 = face_interp_vel(i,cell1);
    RowVectorXd val_f2 = face_interp_vel(i,cell2);

    vel_f.row(i)<<val_f1,val_f2;

  }

  this->velocity_f = vel_f;
  set_c(this->velocity); // To calculate other face parameters

  fin.close();

}


Matrix<double,1,2> Velocity::face_interp_vel(int face_id,int cell_id)
{
  Matrix<double,1,2> mid = mesh.mid(face_id);
  Matrix<double,1,2> centroid = mesh.centroid.row(cell_id);
  int neighb_id = mesh.neighb(face_id,cell_id);
  Matrix<double,1,2> neighb_nodes = (neighb_id!=-1)?mesh.centroid.row(neighb_id):mid;
  Matrix<double,1,2> val_c = velocity.row(cell_id);
  Matrix<double,1,2> val_n;

  if(neighb_id!=-1)
  val_n = velocity.row(neighb_id);

  else
  {
    int ubound_type = mesh.ubound_type(face_id);

    if(ubound_type==1)
    val_n = mesh.ubound_value.row(face_id);

    else if(ubound_type==2)
    val_n = val_c;

    else if(ubound_type==3)
    val_n = mesh.ubound_value.row(face_id);


  }

  Matrix<double,1,2> d_f_n = mid-neighb_nodes;
  Matrix<double,1,2> d_c_n = centroid-neighb_nodes;

  double gc = d_f_n.norm()/d_c_n.norm();
  double gn = 1.0-gc;

  Matrix<double,1,2> val_f_bar = gc*val_c + gn*val_n;

  return val_f_bar;

}



void Velocity::set_c(MatrixXd vel)
{
  velocity = vel;

  MatrixXd vel_f(mesh.n_faces,4);
  for(int i=0;i<mesh.n_faces;i++)
  {
    int cell1 = mesh.face_cell_id(i,0);
    int cell2 = mesh.face_cell_id(i,1);
    if(cell2==-1)
    cell2=cell1;

    Matrix<double,1,2> dummy3 = calc_vel_f(i,cell1);
    Matrix<double,1,2> dummy4 = calc_vel_f(i,cell2);
    vel_f.row(i)<<dummy3,dummy4;

  }

  velocity_f = vel_f;

  MatrixXd vel_grad(mesh.n_cells,4);
  for(int i=0;i<mesh.n_cells;i++)
  {
    Matrix<double,2,2> dummy = calc_grad_vel_c(i);
    vel_grad.row(i) << dummy(0,0),dummy(0,1),dummy(1,0),dummy(1,1);
  }

  grad_vel_c = vel_grad;

  MatrixXd velgrad_f(mesh.n_faces,8);
  for(int i=0;i<mesh.n_faces;i++)
  {
    int cell1 = mesh.face_cell_id(i,0);
    int cell2 = mesh.face_cell_id(i,1);
    if(cell2==-1)
    cell2=cell1;

    Matrix<double,2,2> dummy1 = calc_grad_vel_f(i,cell1);
    Matrix<double,2,2> dummy2 = calc_grad_vel_f(i,cell2);
    velgrad_f.row(i)<<dummy1.row(0),dummy1.row(1),dummy2.row(0),dummy2.row(1);
  }

  grad_vel_f = velgrad_f;
}

Matrix<double,1,2> Velocity::get_f(int face_id,int cell_id)
{
  //RowVectorXd vel_f = calc_vel_f(face_id,cell_id);
  //return vel_f;

  RowVectorXi ids = mesh.face_cell_id.row(face_id);
  Matrix<double,1,2> dummy;

  if(cell_id==ids(0))
  dummy<<velocity_f(face_id,0),velocity_f(face_id,1);
  else
  dummy<<velocity_f(face_id,2),velocity_f(face_id,3);

  return dummy;
}

Matrix<double,2,2> Velocity::grad_c(int cell_id)
{
  //MatrixXd grad_vel_c = calc_grad_vel_c(cell_id);
  //return grad_vel_c;
  Matrix<double,2,2> dummy;
  dummy<< grad_vel_c(cell_id,0),grad_vel_c(cell_id,1),grad_vel_c(cell_id,2),grad_vel_c(cell_id,3);
  return dummy;
}

Matrix<double,2,2> Velocity::grad_f(int face_id,int cell_id)
{
  //MatrixXd grad_vel_f = calc_grad_vel_f(face_id,cell_id);
  //return grad_vel_f;
  RowVectorXi ids = mesh.face_cell_id.row(face_id);
  Matrix<double,2,2> dummy;

  if(cell_id==ids(0))
  dummy<<grad_vel_f(face_id,0),grad_vel_f(face_id,1),grad_vel_f(face_id,2),grad_vel_f(face_id,3);
  else
  dummy<<grad_vel_f(face_id,4),grad_vel_f(face_id,5),grad_vel_f(face_id,6),grad_vel_f(face_id,7);

  return dummy;
}

double Velocity::get_continuity()
{
  double cont = 0.0;
  for(int i=0;i<mesh.n_cells;i++)
  {
    double cell_cont = calc_cell_continuity(i);
    double vol = mesh.volume(i);
    cont+=cell_cont*vol;

  }

  cont /= mesh.volume.sum();

  return cont;
}


Matrix<double,1,2> Velocity::calc_vel_f(int face_id,int cell_id)
{
  Matrix<double,1,2> mid = mesh.mid(face_id);
  Matrix<double,1,2> centroid = mesh.centroid.row(cell_id);
  int neighb_id = mesh.neighb(face_id,cell_id);
  Matrix<double,1,2> neighb_nodes = (neighb_id!=-1)?mesh.centroid.row(neighb_id):mid;

  Matrix<double,1,2> d_f_n = mid-neighb_nodes;
  Matrix<double,1,2> d_c_n = centroid-neighb_nodes;

  double gc = d_f_n.norm()/d_c_n.norm();
  double gn = 1.0-gc;

  Matrix<double,1,2> val_f_bar = face_interp_vel(face_id,cell_id);

  //Rhiew Chow interpolation

  Matrix<double,1,2> del_p = P.grad_f(face_id,cell_id);

  Matrix<double,1,2> pgrad_c = P.grad_c(cell_id);
  Matrix<double,1,2> pgrad_n = (neighb_id!=-1)?P.grad_c(neighb_id):pgrad_c;

  Matrix<double,1,2> del_p_bar = gc*pgrad_c + gn*pgrad_n;

  Matrix<double,1,2> ap_c = Ap.row(cell_id)/mesh.volume(cell_id);
  Matrix<double,1,2> ap_n = (neighb_id!=-1)?Ap.row(neighb_id)/mesh.volume(neighb_id):ap_c;
  Matrix<double,1,2> ap_f = gc*ap_c + gn*ap_n;

  RowVectorXi ids = mesh.face_cell_id.row(face_id);
  Matrix<double,1,2> val_f_old;

  if(cell_id==ids(0))
  val_f_old<<velocity_f(face_id,0),velocity_f(face_id,1);
  else
  val_f_old<<velocity_f(face_id,2),velocity_f(face_id,3);


  Matrix<double,1,2> val_f = val_f_bar - (del_p-del_p_bar).cwiseQuotient(ap_f) + (1-SIMPLE.alpha_u)*(val_f_old-val_f_bar);

  return val_f;
}


Matrix<double,2,2> Velocity::calc_grad_vel_c(int cell_id)
{
  string grad_scheme = schemes.grad_scheme;
  RowVectorXi face_id = mesh.cell_face_id.row(cell_id);
  Matrix<double,1,2> centroid = mesh.centroid.row(cell_id);
  double vol = mesh.volume(cell_id);
  int nd = mesh.nd;

  Matrix<double, 2, 2> grad_val_c;

  if(grad_scheme.compare("LSQ")==0)
  {
    double a11,a12,a21,a22,b1_u,b2_u,b1_v,b2_v;
    a11=0.0;a12=0.0;a21=0.0;a22=0.0;b1_u=0.0;b2_u=0.0;b1_v=0.0;b2_v=0.0;
    for(int i=0;i<nd;i++)
    {
      int id = face_id(i);
      Matrix<double,1,2> mid = mesh.mid(id);
      int neighb_id = mesh.neighb(id,cell_id);
      Matrix<double,1,2> neighb_nodes = (neighb_id!=-1)?mesh.centroid.row(neighb_id):mid;
      Matrix<double,1,2> val_c = velocity.row(cell_id);
      Matrix<double,1,2> val_n = (neighb_id!=-1)?velocity.row(neighb_id):get_f(id,cell_id);

      double delta_x = neighb_nodes(0)-centroid(0);
      double delta_y = neighb_nodes(1)-centroid(1);
      double w = 1/(neighb_nodes-centroid).norm();
      Matrix<double,1,2> delta_val = val_n-val_c;

      a11 = a11 + w*delta_x*delta_x;
      a12 = a12 + w*delta_x*delta_y;
      a21 = a21 + w*delta_y*delta_x;
      a22 = a22 + w*delta_y*delta_y;

      b1_u = b1_u + w*delta_x*delta_val(0);
      b2_u = b2_u + w*delta_y*delta_val(0);

      b1_v = b1_v + w*delta_x*delta_val(1);
      b2_v = b2_v + w*delta_y*delta_val(1);


    }


    grad_val_c(1,0) = (b2_u*a11 - a21*b1_u) / (a22*a11 - a21*a12); //dudy
    grad_val_c(0,0) = (b1_u - a12*grad_val_c(1,0))/a11; //dudx

    grad_val_c(1,1) = (b2_v*a11 - a21*b1_v) / (a22*a11 - a21*a12); //dvdy
    grad_val_c(0,1) = (b1_v - a12*grad_val_c(1,1))/a11; //dvdx

    return grad_val_c;

  }

  if(grad_scheme.compare("GN")==0)
  {
    grad_val_c.setZero();
    for(int i=0;i<nd;i++)
    {
      int fid = face_id(i);
      Matrix<int,1,2> node_id = mesh.face_node_id.row(fid);
      Matrix<double,1,2> val_node_1 = get_n(node_id(0));
      Matrix<double,1,2> val_node_2 = get_n(node_id(1));
      Matrix<double,1,2> val_f = 0.5*(val_node_1 + val_node_2);
      double ds = mesh.area(fid);
      Matrix<double,1,2> n = mesh.normal(fid,cell_id);
      grad_val_c.col(0) += (val_f(0)*n*ds).transpose(); //[dudx;dudy]
      grad_val_c.col(1) += (val_f(1)*n*ds).transpose(); //[dvdx;dvdy]

    }

    grad_val_c /= vol;

    return grad_val_c;

  }
}


Matrix<double,2,2> Velocity::calc_grad_vel_f(int face_id,int cell_id)
{
  Matrix<double,1,2> mid = mesh.mid(face_id);
  Matrix<double,1,2> centroid = mesh.centroid.row(cell_id);
  int neighb_id = mesh.neighb(face_id,cell_id);
  Matrix<double,1,2> neighb_nodes = (neighb_id!=-1)?mesh.centroid.row(neighb_id):mid;
  Matrix<double,2,2> grad_val_c = grad_c(cell_id);
  Matrix<double,2,2> grad_val_n = (neighb_id!=-1)?grad_c(neighb_id):grad_val_c;

  Matrix<double,1,2> d_f_n = mid-neighb_nodes;
  Matrix<double,1,2> d_c_n = neighb_nodes-centroid;
  double gc = d_f_n.norm()/d_c_n.norm();
  double gn = 1.0-gc;

  Matrix<double,2,2> dummy_zero; dummy_zero.setZero();

  Matrix<double,2,2> grad_vel_f_avg = (neighb_id!=-1)?(gc*grad_val_c + gn*grad_val_n):dummy_zero;

  Matrix<double,1,2> n = mesh.normal(face_id,cell_id);
  double ds = mesh.area(face_id);
  Matrix<double,1,2> S = n*ds;
  Matrix<double,1,2> d = (neighb_id!=-1)?d_c_n:(S*S.dot(mid-centroid)/S.squaredNorm());
  Matrix<double,1,2> e = d/d.norm();
  Matrix<double,1,2> val_c = get_c(cell_id);
  Matrix<double,1,2> val_n = (neighb_id!=-1)?velocity.row(neighb_id):get_f(face_id,cell_id);

  Matrix<double,2,2> grad_vel_f_corr;

  Matrix<double,1,2> del_u;del_u<<grad_vel_f_avg(0,0),grad_vel_f_avg(1,0);//[dudx,dudy]
  Matrix<double,1,2> del_v;del_v<<grad_vel_f_avg(0,1),grad_vel_f_avg(1,1);//[dvdx,dvdy]

  grad_vel_f_corr(0,0) = grad_vel_f_avg(0,0) + e(0)*((val_n(0)-val_c(0))/d.norm() - e.dot(del_u));//dudx
  grad_vel_f_corr(1,0) = grad_vel_f_avg(1,0) + e(1)*((val_n(0)-val_c(0))/d.norm() - e.dot(del_u));//dudy
  grad_vel_f_corr(0,1) = grad_vel_f_avg(0,1) + e(0)*((val_n(1)-val_c(1))/d.norm() - e.dot(del_v));//dvdx
  grad_vel_f_corr(1,1) = grad_vel_f_avg(1,1) + e(1)*((val_n(1)-val_c(1))/d.norm() - e.dot(del_v));//dvdy

  return grad_vel_f_corr;

}

double Velocity::calc_cell_continuity(int cell_id)
{
  int nd = mesh.nd;
  RowVectorXi face_id = mesh.cell_face_id.row(cell_id);
  double cell_cont = 0.0;
  for(int i=0;i<nd;i++)
  {
    Matrix<double,1,2> n = mesh.normal(face_id(i),cell_id);
    Matrix<double,1,2> vel = get_f(face_id(i),cell_id);
    double ds = mesh.area(face_id(i));
    cell_cont += vel.dot(n)*ds;

  }

  return abs(cell_cont);

}


Matrix<double,1,2> Velocity::calc_nodal_velocity(int node_id)
{
  vector<int> cell_ids = mesh.node_cell_id[node_id];

  Matrix<double,1,2> numerator;numerator.setZero();
  double denominator = 0.0;

  for(int i=0;i<cell_ids.size();i++)
  {
    Matrix<double,1,2> centroid = mesh.centroid.row(cell_ids[i]);
    Matrix<double,1,2> node = mesh.Nodes.row(node_id);
    double inv_dist = 1/(centroid-node).norm();
    Matrix<double,1,2> vel_c = get_c(cell_ids[i]);

    numerator += vel_c*inv_dist;
    denominator += inv_dist;

  }

  Matrix<double,1,2> nodal_vel = numerator/denominator;

  vector<int> face_ids = mesh.node_face_id[node_id];

  for(int i=0;i<face_ids.size();i++)
  {
    if(mesh.ubound_type(face_ids[i])!=0)
    {
      if(mesh.ubound_type(face_ids[i])==1)
      {nodal_vel = mesh.ubound_value.row(face_ids[i]);break;}

      else if(mesh.ubound_type(face_ids[i])==2)
      {nodal_vel = get_c(mesh.neighb(face_ids[i],-1));break;}

      else if(mesh.ubound_type(face_ids[i])==3)
      {nodal_vel = mesh.ubound_value.row(face_ids[i]);break;}

    }
  }

  return nodal_vel;

}
