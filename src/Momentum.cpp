#include "Momentum.hpp"

Momentum_eq::Momentum_eq(Pressure p,Velocity u)
{
  P = p;
  U = u;
  mesh = p.get_mesh();
  solver_setup = p.get_solver_params();
  schemes = p.get_schemes();
  SIMPLE = p.get_simple_params();
}

eq_coeff Momentum_eq::diffusion_matrix(int face_id,int cell_id)
{
  Matrix<double,1,2> n = mesh.normal(face_id,cell_id);
  double ds = mesh.area(face_id);
  Matrix<double,1,2> S = ds*n;
  Matrix<double,1,2> mid = mesh.mid(face_id);
  int neighb_id = mesh.neighb(face_id,cell_id);
  Matrix<double,1,2> centroid = mesh.centroid.row(cell_id);
  Matrix<double,1,2> neighb = (neighb_id!=-1)?mesh.centroid.row(neighb_id):mid;
  Matrix<double,1,2> d = neighb-centroid;
  if(neighb_id==-1)
  d = S.dot(d)*S/S.squaredNorm();

  Matrix<double,1,2> delta = S.squaredNorm()*d/S.dot(d);
  Matrix<double,1,2> k = S-delta;
  double mu = solver_setup.mu;
  Matrix<double,2,2> grad_vel_f = U.grad_f(face_id,cell_id);
  Matrix<double,1,2> vel_f = U.get_f(face_id,cell_id);
  Matrix<double,1,2> vel_c = U.get_c(cell_id);

  Matrix<double,1,2> ac,af,bc;
  ac << mu*delta.norm()/d.norm(),mu*delta.norm()/d.norm();
  af << -mu*delta.norm()/d.norm(),-mu*delta.norm()/d.norm();
  bc << (mu*grad_vel_f*k.transpose()).transpose() + (mu*grad_vel_f.transpose()*S.transpose()).transpose(); //Make it into a row matrix

  if(neighb_id==-1)
  {
    int ubound_type = mesh.ubound_type(face_id);
    Matrix<double,1,2> ubound = mesh.ubound_value.row(face_id);

    if(ubound_type==1)
    {
      bc = bc + -af.cwiseProduct(ubound);
      af.setZero(); // ac is the same as previously set
    }

    else if(ubound_type==2)
    {
      ac.setZero();
      af.setZero();
      bc.setZero();
    }

    else if(ubound_type==3)
    {
      ac << mu*(1-n(0)*n(0))*S.norm()/d.norm(),mu*(1-n(1)*n(1))*S.norm()/d.norm();
      af.setZero();
      bc(0) = mu*S.norm()/d.norm()*(ubound(0)*(1-n(0)*n(0)) + (vel_c(1)-ubound(1))*n(0)*n(1));
      bc(1) = mu*S.norm()/d.norm()*(ubound(1)*(1-n(1)*n(1)) + (vel_c(0)-ubound(0))*n(0)*n(1));

    }

  }


  eq_coeff diff;
  diff.ac = ac;
  diff.af = af;
  diff.bc = bc;

  return diff;

}

eq_coeff Momentum_eq::convection_matrix(int face_id,int cell_id)
{
  Matrix<double,1,2> n = mesh.normal(face_id,cell_id);
  double ds = mesh.area(face_id);
  Matrix<double,1,2> vel_f = U.get_f(face_id,cell_id);
  double rho = solver_setup.rho;
  double mdot_f = rho*vel_f.dot(n)*ds;
  int neighb_id = mesh.neighb(face_id,cell_id);
  int upwind = (mdot_f>0)?cell_id:neighb_id;

  if(upwind==-1)
  upwind = cell_id;

  Matrix<double,2,2> grad_vel_up = U.grad_c(upwind);
  Matrix<double,2,2> grad_vel_f = U.grad_f(face_id,cell_id);

  Matrix<double,1,2> mid = mesh.mid(face_id);
  Matrix<double,1,2> centroid = mesh.centroid.row(upwind);

  Matrix<double,1,2> dcf = mid-centroid;

  string div_scheme = schemes.div_scheme;

  double c1,c2;

  if(div_scheme.compare("UD")==0)
  {
    c1=0.0;c2=0.0;
  }
  else if(div_scheme.compare("CD")==0)
  {
    c1=0.0;c2=1.0;
  }
  else if(div_scheme.compare("SOU")==0)
  {
    c1=2.0;c2=-1.0;
  }
  else if(div_scheme.compare("FROMM")==0)
  {
    c1=1.0;c2=0.0;
  }
  else //QUICK
  {
    c1=0.5;c2=0.5;
  }

  Matrix<double,1,2> ac,af,bc;
  ac << max(mdot_f,double(0)),max(mdot_f,double(0));
  af << -max(-mdot_f,double(0)),-max(-mdot_f,double(0));

  Matrix<double,1,2> del_u_up;del_u_up << grad_vel_up(0,0),grad_vel_up(1,0);
  Matrix<double,1,2> del_v_up;del_v_up << grad_vel_up(0,1),grad_vel_up(1,1);
  Matrix<double,1,2> del_u_f;del_u_f << grad_vel_f(0,0),grad_vel_f(1,0);
  Matrix<double,1,2> del_v_f;del_v_f << grad_vel_f(0,1),grad_vel_f(1,1);

  bc << -mdot_f*(c1*del_u_up + c2*del_u_f).dot(dcf),-mdot_f*(c1*del_v_up + c2*del_v_f).dot(dcf);

  if(neighb_id==-1)
  {
    int ubound_type = mesh.ubound_type(face_id);
    Matrix<double,1,2> ubound = mesh.ubound_value.row(face_id);

    if(ubound_type==1)
    {
      ac.setZero();
      af.setZero();
      mdot_f = rho*ubound.dot(n)*ds;
      bc = -mdot_f*ubound;
    }

    else if(ubound_type==2)
    {
      ac<<mdot_f,mdot_f;
      af.setZero();
      bc.setZero();
    }

    else if(ubound_type==3)
    {
      ac.setZero();
      af.setZero();
      bc.setZero();
    }


  }

  eq_coeff conv;
  conv.ac = ac;
  conv.af = af;
  conv.bc = bc;

  return conv;

}

eq_coeff Momentum_eq::pressure_gradient(int cell_id)
{
  Matrix<double,1,2> grad_p = P.grad_c(cell_id);
  Matrix<double,1,2> ac,af,bc;
  ac.setZero();
  af.setZero();
  bc = -grad_p*mesh.volume(cell_id);

  eq_coeff pres_grad;
  pres_grad.ac = ac;
  pres_grad.af = af;
  pres_grad.bc = bc;

  return pres_grad;

}

Matrices Momentum_eq::matrix_assemble()
{
  int n_cells = mesh.n_cells;
  int nd = mesh.nd;
  MatrixXd Bc(n_cells,2);Bc.setZero();
  typedef Eigen::Triplet<double> T;
  vector<T> triplet_U,triplet_V;
  triplet_U.reserve(n_cells*2*nd);
  triplet_V.reserve(n_cells*2*nd);

  for(int i=0;i<n_cells;i++)
  {
    RowVectorXi face_id = mesh.cell_face_id.row(i);
    for(int j=0;j<nd;j++)
    {
      int fid = face_id(j);
      int cid = i;
      eq_coeff conv = convection_matrix(fid,cid);
      eq_coeff diff = diffusion_matrix(fid,cid);

      int nid = mesh.neighb(fid,cid);
      if(nid==-1)
      nid=cid;

      Matrix<double,1,2> ac = conv.ac+diff.ac;
      Matrix<double,1,2> af = conv.af+diff.af;
      Matrix<double,1,2> bc = conv.bc+diff.bc;

      //Under relaxing
      double alpha_u = SIMPLE.alpha_u;
      Matrix<double,1,2> vel_c = U.get_c(cid);
      bc = bc + (1-alpha_u)/alpha_u *vel_c.cwiseProduct(ac);
      ac = ac/alpha_u;

      triplet_U.push_back(T(cid,cid,ac(0)));
      triplet_U.push_back(T(cid,nid,af(0)));

      triplet_V.push_back(T(cid,cid,ac(1)));
      triplet_V.push_back(T(cid,nid,af(1)));

      Bc.row(cid) = Bc.row(cid) + bc;

    }
    eq_coeff pgrad = pressure_gradient(i);
    Bc.row(i) = Bc.row(i) + pgrad.bc;

  }

  SparseMatrix<double> U_mat(n_cells,n_cells);
  SparseMatrix<double> V_mat(n_cells,n_cells);
  U_mat.setFromTriplets(triplet_U.begin(), triplet_U.end());
  V_mat.setFromTriplets(triplet_V.begin(), triplet_V.end());

  Matrices matrices;
  matrices.U_mat = U_mat;
  matrices.V_mat = V_mat;
  matrices.Bc = Bc;



  return matrices;

}

Velocity Momentum_eq::matrix_solve()
{//clock_t t;t = clock();
  Matrices mat = matrix_assemble();
  MatrixXd vel = U.get_vel();
  /////////////////////////////// X_momentim equation ///////////////////////////////
  Matrix_solvers x_momentum(mat.U_mat,mat.Bc.col(0),vel.col(0),schemes.u_matrix_solver);
  solution_info x_momentum_results = x_momentum.solve();

  VectorXd Ux = x_momentum_results.x;
  int iter = x_momentum_results.iterations;
  double res = x_momentum_results.residual;

  cout<<"X momentum equation solved. Residual: "<<res<<" Iterations: "<<iter<<endl;

  if(isnan(res))
  {
    cout<<"Divergence detected. Exiting ........"<<endl;
    exit(0);
  }


  ///////////////////////////////// Y_momentum equation ///////////////////////////

  Matrix_solvers y_momentum(mat.V_mat,mat.Bc.col(1),vel.col(1),schemes.v_matrix_solver);
  solution_info y_momentum_results = y_momentum.solve();

  VectorXd Uy = y_momentum_results.x;
  iter = y_momentum_results.iterations;
  res = y_momentum_results.residual;

  cout<<"Y momentum equation solved. Residual: "<<res<<" Iterations: "<<iter<<endl;

  if(isnan(res))
  {
    cout<<"Divergence detected. Exiting ........"<<endl;
    exit(0);
  }

  vel<<Ux,Uy;
  MatrixXd Ap(mesh.n_cells,2);
  Ap<<mat.U_mat.diagonal(),mat.V_mat.diagonal();

  U.Ap = Ap;
  U.set_c(vel);
//t = clock() - t;cout<<"Time elapsed for momentum solution : "<<((double)t)/CLOCKS_PER_SEC<<" s"<<endl;
  return U;
}
