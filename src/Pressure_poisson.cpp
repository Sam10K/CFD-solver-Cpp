#include "Pressure_poisson.hpp"

Pressure_poisson_eq::Pressure_poisson_eq(Pressure p,Velocity u)
{
  Pdash = p;
  U = u;
  mesh = p.get_mesh();
  solver_setup = p.get_solver_params();
  schemes = p.get_schemes();
  SIMPLE = p.get_simple_params();
}

poisson_eq_coeff Pressure_poisson_eq::poisson_matrix(int face_id,int cell_id)
{
  Matrix<double,1,2> n = mesh.normal(face_id,cell_id);
  double ds = mesh.area(face_id);
  Matrix<double,1,2> S = ds*n;
  Matrix<double,1,2> mid = mesh.mid(face_id);
  int neighb_id = mesh.neighb(face_id,cell_id);
  Matrix<double,1,2> centroid = mesh.centroid.row(cell_id);
  Matrix<double,1,2> neighb = (neighb_id!=-1)?mesh.centroid.row(neighb_id):mid;

  RowVectorXd d_f_n = mid-neighb;
  RowVectorXd d_c_n = centroid-neighb;
  double gc = d_f_n.norm()/d_c_n.norm();
  double gn = 1.0-gc;

  RowVectorXd Ap_c = U.Ap.row(cell_id)/mesh.volume(cell_id);
  RowVectorXd Ap_n = (neighb_id!=-1)?U.Ap.row(neighb_id)/mesh.volume(neighb_id):Ap_c;
  RowVectorXd Ap_f = gc*Ap_c + gn*Ap_n;
  S = S.cwiseQuotient(Ap_f); // D = S.*vol/Ap_f
  Matrix<double,1,2> d = neighb-centroid;
  if(neighb_id==-1)
  d = S.dot(d)*S/S.squaredNorm();

  Matrix<double,1,2> delta = S.squaredNorm()*d/S.dot(d);
  Matrix<double,1,2> k = S-delta;

  RowVectorXd grad_p_f = Pdash.grad_f(face_id,cell_id);

  double ac,af,bc;
  ac = -delta.norm()/d.norm();
  af =  delta.norm()/d.norm();
  bc = -k.dot(grad_p_f);

  if(neighb_id==-1)
  {
    int pbound_type = mesh.pbound_type(face_id);

    if(pbound_type==1)
    {
      af=0.0;
      bc = -delta.norm()*Pdash.get_f(face_id,cell_id)/d.norm();
    }

    else if(pbound_type==2)
    {
      ac=0.0;
      af=0.0;
      bc=0.0;
    }

    else if(pbound_type==3)
    {
      ac=0.0;
      af=0.0;
      bc=0.0;
    }

  }

  poisson_eq_coeff poisson;

  poisson.ac = ac;
  poisson.af = af;
  poisson.bc = bc;

  return poisson;
}

poisson_eq_coeff Pressure_poisson_eq::rhs_matrix(int face_id,int cell_id)
{
  Matrix<double,1,2> n = mesh.normal(face_id,cell_id);
  double ds = mesh.area(face_id);
  Matrix<double,1,2> vel_f = U.get_f(face_id,cell_id);
  double mdot_f_no_rho = vel_f.dot(n)*ds;

  poisson_eq_coeff rhs;

  rhs.ac = 0.0;
  rhs.af = 0.0;
  rhs.bc = mdot_f_no_rho;

  int neighb_id = mesh.neighb(face_id,cell_id);

  if(neighb_id==-1)
  {
    int ubound_type = mesh.ubound_type(face_id);
    RowVectorXd ubound = mesh.ubound_value.row(face_id);

    if(ubound_type==1)
    {
      rhs.ac=0.0;
      rhs.af=0.0;
      mdot_f_no_rho = ubound.dot(n)*ds;
      rhs.bc=mdot_f_no_rho;
    }


    if(ubound_type==3)
    {
      rhs.ac=0.0;
      rhs.af=0.0;
      rhs.bc=0.0;
    }

  }

  return rhs;
}

Matrices_poisson Pressure_poisson_eq::matrix_assemble()
{
  int n_cells = mesh.n_cells;
  int nd = mesh.nd;
  VectorXd Bc(n_cells);Bc.setZero();
  typedef Eigen::Triplet<double> T;
  vector<T> triplet_P;
  triplet_P.reserve(n_cells*2*nd);

  for(int i=0;i<n_cells;i++)
  {
    RowVectorXi face_id = mesh.cell_face_id.row(i);
    for(int j=0;j<nd;j++)
    {
      int fid = face_id(j);
      int cid = i;
      poisson_eq_coeff poisson = poisson_matrix(fid,cid);
      poisson_eq_coeff rhs = rhs_matrix(fid,cid);

      int nid = mesh.neighb(fid,cid);
      if(nid==-1)
      nid=cid;

      double ac = poisson.ac+rhs.ac;
      double af = poisson.af+rhs.af;
      double bc = poisson.bc+rhs.bc;

      triplet_P.push_back(T(cid,cid,ac));
      triplet_P.push_back(T(cid,nid,af));

      Bc(i) = Bc(i) + bc;


    }

  }

  SparseMatrix<double> P_mat(n_cells,n_cells);
  P_mat.setFromTriplets(triplet_P.begin(), triplet_P.end());

  Matrices_poisson matrices;
  matrices.P_mat = P_mat;
  matrices.Bc = Bc;

  return matrices;

}

Pressure Pressure_poisson_eq::matrix_solve()
{//clock_t t;t = clock();
  Matrices_poisson mat = matrix_assemble();
  VectorXd pdash_old = Pdash.get_pressure();
  Matrix_solvers poisson(mat.P_mat,mat.Bc,pdash_old,schemes.p_matrix_solver);
  solution_info poisson_results;

  poisson_results = poisson.solve();

  VectorXd pdash_new = poisson_results.x;
  int iter = poisson_results.iterations;
  double res = poisson_results.residual;

  cout<<"Pressure_poisson solved. Residual: "<<res<<" Iterations: "<<iter<<endl;

  if(isnan(res))
  {
    cout<<"Divergence detected. Exiting ........"<<endl;
    exit(0);
  }

  Pdash.set_c(pdash_new);
//t = clock() - t;cout<<"Time elapsed for pressure poisson solution : "<<((double)t)/CLOCKS_PER_SEC<<" s"<<endl;
  return Pdash;
}
