#include "Mesh.hpp"


RowVectorXd Mesh::mid(int face_id)
{
  VectorXi node_id = face_node_id.row(face_id);
  RowVectorXd nodes1 = Nodes.row(node_id(0));
  RowVectorXd nodes2 = Nodes.row(node_id(1));
  RowVectorXd mid = 0.5*(nodes1+nodes2);

  return mid;

}

double Mesh::area(int face_id)
{
  VectorXi node_id = face_node_id.row(face_id);
  RowVectorXd nodes1 = Nodes.row(node_id(0));
  RowVectorXd nodes2 = Nodes.row(node_id(1));
  double area = (nodes1-nodes2).norm();

  return area;

}

RowVectorXd Mesh::normal(int face_id,int owner)
{
  VectorXi node_id = face_node_id.row(face_id);
  RowVectorXd nodes1 = Nodes.row(node_id(0));
  RowVectorXd nodes2 = Nodes.row(node_id(1));

  double dx = nodes1(0)-nodes2(0);
  double dy = nodes1(1)-nodes2(1);

  double nx = -dy/sqrt(pow(dx,2)+pow(dy,2));
  double ny = dx/sqrt(pow(dx,2)+pow(dy,2));

  RowVectorXd cent_nodes = centroid.row(owner);
  RowVectorXd mid_nodes = mid(face_id);

  double d1 = mid_nodes(0)-cent_nodes(0);
  double d2 = mid_nodes(1)-cent_nodes(1);
  double dot_d_n = d1*nx + d2*ny;

  if(dot_d_n<0)
  {
    nx = -nx;
    ny = -ny;
  }


  RowVectorXd normals(1,2);
  normals<<nx,ny;

  return normals;
}

int Mesh::neighb(int face_id,int owner)
{
  int cell1 = face_cell_id(face_id,0);
  int cell2 = face_cell_id(face_id,1);

  int neighb=(cell1!=owner)?cell1:cell2;
  return neighb;
}

void Mesh::Read_Mesh()
{
  ifstream fin;
  fin.open(file_mesh);
  cout<<"Reading Mesh......."<<endl;
  if(!fin)
  {
    cout<<"Mesh file not found. Exiting .........\n";
    exit(0);
  }

  int n1,n2,n3,n4;
  string line;

  while(fin)
  {
    getline(fin,line);
    char buffer[line.length()+1];
    strcpy(buffer,line.c_str());

    if(strcmp(buffer,"(0 \"Node Section\")")==0)
    {
      getline(fin,line);
      char buffer[line.length()+1];
      strcpy(buffer,line.c_str());

      if(sscanf(buffer, "(10 (0 %x %x ", &n1, &n2) && strlen(buffer)!=0)
      break;

    }
  }

  int n_nodes = n2-n1+1;
  double x,y;
  MatrixXd Nodes(n_nodes,2);
  int pos=0;
  int cell_type_id;

  while(fin)
  {
    getline(fin,line);
    char buffer[line.length()+1];
    strcpy(buffer,line.c_str());

    if(sscanf(buffer, "%lf %lf", &x, &y) && strlen(buffer)!=0)
    {
      Nodes(pos,0)=x;
      Nodes(pos,1)=y;
      pos++;
    }

    if(strcmp(buffer,"))")==0)
    break;
  }


  while(fin)
  {
    getline(fin,line);
    char buffer[line.length()+1];
    strcpy(buffer,line.c_str());

    if(sscanf(buffer, "(12 (0 %x %x ", &n1, &n2) && strlen(buffer)!=0)
    {
      getline(fin,line);
      char buffer[line.length()+1];
      strcpy(buffer,line.c_str());
      if(sscanf(buffer, "(12 (%*x %*x %*x %*x %x", &cell_type_id) && strlen(buffer)!=0)
      continue;
    }


    if(sscanf(buffer, "(13 (0 %x %x ", &n3, &n4) && strlen(buffer)!=0)
    break;
  }

  int nd=(cell_type_id==3)?4:3;
  int n_cells = n2-n1+1;
  int n_faces = n4-n3+1;
  int id1,id2,id3,id4;
  MatrixXi face_node_id(n_faces,2);
  MatrixXi face_cell_id(n_faces,2);

  char s[20];

  pos=0;
  Boundary bound_dummy,interior;
  vector <Boundary> bound;

  while(fin)
  {
    getline(fin,line);
    char buffer[line.length()+1];
    strcpy(buffer,line.c_str());

    if(sscanf(buffer, "(0 \"Interior faces of zone %[^\")]",s) && strlen(buffer)!=0)
    {
      interior.name = s;
      getline(fin,line);
      char buffer[line.length()+1];
      strcpy(buffer,line.c_str());

      if(sscanf(buffer, "(13 (%*x %x %x ",&n1,&n2) && strlen(buffer)!=0)
      {
        interior.face_id_start = n1-1;
        interior.face_id_end = n2-1;
        interior.nfaces = interior.face_id_end - interior.face_id_start + 1;
      }


    }

    if(sscanf(buffer, "(0 \"Faces of zone %[^\")]",s) && strlen(buffer)!=0)
    {
      bound_dummy.name = s;
      getline(fin,line);
      char buffer[line.length()+1];
      strcpy(buffer,line.c_str());

      if(sscanf(buffer, "(13 (%*x %x %x ",&n1,&n2) && strlen(buffer)!=0)
      {
        bound_dummy.face_id_start = n1-1;
        bound_dummy.face_id_end = n2-1;
        bound_dummy.nfaces = bound_dummy.face_id_end - bound_dummy.face_id_start + 1;
      }

      bound.push_back(bound_dummy);

    }

    if(sscanf(buffer, "%x %x %x %x", &id1, &id2, &id3, &id4) && strlen(buffer)!=0)
    {
      face_node_id(pos,0)=id1-1;
      face_node_id(pos,1)=id2-1;
      face_cell_id(pos,0)=id3-1;
      face_cell_id(pos,1)=id4-1;
      pos++;
    }

  }

  fin.close();



  MatrixXi cell_face_id,pos_array;
  cell_face_id.setZero(n_cells,nd);
  pos_array.setZero(n_cells,1);

  for(int i=0;i<n_faces;i++)
  {
    int cell1 = face_cell_id(i,0);
    cell_face_id(cell1,pos_array(cell1)) = i;
    pos_array(cell1)++;

    int cell2 = face_cell_id(i,1);
    if(cell2!=-1)
    {
      cell_face_id(cell2,pos_array(cell2)) = i;
      pos_array(cell2)++;
    }

  }

  MatrixXi cell_node_id(n_cells,nd);

  for(int i=0;i<n_cells;i++)
  {
    RowVectorXi face_id = cell_face_id.row(i);
    vector<int> v;

    for(int j=0;j<nd;j++)
    {
      v.push_back(face_node_id(face_id(j),0));
      v.push_back(face_node_id(face_id(j),1));
    }
    sort(v.data(),v.data()+v.size());
    vector<int>::iterator it;
    it = unique(v.begin(), v.end());
    v.resize(distance(v.begin(),it));

    Map<RowVectorXi> v2(v.data(), v.size());

    cell_node_id.row(i) = v2;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Centroid and volume calculation

  MatrixXd centroid;
  VectorXd volume;
  centroid.setZero(n_cells,2);
  volume.setZero(n_cells);


  for(int cell_id=0;cell_id<n_cells;cell_id++)
  {
    VectorXi node_id = cell_node_id.row(cell_id);
    MatrixXd nodes(nd,2);
    for(int i=0;i<nd;i++)
    {
      nodes(i,0)=Nodes(node_id(i),0);
      nodes(i,1)=Nodes(node_id(i),1);
    }

    RowVectorXd bary_centre = nodes.colwise().sum()/nd;
    VectorXd xdiff = nodes.col(0).array()-bary_centre(0);
    VectorXd ydiff = nodes.col(1).array()-bary_centre(1);
    VectorXd angle(nd);

    for(int i=0;i<nd;i++)
    angle(i) = atan2(ydiff(i),xdiff(i));

    vector<int> V(nd);
    int x=0;
    iota(V.begin(),V.end(),x++);
    sort( V.begin(),V.end(), [&](int i,int j){return angle(i)<angle(j);} );

    Map<VectorXi> sorted_index(V.data(), V.size());

    MatrixXd nodes_sorted(nd,2);
    RowVectorXi cell_node_id_sorted(nd);

    for(int i=0;i<nd;i++)
    {
      nodes_sorted.row(i) = nodes.row(sorted_index(i));
      cell_node_id_sorted(i) = cell_node_id(cell_id,sorted_index(i));
    }

    cell_node_id.row(cell_id) = cell_node_id_sorted;

    nodes = nodes_sorted;

    double signedArea = 0.0;
    double x0,y0,x1,y1,a;

    // For all vertices except last
    for (int i=0; i<nd; ++i)
    {

      if(i<nd-1)
      {
        x0 = nodes(i,0);    y0 = nodes(i,1);
        x1 = nodes(i+1,0);  y1 = nodes(i+1,1);

      }

      else
      {
        x0 = nodes(i,0);  y0 = nodes(i,1);
        x1 = nodes(0,0);  y1 = nodes(0,1);

      }

      a = x0*y1 - x1*y0;
      signedArea += a;
      centroid(cell_id,0) += (x0 + x1)*a;
      centroid(cell_id,1) += (y0 + y1)*a;

    }

    signedArea *= 0.5;
    centroid(cell_id,0) /= (6.0*signedArea);
    centroid(cell_id,1) /= (6.0*signedArea);
    volume(cell_id) = abs(signedArea);

  }


  Mesh::n_nodes = n_nodes;
  Mesh::n_faces = n_faces;
  Mesh::n_cells = n_cells;
  Mesh::nd = nd;
  Mesh::Nodes = Nodes;
  Mesh::centroid = centroid;
  Mesh::volume = volume;
  Mesh::face_node_id = face_node_id;
  Mesh::face_cell_id = face_cell_id;
  Mesh::cell_face_id = cell_face_id;
  Mesh::cell_node_id = cell_node_id;
  Mesh::Boundaries = bound;
  Mesh::Interior = interior;

  cout<<"Mesh read and calculations done"<<endl;
  cout<<"It has "<<n_nodes<<" nodes, "<<n_faces<<" faces and "<<n_cells<<" cells"<<endl;
  cout<<"Boundaries info :"<<endl;
  for(int i=0;i<Boundaries.size();i++)
  cout<<"\t"<<Boundaries[i].name<<" has "<<Boundaries[i].nfaces<<" faces"<<endl;

  cout<<endl;

}



void Mesh::Read_P_BC()
{

  ifstream fin;
  fin.open(file_bc_p);
  cout<<"Reading Pressure boundary condition......."<<endl;
  if(!fin)
  {
    cout<<"Pressure boundary file not found. Exiting .........\n";
    exit(0);
  }

  string line;
  int pos=0;
  vector<string> bd_names,types;
  vector<double> values;

  while(fin)
  {
    getline(fin,line);
    char buffer[line.length()+1];
    strcpy(buffer,line.c_str());

    char dummy[20];
    double dummy_f;

    if(sscanf(buffer, "ZONE-%s", dummy) && strlen(buffer)!=0)
    {
      bd_names.push_back(dummy);

      while(fin)
      {
        getline(fin,line);
        char buffer[line.length()+1];
        strcpy(buffer,line.c_str());

        if(sscanf(buffer, "type %[^;]", dummy) && strlen(buffer)!=0)
        types.push_back(dummy);
        if(sscanf(buffer, "value %lf;", &dummy_f) && strlen(buffer)!=0)
        values.push_back(dummy_f);

        if(strcmp(buffer,"}")==0)
        break;


      }

    }

  }

  fin.close();

  vector<Boundary> bound = Mesh::Boundaries;

  if(bound.size()!=bd_names.size())
  {
    cout<<"Improper BC specification. Exiting......."<<endl;
    exit(0);
  }

  for(int i=0;i<bound.size();i++)
  {
    for(int j=0;j<bd_names.size();j++)
    {
      if(bound[i].name.compare(bd_names[j])==0)
      {
        bound[i].pbc.type = types[j];
        bound[i].pbc.value = values[j];

      }
    }
  }

  this->pbound_type.setZero(n_faces,1);
  this->pbound_value.setZero(n_faces);

  for(int i=0;i<bound.size();i++)
  {
    if(bound[i].pbc.type.compare("fixed")==0)
    {
      for(int j=bound[i].face_id_start;j<=bound[i].face_id_end;j++)
      {
        this->pbound_type(j) = 1;
        this->pbound_value(j) = bound[i].pbc.value;
      }
    }

    if(bound[i].pbc.type.compare("zeroGradient")==0)
    {
      for(int j=bound[i].face_id_start;j<=bound[i].face_id_end;j++)
      {
        this->pbound_type(j) = 2;
        this->pbound_value(j) = bound[i].pbc.value;
      }
    }

    if(bound[i].pbc.type.compare("wall")==0)
    {
      for(int j=bound[i].face_id_start;j<=bound[i].face_id_end;j++)
      {
        this->pbound_type(j) = 3;
        this->pbound_value(j) = bound[i].pbc.value;
      }
    }

  }

  Mesh::Boundaries = bound;
  cout<<"Done"<<endl<<endl;



}




void Mesh::Read_U_BC()
{

  ifstream fin;
  fin.open(file_bc_u);
  cout<<"Reading Velocity boundary condition......."<<endl;
  if(!fin)
  {
    cout<<"Velocity boundary file not found. Exiting .........\n";
    exit(0);
  }

  string line;
  int pos=0;
  vector<string> bd_names,types;
  vector<double> values_u,values_v;

  while(fin)
  {
    getline(fin,line);
    char buffer[line.length()+1];
    strcpy(buffer,line.c_str());

    char dummy[20];
    double dummy_f1,dummy_f2;

    if(sscanf(buffer, "ZONE-%s", dummy) && strlen(buffer)!=0)
    {
      bd_names.push_back(dummy);

      while(fin)
      {
        getline(fin,line);
        char buffer[line.length()+1];
        strcpy(buffer,line.c_str());

        if(sscanf(buffer, "type %[^;]", dummy) && strlen(buffer)!=0)
        types.push_back(dummy);
        if(sscanf(buffer, "value (%lf %lf);", &dummy_f1, &dummy_f2) && strlen(buffer)!=0)
        {values_u.push_back(dummy_f1);values_v.push_back(dummy_f2);}

        if(strcmp(buffer,"}")==0)
        break;


      }

    }

  }

  fin.close();

  vector<Boundary> bound = Mesh::Boundaries;

  if(bound.size()!=bd_names.size())
  {
    cout<<"Improper BC specification. Exiting......."<<endl;
    exit(0);
  }

  for(int i=0;i<bound.size();i++)
  {
    for(int j=0;j<bd_names.size();j++)
    {
      if(bound[i].name.compare(bd_names[j])==0)
      {
        bound[i].ubc.type = types[j];
        bound[i].ubc.value(0) = values_u[j];
        bound[i].ubc.value(1) = values_v[j];

      }
    }
  }

  this->ubound_type.setZero(n_faces,1);
  this->ubound_value.setZero(n_faces,2);

  for(int i=0;i<bound.size();i++)
  {
    if(bound[i].ubc.type.compare("fixed")==0)
    {
      for(int j=bound[i].face_id_start;j<=bound[i].face_id_end;j++)
      {
        this->ubound_type(j) = 1;
        this->ubound_value.row(j) = bound[i].ubc.value;
      }
    }

    if(bound[i].ubc.type.compare("zeroGradient")==0)
    {
      for(int j=bound[i].face_id_start;j<=bound[i].face_id_end;j++)
      {
        this->ubound_type(j) = 2;
        this->ubound_value.row(j) = bound[i].ubc.value;
      }
    }

    if(bound[i].ubc.type.compare("wall")==0)
    {
      for(int j=bound[i].face_id_start;j<=bound[i].face_id_end;j++)
      {
        this->ubound_type(j) = 3;
        this->ubound_value.row(j) = bound[i].ubc.value;
      }
    }

  }

  Mesh::Boundaries = bound;

  cout<<"Done"<<endl<<endl;

  cout<<"Boundary\tP_bc_type\tP_bc_value\tU_bc_type\tU_bc_value"<<endl;
  for(int i=0;i<Boundaries.size();i++)
  cout<<Boundaries[i].name<<"\t"<<Boundaries[i].pbc.type<<"\t"<<Boundaries[i].pbc.value<<"\t"<<Boundaries[i].ubc.type<<"\t"<<Boundaries[i].ubc.value<<endl;
  cout<<endl;


}
