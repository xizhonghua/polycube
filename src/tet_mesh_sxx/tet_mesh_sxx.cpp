#include <cassert>
#include <fstream>
#include <string>
#include <zjucad/matrix/io.h>
#include "tet_mesh_sxx.h"
#include "../tetmesh/tetmesh.h"
#include "../tetmesh/util.h"
#include "../common/IO.h"
#include "../common/vtk.h"
#include "../common/util.h"

using namespace zjucad::matrix;
using namespace std;

using namespace sxx;

double tet_mesh::compute_area(const face &f)
{
  return fabs(norm(cross(vertex_vec[f[1]] - vertex_vec[f[0]],
      vertex_vec[f[2]] - vertex_vec[f[0]]))) / 2;
}

tet_mesh::tet_mesh()
{
}

int tet_mesh::split_tet(const vector<size_t> &v_id)
{
  assert(v_id.size() == 4);
  vector<size_t> temp_id = v_id;
  sort(temp_id.begin(), temp_id.end());
  tet2it_type::iterator tet2it_it = tet2it.find(temp_id);
  if(tet2it_it == tet2it.end())
    return -1;
  tet_list_it tet_it = tet2it_it->second ;
  return tet_split(tet_it, NULL, NULL, NULL);
  //return 0;
}

int tet_mesh::split_tet(const size_t & p0, const size_t & p1,
                        const size_t & p2, const size_t & p3,
                        const bool debug_info)
{
  vector<size_t> temp_id(4) ;
  temp_id[0] = p0;
  temp_id[1] = p1;
  temp_id[2] = p2;
  temp_id[3] = p3;
  sort(temp_id.begin(), temp_id.end());
  tet2it_type::iterator tet2it_it = tet2it.find(temp_id);
  if(tet2it_it == tet2it.end())
    {
      if(debug_info)
        cerr << "can not find the tet: "
             << temp_id[0] << "," << temp_id[1]
             << temp_id[2] << "," << temp_id[3] << endl;

      return -1;
    }
  tet_list_it tet_it = tet2it_it->second ;
  return tet_split(tet_it, NULL, NULL, NULL);
  //return 0;
}

int tet_mesh::split_tet(const size_t * tet_p, const bool debug_info)
{
  assert(tet_p);
  return split_tet(tet_p[0], tet_p[1], tet_p[2], tet_p[3], debug_info);
}

int tet_mesh::split_face(const size_t id0, const size_t id1, const size_t id2,
                         const bool debug_info)
{
  vector<size_t> temp_face(3);
  temp_face[0] = id0;
  temp_face[1] = id1;
  temp_face[2] = id2;
  sort(temp_face.begin(), temp_face.end());
  face2it_type::iterator face2it_it = face2it.find(temp_face);
  if(face2it_it == face2it.end())
    {
      if(debug_info)
        cerr << "can not find the face: " << temp_face[0]
             << temp_face[1] << temp_face[2] << endl;

      return -1;
    }
  face_list_it face_it = face2it_it->second;
  face_split_at_edge(face_it, NULL, NULL, NULL);
  return 0;
}

int tet_mesh::del_tet(const size_t & p0, const size_t & p1,
                      const size_t & p2, const size_t & p3,
                      const bool debug_info)
{
  vector<size_t> temp_id(4) ;
  temp_id[0] = p0;
  temp_id[1] = p1;
  temp_id[2] = p2;
  temp_id[3] = p3;
  sort(temp_id.begin(), temp_id.end());
  tet2it_type::iterator tet2it_it = tet2it.find(temp_id);
  if(tet2it_it == tet2it.end())
    {
      if(debug_info)
        {
          cerr << "can not find the tet: "
               << temp_id[0] << "," << temp_id[1] << ","
               << temp_id[2] << "," << temp_id[3] << endl;
        }
      return -1;
    }
  tet_list_it tet_it = tet2it_it->second ;
  delete_tet(tet_it, NULL, NULL, NULL);
  cout << "del tet: <" << temp_id[0] << "," << temp_id[1] << ","
       << temp_id[2] << "," << temp_id[3] << ">" << endl;
  return 0;
}

int tet_mesh::split_face(const size_t * face_p, const bool debug_info)
{
  assert(face_p);
  return split_face(face_p[0], face_p[1], face_p[2]);
}

size_t tet_mesh::split_edge(const edge &e, const bool debug_info)
{
  edge temp_e = e;
  if(temp_e.first > temp_e.second)
    swap(temp_e.first, temp_e.second);
  edge2it_type::iterator edge2it_it;
  edge2it_it = edge2it.find(temp_e);
  if(edge2it_it == edge2it.end()){
      if(debug_info)
        cerr << "can not find the edge: <" << temp_e.first << ","
             << temp_e.second << ">" << endl;
      return -1;
    }
  edge_list_it edge_it = edge2it_it->second;
  edge_split(edge_it, NULL, NULL, NULL);
  size_t new_index = vertex_vec.size() - 1;
  return new_index;
}

int tet_mesh::collapse_edge(const edge &e, const bool debug_info)
{
  edge temp_e = e;
  if(temp_e.first > temp_e.second)
    swap(temp_e.first, temp_e.second);
  edge2it_type::iterator edge2it_it;
  edge2it_it = edge2it.find(temp_e);
  if(edge2it_it == edge2it.end())
    {
      if(debug_info)
        cerr << "can not find the edge" << endl;
      return -1;
    }

  face_it_set_type del_face_set;
  edge_list_it edge_it = edge2it_it->second;
  current_face_it = face_list.begin();

  del_face_index = 0;
  return edge_collapse(edge_it, NULL, NULL, &del_face_set, debug_info);
}

int tet_mesh::add_info_ver2tet(const size_t id, const tet_list_it &it)
{
  assert(it->size() == 4);
  vertex2tet_type::iterator v2t_it;
  v2t_it = vertex2tet.find(id);
  if(v2t_it != vertex2tet.end())       // the vertex has appeared in the vertex2tet map
    {
      v2t_it->second.push_back(it);   //add the tet into the list of the vertex2tet
    }
  else                                //the vertex has not appeared in the vertex2tet map
    {
      list<tet_list_it> tet_it_list;
      tet_it_list.push_back(it);
      vertex2tet.insert(make_pair(id, tet_it_list)); //add the pair into the vertx2tet map
    }
  return 0;
}

int tet_mesh::add_info_edge2tet(const size_t id1, const size_t id2, const tet_list_it &it)
{
  assert(id1 != id2);
  edge e(get_min(id1, id2), get_max(id1, id2));
  assert(e.first < e.second);  //make sure the id of the edge vertex is in ascending order
  edge2it_type::iterator edge2it_it;
  edge2it_it = edge2it.find(e);
  if(edge2it_it != edge2it.end()) //the edge has appeared
    {
      edge2tet[edge2it_it->second].push_back(it); //add the tet into the list of the edge2tet
    }
  else
    {
      edge_list.push_back(e);  //add the edge into the edge list
      list<tet_list_it> tet_it_list;
      tet_it_list.push_back(it);
      edge2it.insert(make_pair(e, --edge_list.end()));//add the pair into the edge2it map
      edge2tet.insert(make_pair(--edge_list.end(),
                                tet_it_list));//add the pair into the edge2tet map
    }
  tet2edge[it].push_back(edge2it[e]); //add the edge into the list of the tet2edge
  return 0;
}

int tet_mesh::sort_face_index(const size_t id1, const size_t id2, const size_t id3,
                              face &f)
{
  assert(f.size() == 3); //sort the id of face vertex
  assert((id1 != id2) && (id1 != id3) && (id2 != id3));
  f[2] = get_max(id1, get_max(id2, id3));
  f[0] = get_min(id1, get_min(id2, id3));
  f[1] = id1 + id2 + id3 - f[0] - f[2];
  return 0;
}

int tet_mesh::add_info_face2tet(const size_t id1, const size_t id2, 
                                const size_t id3, const tet_list_it &it)
{
  face f(3);
  sort_face_index(id1, id2, id3, f);
  assert((f[0] < f[1]) && (f[1] < f[2])); //make sure id of the face index in ascending order
  face2it_type::iterator face2it_it;
  face2it_it = face2it.find(f);
  face_list_it new_face_it;
  if(face2it_it != face2it.end()) //the face has appeared
    {
      new_face_it = face2it_it->second;
      assert(new_face_it->size() == 3);
      face2tet[new_face_it].push_back(it); //add the tet into the list of the face2tet
    }
  else
    {
      face_list.push_back(f);   //add the face into the face list
      list<tet_list_it> tet_it_list;
      tet_it_list.push_back(it);
      new_face_it = --face_list.end();
      face2it.insert(make_pair(f, new_face_it)); //add the pair into the face2it map
      face2tet.insert(make_pair(new_face_it,
                                tet_it_list)); //add the pair into the face2tetlist
      if(face2orginal.find(f) == face2orginal.end())
        face2orginal.insert(make_pair(f, f));
    }
  tet2face[it].push_back(face2it[f]);    //add the face into the list of the tet2facemap
  return 0;
}

int tet_mesh::create_tet(const vector<size_t> &tet_veretx)
{
  assert(tet_veretx.size() == 4);
  vector<size_t> tet_veretx_id = tet_veretx;
  sort(tet_veretx_id.begin(), tet_veretx_id.end());
  tet_list.push_back(tet_veretx_id);  //add the tet into the tet list
  tet_list_it tet_it = --tet_list.end();  //get the iterator in list of the current tet
  tet2it.insert(make_pair(tet_veretx_id, tet_it));

  list<edge_list_it> edge_it_list;
  list<face_list_it> face_it_list;
  if(tet2edge.find(tet_it) == tet2edge.end())
    {
      tet2edge.insert(make_pair(tet_it, edge_it_list)); //add the pair into the tet2edge map
    }
  if(tet2face.find(tet_it) == tet2face.end())
    {
      tet2face.insert(make_pair(tet_it, face_it_list)); //add the pair into the tet2face map
    }

  for(size_t i = 0; i < tet_veretx_id.size(); ++i)
    {
      add_info_ver2tet(tet_veretx_id[i], tet_it);   //add the relation between vertex and tet
      for(size_t j = i+1; j < tet_veretx_id.size(); ++j)
        {
          add_info_edge2tet(tet_veretx_id[i],
                            tet_veretx_id[j], tet_it);//add the relation between edge and tet
          for(size_t k = j+1; k < tet_veretx_id.size(); ++k)
            {
              add_info_face2tet(tet_veretx_id[i], tet_veretx_id[j],
                                tet_veretx_id[k], tet_it);//add the relation between face and tet
            }
        }
    }
  return 0;
}

int tet_mesh::create_tet(const set<size_t> &tet_set)
{
  assert(tet_set.size() == 4);
  set<size_t>::iterator set_it;
  vector<size_t> new_tet(4);
  size_t i = 0;
  for(set_it = tet_set.begin(); set_it != tet_set.end(); ++set_it, ++i)
    new_tet[i] = *set_it;
  return create_tet(new_tet);
}

int tet_mesh::create_tetmesh(const char *path)
{
  zjucad::matrix::matrix<double> node;
  zjucad::matrix::matrix<size_t> tet_matrix;
  jtf::mesh::tet_mesh_read_from_zjumat(path, &node, &tet_matrix);
  for(size_t i = 0; i < node.size(2); ++i)
    vertex_vec.push_back(node(zjucad::matrix::colon(), i));
  const size_t tet_size = 4;
  vector<size_t> v(tet_size, 0);

  set<size_t> vertex_in_order;

  for(size_t i = 0; i < tet_matrix.size(2); ++i)
    {
      assert(tet_matrix(zjucad::matrix::colon(), i).size() == tet_size);

      vertex_in_order.clear();

      for(size_t j = 0; j < tet_size; ++j)
        {
          v[j] = tet_matrix(j, i);

          vertex_in_order.insert(tet_matrix(j,i));
        }
      create_tet(v);

      assert(vertex_in_order.size() == tet_size);
      tet2original.insert(make_pair(vertex_in_order, i));

    }
  //  edge_list_it edge_it;
  //  ofstream ofs("all_edges");
  //  for(edge_it = edge_list.begin(); edge_it != edge_list.end(); ++edge_it)
  //    ofs<< edge_it->first << " " << edge_it->second << endl;
  //  cout << "tet num in the list: " << tet_list.size() << endl;
  //  cout << "tet num in the tet2edge: " << tet2edge.size() << endl;
  //  cout << "tet num in the tet2face: " << tet2face.size() << endl ;
  //  cout << "edge num in the edge list:" << edge_list.size() << endl;
  //  cout << "edge num in the edge2it: " << edge2it.size() << endl;
  //  cout << "edge num in the edge2tet: " << edge2tet.size() << endl ;
  //  cout << "face num in the face list:" << face_list.size() << endl;
  //  cout << "face num in the face2it: " << face2it.size() << endl;
  //  cout << "face num in the face2tet: " << face2tet.size() << endl << endl;
  //  cout << "vertex num in the vertex2tet: " << vertex2tet.size() << endl;

  //    face2orginal_type::iterator face2orginal_it;
  //    double sum = 0.0;
  //    for(face2orginal_it = face2orginal.begin(); face2orginal_it != face2orginal.end(); ++face2orginal_it)
  //      {
  //        sum += compute_area(face2orginal_it->first);
  //      }
  //    cout << "no split area: " << sum << endl;

  //   initial_tet_size = tet_matrix.size(2);
  //   tet2original_type::iterator it;
  //   cout << "orginal tet: ";
  // for(it = tet2original.begin(); it != tet2original.end(); ++it)
  //   cout << it->second << " ";
  // cout << endl;

  //   vector<double> initial_volume(initial_tet_size, 0.0);
  //   zjucad::matrix::matrix<size_t> new_tet;
  //   set<size_t>::iterator set_it;
  //   new_tet.resize(4, initial_tet_size);
  //   size_t m=0, n=0;
  //   for(it = tet2original.begin(); it != tet2original.end(); ++it, ++n)
  //     {
  //       m=0;
  //       for(set_it = it->first.begin(); set_it != it->first.end(); ++set_it, ++m)
  //        new_tet(m,n) = *set_it;
  //     }
  ////   for(size_t i =0; i < initial_tet_size; ++i)
  ////     cout << new_tet(0,i) << " "<< new_tet(1,i) << " "<< new_tet(2,i) << " "
  ////         << new_tet(3,i) << endl;
  //   orient_tet(node, new_tet);
  //   m = 0;
  //   for(it = tet2original.begin(); it != tet2original.end(); ++it,++m)
  //     {
  ////       cout << jtf::mesh::cal_tet_vol(node(colon(), new_tet(colon(), m))) << " ";
  //       initial_volume[it->second] += jtf::mesh::cal_tet_vol(node(colon(), new_tet(colon(), m)));
  //     }
  //   for(size_t i = 0; i < initial_tet_size; ++i)
  //     cout << "the volume of the tet " << i << "  " << initial_volume[i] << endl;
  return 0;
}

int tet_mesh::create_tetmesh(const zjucad::matrix::matrix<double> &node,
                             const zjucad::matrix::matrix<size_t> &tet_matrix)
{
  for(size_t i = 0; i < node.size(2); ++i)
    vertex_vec.push_back(node(zjucad::matrix::colon(), i));
  const size_t tet_size = 4;
  vector<size_t> v(tet_size, 0);

  set<size_t> vertex_in_order;

  for(size_t i = 0; i < tet_matrix.size(2); ++i)
    {
      assert(tet_matrix(zjucad::matrix::colon(), i).size() == tet_size);

      vertex_in_order.clear();

      for(size_t j = 0; j < tet_size; ++j)
        {
          v[j] = tet_matrix(j, i);

          vertex_in_order.insert(tet_matrix(j,i));
        }
      create_tet(v);

      assert(vertex_in_order.size() == tet_size);
      tet2original.insert(make_pair(vertex_in_order, i));
    }
  //  face_list_it face_it;
  //  for(face_it = face_list.begin(); face_it != face_list.end(); ++face_it)
  //    face2orginal.insert(make_pair(*face_it, *face_it));
  return 0;
}

int tet_mesh::write_tetmesh_to_matrix(zjucad::matrix::matrix<double> &node,
                                      zjucad::matrix::matrix<size_t> &tet_matrix,
                                      bool remove_extra_node_flag) const
{
  list<tet>::const_iterator tet_it = tet_list.begin();
  node.resize(3, vertex_vec.size());
  assert(tet_list.size() != 0);
  tet_matrix.resize(4, tet_list.size());
  for(size_t i = 0; i < vertex_vec.size(); ++i)
    node(zjucad::matrix::colon(), i) = vertex_vec[i];
  for(size_t i = 0; tet_it != tet_list.end(); ++tet_it, ++i)
    for(size_t j = 0; j < 4; ++j)
      tet_matrix(j, i) = (*tet_it)[j];
  orient_tet(node, tet_matrix);
  if(remove_extra_node_flag)
    remove_extra_node(tet_matrix, node);
  return 0;
}


int tet_mesh::get_tet2orginal_index(const zjucad::matrix::matrix<size_t> &tet_matrix,
                                    zjucad::matrix::matrix<size_t> &tet_index_map) const
{
  assert(tet_matrix.size(1) == 4);
  assert(tet_matrix.size(2) == tet_list.size());
  tet_index_map.resize(1,tet_matrix.size(2));
  set<size_t> vertex_in_order;
  tet2original_type::const_iterator it;

  for(size_t i = 0; i < tet_matrix.size(2); ++i)
    {
      vertex_in_order.clear();
      for(size_t j = 0; j < 4; ++j)
        vertex_in_order.insert(tet_matrix(j, i));
      it = tet2original.find(vertex_in_order);
      assert(it != tet2original.end());
      tet_index_map[i] = it->second;
    }
  return 0;
}

int tet_mesh::write_tetmesh_to_file(const char *path1)
{
  zjucad::matrix::matrix<double> node;
  zjucad::matrix::matrix<size_t> tet_matrix;
  tet_list_it tet_it = tet_list.begin();
  node.resize(3, vertex_vec.size());
  assert(tet_list.size() != 0);
  tet_matrix.resize(4, tet_list.size());
  for(size_t i = 0; i < vertex_vec.size(); ++i)
    node(zjucad::matrix::colon(), i) = vertex_vec[i];
  for(size_t i = 0; tet_it != tet_list.end(); ++tet_it, ++i)
    for(size_t j = 0; j < 4; ++j)
      tet_matrix(j, i) = (*tet_it)[j];
  orient_tet(node, tet_matrix);
  remove_extra_node(tet_matrix, node);
  jtf::mesh::tet_mesh_write_to_zjumat(path1, &node, &tet_matrix);

  // zjucad::matrix::matrix<size_t> del_tet;
  // del_tet.resize(4, tet_vec.size());
  // set<tet>::iterator it;
  // size_t i ;
  // for(it = tet_vec.begin(), i = 0; it != tet_vec.end(); ++it, ++i)
  //   for(size_t j = 0; j < 4; ++j)
  //     del_tet(j, i) = (*it)[j];
  // jtf::mesh::tet_mesh_write_to_zjumat(path2, &node, &del_tet);
  face_list_it face_it;
  face2tet_type::iterator face2tet_it;
  for(face_it = face_list.begin(); face_it != face_list.end(); ++face_it)
    {
      face2tet_it = face2tet.find(face_it);
      assert(face2tet_it != face2tet.end());
      if(face2tet_it->second.size() >2)
        {
          cout << "exist at least 3 tet share a face"
               << (*face_it)[0] << " " << (*face_it)[1] <<" " << (*face_it)[2] << endl;
        }
    }
  return 0;
}

int tet_mesh::del_info_ver2tet(const size_t id, const tet_list_it &it)
{
  vertex2tet_type::iterator ver2tet_it;
  ver2tet_it = vertex2tet.find(id);
  if(ver2tet_it == vertex2tet.end())  //can not find the vertex
    {
      cout << " The vertex was not in the mesh" << endl;
      return 1;
    }
  else
    {
      assert(! ver2tet_it->second.empty());
      ver2tet_it->second.remove(it);   //delete the tet from the list of vertex2tet
      if(ver2tet_it->second.empty())  //the list is empty
        {
          vertex2tet.erase(id);  //delete the pair from the vertex2tet
        }
    }
  return 0;
}

int tet_mesh::del_info_edge2tet(const edge_list_it &edge_it, const tet_list_it &tet_it,
                                edge_it_set_type *del_edge_set)
{
  assert(edge_it != edge_list.end());
  edge2tet_type::iterator edge2tet_it;
  edge2tet_it = edge2tet.find(edge_it);
  edge_list_it edge_it1 = edge_it;
  if(edge2tet_it == edge2tet.end()) //can not find the edge
    {
      cout << "The edge was not in the mesh" << endl;
      return 1;
    }
  else
    {
      assert(! edge2tet_it->second.empty());
      edge2tet_it->second.remove(tet_it); //delete the tet from the list of edge2tet
      if(edge2tet_it->second.empty())  //the list is empty
        {
          if(del_edge_set != NULL)
            {
              edge_it1 = edge_it;
              ++edge_it1;
              int d = signed_distance(current_edge_it, edge_it, edge_list.end());
              del_edge_set->insert(make_pair(d, edge_it1));
            }
          else
            {
              edge e = *edge_it;
              edge2it.erase(e);        //delete the pair from the edge2it
              edge2tet.erase(edge_it); // delete the  pair from the edge2tet
              edge_list.erase(edge_it);      //delete the edge from the edge list
            }
        }
    }
  return 0;
}

int tet_mesh::del_info_face2tet(const face_list_it &face_it, const tet_list_it &tet_it,
                                face_it_set_type *del_face_set)
{
  face f = *face_it;
  assert(face_it != face_list.end());
  face2tet_type::iterator face2tet_it;
  face2tet_it = face2tet.find(face_it);
  face_list_it face_it1 = face_it;
  if(face2tet_it == face2tet.end()) //can not find the face
    {
      cout << "The face was not in the mesh" << endl;
      return 1;
    }
  else
    {
      assert(! face2tet_it->second.empty());
      face2tet_it->second.remove(tet_it); //delete the tet from the list of the face2tet
      if(face2tet_it->second.empty())   //the list is empty
        {
          if(del_face_set != NULL)
            {
              face_it1 = face_it;
              ++face_it1;
              //        int d = signed_distance(current_face_it, face_it, face_list.end());
              //        del_face_set->insert(make_pair(d, face_it1));
              del_face_set->insert(make_pair(int(++del_face_index), face_it1));
            }
          else
            {
              face2it.erase(f);     //delete the pair from  face2it
              face2tet.erase(face_it); //delete the pair from face2tet
              face_list.erase(face_it);   //delete the pair from face list
            }
        }
    }
  return 0;
}

int tet_mesh::delete_tet(const tet_list_it &tet_it,
                         tet_it_set_type *del_tet_set,
                         edge_it_set_type *del_edge_set,
                         face_it_set_type *del_face_set)
{
  assert(tet_it != tet_list.end() && tet_it->size() == 4);
  tet_list_it tet_it1 = tet_it;
  const tet &ver_id = *tet_it;
  for(size_t i = 0; i < ver_id.size(); ++i) //delete the information between tet and vertex
    del_info_ver2tet(ver_id[i], tet_it);
  list<edge_list_it>::iterator edgelist_it;
  edge_list_it edge_it;
  for(edgelist_it = tet2edge[tet_it].begin(); edgelist_it != tet2edge[tet_it].end();
      ++edgelist_it)   //delete the information between tet and edge
    {
      edge_it = *edgelist_it;
      del_info_edge2tet(edge_it, tet_it, del_edge_set);
    }
  list<face_list_it>::iterator facelist_it;
  face_list_it face_it;
  for(facelist_it = tet2face[tet_it].begin(); facelist_it != tet2face[tet_it].end();
      ++facelist_it)   //delete the information between tet and face
    {
      face_it = *facelist_it;
      del_info_face2tet(face_it, tet_it, del_face_set);
    }
  if(del_tet_set != NULL)
    {
      tet_it1 = tet_it;
      ++tet_it1;
      int d = signed_distance(current_tet_it, tet_it, tet_list.end());
      del_tet_set->insert(make_pair(d, tet_it1));
    }
  else
    {
      tet2it.erase(ver_id);
      tet2face.erase(tet_it);   //delete the information of the tet
      tet2edge.erase(tet_it);
      tet_list.erase(tet_it);
    }
  return 0;
}

int tet_mesh::delete_all_tet()
{
  list<tet>::iterator tet_it, tet_it1;
  for(tet_it = tet_list.begin(); tet_it != tet_list.end();)
    {
      tet_it1 = tet_it++;
      //delete_tet(tet_it1);
    }
  cout << "tet num in the list: " << tet_list.size() << endl;
  cout << "tet num in the tet2edge: " << tet2edge.size() << endl;
  cout << "tet num in the tet2face: " << tet2face.size() << endl ;
  cout << "edge num in the edge list:" << edge_list.size() << endl;
  cout << "edge num in the edge2it: " << edge2it.size() << endl;
  cout << "edge num in the edge2tet: " << edge2tet.size() << endl ;
  cout << "face num in the face list:" << face_list.size() << endl;
  cout << "face num in the face2it: " << face2it.size() << endl;
  cout << "face num in the face2tet: " << face2tet.size() << endl << endl;
  return 0;
}

tet_mesh::~tet_mesh()
{
}

int tet_mesh::tet_split(tet_list_it &tet_it,
                        tet_it_set_type *del_tet_set,
                        edge_it_set_type *del_edge_set,
                        face_it_set_type *del_face_set)
{
  const double rate = 0.25;
  const tet cur_tet = *tet_it;
  assert(cur_tet.size() == 4);
  zjucad::matrix::matrix<double> new_point;
  new_point = rate * ((vertex_vec[cur_tet[1]] - vertex_vec[cur_tet[0]])
      + (vertex_vec[cur_tet[2]] - vertex_vec[cur_tet[0]])
      + (vertex_vec[cur_tet[3]] - vertex_vec[cur_tet[0]]));
  new_point += vertex_vec[cur_tet[0]];
  vertex_vec.push_back(new_point);
  const size_t new_id = vertex_vec.size() - 1;

  size_t original_index;
  set<size_t> vertex_in_order;
  for(size_t i = 0; i < cur_tet.size(); ++i)
    vertex_in_order.insert(cur_tet[i]);
  tet2original_type::iterator tet2original_it = tet2original.find(vertex_in_order);
  assert(tet2original_it != tet2original.end());
  original_index = tet2original_it->second;
  tet2original.erase(vertex_in_order);

  delete_tet(tet_it, del_tet_set, del_edge_set, del_face_set);
  delete_info(del_tet_set, del_edge_set, del_face_set);

  vector<size_t> v_id(4);
  v_id[0] = new_id;



  for(size_t i = 0; i < v_id.size(); ++i)
    for(size_t j = i+1; j < v_id.size(); ++j)
      for(size_t k = j+1; k < v_id.size(); ++k)
        {
          v_id[1] = cur_tet[i];
          v_id[2] = cur_tet[j];
          v_id[3] = cur_tet[k];
          create_tet(v_id);

          vertex_in_order.clear();
          for(size_t m = 0; m < v_id.size(); ++m)
            vertex_in_order.insert(v_id[m]);
          assert(tet2original.find(vertex_in_order) == tet2original.end());
          tet2original.insert(make_pair(vertex_in_order, original_index));

        }
  return new_id;
}

//int tet_mesh::tet_split_at_edges(const tet_list_it &tet_it,
//                       const tet_it_set_type  *del_tet_set,
//                       const edge_it_set_type *del_edge_set,
//                       const face_it_set_type *del_face_set)
//{
//  const tet the_split_tet = *tet_it;
//  boost::unordered_map<edge, size_t> edge2mid;
//  boost::unordered_map<edge, size_t>::iterator edge2mid_it;
//  tet2edge_type tet2edge_it = tet2edge.find(tet_it);
//  assert(tet2edge_it != tet2edge.end());
//  list<edge_list_it>::iterator edgelist_it1, edgelist_it2;
//  edge_list_it edge_it;
//  edge e;
//  size_t mid = 0;
//  for(edgelist_it1 = tet2edge_it->second.begin(); edgelist_it1 != tet2edge_it->second.end();)
//    {
//      edge_it = *(edgelist_it1++);
//      e = *edge_it;
//      mid = edge_split(edge_it, NULL, NULL, NULL, NULL, &tet_it);
//      edge2mid.insert(make_pair(e, mid));
//    }
//  vector<size_t> new_tet(4);
//  vector<size_t> f(3);
//  vector<vector<size_t> > new_tet_vec(4);
//  for(size_t i = 0; i < the_split_tet.size(); ++i)
//    {
//      get_face_except_point(the_split_tet, the_split_tet[i], f);
//      for(size_t j = 0; j < 3; ++j)
//        {
//          e = make_pair(get_min(the_split_tet[i], f[j]),
//                        get_max(the_split_tet[i], f[j]));
//          edge2mid_it = edge2mid.find(e);
//          assert(edge2mid_it != edge2mid.end());
//          new_tet[j] = edge2mid_it->second;
//        }
//      new_tet[3] = the_split_tet[i];
//      create_tet(new_tet);

//      new_tet[3] = -1;
//      new_tet_vec.push_back[i] = new_tet;
//    }

//  for()

//  return 0;
//}

int tet_mesh::face_split(face_list_it &face_it,
                         tet_it_set_type *del_tet_set,
                         edge_it_set_type *del_edge_set,
                         face_it_set_type *del_face_set)
{
  const face f = *face_it;
  face2orginal_type::iterator face2orginal_it = face2orginal.find(f);
  assert(face2orginal_it != face2orginal.end());
  zjucad::matrix::matrix<double> new_vertex = zeros(3, 1);
  for(size_t i = 0; i < f.size(); ++i)
    {
      new_vertex += vertex_vec[f[i]];
    }
  new_vertex /= 3;
  vertex_vec.push_back(new_vertex);
  const size_t new_id = vertex_vec.size() - 1;
  face2tet_type::iterator face2tet_it = face2tet.find(face_it);
  assert(face2tet_it != face2tet.end());
  list<tet_list_it>::iterator tetlist_it1, tetlist_it2;
  tet_list_it tet_it;
  vector<size_t> old_vertex;
  vector<size_t> original_index;
  set<size_t> vertex_in_order;
  tet temp_tet;
  tet2original_type::iterator tet2original_it;
  for(tetlist_it1 = face2tet_it->second.begin(); tetlist_it1 != face2tet_it->second.end();)
    {
      tetlist_it2 = tetlist_it1++;
      tet_it = *tetlist_it2;
      vertex_in_order.clear();
      temp_tet = *tet_it;
      for(size_t i = 0; i < temp_tet.size(); ++i)
        vertex_in_order.insert(temp_tet[i]);
      tet2original_it = tet2original.find(vertex_in_order);
      assert(tet2original_it != tet2original.end());
      original_index.push_back(tet2original_it->second);
      tet2original.erase(tet2original_it);
      old_vertex.push_back(get_id_except_face(face_it, tet_it));
      delete_tet(tet_it, del_tet_set, del_edge_set, del_face_set);
    }
  delete_info(del_tet_set, del_edge_set, del_face_set);

  vector<size_t> new_tet(4);
  vector<size_t> temp_face(3);
  vector<size_t> null_vector(3);
  face2orginal_type::iterator face2orginal_it1;
  new_tet[0] = new_id;
  for(size_t i = 0; i < old_vertex.size(); ++i)
    {
      new_tet[1] = old_vertex[i];

      for(size_t j = 0; j < f.size(); ++j)
        {
          new_tet[2] = f[j];
          new_tet[3] = f[(j + 1) % 3];
          create_tet(new_tet);

          vertex_in_order.clear();
          for(size_t m = 0; m < new_tet.size(); ++m)
            vertex_in_order.insert(new_tet[m]);
          assert(tet2original.find(vertex_in_order) == tet2original.end());
          tet2original.insert(make_pair(vertex_in_order, original_index[i]));

          temp_face[0] = new_id;
          temp_face[1] = old_vertex[i];
          temp_face[2] = f[j];
          sort(temp_face.begin(), temp_face.end());
          face2orginal_it1 = face2orginal.find(temp_face);
          assert(face2orginal_it1 != face2orginal.end());
          face2orginal_it1->second = null_vector;
        }

    }
  for(size_t i = 0; i < f.size(); ++i)
    {
      temp_face[0] = new_id;
      temp_face[1] = f[i];
      temp_face[2] = f[(i+1) % 3];
      sort(temp_face.begin(), temp_face.end());
      face2orginal_it1 = face2orginal.find(temp_face);
      assert(face2orginal_it1 != face2orginal.end());
      face2orginal_it1->second = face2orginal_it->second;
    }
  face2orginal.erase(f);
  return 0;
}

int tet_mesh::face_split_at_edge(face_list_it &face_it,
                                 tet_it_set_type *del_tet_set,
                                 edge_it_set_type *del_edge_set,
                                 face_it_set_type *del_face_set)
{
  const face f = *face_it;
  assert(face2tet.find(face_it) != face2tet.end());
  face2orginal_type::iterator face2orginal_it = face2orginal.find(f);
  assert(face2orginal_it != face2orginal.end());
  face2tet_type::iterator face2tet_it = face2tet.find(face_it);
  assert(face2tet_it != face2tet.end());
  list<tet_list_it>::iterator tetlist_it1, tetlist_it2;
  tet_list_it tet_it;
  edge_list_it edge_it;
  vector<size_t> old_vertex;
  vector<size_t> original_index;
  set<size_t> vertex_in_order;
  tet temp_tet;
  tet2original_type::iterator tet2original_it;
  vector<vector<size_t> > face_adj_tet;
  for(tetlist_it1 = face2tet_it->second.begin(); tetlist_it1 != face2tet_it->second.end();
      ++tetlist_it1)
    {
      tet_it = *tetlist_it1;
      face_adj_tet.push_back(*tet_it);
    }

  edge e;
  edge2it_type::iterator edge2it_it;
  vector<size_t> new_ver_vec(3);
  //  vector<face_it_set_type> edge_del_face(f.size());
  for(size_t i = 0; i < f.size(); ++i)
    {
      if(f[i] < f[(i+1) % 3])
        e = make_pair(f[i], f[(i+1) % 3]);
      else
        e = make_pair(f[(i+1) % 3], f[i]);
      edge2it_it = edge2it.find(e);
      assert(edge2it_it != edge2it.end());
      edge_it = edge2it_it->second;
      new_ver_vec[i] = edge_split(edge_it, NULL, NULL, NULL, &face_adj_tet);
    }
  if(del_face_set != NULL)
    current_face_it = face_it;
  for(tetlist_it1 = face2tet_it->second.begin(); tetlist_it1 != face2tet_it->second.end();)
    {
      tetlist_it2 = tetlist_it1++;
      tet_it = *tetlist_it2;
      vertex_in_order.clear();
      temp_tet = *tet_it;
      for(size_t i = 0; i < temp_tet.size(); ++i)
        vertex_in_order.insert(temp_tet[i]);
      tet2original_it = tet2original.find(vertex_in_order);
      assert(tet2original_it != tet2original.end());
      original_index.push_back(tet2original_it->second);
      tet2original.erase(tet2original_it);
      old_vertex.push_back(get_id_except_face(face_it, tet_it));
      delete_tet(tet_it, del_tet_set, del_edge_set, del_face_set);
    }
  delete_info(del_tet_set, del_edge_set, del_face_set);
  //  face_it_set_type::iterator faceset_it;
  //  for(size_t i = 0; i < edge_del_face.size(); ++i)
  //    {
  //      for(faceset_it = edge_del_face[i].begin(); faceset_it != edge_del_face[i].end();
  //          ++faceset_it)
  //        {
  //          del_face_set->insert(*faceset_it);
  //        }
  //    }
  vector<size_t> new_tet(4);
  vector<size_t> temp_face(3);
  vector<size_t> null_vector(3);
  face2orginal_type::iterator face2orginal_it1, face2orginal_it2;
  for(size_t i = 0; i < old_vertex.size(); ++i)
    {
      new_tet[0] = old_vertex[i];
      for(size_t j = 0; j < f.size(); ++j)
        {
          new_tet[1] = f[j];
          new_tet[2] = new_ver_vec[j];
          new_tet[3] = new_ver_vec[(j - 1 + 3) % 3];
          create_tet(new_tet);
          vertex_in_order.clear();
          for(size_t m = 0; m < new_tet.size(); ++m)
            vertex_in_order.insert(new_tet[m]);
          assert(tet2original.find(vertex_in_order) == tet2original.end());
          tet2original.insert(make_pair(vertex_in_order, original_index[i]));
        }

      for(size_t j = 0; j < new_ver_vec.size(); ++j)
        new_tet[j] = new_ver_vec[j];
      new_tet[3] = old_vertex[i];
      create_tet(new_tet);
      vertex_in_order.clear();
      for(size_t m = 0; m < new_tet.size(); ++m)
        vertex_in_order.insert(new_tet[m]);
      assert(tet2original.find(vertex_in_order) == tet2original.end());
      tet2original.insert(make_pair(vertex_in_order, original_index[i]));
    }

  for(size_t i = 0; i < old_vertex.size(); ++i)
    for(size_t j = 0; j < f.size(); ++j)
      {
        temp_face[0] = old_vertex[i];
        temp_face[1] = new_ver_vec[j];
        temp_face[2] = new_ver_vec[(j - 1 + 3) % 3];
        sort(temp_face.begin(), temp_face.end());
        face2orginal_it1 = face2orginal.find(temp_face);
        assert(face2orginal_it1 != face2orginal.end());
        face2orginal_it1->second = null_vector;

        temp_face[0] = old_vertex[i];
        temp_face[1] = f[j];
        temp_face[2] =f[(j + 1) % 3];
        sort(temp_face.begin(), temp_face.end());
        face2orginal_it1 = face2orginal.find(temp_face);
        assert(face2orginal_it1 != face2orginal.end());

        if(face2orginal_it1->second != null_vector)
          {
            temp_face[0] = old_vertex[i];
            temp_face[1] = new_ver_vec[j];
            temp_face[2] = f[j];
            sort(temp_face.begin(), temp_face.end());
            face2orginal_it2 = face2orginal.find(temp_face);
            assert(face2orginal_it2 != face2orginal.end());
            face2orginal_it2->second = face2orginal_it1->second;

            temp_face[0] = old_vertex[i];
            temp_face[1] = new_ver_vec[j];
            temp_face[2] = f[(j + 1) % 3];
            sort(temp_face.begin(), temp_face.end());
            face2orginal_it2 = face2orginal.find(temp_face);
            assert(face2orginal_it2 != face2orginal.end());
            face2orginal_it2->second = face2orginal_it1->second;

            face2orginal_it1->second = null_vector;
          }
      }

  for(size_t i = 0; i < f.size(); ++i)
    {
      temp_face[0] = f[i];
      temp_face[1] = new_ver_vec[i];
      temp_face[2] = new_ver_vec[(i - 1 + 3) % 3];
      sort(temp_face.begin(), temp_face.end());
      face2orginal_it1 = face2orginal.find(temp_face);
      assert(face2orginal_it1 != face2orginal.end());
      face2orginal_it1->second = face2orginal_it->second;
    }

  sort(new_ver_vec.begin(), new_ver_vec.end());
  face2orginal_it1 = face2orginal.find(new_ver_vec);
  assert(face2orginal_it1 != face2orginal.end());
  face2orginal_it1->second = face2orginal_it->second;

  face2orginal_it->second = null_vector;
  return 0;
}

// { edge e;
//  edge2it_type::iterator edge2it_it;
//  vector<size_t> new_ver_vec(3);
////  vector<face_it_set_type> edge_del_face(f.size());
//  for(size_t i = 0; i < f.size(); ++i)
//    {
//      if(f[i] < f[(i+1) % 3])
//        e = make_pair(f[i], f[(i+1) % 3]);
//      else
//        e = make_pair(f[(i+1) % 3], f[i]);
//      edge2it_it = edge2it.find(e);
//      assert(edge2it_it != edge2it.end());
//      edge_it = edge2it_it->second;
//      new_ver_vec[i] = edge_split(edge_it, NULL, NULL, NULL, &face_adj_tet);
//    }
//  if(del_face_set != NULL)
//    current_face_it = face_it;
//  for(tetlist_it1 = face2tet_it->second.begin(); tetlist_it1 != face2tet_it->second.end();)
//    {
//      tetlist_it2 = tetlist_it1++;
//      tet_it = *tetlist_it2;
//      vertex_in_order.clear();
//      temp_tet = *tet_it;
//      for(size_t i = 0; i < temp_tet.size(); ++i)
//        vertex_in_order.insert(temp_tet[i]);
//      tet2original_it = tet2original.find(vertex_in_order);
//      assert(tet2original_it != tet2original.end());
//      original_index.push_back(tet2original_it->second);
//      tet2original.erase(tet2original_it);
//      old_vertex.push_back(get_id_except_face(face_it, tet_it));
//      delete_tet(tet_it, del_tet_set, del_edge_set, del_face_set);
//    }
//  delete_info(del_tet_set, del_edge_set, del_face_set);
////  face_it_set_type::iterator faceset_it;
////  for(size_t i = 0; i < edge_del_face.size(); ++i)
////    {
////      for(faceset_it = edge_del_face[i].begin(); faceset_it != edge_del_face[i].end();
////          ++faceset_it)
////        {
////          del_face_set->insert(*faceset_it);
////        }
////    }
//  vector<size_t> new_tet(4);
//  vector<size_t> temp_face(3);
//  vector<size_t> null_vector(3);
//  face2orginal_type::iterator face2orginal_it1, face2orginal_it2;
//  for(size_t i = 0; i < old_vertex.size(); ++i)
//  {
//    new_tet[1] = old_vertex[i];

//    for(size_t j = 0; j < f.size(); ++j)
//    {
//      new_tet[0] = old_vertex[i];
//      for(size_t j = 0; j < f.size(); ++j)
//        {
//          new_tet[1] = f[j];
//          new_tet[2] = new_ver_vec[j];
//          new_tet[3] = new_ver_vec[(j - 1 + 3) % 3];
//          create_tet(new_tet);
//      temp_face[0] = new_id;
//      temp_face[1] = old_vertex[i];
//      temp_face[2] = f[j];
//      sort(temp_face.begin(), temp_face.end());
//      face2orginal_it1 = face2orginal.find(temp_face);
//      assert(face2orginal_it1 != face2orginal.end());
//      face2orginal_it1->second = null_vector;
//    }
//        }

//      for(size_t j = 0; j < new_ver_vec.size(); ++j)
//        new_tet[j] = new_ver_vec[j];
//      new_tet[3] = old_vertex[i];
//      create_tet(new_tet);
//      vertex_in_order.clear();
//      for(size_t m = 0; m < new_tet.size(); ++m)
//        vertex_in_order.insert(new_tet[m]);
//      assert(tet2original.find(vertex_in_order) == tet2original.end());
//      tet2original.insert(make_pair(vertex_in_order, original_index[i]));
//    }

//  for(size_t i = 0; i < old_vertex.size(); ++i)
//    for(size_t j = 0; j < f.size(); ++j)
//      {
//        temp_face[0] = old_vertex[i];
//        temp_face[1] = new_ver_vec[j];
//        temp_face[2] = new_ver_vec[(j - 1 + 3) % 3];
//        sort(temp_face.begin(), temp_face.end());
//        face2orginal_it1 = face2orginal.find(temp_face);
//        assert(face2orginal_it1 != face2orginal.end());
//        face2orginal_it1->second = null_vector;

//        temp_face[0] = old_vertex[i];
//        temp_face[1] = f[j];
//        temp_face[2] =f[(j + 1) % 3];
//        sort(temp_face.begin(), temp_face.end());
//        face2orginal_it1 = face2orginal.find(temp_face);
//        assert(face2orginal_it1 != face2orginal.end());

//        if(face2orginal_it1->second != null_vector)
//          {
//            temp_face[0] = old_vertex[i];
//            temp_face[1] = new_ver_vec[j];
//            temp_face[2] = f[j];
//            sort(temp_face.begin(), temp_face.end());
//            face2orginal_it2 = face2orginal.find(temp_face);
//            assert(face2orginal_it2 != face2orginal.end());
//            face2orginal_it2->second = face2orginal_it1->second;

//            temp_face[0] = old_vertex[i];
//            temp_face[1] = new_ver_vec[j];
//            temp_face[2] = f[(j + 1) % 3];
//            sort(temp_face.begin(), temp_face.end());
//            face2orginal_it2 = face2orginal.find(temp_face);
//            assert(face2orginal_it2 != face2orginal.end());
//            face2orginal_it2->second = face2orginal_it1->second;

//            face2orginal_it1->second = null_vector;
//          }
//      }

//  for(size_t i = 0; i < f.size(); ++i)
//  {
//    if(f[i] < f[(i+1) % 3])
//      e = make_pair(f[i], f[(i+1) % 3]);
//    else
//      e = make_pair(f[(i+1) % 3], f[i]);
//    edge2it_it = edge2it.find(e);
//    assert(edge2it_it != edge2it.end());
//    edge_it = edge2it_it->second;
//    new_ver_vec[i] = edge_split(edge_it, NULL, NULL, NULL, &face_adj_tet);
//  }
//  if(del_face_set != NULL)
//    current_face_it = face_it;
//  for(tetlist_it1 = face2tet_it->second.begin(); tetlist_it1 != face2tet_it->second.end();)
//  {
//    tetlist_it2 = tetlist_it1++;
//    tet_it = *tetlist_it2;
//    vertex_in_order.clear();
//    temp_tet = *tet_it;
//    for(size_t i = 0; i < temp_tet.size(); ++i)
//      vertex_in_order.insert(temp_tet[i]);
//    tet2original_it = tet2original.find(vertex_in_order);
//    assert(tet2original_it != tet2original.end());
//    original_index.push_back(tet2original_it->second);
//    tet2original.erase(tet2original_it);
//    old_vertex.push_back(get_id_except_face(face_it, tet_it));
//    delete_tet(tet_it, del_tet_set, del_edge_set, del_face_set);
//  }
//  delete_info(del_tet_set, del_edge_set, del_face_set);
//  //  face_it_set_type::iterator faceset_it;
//  //  for(size_t i = 0; i < edge_del_face.size(); ++i)
//  //    {
//  //      for(faceset_it = edge_del_face[i].begin(); faceset_it != edge_del_face[i].end();
//  //          ++faceset_it)
//  //        {
//  //          del_face_set->insert(*faceset_it);
//  //        }
//  //    }
//  vector<size_t> new_tet(4);
//  vector<size_t> temp_face(3);
//  vector<size_t> null_vector(3);
//  face2orginal_type::iterator face2orginal_it1, face2orginal_it2;
//  for(size_t i = 0; i < old_vertex.size(); ++i)
//  {
//    new_tet[0] = old_vertex[i];
//    for(size_t j = 0; j < f.size(); ++j)
//    {
//      temp_face[0] = f[i];
//      temp_face[1] = new_ver_vec[i];
//      temp_face[2] = new_ver_vec[(i - 1 + 3) % 3];
//      sort(temp_face.begin(), temp_face.end());
//      face2orginal_it1 = face2orginal.find(temp_face);
//      assert(face2orginal_it1 != face2orginal.end());
//      face2orginal_it1->second = null_vector;

//      temp_face[0] = old_vertex[i];
//      temp_face[1] = f[j];
//      temp_face[2] =f[(j + 1) % 3];
//      sort(temp_face.begin(), temp_face.end());
//      face2orginal_it1 = face2orginal.find(temp_face);
//      assert(face2orginal_it1 != face2orginal.end());

//      if(face2orginal_it1->second != null_vector)
//      {
//        temp_face[0] = old_vertex[i];
//        temp_face[1] = new_ver_vec[j];
//        temp_face[2] = f[j];
//        sort(temp_face.begin(), temp_face.end());
//        face2orginal_it2 = face2orginal.find(temp_face);
//        assert(face2orginal_it2 != face2orginal.end());
//        face2orginal_it2->second = face2orginal_it1->second;

//        temp_face[0] = old_vertex[i];
//        temp_face[1] = new_ver_vec[j];
//        temp_face[2] = f[(j + 1) % 3];
//        sort(temp_face.begin(), temp_face.end());
//        face2orginal_it2 = face2orginal.find(temp_face);
//        assert(face2orginal_it2 != face2orginal.end());
//        face2orginal_it2->second = face2orginal_it1->second;

//        face2orginal_it1->second = null_vector;
//      }
//    }

//  sort(new_ver_vec.begin(), new_ver_vec.end());
//  face2orginal_it1 = face2orginal.find(new_ver_vec);
//  assert(face2orginal_it1 != face2orginal.end());
//  face2orginal_it1->second = face2orginal_it->second;

//  face2orginal_it->second = null_vector;
//  return 0;
//}

int tet_mesh::edge_split( edge_list_it &edge_it,
                          tet_it_set_type *del_tet_set,
                          edge_it_set_type *del_edge_set,
                          face_it_set_type *del_face_set,
                          vector<vector<size_t> > *face_adj_tet)
{
  edge2tet_type::iterator edge2tet_it;
  edge2tet_it = edge2tet.find(edge_it);
  list<tet_list_it>::iterator tetlist_it1, tetlist_it2;
  tet_list_it tet_it;
  edge e1, e2;
  e1 = *edge_it;
  vector<size_t> v_id(4);
  vector<edge> edge_vec;
  zjucad::matrix::matrix<double> mid = (vertex_vec[e1.first] + vertex_vec[e1.second]) / 2;
  vertex_vec.push_back(mid);
  size_t new_id = vertex_vec.size() - 1;
  v_id[0] = new_id;

  // vector<size_t> f(3);
  // f[0] = e1.first;
  // f[1] = e1.second;
  // f[2] = new_id;
  // face_vec.push_back(f);
  vector<size_t> original_index;
  set<size_t> vertex_in_order;
  tet temp_tet;
  size_t temp_index;

  for(tetlist_it1 = edge2tet_it->second.begin(); tetlist_it1 != edge2tet_it->second.end();)
    {
      tetlist_it2 = tetlist_it1++;
      tet_it = *tetlist_it2;
      if(face_adj_tet != NULL && is_ele_in_vec(*tet_it, *face_adj_tet, temp_index))
        {
          continue;
        }
      vertex_in_order.clear();
      temp_tet = *tet_it;
      for(size_t i = 0; i < temp_tet.size(); ++i)
        vertex_in_order.insert(temp_tet[i]);
      assert(tet2original.find(vertex_in_order) != tet2original.end());
      original_index.push_back(tet2original[vertex_in_order]);
      tet2original.erase(vertex_in_order);
      if(face_adj_tet != NULL && is_ele_in_vec(*tet_it, *face_adj_tet, temp_index))
        {
          continue;
        }

      e2 = get_opposite_edge(edge_it, tet_it);
      edge_vec.push_back(e2);
      delete_tet(tet_it, del_tet_set, del_edge_set, del_face_set);
    }

  delete_info(del_tet_set, del_edge_set, del_face_set);

  delete_info(del_tet_set, del_edge_set, del_face_set);

  set<vector<size_t> > temp_del_face_set;
  vector<size_t> del_face1, del_face2;
  vector<size_t> add_face;
  face2orginal_type::iterator face2orginal_it0, face2orginal_it1, face2orginal_it2;
  set<vector<size_t> >::iterator set_it;
  vector<size_t> null_vector(3);

  for(size_t i = 0; i < edge_vec.size(); ++i)
    {
      e2 = edge_vec[i];

      v_id[1] = e2.first;
      v_id[2] = e2.second;

      v_id[3] = e1.first;
      create_tet(v_id);

      vertex_in_order.clear();
      for(size_t j = 0; j < v_id.size(); ++j)
        vertex_in_order.insert(v_id[j]);
      assert(tet2original.find(vertex_in_order) == tet2original.end());
      tet2original.insert(make_pair(vertex_in_order, original_index[i]));



      v_id[3] = e1.second;
      create_tet(v_id);

      vertex_in_order.clear();
      for(size_t j = 0; j < v_id.size(); ++j)
        vertex_in_order.insert(v_id[j]);
      assert(tet2original.find(vertex_in_order) == tet2original.end());
      tet2original.insert(make_pair(vertex_in_order, original_index[i]));

      del_face1.clear();
      del_face2.clear();

      del_face1.push_back(e1.first);
      del_face1.push_back(e1.second);
      del_face1.push_back(e2.first);
      sort(del_face1.begin(), del_face1.end());
      temp_del_face_set.insert(del_face1);
      face2orginal_it1 = face2orginal.find(del_face1);
      assert(face2orginal_it1 != face2orginal.end());

      del_face2.push_back(e1.first);
      del_face2.push_back(e1.second);
      del_face2.push_back(e2.second);
      sort(del_face2.begin(), del_face2.end());
      temp_del_face_set.insert(del_face2);
      face2orginal_it2 = face2orginal.find(del_face2);
      assert(face2orginal_it2 != face2orginal.end());

      add_face.clear();
      add_face.push_back(v_id[0]);
      add_face.push_back(e1.first);
      add_face.push_back(e2.first);
      sort(add_face.begin(), add_face.end());
      face2orginal_it0 = face2orginal.find(add_face);
      //      if(face2orginal_it0 == face2orginal.end())
      //        face2orginal.insert(make_pair(add_face, face2orginal_it1->second));
      assert(face2orginal_it0 != face2orginal.end());
      face2orginal[add_face] = face2orginal[del_face1];


      add_face.clear();
      add_face.push_back(v_id[0]);
      add_face.push_back(e1.second);
      add_face.push_back(e2.first);
      sort(add_face.begin(), add_face.end());
      face2orginal_it0 = face2orginal.find(add_face);
      //      if(face2orginal_it0 == face2orginal.end())
      //        face2orginal.insert(make_pair(add_face, face2orginal_it1->second));
      assert(face2orginal_it0 != face2orginal.end());
      face2orginal[add_face] = face2orginal[del_face1];

      add_face.clear();
      add_face.push_back(v_id[0]);
      add_face.push_back(e1.first);
      add_face.push_back(e2.second);
      sort(add_face.begin(), add_face.end());
      face2orginal_it0 = face2orginal.find(add_face);
      //      if(face2orginal_it0 == face2orginal.end())
      //        face2orginal.insert(make_pair(add_face, face2orginal_it2->second));
      assert(face2orginal_it0 != face2orginal.end());
      face2orginal[add_face] = face2orginal[del_face2];

      add_face.clear();
      add_face.push_back(v_id[0]);
      add_face.push_back(e1.second);
      add_face.push_back(e2.second);
      sort(add_face.begin(), add_face.end());
      face2orginal_it0 = face2orginal.find(add_face);
      //      if(face2orginal_it0 == face2orginal.end())
      //        face2orginal.insert(make_pair(add_face, face2orginal_it2->second));
      assert(face2orginal_it0 != face2orginal.end());
      face2orginal[add_face] = face2orginal[del_face2];

      del_face1.clear();
      del_face1.push_back(e2.first);
      del_face1.push_back(e2.second);
      del_face1.push_back(v_id[0]);
      sort(del_face1.begin(), del_face1.end());
      temp_del_face_set.insert(del_face1);
      face2orginal_it1 = face2orginal.find(del_face1);
      assert(face2orginal_it1 != face2orginal.end());
      face2orginal[del_face1] = null_vector;

    }

  for(set_it = temp_del_face_set.begin(); set_it != temp_del_face_set.end(); ++set_it)
    {
      face2orginal[*set_it] = null_vector;
    }


  return new_id;
}

int tet_mesh::edge_collapse( edge_list_it &edge_it,
                             tet_it_set_type *del_tet_set,
                             edge_it_set_type *del_edge_set,
                             face_it_set_type *del_face_set,
                             const bool is_debug)
{
  edge2tet_type::iterator edge2tet_it;
  const edge e = *edge_it;
  edge2tet_it = edge2tet.find(edge_it);
  list<tet_list_it>::iterator tetlist_it1, tetlist_it2, tetlist_it3, tetlist_it4;
  tet_list_it tet_it;

  set<size_t> vertex_in_order;
  tet temp_tet;

  edge opp_edge1, opp_edge2;
  set<edge> split_edge_set;
  split_edge_set.clear();
  set<edge>::iterator splitedge_it;
  edge2it_type::iterator edge2it_it;
  edge_list_it edge_it1;
  edge2tet_type::iterator edge2tet_it1;
  set<size_t> v_adj_edge;
  //  for(tetlist_it1 = edge2tet_it->second.begin(); tetlist_it1 != edge2tet_it->second.end();)
  //  {
  //    tetlist_it2 = tetlist_it1++;
  //    tet_it = *tetlist_it2;
  //    edge2it_it = edge2it.find(e);
  //    assert(edge2it_it != edge2it.end());
  //    opp_edge1 = get_opposite_edge(edge2it_it->second, tet_it);
  //    edge2it_it = edge2it.find(opp_edge1);
  //    edge_it1 = edge2it_it->second;
  //    edge2tet_it1 = edge2tet.find(edge_it1);
  //    if(edge2tet_it1->second.size() == 3)
  //    {
  //      v_adj_edge.clear();
  //      size_t j = 0;
  //      for(tetlist_it3 = edge2tet_it1->second.begin(); tetlist_it3 != edge2tet_it1->second.end();)
  //      {
  //        tetlist_it4 = tetlist_it3++;
  //        if(*(*tetlist_it4) != *tet_it)
  //        {
  //          opp_edge2 = get_opposite_edge(edge_it1, *tetlist_it4);
  //          assert(edge2it.find(opp_edge2) != edge2it.end());
  //          if(opp_edge2.first == e.first || opp_edge2.first == e.second)
  //          {
  //            ++j;
  //            v_adj_edge.insert(opp_edge2.second);
  //          }
  //          else if(opp_edge2.second == e.first || opp_edge2.second == e.second)
  //          {
  //            ++j;
  //            v_adj_edge.insert(opp_edge2.first);
  //          }
  //        }
  //      }
  //      if(v_adj_edge.size() == 1 && j ==2)
  //      {
  //        split_edge_set.insert(make_pair(e.first, *(v_adj_edge.begin())));
  //      }
  //    }
  //  }
  //  for(splitedge_it = split_edge_set.begin(); splitedge_it != split_edge_set.end();
  //      ++splitedge_it)
  //    split_edge(*splitedge_it);


  edge2it_it = edge2it.find(e);
  assert(edge2it_it != edge2it.end());
  edge_it = edge2it_it->second;
  edge2tet_it = edge2tet.find(edge_it);

  vector<set<size_t> > edge_adj_tet_vec;
  vector<set<size_t> > ver_adj_tet_vec;
  size_t temp_index;
  vertex2tet_type::iterator vertex2tet_it;


  // test
  {
    // find the tets adjacent with the edges
    for(tetlist_it1 = edge2tet_it->second.begin(); tetlist_it1 != edge2tet_it->second.end();)
      {
        tetlist_it2 = tetlist_it1++;
        tet_it = *tetlist_it2;

        vertex_in_order.clear();
        temp_tet = *tet_it;
        for(size_t i = 0; i < temp_tet.size(); ++i)
          vertex_in_order.insert(temp_tet[i]);
        assert(tet2original.find(vertex_in_order) != tet2original.end());
        edge_adj_tet_vec.push_back(vertex_in_order);
      }

    // find the tets adjacent with the vertex except the edges
    vertex2tet_it = vertex2tet.find(e.first);
    if(vertex2tet_it != vertex2tet.end())
      for(tetlist_it1 = vertex2tet_it->second.begin(); tetlist_it1 != vertex2tet_it->second.end();)
        {
          tetlist_it2 = tetlist_it1++;
          tet_it = *tetlist_it2;
          temp_tet = *tet_it;
          vertex_in_order.clear();
          for(size_t i = 0; i < temp_tet.size(); ++i)
            vertex_in_order.insert(temp_tet[i]);
          if(!is_ele_in_vec(vertex_in_order, edge_adj_tet_vec, temp_index))
            ver_adj_tet_vec.push_back(vertex_in_order);
        }

    vertex2tet_it = vertex2tet.find(e.second);
    if(vertex2tet_it != vertex2tet.end())
      for(tetlist_it1 = vertex2tet_it->second.begin(); tetlist_it1 != vertex2tet_it->second.end();)
        {
          tetlist_it2 = tetlist_it1++;
          tet_it = *tetlist_it2;
          temp_tet = *tet_it;
          vertex_in_order.clear();
          for(size_t i = 0; i < temp_tet.size(); ++i)
            vertex_in_order.insert(temp_tet[i]);
          if(!is_ele_in_vec(vertex_in_order, edge_adj_tet_vec, temp_index))
            {
              for(set<size_t>::iterator set_it = vertex_in_order.begin(); set_it != vertex_in_order.end();
                  ++ set_it)
                {
                  if(*set_it == e.second)
                    {
                      vertex_in_order.erase(set_it);
                      assert(vertex_in_order.size() == 3);
                      vertex_in_order.insert(e.first);
                      break;
                    }
                }
              assert(vertex_in_order.size() == 4);
              ver_adj_tet_vec.push_back(vertex_in_order);
            }
        }

    // test the local mesh is valid
    zjucad::matrix::matrix<size_t> local_tet(4, ver_adj_tet_vec.size());
    for(size_t i = 0; i < ver_adj_tet_vec.size(); ++i)
      {
        size_t j = 0;
        for(set<size_t>::iterator set_it = ver_adj_tet_vec[i].begin();
            set_it != ver_adj_tet_vec[i].end(); ++set_it, ++j)
          {
            local_tet(j, i) = *set_it;
          }
      }
    if(is_debug && !test(local_tet))
      {
        cerr << "# edge<" << e.first << "," << e.second << "> can not collapsed" << endl;
        return -1;
      }

    ver_adj_tet_vec.clear();
    edge_adj_tet_vec.clear();
  }

  for(tetlist_it1 = edge2tet_it->second.begin(); tetlist_it1 != edge2tet_it->second.end();)
    {
      tetlist_it2 = tetlist_it1++;
      tet_it = *tetlist_it2;

      vertex_in_order.clear();
      temp_tet = *tet_it;
      for(size_t i = 0; i < temp_tet.size(); ++i)
        vertex_in_order.insert(temp_tet[i]);
      assert(tet2original.find(vertex_in_order) != tet2original.end());

      edge_adj_tet_vec.push_back(vertex_in_order);
      //    tet2original.erase(vertex_in_order);
      delete_tet(tet_it, del_tet_set, del_edge_set, del_face_set);
    }

  zjucad::matrix::matrix<double> mid = (vertex_vec[e.first] + vertex_vec[e.second]) / 2;
  vertex_vec.push_back(mid);

  const size_t new_id = vertex_vec.size() - 1;
  //  cout << "new_id from edge collapse: " << new_id << endl;
  vector<vector<size_t> > new_tet_vec;
  vector<size_t> new_tet(4);
  vector<size_t> original_index;
  tet2original_type::iterator tet2original_it;

  vertex2tet_it = vertex2tet.find(e.first);
  if(vertex2tet_it != vertex2tet.end())
    for(tetlist_it1 = vertex2tet_it->second.begin(); tetlist_it1 != vertex2tet_it->second.end();)
      {
        tetlist_it2 = tetlist_it1++;
        tet_it = *tetlist_it2;
        temp_tet = *tet_it;
        vertex_in_order.clear();
        size_t j = 0;
        for(size_t i = 0; i < temp_tet.size(); ++i)
          {
            vertex_in_order.insert(temp_tet[i]);
            if(temp_tet[i] != e.first)
              new_tet[j++] = temp_tet[i];
          }
        assert(j == 3);
        new_tet[3] = new_id;
        new_tet_vec.push_back(new_tet);

        tet2original_it = tet2original.find(vertex_in_order);
        assert(tet2original_it != tet2original.end());
        original_index.push_back(tet2original_it->second);

        ver_adj_tet_vec.push_back(vertex_in_order);
        //tet2original.erase(tet2original_it);
        delete_tet(tet_it, del_tet_set, del_edge_set, del_face_set);
      }

  vertex2tet_it = vertex2tet.find(e.second);
  if(vertex2tet_it != vertex2tet.end())
    for(tetlist_it1 = vertex2tet_it->second.begin(); tetlist_it1 != vertex2tet_it->second.end();)
      {
        tetlist_it2 = tetlist_it1++;
        tet_it = *tetlist_it2;
        temp_tet = *tet_it;
        vertex_in_order.clear();
        size_t j = 0;
        for(size_t i = 0; i < temp_tet.size(); ++i)
          {
            vertex_in_order.insert(temp_tet[i]);
            if(temp_tet[i] != e.second)
              new_tet[j++] = temp_tet[i];
          }
        assert(j == 3);
        new_tet[3] = new_id;
        new_tet_vec.push_back(new_tet);

        tet2original_it = tet2original.find(vertex_in_order);
        assert(tet2original_it != tet2original.end());
        original_index.push_back(tet2original_it->second);

        ver_adj_tet_vec.push_back(vertex_in_order);
        //      tet2original.erase(tet2original_it);
        delete_tet(tet_it, del_tet_set, del_edge_set, del_face_set);
      }

  vector<vector<size_t> > del_face_vec;
  face_it_set_type::iterator faceset_it;
  face_list_it face_it;

  //std::cout << "del_face_set.size: " << del_face_set->size() << std::endl;
  for(faceset_it = del_face_set->begin(); faceset_it != del_face_set->end(); ++faceset_it)
    {
      face_it = faceset_it->second;
      --face_it;
      del_face_vec.push_back(*face_it);
    }

  delete_info(del_tet_set, del_edge_set, del_face_set);

  vector<set<size_t> > new_tet_in_order;
  for(size_t i = 0; i < new_tet_vec.size(); ++i)
    {
      create_tet(new_tet_vec[i]);
      vertex_in_order.clear();
      for(size_t j = 0; j < new_tet_vec[i].size(); ++j)
        vertex_in_order.insert((new_tet_vec[i])[j]);
      assert(tet2original.find(vertex_in_order) == tet2original.end());

      new_tet_in_order.push_back(vertex_in_order);
      //    tet2original.insert(make_pair(vertex_in_order, original_index[i]));
    }

  face2orginal_type::iterator face2orginal_it1, face2orginal_it2;
  size_t index[2];
  vector<size_t> temp_face;
  bool flag[2];

  {
    //if(!test())
    if(false)
      {
        cerr << "# back to the before state" << endl;
        vertex2tet_type::iterator ver2tet_it = vertex2tet.find(new_id);
        assert(ver2tet_it != vertex2tet.end());
        list<tet_list_it> tet_it_list = ver2tet_it->second;
        for(tetlist_it1 = tet_it_list.begin(); tetlist_it1 != tet_it_list.end();
            ++tetlist_it1)
          {
            //          tetlist_it2 = tetlist_it1++;
            tet_it = *tetlist_it1;
            assert(tet_it->size() == 4);
            delete_tet(tet_it, NULL, NULL, NULL);
          }

        for(size_t i = 0; i < ver_adj_tet_vec.size(); ++i)
          {
            create_tet(ver_adj_tet_vec[i]);
          }

        for(size_t i = 0; i < edge_adj_tet_vec.size(); ++i)
          create_tet(edge_adj_tet_vec[i]);

        cerr << "# edge<" << e.first << "," << e.second << "> can not collapsed" << endl;
        return -1;
      }
    else
      {
        for(size_t i = 0; i < edge_adj_tet_vec.size(); ++i)
          tet2original.erase(edge_adj_tet_vec[i]);

        for(size_t i = 0; i < ver_adj_tet_vec.size(); ++i)
          tet2original.erase(ver_adj_tet_vec[i]);

        for(size_t i = 0; i < new_tet_in_order.size(); ++i)
          tet2original.insert(make_pair(new_tet_in_order[i], original_index[i]));

        for(size_t i = 0; i < del_face_vec.size(); ++i)
          {
            temp_face = del_face_vec[i];
            if(face2it.find(temp_face) == face2it.end())
              continue;
            face2orginal_it1 = face2orginal.find(temp_face);

            assert(face2orginal_it1 != face2orginal.end());
            flag[0] = is_ele_in_vec(e.first, temp_face, index[0]);
            flag[1] = is_ele_in_vec(e.second, temp_face, index[1]);
            if(flag[0] && !flag[1])
              {
                temp_face[index[0]] = new_id;
                sort(temp_face.begin(), temp_face.end());
                face2orginal_it2 = face2orginal.find(temp_face);
                assert(face2orginal_it2 != face2orginal.end());
                face2orginal_it2->second = face2orginal_it1->second;
              }
            else if(!flag[0] && flag[1])
              {
                temp_face[index[1]] = new_id;
                sort(temp_face.begin(), temp_face.end());
                face2orginal_it2 = face2orginal.find(temp_face);
                assert(face2orginal_it2 != face2orginal.end());
                face2orginal_it2->second = face2orginal_it1->second;
              }
            //          if(flag[0] || flag[1])
            //            face2orginal.erase(face2orginal_it1);
          }
      }
    return 0;
  }
}

int tet_mesh::edge_to_face( edge_list_it &edge_it, vector<size_t> &v_id,
                            tet_it_set_type *del_tet_set,
                            edge_it_set_type *del_edge_set,
                            face_it_set_type *del_face_set)
{
  edge2tet_type::iterator edge2tet_it;
  list<tet_list_it>::iterator tetlist_it1, tetlist_it2;
  tet_list_it tet_it;
  edge2tet_it = edge2tet.find(edge_it);
  edge e = *edge_it;
  for(tetlist_it1 = edge2tet_it->second.begin(); tetlist_it1 != edge2tet_it->second.end();)
    {
      tetlist_it2 = tetlist_it1++;
      tet_it = *tetlist_it2;
      delete_tet(tet_it, del_tet_set, del_edge_set, del_face_set);
    }

  delete_info(del_tet_set, del_edge_set, del_face_set);

  v_id.push_back(e.first);
  assert(v_id.size() == 4);
  create_tet(v_id);

  v_id[3] = e.second;
  create_tet(v_id);
  return 0;
}

int tet_mesh::face_to_edge( face_list_it &face_it, vector<size_t> &v_id,
                            tet_it_set_type *del_tet_set,
                            edge_it_set_type *del_edge_set,
                            face_it_set_type *del_face_set)
{
  face2tet_type::iterator face2tet_it;
  list<tet_list_it>::iterator tetlist_it1, tetlist_it2;
  tet_list_it tet_it;
  face2tet_it = face2tet.find(face_it);
  assert(face2tet_it->second.size() == 2);
  face f = *face_it;
  cout << "face: " << (*face_it)[0] << " " << (*face_it)[1] << " " << (*face_it)[2] << endl;
  for(tetlist_it1 = face2tet_it->second.begin(); tetlist_it1 != face2tet_it->second.end();)
    {
      tetlist_it2 = tetlist_it1++;
      tet_it = *tetlist_it2;
      cout << "tet vertex " << (*tet_it)[0] << " " << (*tet_it)[1] << " "
           << (*tet_it)[2] << " " << (*tet_it)[3] << endl;
      delete_tet(tet_it, del_tet_set, del_edge_set, del_face_set);
    }

  delete_info(del_tet_set, del_edge_set, del_face_set);

  assert(v_id.size() == 2);
  v_id.push_back(0);
  v_id.push_back(0);

  for(size_t k = 0; k < f.size(); ++k)
    for(size_t j = k+1; j < f.size(); ++j)
      {
        v_id[2] = f[k];
        v_id[3] = f[j];
        create_tet(v_id);
      }
  return 0;
}

int tet_mesh::edge_to_edge( edge_list_it &edge_it, const vector<edge> &swap_edge,
                            size_t selected_edge_index,
                            tet_it_set_type *del_tet_set,
                            edge_it_set_type *del_edge_set,
                            face_it_set_type *del_face_set)
{
  cout << "swap_edge size: " << swap_edge.size() << endl;
  assert(swap_edge.size() == 2);
  assert(selected_edge_index ==0 || selected_edge_index ==1);
  edge2tet_type::iterator edge2tet_it;
  list<tet_list_it>::iterator tetlist_it1, tetlist_it2;
  tet_list_it tet_it;
  edge2tet_it = edge2tet.find(edge_it);
  assert(edge2tet_it->second.size() == 4);
  edge e = *edge_it;
  for(tetlist_it1 = edge2tet_it->second.begin(); tetlist_it1 != edge2tet_it->second.end();)
    {
      tetlist_it2 = tetlist_it1++;
      tet_it = *tetlist_it2;
      delete_tet(tet_it, del_tet_set, del_edge_set, del_face_set);
    }

  delete_info(del_tet_set, del_edge_set, del_face_set);

  vector<size_t> v_id(4);
  v_id[0] = swap_edge[selected_edge_index].first;
  v_id[1] = swap_edge[selected_edge_index].second;
  vector<edge> left_edges;
  left_edges.push_back(e);
  left_edges.push_back(swap_edge[1 - selected_edge_index]);
  assert(left_edges.size() == 2);
  for(size_t i = 0; i < left_edges.size(); ++i)
    {
      if(i == 0)
        v_id[2] = left_edges[0].first;
      else
        v_id[2] = left_edges[0].second;
      v_id[3] = left_edges[1].first;
      create_tet(v_id);

      v_id[3] = left_edges[1].second;
      create_tet(v_id);
    }
  return 0;
}

edge tet_mesh::get_opposite_edge(const edge_list_it &edge_it, const tet_list_it &tet_it)
{
  const tet &v_id = *tet_it;
  assert(v_id.size() == 4);
  vector<size_t> v(2);
  for(size_t i = 0, j = 0; i < v_id.size() && j < v.size(); ++i)
    if(edge_it->first != v_id[i] && edge_it->second != v_id[i])
      {
        v[j++] = v_id[i];
      }
  edge e(get_min(v[0], v[1]), get_max(v[0], v[1]));
  return e;
}

int tet_mesh::get_edges_vertex(const vector<edge> &edge_vec, vector<size_t> &v_id)
{
  assert (v_id.size() == 0);
  set<size_t> vec_set;
  set<size_t>::iterator vec_set_it;
  for(size_t i = 0; i < edge_vec.size(); ++i)
    {
      vec_set.insert(edge_vec[i].first);
      vec_set.insert(edge_vec[i].second);
    }
  for(vec_set_it = vec_set.begin(); vec_set_it != vec_set.end(); ++vec_set_it)
    v_id.push_back(*vec_set_it);
  return 0;
}

size_t tet_mesh:: get_id_except_face(const face_list_it &face_it, const tet_list_it &tet_it)
{
  tet t = *tet_it;
  face f = *face_it;
  size_t sum1 = 0, sum2 = 0;
  for(size_t i = 0; i < t.size(); ++i)
    sum1 += t[i];
  for(size_t i = 0; i < f.size(); ++i)
    sum2 += f[i];
  return (sum1 - sum2);
}

int tet_mesh::get_face_except_point(const tet &temp_tet, const size_t id, face &f)
{
  assert(f.size() == 3 && temp_tet.size() == 4);
  size_t j = -1;
  for(size_t i = 0; i < temp_tet.size(); ++i)
    if(temp_tet[i] != id)
      f[++j] = temp_tet[i];
  assert(j == 2);
  return 0;
}

int tet_mesh::test_topology_operation()
{
  // vector<edge> edge_vec;
  // edge e;
  // size_t trash;
  // ifstream ifs(path, ifstream::binary);
  // ifs >> trash;
  // while(!ifs.eof())
  //   {
  //     ifs >> trash >> e.first >> e.second;
  //     assert(e.first < e.second);
  //     edge_vec.push_back(e);
  //   }
  // edge_list_it edge_it1, edge_it2;
  // cout << edge_vec.size() << endl;
  // for(size_t i = 0; i < edge_vec.size()-1; ++i)
  //   {
  //     //cout << edge_vec[i].first << "  " << edge_vec[i].second <<endl;
  //     edge_it1 = edge2it[edge_vec[i]];
  //     edge_split(edge_it1, NULL, NULL, NULL);
  //   }
  // cout << "face_vec size: " << face_vec.size() <<endl;
  // ofstream ofs("mid_point");
  // for(size_t i = 0; i < face_vec.size(); ++i)
  //   ofs << (face_vec[i])[0] << " " << (face_vec[i])[1] << " "
  // 	<< (face_vec[i])[2] << endl;
  edge_list_it edge_it, edge_it2;
  size_t edge_size = edge_list.size(), i = 0;
  cout << "edge_size:" << edge_size <<endl;
  list<tet_list_it>::iterator tet_it_list;
  tet_list_it tet_it;
  // for(edge_it = edge_list.begin(); edge_it != edge_list.end(); ++edge_it)
  //   {
  //     cout <<"**edge<" << edge_it->first << "," << edge_it->second
  // 	   << ">**: " <<edge2tet[edge_it].size() << endl;
  //   }
  tet_it_set_type del_tet_set;
  edge_it_set_type del_edge_set;
  face_it_set_type del_face_set;

  ofstream ofs("all_edges");
  for(edge_it = edge_list.begin(); edge_it != edge_list.end(); ++edge_it)
    ofs << edge_it->first << " " << edge_it->second << endl;

  //  i = 0;
  //  for(edge_it = edge_list.begin(); edge_it != edge_list.end() && i < edge_size; ++i)
  //    {
  //      current_edge_it = edge_it;
  //      edge_it2 = edge_it++;
  //      // cout <<"***edge<" << edge_it2->first << "," << edge_it2->second
  //      // 	   << ">***" << endl;
  //      //if((edge_it2->first == 8 || edge_it2->second == 8) &&
  //      //  (edge_it2->first + edge_it2->second < 16))
  //        {
  //          // for(tet_it_list = edge2tet[edge_it2].begin();
  //          //     tet_it_list != edge2tet[edge_it2].end(); ++tet_it_list)
  //          //   {
  //          //     tet_it = *tet_it_list;
  //          //     tet_vec.insert(*tet_it);
  //          //   }
  //          // edge_split(edge_it2, NULL, NULL, NULL);
  //          edge_split(edge_it2, NULL, &del_edge_set, NULL);
  //          get_next_it(del_edge_set, edge_it);
  //          del_edge_set.clear();
  //          // break;
  //        }
  //    }
  //    i = 0;
  //    for(edge_it = edge_list.begin(); edge_it != edge_list.end() && i< 1000 ; ++i)
  //       {
  //  //       cout <<"***edge<" << edge_it->first << "," << edge_it->second
  //  //           << ">***" << endl;
  //         current_edge_it = edge_it;
  //         current_face_it = face_list.begin();
  //         edge_it2 = edge_it++;
  //         edge_collapse(edge_it2, NULL, &del_edge_set, &del_face_set, true);
  //         get_next_it(del_edge_set, edge_it);
  //         del_edge_set.clear();
  //         del_face_set.clear();
  //         // edge_collapse(edge_it, NULL, NULL, NULL);
  //         // edge_it = edge_list.begin();

  //       }
  //   cout << "face2original size: " << face2orginal.size() << endl;
  // vector<size_t> v_id;
  // size_t edge_3 = 0;
  // for(edge_it = edge_list.begin(); edge_it != edge_list.end();)
  //   {
  //     current_edge_it = edge_it;
  //     edge_it2 = edge_it++;
  //     if(can_edge_to_face(edge_it2, v_id))
  // 	{
  // 	  ++edge_3;
  // 	  edge_to_face(edge_it2, v_id, NULL, &del_edge_set, NULL);
  // 	  get_next_it(del_edge_set, edge_it);
  // 	  del_edge_set.clear();
  // 	}
  //   }
  // cout << "the num edge to face: " << edge_3 << endl;

  // for(tet_it = tet_list.begin(); tet_it != tet_list.end(); ++tet_it)
  //   cout << "tet vertex " << (*tet_it)[0] << " " << (*tet_it)[1] << " "
  // 	 << (*tet_it)[2] << " " << (*tet_it)[3] << endl;

  //  split_face(0, 2, 3);
  //  list<face>::iterator face_it, face_it2;
  //  face f;
  //  vector<size_t> v_id;
  //  size_t face_size = face_list.size();
  //   i = 0;
  //   for(face_it = face_list.begin(); face_it != face_list.end(); ++face_it)
  //     {
  //       f = *face_it;
  //       cout << "face <" << f[0] <<"," << f[1] <<","
  //           << f[2] << "> :" << face2tet[face_it].size() << endl;
  //     }
  //  cout << "face_size: " << face_size << endl;
  //  for(face_it = face_list.begin(); face_it != face_list.end() && i < face_size; ++i)
  //  {
  //    //       current_face_it = face_it;
  //    face_it2 = face_it++;
  //    ++face_it;
  //    face_split_at_edge(face_it2, NULL, NULL, &del_face_set);
  //    //       split_face(f[0], f[1], f[2]);
  //    //       test();
  //    get_next_it(del_face_set, face_it);
  //    del_face_set.clear();
  //    //       break;
  //  }
  //  cout << "i: " << i << endl;
  // // i = 0;
  // for(face_it = face_list.begin(); face_it != face_list.end() ; )
  //   {
  //     current_face_it = face_it;
  //     face_it2 = face_it++;
  //     if(can_face_to_edge(face_it2, v_id))
  // 	{
  // 	  face_to_edge(face_it2, v_id, NULL, NULL, &del_face_set);
  // 	  get_next_it(del_face_set, face_it);
  // 	  del_face_set.clear();
  // 	  // face_to_edge(face_it2, NULL, NULL, NULL);
  // 	  // break;
  // 	}
  //   }

  // vector<edge> swap_edge;
  // for(edge_it = edge_list.begin(); edge_it != edge_list.end();)
  //   {
  //     current_edge_it = edge_it;
  //     edge_it2 = edge_it++;
  //     cout << edge2tet[edge_it2].size() << endl;
  //     if(can_edge_to_edge(edge_it2, swap_edge))
  // 	{
  // 	  edge_to_edge(edge_it2, swap_edge, 0, NULL, &del_edge_set, NULL);
  // 	  get_next_it(del_edge_set, edge_it);
  // 	  del_edge_set.clear();
  // 	}
  //   }

  //  tet_list_it tet_it1, tet_it2;
  //  const size_t tet_size = tet_list.size();
  //  for(tet_it1 = tet_list.begin(); tet_it1 != tet_list.end() && i < tet_size; ++i)
  //    {
  //      current_tet_it = tet_it1;
  //      tet_it2 = tet_it1++;
  //      tet_split(tet_it2, &del_tet_set, NULL, NULL);
  //      get_next_it(del_tet_set, tet_it1);
  //      del_tet_set.clear();
  //    }
  //  vector<double> quality_vec(tet_list.size());
  //  i = 0;
  //  for(tet_list_it tet_it = tet_list.begin(); tet_it != tet_list.end(); ++tet_it, ++i)
  //    quality_vec[i] = compute_quality(*tet_it);
  //  sort(quality_vec.begin(), quality_vec.end());
  //  for(size_t i = 0; i < 400 && i < quality_vec.size(); ++i)
  //    cout << quality_vec[/*quality_vec.size() - 1 - */i] << " ";
  //  cout << endl;

  //  cout << "tet num in the list: " << tet_list.size() << endl;
  //  cout << "tet num in the tet2edge: " << tet2edge.size() << endl;
  //  cout << "tet num in the tet2face: " << tet2face.size() << endl ;
  //  cout << "edge num in the edge list:" << edge_list.size() << endl;
  //  cout << "edge num in the edge2it: " << edge2it.size() << endl;
  //  cout << "edge num in the edge2tet: " << edge2tet.size() << endl ;
  //  cout << "face num in the face list:" << face_list.size() << endl;
  //  cout << "face num in the face2it: " << face2it.size() << endl;
  //  cout << "face num in the face2tet: " << face2tet.size() << endl << endl;
  //  cout << "vertex num in the vertex2tet: " << vertex2tet.size() << endl;
  // cout << "orginal tet: ";

  //   tet2original_type::iterator it;
  //   cout << tet2original.size() << endl;
  //   for(it = tet2original.begin(); it != tet2original.end(); ++it)
  //     cout << it->second << " ";
  //   cout << endl;
  //   vector<double> initial_volume(initial_tet_size, 0.0);
  //   zjucad::matrix::matrix<size_t> new_tet;
  //   set<size_t>::iterator set_it;
  //   new_tet.resize(4, tet2original.size());
  //   size_t m=0, n=0;
  //   for(it = tet2original.begin(); it != tet2original.end(); ++it, ++n)
  //     {
  //       m=0;
  //       for(set_it = it->first.begin(); set_it != it->first.end(); ++set_it, ++m)
  //        new_tet(m,n) = *set_it;
  //     }
  //   zjucad::matrix::matrix<double> new_node;
  //   new_node.resize(3, vertex_vec.size());
  //   for(size_t l = 0; l < vertex_vec.size(); ++l)
  //     new_node(colon(), l) = vertex_vec[l];
  //   orient_tet(new_node, new_tet);
  //   m = 0;
  //   for(it = tet2original.begin(); it != tet2original.end(); ++it, ++m)
  //     {
  //       cout << jtf::mesh::cal_tet_vol(new_node(colon(), new_tet(colon(), m))) << " ";
  //       initial_volume[it->second] += jtf::mesh::cal_tet_vol(new_node(colon(), new_tet(colon(), m)));
  //     }
  //   for(size_t l = 0; l < initial_tet_size; ++l)
  //     cout << "the volume of the tet "<< l << " " << initial_volume[l] << endl;

  //  boost::unordered_map<face, double> map_aera;
  //  boost::unordered_map<face, double>::iterator map_it;
  face2orginal_type::iterator face2orginal_it;
  //  for(face2orginal_it = face2orginal.begin(); face2orginal_it != face2orginal.end(); ++face2orginal_it)
  //    map_aera.insert(make_pair(face2orginal_it->second, 0.0));
  //  for(face2orginal_it = face2orginal.begin(); face2orginal_it != face2orginal.end(); ++face2orginal_it)
  //    {
  //      map_aera[face2orginal_it->second] += compute_area(face2orginal_it->first);
  //    }
  //  for(map_it = map_aera.begin(); map_it != map_aera.end(); ++map_it)

  //  double sum = 0.0;
  //  ofstream ofs("face.txt");
  //  vector<size_t> temp_f;
  //  vector<size_t> null_vector(3);
  //  for(face2orginal_it = face2orginal.begin(); face2orginal_it != face2orginal.end(); ++face2orginal_it)
  //  {
  //    if(face2orginal_it ->second != null_vector)
  //    {
  //      sum += compute_area(face2orginal_it->first);
  //    }
  //  }
  //  cout << "after split area: " << sum << endl;
  //  cout << "face size: " << face_list.size() << endl;
  //  cout << "face2orginal size: " << face2orginal.size() << endl;

  //  sum = 0.0;
  //  face_list_it face_it;
  //  for(face_it = face_list.begin(); face_it != face_list.end(); ++face_it)
  //    sum += compute_area(*face_it);
  //  cout << "after split area: " << sum << endl;

  return 0;
}

//v_id store the vertex id of the face 
bool tet_mesh::can_edge_to_face(const edge_list_it &edge_it, vector<size_t> &v_id)
{
  edge2tet_type::iterator edge2tet_it;
  list<tet_list_it>::iterator tetlist_it1, tetlist_it2;
  tet_list_it tet_it;
  v_id.clear();
  edge2tet_it = edge2tet.find(edge_it);
  if(edge2tet_it->second.size() != 3)
    return false;
  else
    {
      vector<edge> edge_vec;
      for(tetlist_it1 = edge2tet_it->second.begin();
          tetlist_it1 != edge2tet_it->second.end();)
        {
          tetlist_it2 = tetlist_it1++;
          tet_it = *tetlist_it2;
          edge_vec.push_back(get_opposite_edge(edge_it, tet_it));
        }
      get_edges_vertex(edge_vec, v_id);
      if(v_id.size() != 3)
        return false;
      if(intersect_triangle(edge_it->first, edge_it->second, v_id))
        return true;
      else
        return false;
    }
}

bool tet_mesh::can_face_to_edge(const face_list_it &face_it, vector<size_t> &v_id)
{
  face2tet_type::iterator face2tet_it;
  list<tet_list_it>::iterator tetlist_it1, tetlist_it2;
  tet_list_it tet_it;
  face2tet_it = face2tet.find(face_it);
  if(face2tet_it->second.size() != 2)
    return false;
  v_id.clear();
  for(tetlist_it1 = face2tet_it->second.begin(); tetlist_it1 != face2tet_it->second.end();)
    {
      tetlist_it2 = tetlist_it1++;
      tet_it = *tetlist_it2;
      v_id.push_back(get_id_except_face(face_it, tet_it));
    }
  if(!intersect_triangle(v_id[1], v_id[0], *face_it))
    return false;
  else
    return true;
}

bool tet_mesh::can_edge_to_edge(const edge_list_it &edge_it, vector<edge> &swap_edge)
{
  edge2tet_type::iterator edge2tet_it;
  list<tet_list_it>::iterator tetlist_it1, tetlist_it2;
  tet_list_it tet_it;
  vector<size_t> v_id;
  edge2tet_it = edge2tet.find(edge_it);
  swap_edge.clear();
  if(edge2tet_it->second.size() != 4)
    return false;
  else
    {
      cout << "edge: " << edge_it->first <<" " << edge_it->second<< endl;
      vector<edge> edge_vec;
      for(tetlist_it1 = edge2tet_it->second.begin();
          tetlist_it1 != edge2tet_it->second.end();)
        {
          tetlist_it2 = tetlist_it1++;
          tet_it = *tetlist_it2;
          cout << "tet vertex " << (*tet_it)[0] << " " << (*tet_it)[1] << " "
               << (*tet_it)[2] << " " << (*tet_it)[3] << endl;
          edge_vec.push_back(get_opposite_edge(edge_it, tet_it));
        }
      assert(edge_vec.size() == 4);
      get_edges_vertex(edge_vec, v_id);
      if(v_id.size() != 4)
        return false;
      cout << "edge_vec: ";
      for(size_t i = 0; i < edge_vec.size(); ++i)
        cout << "<" << edge_vec[i].first << ", " << edge_vec[i].second << ">" << " ";
      cout << endl;
      cout << "v_id: ";
      for(size_t i = 0 ; i < v_id.size(); ++i)
        cout << v_id[i] << " ";
      cout << endl;
      // else
      {
        edge e;
        bool flag;
        for(size_t i = 0; i < v_id.size(); ++i)
          for(size_t j = i+1; j < v_id.size(); ++j)
            {
              flag = true;
              e = make_pair(get_min(v_id[i], v_id[j]), get_max(v_id[i], v_id[j]));
              for(size_t k = 0; k < edge_vec.size(); ++k)
                {
                  if(edge_vec[k] == e)
                    {
                      flag = false;
                      break;
                    }
                }
              if(flag)
                swap_edge.push_back(e);
            }
        cout << "swap_edge: ";
        for(size_t i = 0; i < swap_edge.size(); ++i)
          cout << "<" << swap_edge[i].first << "," << swap_edge[i].second << ">" << " ";
        cout << endl;
        return true;
      }

    }
}

template <typename T>
int tet_mesh::get_next_it(const std::set<pair<int, T>,cmp<T> > &del_set, T &it)
{
  typename std::set<pair<int, T>,cmp<T> >::iterator s_it1, s_it2;
  if(del_set.size() == 0)
    return 1;
  s_it1 = del_set.begin();
  for(s_it1 = del_set.begin(); s_it1 != del_set.end(); ++s_it1)
    if(s_it1->first == 1)
      break;
  if(s_it1 == del_set.end())
    {
      for(s_it2 = del_set.begin(); s_it2 != del_set.end(); ++s_it2)
        if(s_it2->first == 0)
          it = s_it2->second;
      return 2;
    }
  else
    {
      s_it2 = s_it1;
      ++s_it2;
      for(; s_it2 != del_set.end(); ++s_it2)
        {
          if(s_it2->first - s_it1->first != 1)
            break;
          else
            {
              s_it1 = s_it2;
            }
        }
      it = s_it1->second;
    }
  return 0;
}

int tet_mesh::delete_info(const tet_it_set_type *del_tet_set, 
                          const edge_it_set_type *del_edge_set,
                          const face_it_set_type *del_face_set)
{
  tet_it_set_type::iterator tetset_it;
  edge_it_set_type::iterator edgeset_it;
  face_it_set_type::iterator faceset_it;
  tet_list_it tet_it;
  face_list_it face_it;
  edge_list_it edge_it;
  edge e;
  face f;
  if(del_tet_set != NULL)
    {
      for(tetset_it = del_tet_set->begin(); tetset_it != del_tet_set->end(); ++tetset_it)
        {
          tet_it = tetset_it->second;
          --tet_it;
          tet2it.erase(*tet_it);
          tet2face.erase(tet_it);   //delete the information of the tet
          tet2edge.erase(tet_it);
          tet_list.erase(tet_it);
        }
    }
  if(del_edge_set != NULL)
    {
      for(edgeset_it = del_edge_set->begin(); edgeset_it != del_edge_set->end();
          ++edgeset_it)
        {
          edge_it = edgeset_it->second;
          --edge_it;
          e = *edge_it;
          edge2it.erase(e);        //delete the pair from the edge2it
          edge2tet.erase(edge_it); // delete the  pair from the edge2tet
          edge_list.erase(edge_it);      //delete the edge from the edge list
        }
    }
  if(del_face_set != NULL)
    {
      for(faceset_it = del_face_set->begin(); faceset_it != del_face_set->end();
          ++faceset_it)
        {
          face_it = faceset_it->second;
          --face_it;
          f = *face_it;
          face2it.erase(f);     //delete the pair from  face2it
          face2tet.erase(face_it); //delete the pair from face2tet
          face_list.erase(face_it);   //delete the pair from face list
        }
    }
  return 0;
}

// the algorithm was from the directX SDK
bool tet_mesh::intersect_triangle(const size_t id1, const size_t id2, const face &f)
{
  assert(f.size() == 3);
  zjucad::matrix::matrix<double> D, E1, E2, T, P, Q;
  double det, u ,v, t, inv_det;
  cout << "id1: " << id1 << " " << "id2: " << id2 << " "
       << "face: " << f[0] << " " << f[1] << " " << f[2] << endl;
  D = vertex_vec[id2] - vertex_vec[id1];
  E1 = vertex_vec[f[1]] - vertex_vec[f[0]];
  E2 = vertex_vec[f[2]] - vertex_vec[f[0]];
  P = cross(D, E2);
  det = dot(E1, P);
  if(det > 0)
    {
      T = vertex_vec[id1] - vertex_vec[f[0]];
    }
  else
    {
      T = vertex_vec[f[0]] - vertex_vec[id1];
      det = -det;
    }
  cout << "det: " << det << endl;
  Q = cross(T, E1);
  if(det < 1e-6)
    return false;
  inv_det = 1.0 / det;
  u = dot(T, P);
  u *= inv_det;
  cout << "u: " << u << endl;
  if(!(u > 1e-4 && u < 1)) // make sure not intersect with triangle edges
    return false;
  v = dot(D, Q);
  v *= 1.0 / det;
  cout << "v: " << v << endl;
  if(!(v > 1e-4 && u + v < 1))
    return false;
  t = dot(E2, Q);
  t *= inv_det;
  cout << "t: " << t << endl;
  if(!(t > 1e-4 && t < 1))
    return false;
  return true;
}

int tet_mesh::create_tet_bin(const char *path1, const char *path2)
{
  size_t tet_size, node_size;
  zjucad::matrix::matrix<size_t> tet;
  zjucad::matrix::matrix<double> node;
  ifstream ifs(path1);
  ifs >> node_size >> tet_size;
  tet.resize(4, tet_size);
  node.resize(3, node_size);
  for(size_t i = 0; i < node_size; ++i)
    ifs >> node[3*i] >> node[3*i+1] >> node[3*i+2];
  for(size_t i = 0; i < tet_size; ++i)
    ifs >> tet[4*i] >> tet[4*i+1] >> tet[4*i+2] >> tet[4*i+3];
  orient_tet(node, tet);
  cout << node << endl;
  cout << tet<< endl;
  jtf::mesh::tet_mesh_write_to_zjumat(path2, &node, &tet);
  return 0;
}

int tet_mesh::dump_edge_adj_tet(const char *path, const edge &e)
{
  zjucad::matrix::matrix<size_t> adj_tet;
  adj_tet.resize(4, vertex2tet[e.first].size() +
      vertex2tet[e.second].size());
  list<tet_list_it>::iterator tetlist_it;
  size_t i = 0;
  for(tetlist_it = vertex2tet[e.first].begin();
      tetlist_it != vertex2tet[e.first].end(); ++tetlist_it, ++i)
    {
      for(size_t j = 0; j < 4; ++j)
        adj_tet(j, i) = (*(*tetlist_it))[j];
    }
  for(tetlist_it = vertex2tet[e.second].begin();
      tetlist_it != vertex2tet[e.second].end(); ++tetlist_it, ++i)
    {
      for(size_t j = 0; j < 4; ++j)
        adj_tet(j, i) = (*(*tetlist_it))[j];
    }
  zjucad::matrix::matrix<double> node;
  node.resize(3, vertex_vec.size());
  for(size_t i = 0; i < vertex_vec.size(); ++i)
    node(zjucad::matrix::colon(), i) = vertex_vec[i];
  orient_tet(node, adj_tet);
  jtf::mesh::tet_mesh_write_to_zjumat(path, &node, &adj_tet);
}


int tet_mesh::dfs(const size_t index, vector<bool> &is_visited,
                  const zjucad::matrix::matrix<size_t> &face,
                  const jtf::mesh::edge2cell_adjacent &edge_adj)
{
  size_t face_index, edge_index;
  vector<size_t> edge2tri;
  if(!is_visited[index])
    {
      is_visited[index] = true;
      for(size_t i = 0; i < 3; ++i)
        for(size_t j = i+1; j < 3; ++j)
          {
            edge_index = edge_adj.get_edge_idx(face(i, index), face(j, index));
            edge2tri.clear();
            edge2tri.push_back(edge_adj.edge2cell_[edge_index].first);
            edge2tri.push_back(edge_adj.edge2cell_[edge_index].second);
            if(edge2tri.size() != 2)
              continue;
            face_index = edge2tri[0] + edge2tri[1] - index;
            dfs(face_index, is_visited, face, edge_adj);
          }
    }
  return 0;
}

bool tet_mesh::test()
{
  zjucad::matrix::matrix<double> node;
  zjucad::matrix::matrix<size_t> tet;
  write_tetmesh_to_matrix(node, tet);
  return test(tet);
}

bool tet_mesh::test(const zjucad::matrix::matrix<size_t> &tet_ver)
{
  zjucad::matrix::matrix<size_t> face;
  unique_ptr<jtf::mesh::face2tet_adjacent> face_adj(face_adj->create(tet_ver));
  if(!face_adj.get())
    return false;
  jtf::mesh::get_outside_face(*face_adj, face);

  unique_ptr<jtf::mesh::edge2cell_adjacent> edge_adj(edge_adj->create(face));
  if(!edge_adj.get())
    return false;

  vector<bool> is_visited(face.size(2), false);
  dfs(0, is_visited, face, *edge_adj);
  for(size_t i = 0; i < is_visited.size(); ++i)
    if(!is_visited[i])
      {
        cout << "exist a hole" << endl;
        return false;
      }
  return true;
}

double tet_mesh::compute_quality(const tet &temp_tet)
{
  double ave_length = 0.0, sum_face = 0.0, vol = 0.0;
  const double radius =  0.20412;
  vol = fabs(dot(cross(vertex_vec[temp_tet[1]] - vertex_vec[temp_tet[0]],
      vertex_vec[temp_tet[2]] - vertex_vec[temp_tet[0]]),
      vertex_vec[temp_tet[3]] - vertex_vec[temp_tet[0]])) / 6;
  vector<double> side_length(6);
  tet2it_type::iterator tet2it_it = tet2it.find(temp_tet);
  assert(tet2it_it != tet2it.end());
  tet_list_it tet_it= tet2it_it->second;
  tet2edge_type::iterator tet2edge_it = tet2edge.find(tet_it);
  assert(tet2edge_it != tet2edge.end());
  edge temp_e;
  assert(tet2edge_it->second.size() == 6);
  size_t i = 0;
  double length;
  for(list<edge_list_it>::iterator edgelist_it = tet2edge_it->second.begin();
      edgelist_it != tet2edge_it->second.end(); ++ edgelist_it)
    {
      temp_e = *(*(edgelist_it));
      length = norm(vertex_vec[temp_e.first] - vertex_vec[temp_e.second]);
      side_length[i++] = length;
      ave_length += length;
    }
  ave_length = ave_length / 6;

  tet2face_type::iterator tet2face_it = tet2face.find(tet_it);
  assert(tet2face_it != tet2face.end());
  face temp_face;
  assert(tet2face_it->second.size() == 4);
  for(list<face_list_it>::iterator facelist_it = tet2face_it->second.begin();
      facelist_it != tet2face_it->second.end(); ++facelist_it)
    {
      temp_face = *(*(facelist_it));
      sum_face += compute_area(temp_face);
    }
  assert(sum_face > 1e-6);
  double inscribed_sphere_radius = (3 * vol) / sum_face;
  double first_term = 0.0;
  for(size_t i = 0; i < side_length.size(); ++i)
    first_term += (side_length[i] / ave_length - 1) * (side_length[i] / ave_length - 1);
  first_term /= 2;
  double second_term = ((radius * ave_length) / inscribed_sphere_radius - 1) *
      ((radius * ave_length) / inscribed_sphere_radius - 1);
  return first_term + second_term;
}

int tet_mesh::output_all_edges(const char *path)
{
  std::ofstream ofs(path);
  for(edge2it_type::iterator edge2it_it = edge2it.begin();
      edge2it_it != edge2it.end(); ++edge2it_it)
    ofs << (edge2it_it->first).first << " "
        << (edge2it_it->first).second << std::endl;
  return 0;
}
