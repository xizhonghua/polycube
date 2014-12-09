#include "split_for_black_edge_test.h"
#include "../common/transition.h"
#include "../common/vtk.h"
#include "../tetmesh/util.h"

using namespace std;
using namespace zjucad::matrix;
using namespace jtf::tetmesh;

void split_for_black_edge_test::setUp()
{
  tet_.resize(4,8);
  const size_t tet_array[] = {1,0,2,3,
                              1,0,3,4,
                              1,0,4,5,
                              1,0,5,6,
                              1,0,6,2,
                              5,6,0,7,
                              5,6,7,8,
                              5,6,1,8};
  copy(tet_array,tet_array + 4 * 8,tet_.begin());
  node_.resize(3,9);
  const double node_array[] = {0,0,1,
                               0,0,0,
                               -1,0,0,
                               -1,-1,0,
                               1,-1,0,
                               1,0,0,
                               0,1,0,
                               1,1,0,
                               0.5,0.5,-1};
  copy(node_array,node_array + 3 * 9, node_.begin());
  orient_tet(node_,tet_);
  fa.reset(face2tet_adjacent::create(tet_));
  for(size_t t = 0; t < tet_.size(2); ++t)
    ortae.add_tet(tet_(colon(),t),*fa);
  ortae.sort_into_loop(tet_,node_);
  inner_face_jump_type[make_pair(0,1)] = 0;
  inner_face_jump_type[make_pair(1,0)] = type_transition1(trans(type_transition2(0)));

  inner_face_jump_type[make_pair(1,2)] = 3;
  inner_face_jump_type[make_pair(2,1)] = type_transition1(trans(type_transition2(3)));

  ofstream ofs("original_tet.vtk");
  tet2vtk(ofs,&node_[0],node_.size(2),&tet_[0],tet_.size(2));
}

void split_for_black_edge_test::split_for_black_edge_test_with_trivial_model()
{
  matrix<matrixd > frame(tet_.size(2));
  for(size_t t = 0; t < tet_.size(2); ++t) frame[t] = eye<double>(3);

  vector<deque<pair<size_t,size_t> > > chain_list;
  deque<pair<size_t,size_t> > temp;
  temp.push_back(make_pair(1,0));
  chain_list.push_back(temp);
  vector<deque<size_t> > singularities_type_;
  deque<size_t> temp_type;
  temp_type.push_back(type_transition1( type_transition2(0) * type_transition2(3)));
  singularities_type_.push_back(temp_type);
  vector<size_t> tet_rot_type;
  relabel_singularity_chain_by_splitting(tet_,node_,frame,ortae,chain_list,singularities_type_,inner_face_jump_type,*fa,pt,tet_rot_type);
  ofstream ofs("split_tet.vtk");
  tet2vtk(ofs,&node_[0],node_.size(2),&tet_[0],tet_.size(2));

  for(one_ring_tet_at_edge::e2tet_type::const_iterator oecit = ortae.e2t_.begin();
      oecit != ortae.e2t_.end(); ++oecit){
    cerr << " edge " << oecit->first.first << " --> " << oecit->first.second << endl;
    const vector<size_t> &vec = oecit->second;
    for(size_t t = 0; t < vec.size(); ++t){
      cerr << vec[t] << endl;
    }
  }

  set<vector<size_t> > face_set;
  vector<size_t> shared_face(3);
  for(boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator mcit = inner_face_jump_type.begin();
      mcit != inner_face_jump_type.end(); ++mcit){
    if(mcit->second == 9)  continue;
    const pair<size_t,size_t> & tet_pair = mcit->first;

    if(get_shared_face(tet_pair,tet_,&shared_face[0])) {
      cerr << "# [error] can not find shared face." << endl;
      return ;
    }
    sort(shared_face.begin(),shared_face.end());
    face_set.insert(shared_face);
  }
  vector<size_t> face_array;
  face_array.reserve(3 * face_set.size());
  for(set<vector<size_t> >::const_iterator scit = face_set.begin();
      scit != face_set.end(); ++scit){
    const vector<size_t> &face = *scit;
    face_array.insert(face_array.end(),face.begin(),face.end());
  }
  ofstream ofs_tri("inner_face.vtk");
  tri2vtk(ofs_tri,&node_[0],node_.size(2),&face_array[0],face_set.size());
}
