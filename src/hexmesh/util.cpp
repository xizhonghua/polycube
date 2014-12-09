#include "util.h"
#include <iostream>
#include <numeric>
#include <zjucad/matrix/itr_matrix.h>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <jtflib/mesh/util.h>
#include <jtflib/mesh/mesh.h>
#include "../common/util.h"

#define USE_MESQUITE 1

#if USE_MESQUITE
#include <mesquite/Mesquite_all_headers.hpp>
#endif

//#include <vector>
using namespace std;
using namespace zjucad::matrix;

namespace jtf{
  namespace hexmesh{

    size_t get_edge_idx(const size_t & p0, const size_t & p1,
                        const vector<pair<size_t,size_t> > & edge_vec)
    {
      pair<size_t, size_t> edge(p0,p1) ;
      if(edge.first > edge.second)
        swap(edge.first, edge.second);

      vector<pair<size_t,size_t> >::const_iterator cit =
          find(edge_vec.begin(), edge_vec.end(), edge);
      if(cit == edge_vec.end()) return -1;
      else
        return cit-edge_vec.begin();
    }

    int subdivide_hexmesh(matrixst & hex, matrixd & node)
    {
      unique_ptr<jtf::mesh::face2hex_adjacent> fa(jtf::mesh::face2hex_adjacent::create(hex));
      if(!fa.get()){
          cerr << "# [error] can not build face2hex_adjacent." << endl;
          return __LINE__;
        }

      jtf::mesh::one_ring_hex_at_edge orhae;
      for(size_t hi = 0; hi < hex.size(2); ++hi)
        orhae.add_hex(hex(colon(),hi), node, *fa);

      const size_t new_node_num = node.size(2) + hex.size(2) + fa->faces_.size()
          + orhae.e2h_.size();
      matrixd new_node = zeros<double>(3, new_node_num);

      new_node(colon(), colon(0,node.size(2)-1)) = node;

      // add new nodes
      matrixd average_node;
      for(size_t hi = 0; hi < hex.size(2); ++hi){
          cal_average_node(node(colon(),hex(colon(),hi)), average_node);
          new_node(colon(), hi + node.size(2)) = average_node;
        }

      for(size_t fi = 0; fi < fa->faces_.size(); ++fi){
          const vector<size_t> & one_face = fa->faces_[fi];
          itr_matrix<const size_t*> one_face_mat(4,1, &one_face[0]);
          cal_average_node(node(colon(), one_face_mat), average_node);
          new_node(colon(), node.size(2) + hex.size(2) + fi) = average_node;
        }

      size_t ei = node.size(2) + hex.size(2) + fa->faces_.size();
      vector<pair<size_t,size_t> > edge_vec;
      for(map<pair<size_t,size_t>,vector<size_t> >::const_iterator cit =
          orhae.e2h_.begin(); cit != orhae.e2h_.end(); ++cit, ++ei){
          const pair<size_t,size_t> & one_edge = cit->first;
          if(one_edge.first > one_edge.second)
            edge_vec.push_back(make_pair(one_edge.second, one_edge.first));
          else
            edge_vec.push_back(one_edge);
          new_node(colon(), ei) =
              (node(colon(),one_edge.first) + node(colon(), one_edge.second))/2.0;
        }

      // construct the connections
      matrixst new_hex = zeros<size_t>(8, 8 * hex.size(2));
      for(size_t hi = 0; hi < hex.size(2); ++hi){
          const matrix<size_t> &P = hex(colon(),hi);
          // edge node
          const size_t P01 = get_edge_idx(P[0],P[1], edge_vec) + node.size(2) + hex.size(2) + fa->faces_.size();
          const size_t P02 = get_edge_idx(P[0],P[2], edge_vec) + node.size(2) + hex.size(2) + fa->faces_.size();
          const size_t P04 = get_edge_idx(P[0],P[4], edge_vec) + node.size(2) + hex.size(2) + fa->faces_.size();
          const size_t P13 = get_edge_idx(P[1],P[3], edge_vec) + node.size(2) + hex.size(2) + fa->faces_.size();
          const size_t P15 = get_edge_idx(P[1],P[5], edge_vec) + node.size(2) + hex.size(2) + fa->faces_.size();
          const size_t P26 = get_edge_idx(P[2],P[6], edge_vec) + node.size(2) + hex.size(2) + fa->faces_.size();
          const size_t P23 = get_edge_idx(P[2],P[3], edge_vec) + node.size(2) + hex.size(2) + fa->faces_.size();
          const size_t P37 = get_edge_idx(P[3],P[7], edge_vec) + node.size(2) + hex.size(2) + fa->faces_.size();
          const size_t P45 = get_edge_idx(P[4],P[5], edge_vec) + node.size(2) + hex.size(2) + fa->faces_.size();
          const size_t P46 = get_edge_idx(P[4],P[6], edge_vec) + node.size(2) + hex.size(2) + fa->faces_.size();
          const size_t P57 = get_edge_idx(P[5],P[7], edge_vec) + node.size(2) + hex.size(2) + fa->faces_.size();
          const size_t P67 = get_edge_idx(P[6],P[7], edge_vec) + node.size(2) + hex.size(2) + fa->faces_.size();

          // face node
          const size_t Pf0 = fa->get_face_idx(P[0],P[2],P[3],P[1]) + node.size(2) + hex.size(2);
          const size_t Pf1 = fa->get_face_idx(P[7],P[6],P[4],P[5]) + node.size(2) + hex.size(2);
          const size_t Pf2 = fa->get_face_idx(P[0],P[1],P[5],P[4]) + node.size(2) + hex.size(2);
          const size_t Pf3 = fa->get_face_idx(P[7],P[3],P[2],P[6]) + node.size(2) + hex.size(2);
          const size_t Pf4 = fa->get_face_idx(P[0],P[4],P[6],P[2]) + node.size(2) + hex.size(2);
          const size_t Pf5 = fa->get_face_idx(P[7],P[5],P[1],P[3]) + node.size(2) + hex.size(2);

          const size_t Ph = node.size(2) + hi;

          const size_t hex_[] = {
            P[0], P01, P02, Pf0, P04, Pf2, Pf4, Ph,
            P02, Pf0, P[2], P23, Pf4, Ph, P26, Pf3,
            P04, Pf2, Pf4, Ph, P[4], P45, P46, Pf1,
            Pf4, Ph, P26, Pf3, P46, Pf1, P[6], P67,
            P01, P[1], Pf0, P13, Pf2, P15, Ph, Pf5,
            Pf0, P13, P23, P[3], Ph, Pf5, Pf3, P37,
            Pf2, P15, Ph, Pf5, P45, P[5], Pf1, P57,
            Ph, Pf5, Pf3, P37, Pf1, P57, P67, P[7]
          };
          copy(hex_, hex_ + 64, new_hex.begin() + 64 * hi);
        }

      hex = new_hex;
      node = new_node;
      return 0;
    }

    void orient_hex(matrix<size_t> & hex,
                   const matrix<double> & node)
    {
      matrix<double> edge(3,3);
      for(size_t hi = 0; hi < hex.size(2); ++hi){
          edge(colon(),0) = node(colon(), hex(1,hi)) - node(colon(), hex(0,hi));
          edge(colon(),1) = node(colon(), hex(2,hi)) - node(colon(), hex(1,hi));
          edge(colon(),2) = node(colon(), hex(4,hi)) - node(colon(), hex(0,hi));

          // dot(cross(edge0,edge1), edge2) should be >0
          if(dot(cross(edge(colon(),0), edge(colon(),1)), edge(colon(),2)) < 0){
              for(size_t i = 0; i < 4; ++i)
                swap(hex(i*2+0,hi), hex(i*2+1,hi));
            }
        }
    }

    int hexmesh_quality_improver::init()
    {
      fa_.reset(jtf::mesh::face2hex_adjacent::create(hex_));
      if(!fa_.get()){
          cerr << "# [error] can not build face2hex_adjacent." << endl;
          return __LINE__;
        }

      jtf::mesh::get_outside_face(*fa_, outside_face_);
      jtf::mesh::get_outside_face_idx(*fa_, outside_face_idx_);
      return 0;
    }

    int hexmesh_quality_improver::improve_suface_by_one_ring()
    {
      return 0;
      matrixd face_normal(3, outside_face_.size(2));
      jtf::mesh::cal_face_normal(outside_face_, node_, face_normal);
      //  orient_face_normal_outside_tetmesh(
      //        hex_, node_, outside_face_, outside_face_idx_,
      //        *fa_, face_normal);

      boost::unordered_map<size_t, boost::unordered_set<size_t> > p2f, p2adj_p;
      for(size_t fi = 0; fi < outside_face_.size(2); ++fi){
          for(size_t pi = 0; pi < outside_face_.size(1); ++pi){
              p2f[outside_face_(pi,fi)].insert(fi);
              p2adj_p[outside_face_(pi,fi)].insert(
                    outside_face_((pi+1)%outside_face_.size(1),fi));
              p2adj_p[outside_face_(pi,fi)].insert(
                    outside_face_((pi+2)%outside_face_.size(1),fi));
            }
        }
      assert(p2f.size() == p2adj_p.size());

      matrixd normal = zeros<double>(3,1), average_node = zeros<double>(3,1);
      for(boost::unordered_map<size_t, boost::unordered_set<size_t> >::const_iterator
          cit = p2f.begin(); cit != p2f.end(); ++cit){
          normal = zeros<double>(3,1);
          average_node = zeros<double>(3,1);
          for(boost::unordered_set<size_t>::const_iterator ccit = cit->second.begin();
              ccit != cit->second.end(); ++ccit){
              normal += face_normal(colon(), *ccit);
            }
          const double len = norm(normal);
          if(len > 1e-6)
            normal /= len;
          const boost::unordered_set<size_t> & adj_points = p2adj_p[cit->first];
          if(adj_points.empty()) continue;
          for(boost::unordered_set<size_t>::const_iterator ccit = adj_points.begin();
              ccit != adj_points.end(); ++ccit){
              average_node += node_(colon(), *ccit);
            }
          average_node /= adj_points.size();

          // to calculate the project point on normal plane
          // for a point p with normal n, and another point p'
          // the project point from p' to the normal plane will be P = p' + d*n
          // then (P-p)n=0 => d = (p-p')n

          const double d = dot(node_(colon(),cit->first) - average_node, normal);
          node_(colon(), cit->first) = average_node + d * normal;
        }
      return 0;
    }


    int hexmesh_quality_improver::improve(const string strategy)
    {
#if USE_MESQUITE
      using namespace Mesquite;

      //improve_suface_by_one_ring();
      bool fix_surface = true;
      vector<int> fix_node(node_.size(2), 0);

      //  if(surface_strategy == "one_ring"){
      //    //improve_suface_by_one_ring();
      //    fix_surface = true;
      //  }

      if(fix_surface){
          for(size_t pi = 0; pi < outside_face_.size(); ++pi)
            fix_node[outside_face_[pi]] = 1;
        }

      MsqError err;

      zjucad::matrix::matrix<size_t> msqhex_ = hex_;

      const size_t new_map_indice[8] = {7,5,4,6,3,1,0,2};
      const zjucad::matrix::itr_matrix<const size_t*> node_mapping(8,1,new_map_indice);

      for(size_t hi = 0; hi < hex_.size(2); ++hi){
          msqhex_(colon(),hi) = hex_(node_mapping, hi);
        }

      ArrayMesh mesh(3, node_.size(2), &node_[0], &fix_node[0], msqhex_.size(2),
          HEXAHEDRON, &msqhex_[0]);

      //Mesquite::ShapeImprover optimizer;
      //SizeAdaptShapeWrapper optimizer(1e-3);
      if(strategy == "shape"){
          Mesquite::ShapeImprover optimizer;
          optimizer.run_instructions(&mesh, err);
        }else if(strategy == "lap"){
          LaplaceWrapper optimizer;
          optimizer.run_instructions(&mesh, err);
        }else if(strategy == "size"){
          SizeAdaptShapeWrapper optimizer(1e-3);
          optimizer.run_instructions(&mesh, err);
        }

      if (err) {
          cerr << "# [error] hexmesh quality improve fail " << err << endl;
          return __LINE__;
        }
#else
      cerr << "# [info] No Mesquite !!!" << endl;
#endif
      return 0;
    }

		void extend_hexmesh_based_on_faces(const matrix<size_t> & select_faces,
																			 const matrix<double> & select_normal,
																			 const jtf::mesh::meshes & orig_hexmesh,
																			 jtf::mesh::meshes & new_hexmesh,
																			 const double len_percent)
		{
			  double avg_edge_len = 0.0;

				set<size_t> surface_node(select_faces.begin(), select_faces.end());
				vector<size_t> surface_node_vec(surface_node.begin(), surface_node.end());
				
				map<size_t, set<size_t> > p2f;
				for(size_t fi = 0; fi < select_faces.size(2); ++fi){
					for(size_t pi = 0; pi < select_faces.size(1); ++pi){
						p2f[select_faces(pi,fi)].insert(fi);
						avg_edge_len += norm(orig_hexmesh.node_(colon(), select_faces(pi,fi))
																 - orig_hexmesh.node_(colon(), select_faces((pi+1)%select_faces.size(1))));
					}
				}

				avg_edge_len /= 4 * select_faces.size(2);
				
				matrix<double> new_surface_node_= zeros<double>(3, surface_node_vec.size());
				
				matrix<double> normal = zeros<double>(3,1);
				for(size_t pi = 0; pi < new_surface_node_.size(2); ++pi){
					normal = zeros<double>(3,1);
					const set<size_t> & adj_faces = p2f[surface_node_vec[pi]];
					
					for(set<size_t>::const_iterator cit = adj_faces.begin();
							cit != adj_faces.end();  ++cit){
						normal += select_normal(colon(), *cit);
					}
					normal /= adj_faces.size();
					
					new_surface_node_(colon(),pi) = orig_hexmesh.node_(colon(), surface_node_vec[pi]) +
						len_percent * avg_edge_len * normal;
				}

				map<size_t,size_t> p2new_p;
				new_hexmesh.node_ = zeros<double>(3, orig_hexmesh.node_.size(2) + surface_node_vec.size());
				new_hexmesh.node_(colon(), colon(0, orig_hexmesh.node_.size(2)-1)) = orig_hexmesh.node_;
				new_hexmesh.node_(colon(), colon(orig_hexmesh.node_.size(2), new_hexmesh.node_.size(2)-1)) = new_surface_node_;
				
				for(size_t pi = 0; pi < surface_node_vec.size(); ++pi)
					p2new_p[surface_node_vec[pi]] = pi + orig_hexmesh.node_.size(2);
				
				matrix<size_t> new_hex = zeros<size_t>(8, select_faces.size(2));
				
				for(size_t fi = 0; fi < select_faces.size(2); ++fi){
					new_hex(0,fi) = select_faces(0,fi);
					new_hex(1,fi) = select_faces(1,fi);
					new_hex(2,fi) = select_faces(3,fi);
					new_hex(3,fi) = select_faces(2,fi);
					new_hex(4,fi) = p2new_p[select_faces(0,fi)];
					new_hex(5,fi) = p2new_p[select_faces(1,fi)];
					new_hex(6,fi) = p2new_p[select_faces(3,fi)];
					new_hex(7,fi) = p2new_p[select_faces(2,fi)];
				}
				
				new_hexmesh.mesh_ = zeros<size_t>(8, orig_hexmesh.mesh_.size(2) + select_faces.size(2));
				
				new_hexmesh.mesh_(colon(), colon(0, orig_hexmesh.mesh_.size(2)-1)) = orig_hexmesh.mesh_;
				new_hexmesh.mesh_(colon(), colon(orig_hexmesh.mesh_.size(2), new_hexmesh.mesh_.size(2)-1)) = new_hex;
		}
  }
}

