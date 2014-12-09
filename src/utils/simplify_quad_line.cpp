#include <vector>
#include <deque>
#include <fstream>
#include <jtflib/mesh/io.h>

#include "../common/visualize_tool.h"
#include "../common/def.h"

#include "../common/util.h"
#include "../quadmesh/util.h"
using namespace std;

int simplify_quad_line(vector<deque<pair<size_t,size_t> > > &feature_line,
                       const size_t point_link_idx = 4)
{
  vector<bool> chains_valid(feature_line.size(), true);
  vector<pair<size_t,size_t> > feature_ends(feature_line.size());
  map<size_t,vector<size_t> > end2chains;
  for(size_t ci = 0; ci < feature_line.size(); ++ci){
    feature_ends[ci] =
        make_pair(feature_line[ci].front().first, feature_line[ci].back().second);
    if(feature_ends[ci].first == feature_ends[ci].second) {
      chains_valid[ci] = false;
      continue;
    }
    end2chains[feature_ends[ci].first].push_back(ci);
    end2chains[feature_ends[ci].second].push_back(ci);
  }

  while(1){
    size_t singularity_num = 0;
    for(map<size_t,vector<size_t> >::iterator cit = end2chains.begin();
        cit != end2chains.end(); ++cit){
      if(cit->second.size() !=  point_link_idx){ // singularity point
        ++singularity_num;
        const vector<size_t> & linked_chains = cit->second;
        for(size_t i = 0; i < linked_chains.size(); ++i)
          chains_valid[linked_chains[i]] = false;
        end2chains.erase(cit);
        break;
      }
    }
    if(singularity_num == 0)
      break;
  }

  vector<deque<pair<size_t,size_t> > > feature_line_temp;
  for(size_t ci = 0; ci < feature_line.size();++ci){
    if(chains_valid[ci]){
      feature_line_temp.push_back(feature_line[ci]);
    }
  }
  feature_line = feature_line_temp;
  return 0;
}

int reassemble_feature_line(vector<deque<pair<size_t,size_t> > > &feature_lines)
{
  vector<pair<size_t,size_t> > edges;
  for(size_t ci = 0; ci < feature_lines.size(); ++ci){
    const deque<pair<size_t,size_t> >& one_chain = feature_lines[ci];
    for(size_t ei = 0; ei < one_chain.size(); ++ei){
      edges.push_back(one_chain[ei]);
    }
  }

  jtf::util::extract_chain_from_edges(edges, feature_lines);
  return 0;
}

int simplify_quad_line(int argc, char * argv[])
{
  if(argc != 4){
    cerr << "# [info] simplify_quad_line quad line output_line" << endl;
    return __LINE__;
  }

  matrixst mesh;
  matrixd node;
  vector<deque<pair<size_t,size_t> > > feature_lines;

  if(jtf::mesh::load_obj(argv[1], mesh, node))
    return __LINE__;

  if(jtf::mesh::load_feature_line(argv[2], feature_lines))
    return __LINE__;

  simplify_quad_line(feature_lines,4);
  reassemble_feature_line(feature_lines);
  //simplify_quad_line(feature_lines,1);

  dump_singularity_to_vtk("new_fl.vtk", node, feature_lines);
  return 0;
}

int extract_simp_quad_line(int argc, char * argv[])
{
  if(argc != 3){
    cerr << "# [info] simplify_quad_line quad output_line" << endl;
    return __LINE__;
  }

  jtf::quadmesh::extract_quadmesh_singularity_line(argv[1], argv[2]);

  return 0;
}
