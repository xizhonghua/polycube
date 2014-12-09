#include <boost/property_tree/ptree.hpp>
#include <zjucad/ptree/ptree.h>
#include <string>
#include <iostream>
using namespace std;
using boost::property_tree::ptree;

int main(int argc, char *argv[])
{
  ptree pt;

  try {
    zjucad::read_cmdline(argc,argv,pt);
    pt.put("prog.desc","<frame, fram_inner, frame_polycube, frame2vtk, cut_tet,cut_tet_inner, int_param, find_singularities,find_singularities2, find_singularities3,int_pts, init_zyz, init_zyz_inner, init_zyz_polycube, label_polycube, polycube_param, remove_black_lines, find_singularities_type_consistency, map_hex_to_tet, cut_hex_to_tet>");
    pt.put("tet.desc","file name to the input tet model");

#define CALL_SUB_PROG(prog)						\
  int prog(ptree &pt);						  \
  if(pt.get<string>("prog.value") == #prog)				\
  return prog(pt);

//    CALL_SUB_PROG(frame_inner);
    CALL_SUB_PROG(frame_inner_with_surface);
    CALL_SUB_PROG(frame_inner_with_type);
    CALL_SUB_PROG(frame_inner_L1);
    CALL_SUB_PROG(cut_tet_inner);
    CALL_SUB_PROG(find_singularities2_inner);
    CALL_SUB_PROG(find_singularities_use_face_type);
    CALL_SUB_PROG(draw_inner_residual);
    CALL_SUB_PROG(draw_sh_smooth);
    CALL_SUB_PROG(init_zyz_inner_with_surface2);
//    CALL_SUB_PROG(init_zyz_inner_with_surface3);
//    CALL_SUB_PROG(init_zyz_inner_with_surface4);
    CALL_SUB_PROG(init_zyz_inner_with_surface_tri);
    CALL_SUB_PROG(init_zyz_inner_with_surface_tri2);
    CALL_SUB_PROG(label_polycube);
    CALL_SUB_PROG(surf2inner_dijkstra);
    CALL_SUB_PROG(refine_frame_field_after_aligned);
    CALL_SUB_PROG(split_tet_at_zigzag);
    CALL_SUB_PROG(check_jump_type_and_folding);
    CALL_SUB_PROG(pare_hexmesh_from_surface);
    CALL_SUB_PROG(extend_tetmesh);
    CALL_SUB_PROG(pare_hex_outside_tet_surface);
    CALL_SUB_PROG(relax_singularities);
    CALL_SUB_PROG(calc_tet_rotation);
    CALL_SUB_PROG(arap_to_concentrate_points);
    CALL_SUB_PROG(remove_degenerated_edges);
    CALL_SUB_PROG(check_bad_rounding);
    CALL_SUB_PROG(tet_mesh_inflation_for_ball);
    CALL_SUB_PROG(aggressive_extract_jump_type);
    CALL_SUB_PROG(deform_hex_to_original_tet);
    CALL_SUB_PROG(split_tet_around_singularities);
    CALL_SUB_PROG(remove_near_miss_by_minimal_cut);
    CALL_SUB_PROG(remove_surface_wedge);
    CALL_SUB_PROG(remove_surface_degenerated_patch_in_cut);
    CALL_SUB_PROG(remove_surface_degenerated_patch_by_collapsing);
    CALL_SUB_PROG(glue_cut_tet_according_to_cut_face);
    CALL_SUB_PROG(check_polycube_validation_with_equation_graph);
    CALL_SUB_PROG(extract_relax_surface);
    CALL_SUB_PROG(visualize_normal);
    CALL_SUB_PROG(equation_graph_check);
    CALL_SUB_PROG(get_signed_surface_type);
    CALL_SUB_PROG(label_surface_type_after_param);
    CALL_SUB_PROG(revert_polycube_rotation);
    CALL_SUB_PROG(fairing_obj);
    CALL_SUB_PROG(fairing_tet);
    CALL_SUB_PROG(adaptive_frame_field_generation_with_polycube);
    CALL_SUB_PROG(inteploate_frame_field);
    CALL_SUB_PROG(map_tets);
    CALL_SUB_PROG(map_tris);
    CALL_SUB_PROG(meta_mesh_extraction);
    CALL_SUB_PROG(extract_sheets);
    CALL_SUB_PROG(shell_for_polycube_hex);
    CALL_SUB_PROG(frame_field_check);
  } catch(std::exception& e) {
    cerr << endl;
    cerr << "[error] "<<e.what() << endl;
    cerr << "Usage:" << endl;
    zjucad::show_usage_info(std::cerr, pt);
  }
  return 0;
}
