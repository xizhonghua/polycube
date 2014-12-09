
#include <iostream>
#include <string.h>

using namespace std;

#define CALL_SUB_PROG(name, prog) \
  int prog(int argc, char *argv[]); \
  if(!strcmp(argv[1], name)) \
  return prog(argc-1, argv+1);

int main(int argc, char *argv[])
{
  if(argc < 2) {
      cerr << "Usage: utils sub-program" << endl;
      return -1;
    }

  CALL_SUB_PROG("tet2vtk", tet2vtk);
  CALL_SUB_PROG("hex2vtk", hex2vtk);
  CALL_SUB_PROG("tet_surf2obj", tet_surf2obj);
  CALL_SUB_PROG("vtk_tet_surf2obj", vtk_tet_surf2obj);
  CALL_SUB_PROG("hex_surf2obj", hex_surf2obj);
  CALL_SUB_PROG("face_in_tet", face_in_tet);
  CALL_SUB_PROG("interp", interp);
  CALL_SUB_PROG("show_sh_func", show_sh_func);
  CALL_SUB_PROG("zyz2color", zyz2color);
  CALL_SUB_PROG("box_frame", box_frame);
  CALL_SUB_PROG("fast_marching", fast_marching);
  CALL_SUB_PROG("depth2stiff", depth2stiff);
  CALL_SUB_PROG("init_prism", init_prism);
  CALL_SUB_PROG("dump_surf_zyz2fv", dump_surf_zyz2fv);
  CALL_SUB_PROG("dump_tet_surf_normal", dump_tet_surf_normal);
  CALL_SUB_PROG("dump_fv_to_cross", dump_fv_to_cross);
  CALL_SUB_PROG("dump_pv_to_cross", dump_pv_to_cross);
  CALL_SUB_PROG("dump_fv_to_N_field", dump_fv_to_N_field);
  CALL_SUB_PROG("draw_surface_zyz", draw_surface_zyz);
  CALL_SUB_PROG("pull_hexmesh",pull_hexmesh);
  CALL_SUB_PROG("interp_tet",interp_tet);
  CALL_SUB_PROG("zyz_mapping", zyz_mapping);
  CALL_SUB_PROG("dump_jump_type",dump_jump_type);
  CALL_SUB_PROG("merge_hex_vertex",merge_hex_vertex);
  CALL_SUB_PROG("merge_obj",merge_obj);
  CALL_SUB_PROG("test", test);
  CALL_SUB_PROG("dump_surface_type",dump_surface_type);
  CALL_SUB_PROG("dump_surface_normal_type",dump_surface_normal_type);
  CALL_SUB_PROG("dump_jump_face_to_vtk",dump_jump_face_to_vtk);
  CALL_SUB_PROG("sc2vtk", dump_sc_to_vtk);
  CALL_SUB_PROG("mesh2tet", mesh2tet);
  CALL_SUB_PROG("mesh2hex", mesh2hex);
  CALL_SUB_PROG("mesh2obj", mesh2obj);
  CALL_SUB_PROG("hex2obj", hex2obj);
  CALL_SUB_PROG("transform_cell", transform_cell);
  CALL_SUB_PROG("find_hex_singularity",find_hex_singularity);
  CALL_SUB_PROG("find_quad_singularity",find_quad_singularity);
  CALL_SUB_PROG("check_type_setting_files", check_type_setting_files);
  CALL_SUB_PROG("frame2surfv_angle", frame2surfv_angle);
  CALL_SUB_PROG("extract_tets_inside_box", extract_tets_inside_box);
  CALL_SUB_PROG("split_not_aligned_face", split_not_aligned_face);
  CALL_SUB_PROG("check_polycube_surface", check_polycube_surface);
  CALL_SUB_PROG("check_polycube_surface_with_rotation", check_polycube_surface_with_rotation);
  CALL_SUB_PROG("split_surface_tets", split_surface_tets);
  CALL_SUB_PROG("cut_surface_type2orig_type", cut_surface_type2orig_type);
  CALL_SUB_PROG("remove_degenerated_tet", remove_degenerated_tet);
  CALL_SUB_PROG("check_flipped_tet", check_flipped_tet);
  CALL_SUB_PROG("translate_surface_type_from_lq", translate_surface_type_from_lq);
  CALL_SUB_PROG("remove_extra_nodes", remove_extra_nodes);
  CALL_SUB_PROG("improve_tet", improve_tet);
  CALL_SUB_PROG("improve_hex", improve_hex);
  CALL_SUB_PROG("improve_quad",improve_quad);
  CALL_SUB_PROG("check_flipped_surface", check_flipped_surface);
  CALL_SUB_PROG("compress_obj", compress_obj);
  CALL_SUB_PROG("subdivide_hexmesh", subdivide_hexmesh);
  CALL_SUB_PROG("split_edge_on_polycube", split_edge_on_polycube);
  CALL_SUB_PROG("validate_tetmesh", validate_tetmesh);
  CALL_SUB_PROG("vis_arap_distortion", vis_arap_distortion);
  CALL_SUB_PROG("vis_angle_area_distortion", vis_angle_area_distortion);
  CALL_SUB_PROG("vis_slim_part", vis_slim_part);
  CALL_SUB_PROG("extract_polycube_edge", extract_polycube_edge);
  CALL_SUB_PROG("diffustion_weight", diffustion_weight);
  CALL_SUB_PROG("extend_hexmesh2", extend_hexmesh2);
  CALL_SUB_PROG("extend_hexmesh", extend_hexmesh);
  CALL_SUB_PROG("extend_hexmesh3", extend_hexmesh3);
  CALL_SUB_PROG("dump_surface_patch_for_szy", dump_surface_patch_for_szy);
  CALL_SUB_PROG("tet_jac", tet_jacobian);
  CALL_SUB_PROG("hex_jac", hex_jacobian);
  CALL_SUB_PROG("find_polycube_edge", find_polycube_edge);
  CALL_SUB_PROG("obj2faketet",obj2faketet);
  CALL_SUB_PROG("polycube_concave_detection", polycube_concave_detection);
  CALL_SUB_PROG("draw_quad_type", draw_quad_type);
  CALL_SUB_PROG("obj2vtk", obj2vtk);
  CALL_SUB_PROG("obj2obj", obj2obj);
  CALL_SUB_PROG("vtk2tet", vtk2tet);
  CALL_SUB_PROG("vtk2hex", vtk2hex);
  CALL_SUB_PROG("yf2jtf",yf2jtf);
  CALL_SUB_PROG("tf2yf", tf2yf);
  CALL_SUB_PROG("vol2tet", vol2tet);
  CALL_SUB_PROG("yfsingularity2vtk",yfsingularity2vtk);
  CALL_SUB_PROG("draw_singularity_yfl",draw_singularity_yfl);
  CALL_SUB_PROG("draw_feature_line", draw_feature_line);
  CALL_SUB_PROG("rigid_shape_error", rigid_shape_error);
  CALL_SUB_PROG("extract_quad_polyedge", extract_quad_polyedge);
  CALL_SUB_PROG("simplify_quad_line", simplify_quad_line);
  CALL_SUB_PROG("orient_tet", orient_tet);
  CALL_SUB_PROG("dump_surface_feature", dump_surface_feature);
  CALL_SUB_PROG("fit_position", fit_position);
  CALL_SUB_PROG("dump_surface_patch", dump_surface_patch);
  CALL_SUB_PROG("extract_simp_quad_line",extract_simp_quad_line);
  CALL_SUB_PROG("remove_degenerated_tri",remove_degenerated_tri);
  CALL_SUB_PROG("find_shortest_path",find_shortest_path);
  CALL_SUB_PROG("lin_para", lin_para);
  CALL_SUB_PROG("normal_map", normal_map);
  CALL_SUB_PROG("smooth_cell", smooth_cell);
  CALL_SUB_PROG("sgp_deform", sgp_deform);
  CALL_SUB_PROG("tri_size2quad", tri_size2quad);
  CALL_SUB_PROG("heying2obj",heying2obj);
  CALL_SUB_PROG("off2cell", off2cell);
  CALL_SUB_PROG("zyz2frame", zyz2frame);
  CALL_SUB_PROG("fix_type", fix_type);
  CALL_SUB_PROG("optimize_polycube_surface",optimize_polycube_surface);
  CALL_SUB_PROG("tet2ascii", tet2ascii);
  CALL_SUB_PROG("draw_zyz_normal_error", draw_zyz_normal_error);
  CALL_SUB_PROG("map_normal_type_to_cut_tet", map_normal_type_to_cut_tet);
  CALL_SUB_PROG("dump_cut_face_pair2vtk", dump_cut_face_pair2vtk);
  CALL_SUB_PROG("check_redundant_node", check_redundant_node);
  CALL_SUB_PROG("simply_cut_tet", simply_cut_tet);
  CALL_SUB_PROG("dump_surface_zyz_type", dump_surface_zyz_type);
  CALL_SUB_PROG("dump_wireframe", dump_wireframe);
  CALL_SUB_PROG("fit_centroid", fit_centroid);
  CALL_SUB_PROG("sphere_obj",sphere_obj);
  CALL_SUB_PROG("sphere_obj_with_boundary",sphere_obj_with_boundary);
  CALL_SUB_PROG("sphere_tet",sphere_tet);
  CALL_SUB_PROG("extract_4sys_field_singularity",extract_4sys_field_singularity);
  CALL_SUB_PROG("transplant_feature_line", transplant_feature_line);
  CALL_SUB_PROG("map_zyz",map_zyz);
  CALL_SUB_PROG("select_edges", select_edges);
  CALL_SUB_PROG("split_tets_on_edges", split_tets_on_edges);
  CALL_SUB_PROG("collapse_tet_edges", collapse_tet_edges);
  CALL_SUB_PROG("draw_frame_cube", draw_frame_cube);
  CALL_SUB_PROG("dump_singularity_vtk_to_netgen_size", dump_singularity_vtk_to_netgen_size);
  CALL_SUB_PROG("draw_frame_axis", draw_frame_axis);
  CALL_SUB_PROG("draw_frame_axis_on_surface", draw_frame_axis_on_surface);
  CALL_SUB_PROG("obj2hex", obj2hex);
  CALL_SUB_PROG("deform_gradient_2d",deform_gradient_2d);
  CALL_SUB_PROG("test4",test4);
  CALL_SUB_PROG("test5",test5);
  CALL_SUB_PROG("mesh_curvature", mesh_curvature);
  CALL_SUB_PROG("hj_boundary",hj_boundary);
  CALL_SUB_PROG("map_obj2hex", map_obj2hex);
  cerr << "no such sub-program." << endl;
  return -2;
}
