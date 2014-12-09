HEX_UI=/home/tfjiang/Projects/hex/build/bin/hex_ui

TET_FILE=$(wildcard ../tet/*.tet)
ZYZ_FILE=$(wildcard *.inner.with_surface_normal.zyz)
TEMP_ZYZ=$(patsubst %.inner.with_surface_normal.zyz,%.inner.with_surface_normal.temp.zyz, $(ZYZ_FILE))
NEW_TET_FILE=$(patsubst %.inner.with_surface_normal.zyz,%.inner.with_surface_normal.new.tet, $(ZYZ_FILE))
INNER_FACE_TYPE=$(patsubst %.inner.with_surface_normal.zyz,%.inner.with_surface_normal.zyz.inner_face_jump_type, $(ZYZ_FILE))
SURFACE_FACE_TYPE=$(patsubst %.inner.with_surface_normal.zyz,%.inner.with_surface_normal.zyz.surface_type, $(ZYZ_FILE))

MODIFIED_TOPOLOGY=$(TEMP_ZYZ) $(NEW_TET_FILE) $(INNER_FACE_TYPE) $(SURFACE_FACE_TYPE)

default: $(TET_FILE) $(HEX_UI) $(MODIFIED_TOPOLOGY)
	$(HEX_UI) prog=aggressive_extract_jump_type tet=$(NEW_TET_FILE) zyz=$(TEMP_ZYZ) inner_face_jump_type=$(INNER_FACE_TYPE) surface_type=$(SURFACE_FACE_TYPE)
#	$(HEX_UI) prog=remove_near_miss_by_minimal_cut cut_tet=cut_tet.tet  tet=$(NEW_TET_FILE) inner_face_jump_type=inner_face_jump_type.debug surface_type=surface_type.debug g_unknown_face=g_unknown_jump_face fnode_group=group_file
#	$(HEX_UI) prog=aggressive_extract_jump_type tet=$(NEW_TET_FILE) zyz=$(TEMP_ZYZ) inner_face_jump_type=inner_face_jump_type_after_removing_near_miss surface_type=$(SURFACE_FACE_TYPE)
#	$(HEX_UI) prog=remove_near_miss_by_minimal_cut cut_tet=cut_tet.tet  tet=$(NEW_TET_FILE) inner_face_jump_type=inner_face_jump_type.debug surface_type=surface_type.debug g_unknown_face=g_unknown_jump_face fnode_group=group_file
#	$(HEX_UI) prog=remove_surface_wedge tet=$(NEW_TET_FILE) cut_tet=cut_tet.tet inner_face_jump_type=inner_face_jump_type_after_removing_near_miss surface_type=surface_type.debug loop_points=loop_points	g_unknown_face=g_unknown_jump_face
#	$(HEX_UI) prog=aggressive_extract_jump_type tet=$(NEW_TET_FILE) zyz=$(TEMP_ZYZ) inner_face_jump_type=inner_face_jump_type_after_removing_near_miss surface_type=surface_type_original_after_remove_fractional_patch

$(MODIFIED_TOPOLOGY): $(TET_FILE) $(HEX_UI)
	$(HEX_UI) prog=refine_frame_field_after_aligned tet=$(TET_FILE) zyz=$(ZYZ_FILE)  surface_align_w=1 jump_smooth_w=1 dump_out_inner_face_jump_type=$(INNER_FACE_TYPE) dump_out_surface_type=$(SURFACE_FACE_TYPE) new_tet=$(NEW_TET_FILE) temp_zyz=$(TEMP_ZYZ)

$(TET_FILE) $(HEX_UI) $(ZYZ_FILE):

clean:
	-rm $(MODIFIED_TOPOLOGY)
