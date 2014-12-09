#ifndef HJ_MESH_ROUTINE_CONF_H_
#define HJ_MESH_ROUTINE_CONF_H_

#ifdef WIN32
#  ifdef hj_mesh_EXPORTS
	#define HJ_MESH_API __declspec(dllexport)
#  endif
#endif

#ifndef HJ_MESH_API 
#  define HJ_MESH_API 
#endif

#endif
