#ifndef LOAD_CONFIGURE_H
#define LOAD_CONFIGURE_H
#include <string>

namespace lq {

  class load_configure {

 public:

    static load_configure *New()
    {
      return new load_configure;
    }
    void get_background_color(double color[3]);
    void get_surface_edge_color(double color[3]);
    double get_loop_edge_width();
    double get_loop_point_size();
    double get_surface_edge_width();
    void load_file(const char* path);
    
 private:
    std::string file_path;
    double backgroubd_color[3];
    double surface_edge_color[3];
    double loop_edge_width;
    double loop_point_size;
    double surface_edge_width;
  };

}

#endif
