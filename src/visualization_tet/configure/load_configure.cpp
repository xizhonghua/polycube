#include "load_configure.h"
#include <fstream>
#include <iostream>
#include <string.h>
#include <cstdlib>
namespace lq {

using namespace std;

void load_configure::get_background_color(double color[3])
{
  for(size_t i = 0 ;i < 3; ++i)
    color[i] = this->backgroubd_color[i];
}

void load_configure::get_surface_edge_color(double color[3])
{
  for(size_t i = 0; i < 3; ++i)
    color[i] = this->surface_edge_color[i];
}

double load_configure::get_loop_edge_width()
{
  return this->loop_edge_width;
}

double load_configure::get_loop_point_size()
{
  return this->loop_point_size;
}

double load_configure::get_surface_edge_width()
{
  return this->surface_edge_width;
}

void load_configure::load_file(const char* path)
{
  ifstream fin;
  // this->file_path = path;
  fin.open(path);
  string input, type_flag, data_str;
  if(fin.fail())
  {
    cerr << "fail to open : " << path << endl;
    cerr << "error happen in : " << __FILE__ << " line : "
         << __LINE__ <<endl;
    return ;
  }
  while(getline(fin, input))
  {
    size_t space_flag;
    space_flag = input.find(" ");
    type_flag = input.substr(0, space_flag);
    data_str = input.substr(space_flag +1);
    if(input[0] == '#')
    {
      input.clear();
      continue;
    }
    if(strcmp(type_flag.c_str(), "background_color") == 0)
      sscanf(data_str.c_str(), "%lf %lf %lf", &backgroubd_color[0],
             &backgroubd_color[1], &backgroubd_color[2]);
    else if(strcmp(type_flag.c_str(), "surface_edge_color") == 0)
      sscanf(data_str.c_str(), "%lf %lf %lf", &surface_edge_color[0],
             &surface_edge_color[1], &surface_edge_color[2]);
    else if(strcmp(type_flag.c_str(), "loop_edge_width") == 0)
      loop_edge_width = atof(data_str.c_str());
    else if(strcmp(type_flag.c_str(), "loop_point_size") == 0)
      loop_point_size = atof(data_str.c_str());
    else if(strcmp(type_flag.c_str(), "surface_edge_width") == 0)
      surface_edge_width = atof(data_str.c_str());
    input.clear();
  }
}
}
