#ifndef LQ_VTK_IO_H
#define LQ_VTK_IO_H

#include "../configure/data_type.h"
#include "../configure/load_configure.h"
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkDataSetMapper.h>
#include <vtkUnsignedCharArray.h>
#include <string>
#include <vector>
#include <fstream>

namespace lq {


  
  class vtk_io {

 public:
    
    vtk_io();
    ~vtk_io();
    std::string get_data_type();
    std::string get_data_set();
    std::string get_point_type();
    void set_configure(const load_configure &config);
    size_t get_cell_type();
    void set_cell_index(std::vector<size_t> &index);
    int load_vtk(const char *path, vtkSmartPointer<vtkPolyData> &poly,
                 vtkSmartPointer<vtkPolyDataMapper> &mapper,
                 vtkSmartPointer<vtkDataSetMapper> &line_mapper,
                 vtkSmartPointer<vtkDataSetMapper> &point_mapper,
                 vtkSmartPointer<vtkUnsignedCharArray> &cell_color);

    void reset_color(const std::vector<size_t> &cell_index,
                     const size_t flag);
    void output_cell();
    void clear();
    
 private:
    
    load_configure configure;
    enum CELL_TPYE{POINT = 1, LINE, TRIANGLE,ERROR};
    size_t min_type;
    size_t max_type;
    size_t cell_type;
    std::string data_type;
    std::string data_set;
    std::string point_type;
    std::vector<size_t> cell_look_table;
    std::vector<point> vertex_list;
    std::vector<triangle> triangle_list;
    std::vector<size_t> point_cell_list;
    std::vector<size_t> cell_index;
    std::string file_name;
    int read_cell(string &input, const size_t index, ifstream &fin);
    void read_cell_data(ifstream &fin, const size_t index,
                        std::string &input);
    void color_map(const std::vector<size_t> &cell_type,
                   vtkSmartPointer<vtkPolyData> &poly,
                   vtkSmartPointer<vtkUnsignedCharArray> &cell_color);
    
    
  };

}

#endif
