/**
 *@class vtk_operation
 *@brief support some operation about vtk library

 *@detail
 *(1).initialition vtk render
 *(2).load the data of a tet model
 *(3).add box and plane
 *(4).get tet index of being included by box
 *(5).get tet index of being crossed bu plane
 *(6).rotate the box around one axis which set by user
 *(7).save box position into the file
 *(8).load the box position file and add its in render
 *
 *@author Li Quan
 *@data 2012-11-20
 *@version 1.0
 */


#ifndef VTK_OPERATION_H
#define VTK_OPERATION_H

#include "../configure/data_type.h"
#include "../widget/box_widget.h"
#include "../call_back/box_callback.h"
#include "QVTKWidget.h"
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkCellData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkPoints.h>
#include <vtkLineWidget2.h>
#include <vtkLineRepresentation.h>
#include <vtkDataSetMapper.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkStructuredGridReader.h>
#include <vector>


namespace lq {
  class plane_callback;
  class plane_widget;
  class vtk_operation {
 public:

    vtk_operation();
    void set_points(const vtkSmartPointer<vtkPoints> &point);
    void set_qvtk(QVTKWidget *qvtk_widget);
    // void set_render(const vtkSmartPointer<vtkRenderer> &r);
    void init_tet_render(const char *path, matrixd &point,
                         matrixst &tet, matrixst &tri,
                         vtkSmartPointer<vtkRenderer> &render);
    void init_vtk_render(const char *path, const size_t flag_win,
                         vtkSmartPointer<vtkRenderer> &render);
    void init_obj_render(QVTKWidget *qvtk_widget, const char *path,
                         vtkSmartPointer<vtkRenderer> &render);

    void reset_render(size_t flag);
    void add_box_widget(const matrixd &point,
                        QVTKWidget *qvtk_widget,
                        std::vector<vtkSmartPointer<box_widget> > &boxwidget,
                        const vtkSmartPointer<box_callback> &callback);

    void add_plane_widget(const matrixd &point,
                          QVTKWidget *qvtk_widget,
                          std::vector<vtkSmartPointer<plane_widget> > &planewidget,
                          const vtkSmartPointer<plane_callback> &callback);

    void get_tet_data(std::vector<vtkSmartPointer<box_widget> > &boxwidget,
                      const matrixd &tet);
    
    int  get_plane_tet(const std::vector<vtkSmartPointer<plane_widget> > &widget,
                        const matrixd &tet, const matrixd &vertex);
    void rotate_box(const int zoom, const int flag,
                    const double axis[3], vtkSmartPointer<box_widget> &box);

    void output_box_position(const std::vector<vtkSmartPointer<box_widget> > &boxwidget,
                             ofstream &fout);

    void import_box_position(ifstream &fin, QVTKWidget *qvtk_widget,
                             std::vector<vtkSmartPointer<box_widget> > &boxwidget,
                             const vtkSmartPointer<box_callback> &callback);

    void add_line_widget(QVTKWidget *qvtk_widget, const vtkSmartPointer<vtkRenderer> &render,
                         std::vector<vtkSmartPointer<vtkLineWidget2> > &line_list);
    void add_axes(vtkSmartPointer<vtkRenderer> &render);
    void add_area_picker(QVTKWidget *qvtk_widget, const size_t flag_win,
                         vtkSmartPointer<vtkRenderer> &render);
   
    void change_color(const size_t flag, vtkSmartPointer<vtkRenderer> &render,
                      const size_t win_flag);
    void display_surface_edge_mesh(const size_t flag_win);

    void output_cell_data(const size_t flag);
    void switch_pick_type(size_t flag, size_t win_flag);

    
    
 private:
   
    matrixd point;
    matrixst tet;
    matrixst tri;
    QVTKWidget *qvtk;
    vtkSmartPointer<vtkPoints> points;
    vtkSmartPointer<vtkPolyData> polydata;
    vtkSmartPointer<vtkActor> actor;
    vtkSmartPointer<vtkActor> actor_vtk;
    vtkSmartPointer<vtkUnstructuredGridReader> reader;
    vtkSmartPointer<vtkStructuredGridReader> reader_s;
    std::vector<size_t> cell_index;
   
  };
}




template <class T> const T& max ( const T& a, const T& b )
{
  return (a<b)?b:a;
}

template <class T> const T& min ( const T& a, const T& b )
{
  return !(b<a)?a:b;   
}

template <typename T>
inline void get_bounding(T point[4][3], T *bounding)
{
  T min_x, min_y, min_z;
  T max_x, max_y, max_z;
  min_x = max_x = point[0][0];
  min_y = max_y = point[0][1];
  min_z = max_z = point[0][2];
  for(size_t i = 1; i < 4; ++i)
  {
    min_x = min(min_x, point[i][0]);
    max_x = max(max_x, point[i][0]);
    min_y = min(min_y, point[i][1]);
    max_y = max(max_y, point[i][1]);
    min_z = min(min_z, point[i][2]);
    max_z = max(max_z, point[i][2]);
  }
  bounding[0] = min_x;
  bounding[1] = max_x;
  bounding[2] = min_y;
  bounding[3] = max_y;
  bounding[4] = min_z;
  bounding[5] = max_z;
}
#endif
