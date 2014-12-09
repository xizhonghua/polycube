#ifndef AREA_SELECT_CALLBACK_H
#define AREA_SELECT_CALLBACK_H
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkActor.h>
#include <vtkCommand.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkDataSetMapper.h>
#include <vtkDataSetSurfaceFilter.h>
#include "../configure/data_type.h"

namespace lq {

  class area_select_callback : public vtkCommand
  {

 public :
    area_select_callback()
    {
      select_actor = vtkSmartPointer<vtkActor>::New();
      select_render = vtkSmartPointer<vtkRenderer>::New();
      unstructured = vtkSmartPointer<vtkUnstructuredGrid>::New();
      sphere = vtkSmartPointer<vtkSphereSource>::New();
      select_poly = vtkSmartPointer<vtkPolyData>::New();
      select_mapper = vtkSmartPointer<vtkDataSetMapper>::New();
      surface_filter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
      
    }
    static area_select_callback *New()
    {
      return new area_select_callback;
    }

    void get_cell_list(std::vector<triangle> &cell);
    void get_cell_index(std::vector<size_t> &index);
    void get_mapper(vtkSmartPointer<vtkDataSetMapper> &mapper);
    size_t get_mapper_num();
    void get_mapper(size_t index, vtkSmartPointer<vtkDataSetMapper> &mapper);
    void set_actor(vtkSmartPointer<vtkActor> &actor);
    void set_render(vtkSmartPointer<vtkRenderer> &render);
    void set_data(vtkSmartPointer<vtkUnstructuredGrid> &data);
    void set_sphere(vtkSmartPointer<vtkSphereSource> &s);
    void set_poly(vtkSmartPointer<vtkPolyData> &poly);
    void set_mapper(vtkSmartPointer<vtkDataSetMapper> &mapper);
    void set_surface_filter(const vtkSmartPointer<vtkDataSetSurfaceFilter> &filter);
    void set_cell_index(std::vector<size_t> &index);
    void set_pick_type(const size_t flag);
    void search_point_cell(size_t point_index,
                           vtkSmartPointer<vtkIdList> index_list);
    virtual void Execute(vtkObject *obj, unsigned long, void *cell_data);
    
 private:
    
    vtkSmartPointer<vtkActor> select_actor;
    vtkSmartPointer<vtkRenderer> select_render;
    vtkSmartPointer<vtkUnstructuredGrid> unstructured;
    vtkSmartPointer<vtkSphereSource> sphere;
    vtkSmartPointer<vtkPolyData> select_poly;
    vtkSmartPointer<vtkDataSetMapper> select_mapper;
    std::vector<vtkSmartPointer<vtkDataSetMapper> > select_mapper_list;
    vtkSmartPointer<vtkDataSetSurfaceFilter> surface_filter;
    std::vector<triangle> cell_list;
    std::vector<size_t> cell_index;
    size_t pick_type;
    int get_cell_index(const std::vector<triangle> &tri_list);
  };
}
#endif
