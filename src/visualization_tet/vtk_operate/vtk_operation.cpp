
/**
 * @file  vet_operation.cpp
 * @author Li Quan
 * @version 1.0
 * @date 2012-11-20
 * 
 * @ingroup visualization_tet
 * @brief support some operation about vtk library
 * @details
 *-initialition vtk render
 *-load the data of a tet model
 *-add box and plane
 *-get tet index of being included by box
 *-get tet index of being crossed bu plane
 *-rotate the box around one axis which set by user
 *-save box position into the file
 *-load the box position file and add its in render
 * 
 * @history
 *
 */

#include "vtk_operation.h"
#include "../../tetmesh/tetmesh.h"
#include "../../common/util.h"
#include "../widget/plane_widget.h"
#include "../call_back/plane_callback.h"
#include "../call_back/select_interactor_style.h"
#include "../call_back/area_select_callback.h"
#include "../io/vtk_io.h"

#include <jtflib/mesh/io.h>

//#include "load_configure.h"
#include <vtkTetra.h>
#include <vtkCellType.h>
#include <vtkCellArray.h>
#include <vtkTriangle.h>
#include <vtkCellData.h>
#include <vtkClipPolyData.h>
#include <vtkStructuredGridGeometryFilter.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPointData.h>
#include <vtkCommand.h>
#include <vtkBoxRepresentation.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkTransform.h>
#include <vtkMatrix4x4.h>
#include <vtkProperty.h>

#include <vtkOrientationMarkerWidget.h>
#include <vtkPropAssembly.h>
#include <vtkAxesActor.h>
#include <vtkInteractorStyleRubberBandPick.h>
#include <vtkCallbackCommand.h>
#include <vtkHardwareSelector.h>
#include <vtkRenderedAreaPicker.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkSelectionNode.h>
#include <vtkSelection.h>
#include <vtkExtractSelectedPolyDataIds.h>
#include <vtkTriangle.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkIdList.h>
#include <cmath>


#include <vtkRenderWindow.h>
#include <vtkAreaPicker.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkPolyDataMapper.h>
#include <vtkLookupTable.h>
#include <vtkUnsignedCharArray.h>
//#define TEST_ROTATE_BOX
//#define TEST_MATRIX4X4
namespace lq {
using namespace std;
using namespace zjucad::matrix;

#define CREATE_NEW(class, variable)                              \
  vtkSmartPointer<class> variable = vtkSmartPointer<class>::New();


enum TYPE{POINT = 1, LINE, TRIANGLE, ERROR};
const double box_zoom = 10.0;
const double pos_angle = 10.0;
const double neg_angle = -10.0;
const double pos_min_angle = 1.0;
const double neg_min_angle = -1.0;
const double PI = 3.14159265;
const double WHEAT[3] = {0.9608, 0.8706, 0.7020};
int num = 0;
bool flag = false;
size_t win_selected = 0;
CREATE_NEW(vtkActor, select_actor);
CREATE_NEW(vtkRenderer, render_call);
CREATE_NEW(vtkUnstructuredGrid, unstructured);
CREATE_NEW(vtkDataSetMapper, mapper);
CREATE_NEW(vtkDataSetMapper, vtk_point_mapper_left);
CREATE_NEW(vtkDataSetMapper, vtk_line_mapper_left);
CREATE_NEW(vtkDataSetMapper, vtk_point_mapper_right)
CREATE_NEW(vtkDataSetMapper, vtk_line_mapper_right)
CREATE_NEW(vtkPolyDataMapper, poly_mapper);
CREATE_NEW(vtkPolyDataMapper, poly_mapper_left);
CREATE_NEW(vtkPolyDataMapper, poly_mapper_right);
CREATE_NEW(vtkPolyData, tet_poly);
CREATE_NEW(vtkPolyData, obj_poly);
CREATE_NEW(vtkPolyData, vtk_poly_left);
CREATE_NEW(vtkPolyData, vtk_poly_right);
CREATE_NEW(vtkUnsignedCharArray, cell_color_left);
CREATE_NEW(vtkUnsignedCharArray, cell_color_right);
CREATE_NEW(vtkActor, vtk_point_actor_left);
CREATE_NEW(vtkActor, vtk_line_actor_left);
CREATE_NEW(vtkActor, vtk_point_actor_right);
CREATE_NEW(vtkActor, vtk_line_actor_right);
CREATE_NEW(vtkActor, actor_left);
CREATE_NEW(vtkActor, actor_right);
CREATE_NEW(vtkDataSetSurfaceFilter, surface_filter);
CREATE_NEW(area_select_callback, call_back_left);
CREATE_NEW(area_select_callback, call_back_right);
vtk_io vtk_file_op;
vtk_io vtk_file_op_left, vtk_file_op_right;
load_configure configure;
size_t cnt_surface_edge_left = 0;
size_t cnt_surface_edge_right = 0;
vector<point> vertex_list;
vector<triangle> tri_list;
//vtkSmartPointer<vtkRenderer> render_call; 
//vtkSmartPointer<vtkDataSetMapper> mapper;
//vtkSmartPointer<vtkUnstructuredGrid> unstructured;
vtk_operation::vtk_operation()
{
  // points = vtkSmartPointer<vtkPoints>::New();
  polydata = vtkSmartPointer<vtkPolyData>::New();
  //mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  actor = vtkSmartPointer<vtkActor>::New();
  actor_vtk = vtkSmartPointer<vtkActor>::New();
  //  unstructured = vtkSmartPointer<vtkUnstructuredGrid>::New();
  reader =  vtkSmartPointer<vtkUnstructuredGridReader>::New();
  reader_s =  vtkSmartPointer<vtkStructuredGridReader>::New();
  configure.load_file("./configure.txt");
  //  mapper = vtkSmartPointer<vtkDataSetMapper>::New();
}

void vtk_operation::set_points(const vtkSmartPointer<vtkPoints> &point)
{
  this->points = point;
}

void vtk_operation::set_qvtk(QVTKWidget *qvtk_widget)
{
  this->qvtk = qvtk_widget;
}
int get_obj(const char *path,
            std::vector<point> &vertex_list,
            std::vector<triangle> &triangle_list)
{
  std::ifstream fin;
  //std::ofstream fout("test.txt");
  fin.open(path);
  if(fin.fail())
  {
    std::cerr << "error happen at:" << std::endl;
    std::cerr << "file: " << __FILE__ << std::endl;
    std::cerr << "line: " << __LINE__ << std::endl;
    return 1;
  }
  else
    std::cout << "file " << path << " open successful!" << std::endl;
  std::string input;
  std::string sub_str;
  point tmp_vertex;
  int index1, index2, index3;
  while(getline(fin, input))
  {
    if(input[0] == 'v' && input[1] == 't')
      continue;
    else if(input[0] == 'v' && input[1] == 'n')
      continue;
    else if(input[0] == 'v')
    {
      sub_str = input.substr(1);
      sscanf(sub_str.c_str(), "%lf %lf %lf", &tmp_vertex.x, &tmp_vertex.y, &tmp_vertex.z);
      vertex_list.push_back(tmp_vertex);
    }
    else if(input[0] == 'f')
    {
      triangle tmp_triangle;
      // triangle tmp_triangle;
      int tmp_index, texture, begin = 0, end = 0;
      input = input.substr(2);
      end = input.find(" ");
      std::string tmp;
      tmp = input.substr(begin, end);
      int flag;
      flag = sscanf(tmp.c_str(), "%d/%d", &tmp_index, &texture);
      if(flag == 2)
      {
        while(end != std::string::npos)
        {
          sub_str = input.substr(begin, end);
          sscanf(sub_str.c_str(), "%d/%d", &tmp_index, &texture);
          tmp_triangle.vertex.push_back(tmp_index);
          input = input.substr(end + 1, input.length() - 1);
          end = input.find(" ");
        }
        sscanf(input.c_str(), "%d/%d", &tmp_index, &texture);
        tmp_triangle.vertex.push_back(tmp_index);
        triangle_list.push_back(tmp_triangle);
        tmp_triangle.vertex.clear();
      }
      else if(flag == 1)
      {
        sscanf(input.c_str(), "%d %d %d", &index1, &index2, &index3);
        tmp_triangle.vertex.push_back(index1);
        tmp_triangle.vertex.push_back(index2);
        tmp_triangle.vertex.push_back(index3);
        triangle_list.push_back(tmp_triangle);
      }
    }
    input.clear();
    sub_str.clear();
  }
  std::cout << "right tri_size " << triangle_list.size() << std::endl;
  return 0;
}

// void vtk_operation::set_render(const vtkSmartPointer<vtkRenderer> &r)
// {
//   this->renderer = r;
// }
/**
 *@brief load tet model's data into point, tet, tri, and initialization render ready for
 *display the model
 *
 *@param[in]
 *-path tet model's path
 *-point save point set of tet model
 *-tet save tet data set
 *-tri
 *-render vtkRenderer
 *
 *@return void
 */
void vtk_operation::init_tet_render(const char *path, matrixd &point,
                                    matrixst &tet, matrixst &tri,
                                    vtkSmartPointer<vtkRenderer> &render)
{
  // actor = vtkSmartPointer<vtkActor>::New();
  // mapper = vtkSmartPointer<vtkDataSetMapper>::New();
  if(flag)
  {
    render->RemoveActor(actor_vtk);
    flag = false;
  }
  point.resize(0, 0);
  tet.resize(0, 0);
  tri.resize(0, 0);
  jtf::mesh::tet_mesh_read_from_zjumat(path, &point, &tet, &tri);
  points->Reset();
  for(size_t i = 0; i < point.size(2); ++i)
    points->InsertNextPoint(point(0, i), point(1, i), point(2, i));
  unstructured->Initialize();
  unstructured->SetPoints(points);
  vtkSmartPointer<vtkCellArray> cell_array =
      vtkSmartPointer<vtkCellArray>::New();
  CREATE_NEW(vtkCellArray, tri_array);
  CREATE_NEW(vtkTriangle, triangle);
  for(size_t i = 0; i < tet.size(2); ++i)
  {
    vtkSmartPointer<vtkTetra> tetra =
        vtkSmartPointer<vtkTetra>::New();
    tetra->GetPointIds()->SetId(0, tet(0, i));
    tetra->GetPointIds()->SetId(1, tet(1, i));
    tetra->GetPointIds()->SetId(2, tet(2, i));
    tetra->GetPointIds()->SetId(3, tet(3, i));
    for(size_t i = 0 ; i < tetra->GetNumberOfFaces(); ++i)
    {
      vtkSmartPointer<vtkCell> tmp_cell;
      tmp_cell = tetra->GetFace(i);
      for(size_t j = 0; j < 3; ++j)
        triangle->GetPointIds()->SetId(j, tmp_cell->GetPointId(j));
      tri_array->InsertNextCell(triangle);
    }
    cell_array->InsertNextCell(tetra);
  }
  tet_poly->SetPoints(points);
  tet_poly->SetPolys(tri_array);
  vtkSmartPointer<vtkUnsignedCharArray> cell_color =
      vtkSmartPointer<vtkUnsignedCharArray>::New();
  cell_color->SetNumberOfComponents(3);
  unsigned char tmp_color[3] = {70, 93, 207};
  for(size_t i = 0; i < tet.size(2); ++i)
    cell_color->InsertNextTupleValue(tmp_color);
  // for(size_t i = 0; i < tri_array->GetNumberOfCells(); ++i)
  //   cell_color->InsertNextTupleValue(tmp_color);
  tet_poly->GetCellData()->SetScalars(cell_color);
  unstructured->SetCells(VTK_TETRA, cell_array);
  unstructured->GetCellData()->SetScalars(cell_color);
  surface_filter->SetInputConnection(unstructured->GetProducerPort());
  surface_filter->PassThroughCellIdsOff();
  surface_filter->Update();
#if VTK_MAJOR_VERSION <= 5
  poly_mapper->SetInputConnection(surface_filter->GetOutputPort());
#else
  poly_mapper->SetInputData(surface_filter);
#endif
    
// #if VTK_MAJOR_VERSION <= 5
//   poly_mapper->SetInputConnection(tet_poly->GetProducerPort());
// #else
//   poly_mapper->SetInputData(tet_poly);
// #endif
  actor->SetMapper(poly_mapper);
  actor->PickableOn();
  render->AddActor(actor);
 
  return;
}

/**
 *@brief load vtk model
 *
 *@param[in]
 *-path tet model's path
 *-render vtkRenderer
 *
 *@return void
 */

void vtk_operation::init_vtk_render(const char *path, const size_t flag_win,
                                    vtkSmartPointer<vtkRenderer> &render)
{
  flag = true;
  reader->SetFileName(path);
  win_selected = flag_win;
  size_t flag;
  if(reader->IsFileUnstructuredGrid() == 1)
  {
    cout << "Unstructed vtk file loaded!" <<endl;
    if(flag_win == LEFT)
    {
      vtk_file_op_left.clear();
      cell_index.clear();
      flag =  vtk_file_op_left.load_vtk(path, vtk_poly_left, poly_mapper_left,
                                        vtk_line_mapper_left,
                                        vtk_point_mapper_left, cell_color_left);
      vtk_file_op_left.set_cell_index(cell_index);
    }
    else if(flag_win == RIGHT)
    {
      vtk_file_op_right.clear();
      cell_index.clear();
      flag = vtk_file_op_right.load_vtk(path, vtk_poly_right, poly_mapper_right,
                                        vtk_line_mapper_right,
                                        vtk_point_mapper_right, cell_color_right);
      vtk_file_op.set_cell_index(cell_index);
    }
  }
  else if(reader->IsFileStructuredGrid() == 1)
  {
    cout << "vtk structed loaded!" <<endl;
    vtkSmartPointer<vtkDataSetMapper> mapper =
        vtkSmartPointer<vtkDataSetMapper>::New();
    reader_s->SetFileName(path);
    reader_s->Update();
    vtkSmartPointer<vtkStructuredGridGeometryFilter> geometryFilter =
        vtkSmartPointer<vtkStructuredGridGeometryFilter>::New();
    geometryFilter->SetInputConnection(reader->GetOutputPort());
    geometryFilter->Update();
    mapper->SetInputConnection(geometryFilter->GetOutputPort());
  }
  if(flag == POINT)
  {
    double size;
    size = configure.get_loop_point_size();
    if(flag_win == LEFT)
    {
      vtk_point_actor_left->SetMapper(vtk_point_mapper_left);
      vtk_point_actor_left->GetProperty()->SetPointSize(size);
      render->AddActor(vtk_point_actor_left);
    }
    else if(flag_win == RIGHT)
    {
      vtk_point_actor_right->SetMapper(vtk_point_mapper_right);
      vtk_point_actor_right->GetProperty()->SetPointSize(size);
      render->AddActor(vtk_point_actor_right);
    }
  }
  else if(flag == LINE)
  {
    double width;
    width = configure.get_loop_edge_width();
    if(flag_win == LEFT)
    {
      vtk_line_actor_left->SetMapper(vtk_line_mapper_left);
      vtk_line_actor_left->GetProperty()->EdgeVisibilityOn();
      vtk_line_actor_left->GetProperty()->SetLineWidth(width);
      render->AddActor(vtk_line_actor_left);
    }
    else if(flag_win == RIGHT)
    {
      vtk_line_actor_right->SetMapper(vtk_line_mapper_right);
      vtk_line_actor_right->GetProperty()->EdgeVisibilityOn();
      vtk_line_actor_right->GetProperty()->SetLineWidth(width);
      render->AddActor(vtk_line_actor_right);
    }
  }
  else if(flag == TRIANGLE)
  {
    if(flag_win == LEFT)
    {
      actor_left->SetMapper(poly_mapper_left);
      render->AddActor(actor_left);
    }
    else if(flag_win == RIGHT)
    {
      actor_right->SetMapper(poly_mapper_right);
      render->AddActor(actor_right);
    }
    // actor_vtk->SetMapper(poly_mapper);
    // render->AddActor(actor_vtk);
  }
  double backgroubd_color[3];
  configure.get_background_color(backgroubd_color);
  render->SetBackground(backgroubd_color);
  return;
}


void vtk_operation::init_obj_render(QVTKWidget *qvtk_widget, const char *path,
                                    vtkSmartPointer<vtkRenderer> &render)
{
  get_obj(path, vertex_list, tri_list);
  points = vtkSmartPointer<vtkPoints>::New();
  cout << " vertex_size " << vertex_list.size() << endl;
  for(size_t i = 0; i < vertex_list.size(); ++i)
    points->InsertNextPoint(vertex_list[i].x, vertex_list[i].y,
                            vertex_list[i].z);
  vtkSmartPointer<vtkTriangle> tmp_triangle =
      vtkSmartPointer<vtkTriangle>::New();
   vtkSmartPointer<vtkCellArray> triangle_cell =
       vtkSmartPointer<vtkCellArray>::New();

   for(size_t i = 0; i < tri_list.size(); ++i)
   {
     tmp_triangle->GetPointIds()->SetId(0, tri_list[i].vertex[0] - 1);
     tmp_triangle->GetPointIds()->SetId(1, tri_list[i].vertex[1] - 1);
     tmp_triangle->GetPointIds()->SetId(2, tri_list[i].vertex[2] - 1);
     triangle_cell->InsertNextCell(tmp_triangle);
   }
   obj_poly->SetPoints(points);
   obj_poly->SetPolys(triangle_cell);
   // #if VTK_MAJOR_VERSION <= 5
   //   mapper->SetInputConnection(obj_poly->GetProducerPort());
   // #else
   //   mapper->SetInputData(obj_poly);
   // #endif
#if VTK_MAJOR_VERSION <= 5
   poly_mapper->SetInputConnection(obj_poly->GetProducerPort());
#else
   poly_mapper->SetInputData(obj_poly);
#endif
   // CREATE_NEW(vtkActor, actor);
   actor->SetMapper(poly_mapper);
   actor->PickableOn();
   render->AddActor(actor);
}

void vtk_operation::reset_render(size_t flag)
{
  if(flag == LEFT)
  {
    vtk_poly_left->Initialize();
    vtk_file_op_left.clear();
  }
  else if(flag == RIGHT)
  {
    vtk_poly_right->Initialize();
    vtk_file_op_right.clear();
  }
}

/**
 *@brief add a contral box into the QVTKWidget
 *
 *@param[in]
 *-point point set of tet model
 *-qvtk_widget 
 *-boxwidget box widget set
 *-callback a smart pointer to a instance of callback class which get box interaction
 *
 *@return void
 */

void vtk_operation::add_box_widget(const matrixd &point,
                                   QVTKWidget *qvtk_widget,
                                   std::vector<vtkSmartPointer<box_widget> > &boxwidget,
                                   const vtkSmartPointer<box_callback> &callback)
{
  vtkSmartPointer<box_widget> box = vtkSmartPointer<box_widget>::New();
  vtkSmartPointer<vtkBoxRepresentation> boxRepresentation = 
      vtkSmartPointer<vtkBoxRepresentation>::New();
  vtkSmartPointer<vtkPolyData> box_boundary =
       vtkSmartPointer<vtkPolyData>::New();
  box->SetRepresentation(boxRepresentation);
#if TEST_ADD_BOX
  double *vertices;
  vertices =  boxRepresentation->GetBounds();
  for(size_t i = 0; i < 6; ++i)
    cout << "points" << vertices[i] << endl;
  vtkSmartPointer<vtkTransform> transform =
      vtkSmartPointer<vtkTransform>::New();
  boxRepresentation->GetTransform(transform);
  double orient[3];
  transform->GetOrientation(orient);
  cout << "before" << orient[0] << " " << orient[1] << " " << orient[2] <<endl;
#endif
  double diameter;
  diameter = calc_bounding_sphere_size(point) / box_zoom;
  
  double bounding[6] = {-diameter, diameter, -diameter,
                        diameter, -diameter, diameter};
  boxRepresentation->PlaceWidget(bounding);
  double *bound;
  bound = boxRepresentation->GetBounds();
  box->SetInteractor(qvtk_widget->GetInteractor());
  boxwidget.push_back(box);
  box->AddObserver(vtkCommand::EndInteractionEvent, callback);
  box->On();
  cout << "Add a box_widget..." << endl;
}

/**
 *@brief add a contral box into the QVTKWidget
 *
 *@param[in]
 *-point point set of tet model
 *-qvtk_widget 
 *-planwidget plane widget set
 *-callback a smart pointer to a instance of callback class which get box interaction
 *
 *@return void
 */

void vtk_operation::add_plane_widget(const matrixd &point,
                                     QVTKWidget *qvtk_widget,
                                     std::vector<vtkSmartPointer<plane_widget> > &planewidget,
                                     const vtkSmartPointer<plane_callback> &callback)
{
  vtkSmartPointer<plane_widget> plane =
      vtkSmartPointer<plane_widget>::New();
  plane->SetInteractor(qvtk_widget->GetInteractor());
  plane->AddObserver(vtkCommand::EndInteractionEvent, callback);
  plane->SetResolution(1);
  planewidget.push_back(plane);
  plane->On();
#ifdef TEST_PLANE
  double p1[3], p2[3];
  double p3[3];
  plane->GetPoint1(p1);
  plane->GetPoint2(p2);
  plane->GetOrigin(p3);
  cout << "p1 : " << p1[0] << " " << p1[1] << " " << p1[2] << endl;
  cout << "p2 : " << p2[0] << " " << p2[1] << " " << p2[2] << endl;
  cout << "p3 : " << p3[0] << " " << p3[1] << " " << p3[2] << endl;
  vtkSmartPointer<vtkPolyData> poly =
      vtkSmartPointer<vtkPolyData>::New();
  plane->GetPolyData(poly);
  size_t num = poly->GetNumberOfPoints();
  cout << "point num " << poly->GetNumberOfPoints() << endl;
  double p[3];
  for(size_t i = 0; i < num; ++i)
  {
    poly->GetPoint(i, p);
    cout << "point " << i << ":" << p[0] << " " << p[1] << " " << p[2] << endl;
  }
  double normal[3];
  plane->GetNormal(normal);
  cout << "normal :" << normal[0] << " " << normal[1] << " " << normal[2] << endl;
#endif
}
class vtkLineCallback : public vtkCommand
{
  public:
    static vtkLineCallback *New()
    {
      return new vtkLineCallback;
    }
 
    virtual void Execute(vtkObject *caller, unsigned long, void*)
    {
 
      vtkLineWidget2 *lineWidget = 
          reinterpret_cast<vtkLineWidget2*>(caller);
 
      // Get the actual box coordinates of the line
      vtkSmartPointer<vtkPolyData> polydata = 
          vtkSmartPointer<vtkPolyData>::New();
      static_cast<vtkLineRepresentation*>(lineWidget->GetRepresentation())->GetPolyData (polydata);
 
      // Display one of the points, just so we know it's working
      double p[3];
      polydata->GetPoint(0,p);
      std::cout << "P: " << p[0] << " " << p[1] << " " << p[2] << std::endl;
    }
    vtkLineCallback(){}
 
};
void vtk_operation::add_line_widget(QVTKWidget *qvtk_widget, const vtkSmartPointer<vtkRenderer> &render,
                                    std::vector<vtkSmartPointer<vtkLineWidget2> > &line_list)
{
  vtkSmartPointer<vtkLineWidget2> lineWidget = 
      vtkSmartPointer<vtkLineWidget2>::New();
  lineWidget->SetInteractor(qvtk_widget->GetInteractor());
  lineWidget->CreateDefaultRepresentation();
  vtkSmartPointer<vtkLineCallback> lineCallback = 
    vtkSmartPointer<vtkLineCallback>::New();
  cout << "------" << endl;
  lineWidget->AddObserver(vtkCommand::InteractionEvent,lineCallback);
  lineWidget->SetEnabled(1);

  lineWidget->On();
}

/**
 *@brief ouput tet index which being included by contral box 
 *
 *@param[in]
 *-select use to get included tet
 *-box control box widget
 *-box_boundary the boundary of control box widget
 *-fout the output stream
 *-tet the tet set of tet model
 *
 *@return void
 */

void output_tet(const vtkSmartPointer<vtkSelectEnclosedPoints> &select,
                const vtkSmartPointer<box_widget> &box,
                const vtkSmartPointer<vtkPolyData> &box_boundary,
                ofstream &fout, const matrixd &tet)
{
  matrixd x_vec, y_vec, z_vec;
  double tmp[6][3];
  vector<int> index;
  for(size_t i = 0; i < 6; ++i)
    box_boundary->GetPoint(i + 8, tmp[i]);
#if TEST_OUTPUT_TET
  for(size_t i = 0; i < 6; ++i)
  {
    for(size_t j = 0; j < 3; ++j)
      cout << tmp[i][j] << " " ;
    cout << endl;
  }
#endif
  x_vec = zeros<double>(1, 3);
  y_vec = zeros<double>(1, 3);
  z_vec = zeros<double>(1, 3);
  for(size_t i = 0; i < 3; ++i)
    x_vec(0, i) = tmp[1][i] - tmp[0][i];
  for(size_t i = 0; i < 3; ++i)
    y_vec(0, i) = tmp[3][i] - tmp[2][i];
  for(size_t i = 0; i < 3; ++i)
    z_vec(0, i) = tmp[5][i] - tmp[4][i];
  x_vec /= norm(x_vec);
  y_vec /= norm(y_vec);
  z_vec /= norm(z_vec);
  for(size_t i = 0; i < 3; ++i)
    fout << x_vec(0, i) << " ";
  for(size_t i = 0; i < 3; ++i)
    fout << y_vec(0, i) << " ";
  for(size_t i = 0; i < 3; ++i)
    fout << z_vec(0, i) << " ";
  fout << endl;
  for(size_t i = 0; i < tet.size(2); ++i)
  {
    if(select->IsInside(tet(0, i)) && select->IsInside(tet(1, i)) &&
       select->IsInside(tet(2, i)) && select->IsInside(tet(3, i)))
      index.push_back(i);
  }
  fout << index.size() << endl;
  for(size_t i = 0; i < index.size(); ++i)
    fout << index[i] << " ";
  fout << endl;
}


/**
 *@brief get tet index which included by control box widget
 *
 *@param[in]
 *-boxwidget control box widget set which include all exist boxs
 *-tet the tet set of tet model
 *
 *@return void
 */

void vtk_operation::get_tet_data(std::vector<vtkSmartPointer<box_widget> > &boxwidget,
                                 const matrixd &tet)
{
  ofstream fout;
  fout.open("tet_index.txt");
  if(fout.fail())
  {
    cerr << "fail to open tet_index.txt" << endl;
    cerr << "file: " << __FILE__ << endl;
    cerr << "line: " << __LINE__ << endl;
    return;
  }
  std::vector<vtkSmartPointer<box_widget> >::iterator box_it;
  fout << boxwidget.size() << endl;
  for(box_it = boxwidget.begin(); box_it != boxwidget.end(); ++box_it)
  {
    vtkSmartPointer<vtkSelectEnclosedPoints> select =
        vtkSmartPointer<vtkSelectEnclosedPoints>::New();
    vtkSmartPointer<vtkPolyData> box_boundary =
        vtkSmartPointer<vtkPolyData>::New();
    static_cast<vtkBoxRepresentation*>((*box_it)->GetRepresentation())->GetPolyData(box_boundary);
    vtkSmartPointer<vtkPolyData> pointsPolydata = 
        vtkSmartPointer<vtkPolyData>::New();
    pointsPolydata->SetPoints(points);
      
#if TEST_GET_TET
    cout << "line " << box_boundary->GetNumberOfLines() << endl;
    cout << "poly " << box_boundary->GetNumberOfPolys() << endl;
    cout << "verts " << box_boundary->GetNumberOfVerts() << endl;
    cout << "strip " << box_boundary->GetNumberOfStrips() << endl;
    cout << "cell " << box_boundary->GetNumberOfCells() << endl;
    cout << "points " << box_boundary->GetNumberOfPoints() << endl;
#endif
      
   
#if VTK_MAJOR_VERSION <= 5
    select->SetInput(pointsPolydata);
#else
    select->SetInputData(pointsPolydata);
#endif
#if VTK_MAJOR_VERSION <= 5
    select->SetSurface(box_boundary);
#else
    select->SetSurfaceData(box_boundary);
#endif
    select->Update();
    output_tet(select, (*box_it), box_boundary, fout, tet);
  }
  fout.close();

}

/**
 *@brief compute cross point of line and plane
 *
 *@param[in]
 *-vec_line the direction of line
 *-point_line a point on line
 *-vec_plane the normal vector of plane
 *-point_plane a point on plane
 *-flag indicate which point
 *-cross_point store cross point
 *
 *@return void
 */

inline void line_plane_cross_point(const matrixd &vec_line,
                                   const matrixd &point_line,
                                   const matrixd &vec_plane,
                                   const matrixd &point_plane,
                                   const size_t flag,
                                   matrixd &cross_point)
{
  double t;
  matrixd tmp;
  tmp = zeros(1, 3);
  tmp = point_plane - point_line;
  t = dot(vec_plane, tmp) / dot(vec_plane, vec_line);
  cross_point(flag, 0) = vec_line[0] * t + point_line[0];
  cross_point(flag, 1) = vec_line[1] * t + point_line[1];
  cross_point(flag, 2) = vec_line[2] * t + point_line[2];
}

bool is_same_direction(const matrixd &vec_list)
{
  size_t cnt = 0;
  matrixd sum_vec = zeros(1, 3);
  for(size_t i = 0; i < vec_list.size(1); ++i)
  {
    for(size_t j = 0; j < 3; ++j)
      sum_vec(0, j) += vec_list(i, j);
    if(is_zero(vec_list(i, 0)) && is_zero(vec_list(i, 1)) && is_zero(vec_list(i, 2)))
      cnt++;
  }
  matrixd sum = zeros(1, 3);
  sum = (vec_list.size(1) - cnt) * vec_list(0, colon());
  for(size_t i = 0; i < 3; ++i)
  {
    if(!equal(sum(0, i), sum_vec(0, i)))
      return false;
  }
  return true;
    
}

template <typename T>
inline bool is_between(const T &min, const T &max, const T &value)
{
  if(is_larger_equal(value, min) && is_less_equal(value, max))
    return true;
  else
    return false;
}

inline void copy_to_matrix(const double p[], matrixd &mat)
{
  mat(0, 0) = p[0];
  mat(0, 1) = p[1];
  mat(0, 2) = p[2];
  //cout << "mat point" << endl;
  // cout << mat << endl;
}

/**
 *@brief  judge a point whether inside the rectangular
 *
 *@param[in]
 *-p the point which will be judged
 *-rect_point: four vertex of the plane
 *
 *@return void
 */

inline bool is_inside_rect(matrixd &p, const double rect_point[4][3])
{
  matrixd tmp1 = zeros(1, 3);
  matrixd tmp2 = zeros(1, 3);
  double area_tri, area_rect;
  area_tri = 0.0;
  for(size_t i = 0; i < 4; ++i)
  {
    copy_to_matrix(rect_point[i], tmp1);
    copy_to_matrix(rect_point[(i + 1) % 4], tmp2);
    area_tri =  area_tri + tri_area(p, tmp1, tmp2);
  }
  matrixd tmp3 = zeros(1, 3);
  copy_to_matrix(rect_point[0], tmp1);
  copy_to_matrix(rect_point[1], tmp2);
  copy_to_matrix(rect_point[2], tmp3);
  area_rect = 2 * tri_area(tmp1, tmp2, tmp3);
  if(equal(area_tri, area_rect))
    return true;
  else
    return false;
}

/**
 *@brief output the tet index which crossed by plane
 *
 *@param[in]
 *-tet data set of tet
 *-vertex point set of tet model
 *-rect_point four vertex of the plane
 *-fout output stream
 *-widget control plane widget
 *
 *@return void
 */

void output_plane_tet(const matrixd &tet, const matrixd &vertex,
                      const double rect_point[4][3], ofstream &fout,
                      const vtkSmartPointer<plane_widget> &widget)
{
  double vec[3], point[3];
  vector<int> index;
  matrixd vec_line, vec_plane, point_line, point_plane;
  vec_line = zeros(1, 3);
  vec_plane = zeros(1, 3);
  point_line = zeros(1, 3);
  point_plane = zeros(1, 3);
  widget->GetCenter(point);
  for(size_t j = 0; j < 3; ++j)
    point_plane(0, j) = point[j];
  widget->GetNormal(vec);
   
  for(size_t j = 0; j < 3; ++j)
    vec_plane(0, j) = vec[j];
  for(size_t i = 0; i < 3; ++i)
    fout << vec_plane(0, i) << " " ;
  fout << endl;
  matrixd cross_point = zeros(4, 3);
  matrixd vec_list = zeros(4, 3);
  for(size_t i = 0; i < tet.size(2); ++i)
  {
    bool flag = true;
    for(size_t j = 0; j < 4; ++j)
    {
      point_line = trans(vertex(colon(), tet(j, i)));
      
      line_plane_cross_point(vec_plane, point_line,
                             vec_plane, point_plane, j,
                             cross_point);
      vec_list(j, colon()) = (point_line - cross_point(j, colon()))/ norm(point_line - cross_point(j, colon()));
      
    }
    if(is_same_direction(vec_list))
      continue;
    matrixd tmp = zeros(1, 3);
    for(size_t j = 0; j < 4; ++j)
    {
      tmp = cross_point(j, colon());
      if(!is_inside_rect(tmp, rect_point))
      {
        flag = false;
        break;
      }
    }
    if(flag)
      index.push_back(i);
  }
  fout << index.size() << endl;
  for(size_t i = 0; i < index.size(); ++i)
    fout << index[i] << " ";
  fout << endl;
}

/**
 *@brief get tet index which crossed by control plane
 *
 *@param[in]
 *-widget control plane widget list which include all exist int qvtk_widget now
 *-tet tet set of tet model
 *-vertex point set of tet model
 *
 *@return void
 */


int vtk_operation::get_plane_tet(const std::vector<vtkSmartPointer<plane_widget> > &widget,
                                  const matrixd &tet, const matrixd &vertex)
{
  double tmp[4][3];
  double bounding[6];
  ofstream fout;
  fout.open("plane_tet_index.txt");
  if(fout.fail())
  {
    cerr << "fail to open plane_tet_index.txt." << endl;
    cerr << "file:" << __FILE__ << endl;
    cerr << "line:" << __LINE__ << endl;
    return 1;
  }
  fout << widget.size() << endl;
  vtkSmartPointer<vtkPolyData> poly =
      vtkSmartPointer<vtkPolyData>::New();
  for(size_t i = 0; i < widget.size(); ++i)
  {
    widget[i]->GetPolyData(poly);
    for(size_t j = 0; j < 4; ++j)
      poly->GetPoint(j, tmp[j]);
    swap(tmp[2], tmp[3]);
    // get_bounding(tmp, bounding);
    output_plane_tet(tet, vertex, tmp, fout, widget[i]);
  }
  fout.close();
}

/**
 *@brief compute rotation matrix and save into mat
 *
 *@param[in]
 *-angle: the degree of rotate 
 *-mat: save rotation matrix
 *-vec: the axis which widget around it rotate
 *
 *@return void
 */

void set_mutil_matrix(double angle, matrixd &mat, const matrixd &vec)
{
  angle = angle * PI / 180;
  double cos_angle = cos(angle/2);
  double sin_angle = sin(angle/2);
  mat = zeros(3, 3);
  double e0 = cos_angle;
  double e1 = vec(0, 0) * sin_angle;
  double e2 = vec(0, 1) * sin_angle;
  double e3 = vec(0, 2) * sin_angle;
  mat(0, 0) = 1 - 2 * e2 * e2 - 2 * e3 * e3;
  mat(0, 1) = 2 * e1 * e2 + 2 *  e0 * e3;
  mat(0, 2) = 2 * e1 * e3 - 2 * e0 * e2;
  mat(1, 0) = 2 * e1 * e2 - 2 * e0 * e3;
  mat(1, 1) = 1 - 2 * e1 * e1 - 2 * e3 * e3;
  mat(1, 2) = 2 * e2 * e3 + 2 * e0 * e1;
  mat(2, 0) = 2 * e1 * e3 + 2 * e0 * e2;
  mat(2, 1) = 2 * e2 * e3 - 2 * e0 * e1;
  mat(2, 2) = 1 - 2 * e1 * e1 - 2 * e2 *e2;
#ifdef TEST_ROTATE_MATRIX
  double a1;
  a1 = mat[0, 0] * (mat(1, 1) * mat(2, 2) - mat(1, 2) * mat(2, 1));
  double a2;
  a2 = mat(0, 1) * (mat(1, 0) * mat(2, 2) - mat(1, 2) * mat(2, 0));
  double a3;
  a2 = mat(0, 2) * (mat(1, 0) * mat(2, 1) - mat(1, 1) * mat(2, 0));
  cout << "mat size " << a1 -a2 + a3 << endl;
#endif
}


/**
 *@brief rotate the control box widget
 *
 *@param[in]
 *-zoom: degree of rotate
 *   -1  rotate 10
 *   -2  rotate 1
 *-flag clockwise or counterclockwise
 *   -1  clockwise
 *   -2  counterclockwise
 *-axis: the axis which widget around it rotate
 *-box: control box widget which wil be rotated
 *
 *@return void
 */

void vtk_operation::rotate_box(const int zoom, const int flag,
                               const double axis[3], vtkSmartPointer<box_widget> &box)
{
  
  vtkSmartPointer<vtkBoxRepresentation> boxRepresentation = 
      vtkSmartPointer<vtkBoxRepresentation>::New();
  boxRepresentation = static_cast<vtkBoxRepresentation*>(box->GetRepresentation());
  vtkSmartPointer<vtkPolyData> box_boundary =
      vtkSmartPointer<vtkPolyData>::New();
  static_cast<vtkBoxRepresentation*>(box->GetRepresentation())->GetPolyData(box_boundary);
  vtkSmartPointer<vtkTransform> transform =
      vtkSmartPointer<vtkTransform>::New();
  boxRepresentation->GetTransform(transform);
  vtkSmartPointer<vtkMatrix4x4> pre_mat =
      vtkSmartPointer<vtkMatrix4x4>::New();
  double angle;
  if(zoom == 1)
  {
    if(flag == 1)
      angle = pos_angle;
    else if(flag == 2)
      angle = pos_min_angle;
  }
  else if(zoom == 2)
  {
    if(flag == 1)
      angle = neg_angle;
    else if(flag == 2)
      angle = neg_min_angle;
  }
  matrixd rotate_mat, vec;
  double neg_x[3], pos_x[3];
  if(axis[0] == 1)
  {
    box_boundary->GetPoint(8, neg_x);
    box_boundary->GetPoint(9, pos_x);
  }
  else if(axis[1] == 1)
  {
    box_boundary->GetPoint(10, neg_x);
    box_boundary->GetPoint(11, pos_x);
  }
  else if(axis[2] == 1)
  {
    box_boundary->GetPoint(12, neg_x);
    box_boundary->GetPoint(13, pos_x);
  }
  //set rotation axis
  vec = zeros(1, 3);
  vec(0, 0) = pos_x[0] - neg_x[0];
  vec(0, 1) = pos_x[1] - neg_x[1];
  vec(0, 2) = pos_x[2] - neg_x[2];
  vec /= norm(vec);
  //set totation matrix
  set_mutil_matrix(angle, rotate_mat, vec);
  matrixd mat = zeros(3, 3);
  //get the matrix 4x4 of box
  transform->GetMatrix(pre_mat);
  for(size_t i = 0; i < 3; ++i)
    for(size_t j = 0; j < 3; ++j)
      mat(i, j) = pre_mat->Element[i][j];
  //rotation
  mat =  temp(rotate_mat * mat);
  for(size_t i = 0; i < 3; ++i)
    for(size_t j = 0; j < 3; ++j)
      pre_mat->Element[i][j] = mat(i, j);
  transform->SetMatrix(pre_mat);
  boxRepresentation->SetTransform(transform);
  box->SetRepresentation(boxRepresentation);
}


void vtk_operation::output_box_position(const std::vector<vtkSmartPointer<box_widget> > &boxwidget,
                                        ofstream &fout)
{
  vtkSmartPointer<vtkBoxRepresentation> representation = 
      vtkSmartPointer<vtkBoxRepresentation>::New();
  double *bound;
  fout << boxwidget.size() << endl;
  for(size_t i = 0; i < boxwidget.size(); ++i)
  {
    representation = static_cast<vtkBoxRepresentation*>(boxwidget[i]->GetRepresentation());
    bound = representation->GetBounds();
    for(size_t j = 0; j < 6; ++j)
      fout << bound[j] << endl;
    vtkSmartPointer<vtkTransform> transform =
        vtkSmartPointer<vtkTransform>::New();
    vtkSmartPointer<vtkMatrix4x4> pre_mat =
        vtkSmartPointer<vtkMatrix4x4>::New();
    representation->GetTransform(transform);
    transform->GetMatrix(pre_mat);
    for(size_t j = 0; j < 4; ++j)
      for(size_t k = 0; k < 4; ++k)
        fout << pre_mat->Element[j][k] << " ";
    fout << endl;
  }
  
}

void vtk_operation::import_box_position(ifstream &fin, QVTKWidget * qvtk_widget,
                                        std::vector<vtkSmartPointer<box_widget> > &boxwidget,
                                        const vtkSmartPointer<box_callback> &callback)
{
  double bound[6];
  size_t num;
  fin >> num;
  cout << "num " << num << endl;
  vtkSmartPointer<box_widget> box = vtkSmartPointer<box_widget>::New();
  for(size_t i = 0; i < num; ++i)
  {
    for(size_t j = 0; j < 6; ++j)
    {
      fin >> bound[j];
      //bound[j] *= 2;
      cout << "bound " << bound[j] << endl;
    }
   
    vtkSmartPointer<vtkBoxRepresentation> representation = 
        vtkSmartPointer<vtkBoxRepresentation>::New();
    vtkSmartPointer<vtkTransform> transform =
        vtkSmartPointer<vtkTransform>::New();
    vtkSmartPointer<vtkMatrix4x4> pre_mat =
        vtkSmartPointer<vtkMatrix4x4>::New();
    representation->GetTransform(transform);
    transform->GetMatrix(pre_mat); 
    for(size_t j = 0; j < 4; ++j)
      for(size_t k = 0; k < 4; ++k)
        fin >> pre_mat->Element[j][k];
    transform->SetMatrix(pre_mat);
   
    representation->SetTransform(transform); 
    representation->PlaceWidget(bound);
    box->SetRepresentation(representation);
    box->SetInteractor(qvtk_widget->GetInteractor());
    boxwidget.push_back(box);
    box->AddObserver(vtkCommand::EndInteractionEvent, callback);
    box->On();
  }
 
}

void vtk_operation::add_axes(vtkSmartPointer<vtkRenderer> &render)
{
   vtkSmartPointer<vtkTransform> transform =
    vtkSmartPointer<vtkTransform>::New();
  transform->Translate(1.0, 0.0, 0.0);
 
  vtkSmartPointer<vtkAxesActor> axes =
    vtkSmartPointer<vtkAxesActor>::New();
 
  // The axes are positioned with a user transform
  axes->SetUserTransform(transform);
  render->AddActor(axes);
}

/**
 *@brief this function add an area picker
 *
 */
void vtk_operation::add_area_picker(QVTKWidget *qvtk_widget, const size_t flag_win,
                                    vtkSmartPointer<vtkRenderer> &render)
{
  //  vtkSmartPointer<vtkIdFilter> id_filter =
  //     vtkSmartPointer<vtkIdFilter>::New();
  // id_filter->SetInputConnection(unstructured->GetProducerPort());
  // id_filter->SetIdsArrayName("OriginalIds");
  // id_filter->Update();
  // vtkSmartPointer<vtkDataSetSurfaceFilter> surface_filter =
  //     vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
  // surface_filter->SetInputConnection(id_filter->GetOutputPort());
  // surface_filter->PassThroughCellIdsOff();
  // surface_filter->Update();
  // CREATE_NEW(vtkPolyData, input);
  // input = surface_filter->GetOutput();
  
  CREATE_NEW(vtkInteractorStyleRubberBandPick, interactor_pick);
  // it must be use vtkRenderAreaPicker
  CREATE_NEW(vtkRenderedAreaPicker, area_picker);
  qvtk_widget->GetRenderWindow()->GetInteractor()->SetInteractorStyle(interactor_pick);
  qvtk_widget->GetRenderWindow()->GetInteractor()->SetPicker(area_picker);
  if(flag_win == LEFT)
  {
    call_back_left->set_render(render);
    call_back_left->set_pick_type(1);
    //call_back_left->set_data(unstructured);
    //call_back_left->set_actor(select_actor);

    //call_back_left->set_poly(surface_filter->GetOutput());
    call_back_left->set_poly(vtk_poly_left);
    call_back_left->set_cell_index(cell_index);
    qvtk_widget->GetRenderWindow()->GetInteractor()->AddObserver(vtkCommand::EndPickEvent,
                                                                 call_back_left);
  }
  else if(flag_win == RIGHT)
  {
     call_back_right->set_render(render);
     call_back_right->set_pick_type(1);
     //call_back_right->set_data(unstructured);
     // call_back_right->set_actor(select_actor);
     call_back_right->set_poly(vtk_poly_right);
     call_back_right->set_cell_index(cell_index);
     qvtk_widget->GetRenderWindow()->GetInteractor()->AddObserver(vtkCommand::EndPickEvent,
                                                                  call_back_right);
  }
}


/**
 *@brief change the cells color of area picker selected. 
 *
 *@param[in]
 *-flag: the color type that want to convert
 *   -0 blue
 *   -1 white
 *   -2 red
 *
 *-render: 
 *
 *@return void
 */

void vtk_operation::change_color(const size_t flag, vtkSmartPointer<vtkRenderer> &render,
                                 const size_t win_flag)
{
 
  CREATE_NEW(vtkDataSetMapper, tmp_mapper);
  cell_index.clear();
  size_t num ;
  if(win_flag == LEFT)
    num = call_back_left->get_mapper_num();
  else
    num = call_back_right->get_mapper_num();
  for(size_t i = 0; i < num; ++i)
  {
    if(win_flag == LEFT)
    {
      call_back_left->get_cell_index(cell_index);
    
      if(cell_index.size() == 0)
        return;
      // call_back_left->get_mapper(tmp_mapper);
      call_back_left->get_mapper(i, tmp_mapper);
    }
    else if(win_flag == RIGHT)
    {
      call_back_right->get_cell_index(cell_index);
      if(cell_index.size() == 0)
        return ;
      //call_back_right->get_mapper(tmp_mapper);
      call_back_right->get_mapper(i, tmp_mapper);
    }
    //  call_back->get_mapper(tmp_mapper);
    CREATE_NEW(vtkActor, tmp_actor);
    tmp_actor->SetMapper(tmp_mapper);
    if(flag == RED)
      tmp_actor->GetProperty()->SetColor(1, 0, 0);
    else if(flag == BLUE)
      tmp_actor->GetProperty()->SetColor(0, 0, 1);
    else if(flag == WHITE)
      tmp_actor->GetProperty()->SetColor(1, 1, 1);
    render->AddActor(tmp_actor);
    if(win_flag == LEFT)
      vtk_file_op_left.reset_color(cell_index, flag);
    else if(win_flag == RIGHT)
      vtk_file_op_right.reset_color(cell_index, flag);
  }
}

void surface_edge_mesh(size_t &cnt, vtkSmartPointer<vtkActor> &actor,
                       double surface_edge_color[])
{
  if(cnt & 1)
    actor->GetProperty()->EdgeVisibilityOff();
  else
  {
    double width;
    width = configure.get_surface_edge_width();
    actor->GetProperty()->EdgeVisibilityOn();
    actor->GetProperty()->SetEdgeColor(surface_edge_color);
    actor->GetProperty()->SetLineWidth(width);
  }
  cnt++;
}
/**
 *@brief display surface with edge mesh
 *
 *@param[in]
 *-void
 *
 *@return void
 */
void vtk_operation::display_surface_edge_mesh(const size_t flag_win)
{
  double surface_edge_color[3];
  configure.get_surface_edge_color(surface_edge_color);
  if(flag_win == LEFT)
    surface_edge_mesh(cnt_surface_edge_left, actor_left,
                      surface_edge_color);
  else if(flag_win == RIGHT)
    surface_edge_mesh(cnt_surface_edge_right, actor_right,
                      surface_edge_color);
}

void vtk_operation::output_cell_data(const size_t flag)
{
  if(flag == LEFT)
    vtk_file_op_left.output_cell();
  else if(flag == RIGHT)
    vtk_file_op_right.output_cell();
}

void vtk_operation::switch_pick_type(size_t flag, size_t win_flag)
{
  if(win_flag == LEFT)
    call_back_left->set_pick_type(flag);
  else if(win_flag == RIGHT)
    call_back_right->set_pick_type(flag);
  else
  {
    cerr << "win_flag is error, error happen in file : " << endl;
    cerr <<  __FILE__ << " and in line  " << __LINE__ << endl;
    return ;
  }
}

}
