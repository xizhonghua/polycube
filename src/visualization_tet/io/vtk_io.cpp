#include "vtk_io.h"
#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <string.h>
#include <vtkTriangle.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>

#include <vtkCellData.h>

#include <vtkUnstructuredGridReader.h>
namespace lq {

using namespace  std;
size_t mid_type;

vtk_io::vtk_io()
{

  
}

vtk_io::~vtk_io()

{

  
}

string vtk_io::get_data_type()
{
  return this->data_type;
}

string vtk_io::get_data_set()
{
  return this->data_set;
}

string vtk_io::get_point_type()
{
  return this->point_type;
}

size_t vtk_io::get_cell_type()
{
  return this->cell_type;
}

void vtk_io::set_cell_index(std::vector<size_t> &index)
{
  this->cell_index = index;
}

void vtk_io::set_configure(const load_configure &config)
{
  this->configure = config;
}

void vtk_io::clear()
{
  cell_index.clear();
  cell_look_table.clear();
  vertex_list.clear();
  triangle_list.clear();
  point_cell_list.clear();
  file_name.clear();
}

int vtk_io::read_cell(string &input, const size_t index, ifstream &fin)
{
  string tmp;
  tmp = input.substr(index + 1);
  size_t num = 0;
  size_t sum = 0;
  size_t tmp_num;
  sscanf(tmp.c_str(), "%d %d", &num, &sum);
  string sub_cell;
  size_t cell_num;
  if(num == 0 )
    return ERROR;
  cell_num =sum / num;
  cell_type = cell_num - 1;
  triangle tmp_triangle;
  size_t index1, index2, index3;
  vtkSmartPointer<vtkCellArray> vertices =
      vtkSmartPointer<vtkCellArray>::New();
  for(size_t i = 0; i < num; ++i)
  {
    getline(fin, sub_cell);
    if(cell_type == POINT)
      return POINT;
    else if(cell_type == LINE)
      return LINE;
    else if(cell_type == TRIANGLE)
    {
      sscanf(sub_cell.c_str(), "%d %d %d %d", &tmp_num, &index1, &index2, &index3);
      tmp_triangle.vertex.push_back(index1);
      tmp_triangle.vertex.push_back(index2);
      tmp_triangle.vertex.push_back(index3);
      triangle_list.push_back(tmp_triangle);
      tmp_triangle.vertex.clear();
    }
  }
 
}

void set_default_color(const size_t cell_num,
                       vtkSmartPointer<vtkPolyData> &poly,
                       vtkSmartPointer<vtkUnsignedCharArray> &cell_color)
{
  cell_color->Initialize();
  // cell_color->Reset();
  cell_color->SetNumberOfComponents(3);
  cell_color->SetName("cell_color");
  for(size_t i = 0; i < cell_num; ++i)
    cell_color->InsertNextTupleValue(white);
  poly->GetCellData()->SetScalars(cell_color);
}

void linear_inter(unsigned char color[], const size_t max,
                  const size_t min, const size_t type,
                  const size_t index_min, const size_t index_max)
{
  double coeff1, coeff2;
  coeff1 = (max - type) / (max - min);
  coeff2 = (type - min) / (max - min);
  for(size_t i = 0; i < 3; ++i)
  {
    color[i] += static_cast<unsigned char>(coeff1 * (static_cast<int>(COLOR[index_min][i])));
    color[i] += static_cast<unsigned char>(coeff2 * (static_cast<int>(COLOR[index_max][i])));
  }
}

void vtk_io::color_map(const std::vector<size_t> &cell_look_table,
                       vtkSmartPointer<vtkPolyData> &poly,
                       vtkSmartPointer<vtkUnsignedCharArray> &cell_color)
{
  cell_color->Reset();
  cell_color->SetNumberOfComponents(3);
  cell_color->SetName("cell_color");
  size_t type;
  for(size_t i = 0; i < cell_look_table.size(); ++i)
  {
    type = cell_look_table[i];
    unsigned char color[3] = {0 , 0, 0};
    if(type == mid_type)
      cell_color->InsertNextTupleValue(white);
    else if( type >= min_type && type < mid_type)
    {
      linear_inter(color, mid_type, min_type, type, BLUE, WHITE);
      cell_color->InsertNextTupleValue(color);
    }
    else if( type > mid_type && type <= max_type)
    {
      linear_inter(color, max_type, mid_type, type, WHITE, RED);
      cell_color->InsertNextTupleValue(color);
    }
  }
  poly->GetCellData()->SetScalars(cell_color);
}

void set_color(std::vector<size_t> &cell_look_table,
               const size_t cell_num,
               vtkSmartPointer<vtkPolyData> &poly,
               vtkSmartPointer<vtkUnsignedCharArray> &cell_color)
{
  //set_default_color(cell_num, poly);
  //cell_color->Initialize();
  cell_color->Reset();
  cell_color->SetNumberOfComponents(3);
  cell_color->SetName("cell_color");
  for(size_t i = 0; i < cell_num; ++i)
  {
    if(cell_look_table[i] == RED)
      cell_color->InsertNextTupleValue(red);
    else if(cell_look_table[i] == BLUE)
      cell_color->InsertNextTupleValue(blue);
    else if(cell_look_table[i] == WHITE)
      cell_color->InsertNextTupleValue(white);
  }
  poly->GetCellData()->SetScalars(cell_color);
}

void update_color(const std::vector<size_t> &cell_index,
                  vtkSmartPointer<vtkPolyData> &poly,
                  const size_t flag,
                  vtkSmartPointer<vtkUnsignedCharArray> &cell_color)
{
  size_t index;
  for(size_t i = 0; i < cell_index.size(); ++i)
  {
    index = cell_index[i];
    if(flag == RED)
      cell_color->SetTupleValue(index, red);
    else if(flag == BLUE)
      cell_color->SetTupleValue(index, blue);
    else if(flag == WHITE)
      cell_color->SetTupleValue(index, white);
  }
  poly->GetCellData()->SetScalars(cell_color);
}

void vtk_io::read_cell_data(ifstream &fin, const size_t index,
                            string &input)
{
  
  size_t num = 0;
  size_t space_index;
  int min = 256;
  int max = -1;
  num = atoi((input.substr(index + 1)).c_str());
  string look_table, flag_str;
  getline(fin, look_table);
  space_index = look_table.find(" ");
  flag_str = look_table.substr(0, space_index);
  if(strcmp(flag_str.c_str(), "SCALARS") == 0)
  {
    getline(fin, look_table);
    look_table.clear();
    for(size_t i = 0; i < num; ++i)
    {
      getline(fin, look_table);
      int tmp;
      tmp = atoi(look_table.c_str());
      if(tmp < min)
        min = tmp;
      else if (tmp > max)
        max = tmp;
      cell_look_table.push_back(tmp);
      look_table.clear();
    }
    min_type = min;
    max_type = max;
    mid_type = (min + max) /2;
  }
  // else if(flag_str.c_str(), "")
}


int vtk_io::load_vtk(const char *path, vtkSmartPointer<vtkPolyData> &poly,
                     vtkSmartPointer<vtkPolyDataMapper> &mapper,
                     vtkSmartPointer<vtkDataSetMapper> &line_mapper,
                     vtkSmartPointer<vtkDataSetMapper> &point_mapper,
                     vtkSmartPointer<vtkUnsignedCharArray> &cell_color)
{
  triangle_list.clear();
  vertex_list.clear();
  cell_look_table.clear();
  file_name = path;
  cout << "file_name " << file_name.c_str() << endl;
  ifstream fin;
  fin.open(path);
  if(fin.fail())
  {
    cerr << "Fail to open : " << path << endl;
    cerr << "File : " << __FILE__ << endl;
    cerr << "Line :" << __LINE__ << endl;
    return 1;
  }
  string input, sub_str, tmp;
  char type[10];
  while(getline(fin, input))
  {
    if(input[0] == '#')
    {
      input.clear();
      continue;
    }
    else if(strcmp(input.c_str(), "ASCII") == 0 ||
            strcmp(input.c_str(), "BINARY") == 0)
    {
      data_type = input;
      input.clear();
      continue;
    }
    else
    {
      size_t flag;
      flag = input.find(" ");
      sub_str = input.substr(0, flag);
      if(strcmp(sub_str.c_str(), "DATASET") == 0)
      {
        data_set = input.substr(flag);
        input.clear();
        continue;
      }
      else if(strcmp(sub_str.c_str(), "POINTS") == 0)
      {
        tmp = input.substr(flag + 1);
        size_t num = 0;
        sscanf(tmp.c_str(), "%d %s", &num, &type);
       
        string point_str;
        double x, y ,z;
        point tmp_point;
        for(size_t i = 0; i < num; ++i)
        {
          getline(fin, point_str);
          sscanf(point_str.c_str(), "%lf %lf %lf", &tmp_point.x, &tmp_point.y, &tmp_point.z);
         
          vertex_list.push_back(tmp_point);
          point_str.clear();
        }
      }
      else if(strcmp(sub_str.c_str(), "CELLS") == 0)
      {
        vtkSmartPointer<vtkUnstructuredGridReader> reader =
            vtkSmartPointer<vtkUnstructuredGridReader>::New();
        if(read_cell(input, flag, fin) == POINT)
        {
          reader->SetFileName(path);
          point_mapper->SetInputConnection(reader->GetOutputPort());
          return POINT;
        }
        else  if(read_cell(input, flag, fin)  == LINE)
        {
          reader->SetFileName(path);
          line_mapper->SetInputConnection(reader->GetOutputPort());
          return LINE;
        }
      }
      else if(strcmp(sub_str.c_str(), "CELL_TYPES") == 0)
      {
        size_t num = 0;
        num = atoi((input.substr(flag + 1)).c_str());
        string sub_cell;
        for(size_t i = 0; i < num; ++i)
        {
          getline(fin, sub_cell);
          sub_cell.clear();
        }
      }
      else if(strcmp(sub_str.c_str(), "CELL_DATA") == 0)
        read_cell_data(fin, flag, input);
      input.clear();
      sub_str.clear();
      tmp.clear();
    }
  }
  vtkSmartPointer<vtkPoints> points =
      vtkSmartPointer<vtkPoints>::New();
  for(size_t i =0 ; i < vertex_list.size(); ++i)
    points->InsertNextPoint(vertex_list[i].x, vertex_list[i].y,
                            vertex_list[i].z);
  if(cell_type == TRIANGLE)
  {
    vtkSmartPointer<vtkTriangle> tmp_triangle =
        vtkSmartPointer<vtkTriangle>::New();
    vtkSmartPointer<vtkCellArray> triangle_cell =
        vtkSmartPointer<vtkCellArray>::New();
    for(size_t i = 0; i < triangle_list.size(); ++i)
    {
      tmp_triangle->GetPointIds()->SetId(0, triangle_list[i].vertex[0]);
      tmp_triangle->GetPointIds()->SetId(1, triangle_list[i].vertex[1]);
      tmp_triangle->GetPointIds()->SetId(2, triangle_list[i].vertex[2]);
      triangle_cell->InsertNextCell(tmp_triangle);
    }
    poly->SetPoints(points);
    poly->SetPolys(triangle_cell);
    //set_color(cell_look_table, triangle_list.size(), poly, cell_color);
    color_map(cell_look_table, poly, cell_color);
    //set_default_color(triangle_list.size(), poly);
    mapper->SetInputConnection(poly->GetProducerPort());
    return TRIANGLE;
  }
}

void vtk_io::reset_color(const std::vector<size_t> &cell_index,
                         const size_t flag)
{
  std::cout << "[#info]Reset cell edge color." << std::endl;
  std::cout << "[#info]cell_index size : " << cell_index.size() << std::endl;
  std::cout << "[#info]the color after fix : " << flag << std::endl;
  for(size_t i = 0; i < cell_index.size(); ++i)
  {
    cell_look_table[cell_index[i]] = flag;
    //std::cout << "index : " << cell_index[i] << std::endl;
  }
}

void vtk_io::output_cell()
{
  
  size_t flag;
  flag = file_name.find_last_of("/");
  string sub_str;
  sub_str = file_name.substr(flag + 1);
  flag = sub_str.find_last_of(".");
  if(flag != string::npos)
    sub_str.insert(flag, "-lq");
  cout << "out put file " << sub_str << endl;
  ofstream fout;
  fout.open(sub_str.c_str());
  size_t index1, index2, index3;
  cout << "cell_look_table size" << cell_look_table.size() << endl;
  for(size_t i = 0; i < cell_look_table.size(); ++i)
  {
    index1 = triangle_list[i].vertex[0];
    index2 = triangle_list[i].vertex[1];
    index3 = triangle_list[i].vertex[2];
    fout << index1 << " " << index2 << " " << index3 << " " << cell_look_table[i] << endl;
  }
  fout.close();
}

}







