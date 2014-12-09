#include "area_select_callback.h"
#include <vtkHardwareSelector.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkExtractSelectedPolyDataIds.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkDataSetMapper.h>
#include <vtkCellType.h>
#include <vtkCellArray.h>
#include <vtkProperty.h>
#include <vtkExtractSelection.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkRendererCollection.h>
#include <vtkIdList.h>
#include <vtkIdFilter.h>
#include <vtkIdTypeArray.h>
#include <vtkCellData.h>
#include <string.h>
#include <fstream>
#include <map>

#define CREATE_NEW(class, variable)                              \
  vtkSmartPointer<class> variable = vtkSmartPointer<class>::New();
int num = 0;
namespace lq {

  void area_select_callback::get_cell_list(std::vector<triangle> &cell)
  {
    cell = this->cell_list;
  }

  void area_select_callback::get_cell_index(std::vector<size_t> &index)
  {
    index = this->cell_index;
  }

  void area_select_callback::get_mapper(vtkSmartPointer<vtkDataSetMapper> &mapper)
  {
    mapper = this->select_mapper;
  }

  size_t area_select_callback::get_mapper_num()
  {
    return select_mapper_list.size();
  }

  void area_select_callback::get_mapper(size_t index, vtkSmartPointer<vtkDataSetMapper> &mapper)
  {
    if(index >= select_mapper_list.size())
      {
        cerr << "the index is over the number of mapper list." << endl;
        cerr << "error happen in file : " << __FILE__ << endl;
        cerr << "Line : " << __LINE__ << endl;
        return ;
      }
    mapper = select_mapper_list[index];
  }

  void area_select_callback::set_actor(vtkSmartPointer<vtkActor> &actor)
  {
    this->select_actor = actor;
  }

  void area_select_callback::set_render(vtkSmartPointer<vtkRenderer> &render)
  {
    this->select_render = render;

  }

  void area_select_callback::set_data(vtkSmartPointer<vtkUnstructuredGrid> &data)
  {
    this->unstructured = data;
  }

  void area_select_callback::set_sphere(vtkSmartPointer<vtkSphereSource> &s)
  {
    this->sphere = s;
  }

  void area_select_callback::set_poly(vtkSmartPointer<vtkPolyData> &poly)
  {
    this->select_poly = poly;
  }

  void area_select_callback::set_mapper(vtkSmartPointer<vtkDataSetMapper> &mapper)
  {
    this->select_mapper = mapper;
  }

  void area_select_callback::set_surface_filter(const vtkSmartPointer<vtkDataSetSurfaceFilter> &filter)
  {
    this->surface_filter = filter;
  }
  void area_select_callback::set_cell_index(std::vector<size_t> &index)
  {
    this->cell_index = index;
  }

  void area_select_callback::set_pick_type(const size_t flag)
  {
    this->pick_type = flag;
  }


  int area_select_callback::get_cell_index(const std::vector<triangle> &tri_list)
  {
    std::map<size_t, size_t> filter;
    std::map<size_t, size_t>::iterator itor;
    for(size_t i = 0; i < tri_list.size(); ++i)
      for(size_t j = 0; j < tri_list[i].vertex.size(); ++j)
        {
          if(filter.find(tri_list[i].vertex[j]) == filter.end())
            filter[tri_list[i].vertex[j]] = 1;
          else
            filter[tri_list[i].vertex[j]]++;
        }
    for(itor = filter.begin(); itor != filter.end(); ++itor)
      {
        if((*itor).second == tri_list.size())
          return (*itor).first;
      }
  }

  vtkSmartPointer<vtkPolyData> emptyPD = vtkSmartPointer<vtkPolyData>::New();
  void area_select_callback::Execute(vtkObject *obj, unsigned long, void *cell_data)
  {
    // CREATE_NEW(vtkRenderWindowInteractor, itor);
    // CREATE_NEW(vtkRenderWindow, render_win);
    // itor = reinterpret_cast<vtkRenderWindowInteractor*>(obj);
    // render_win = itor->GetRenderWindow ();
    std::cout << "execute pick" << std::endl;
    CREATE_NEW(vtkDataSetMapper, sMap);

    vtkSmartPointer<vtkIdFilter> idFilter =
        vtkSmartPointer<vtkIdFilter>::New();
    idFilter->SetInput(select_poly);
    idFilter->SetIdsArrayName("OriginalIds");
    idFilter->Update();

    vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter =
        vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surfaceFilter->SetInputConnection(idFilter->GetOutputPort());
    surfaceFilter->PassThroughCellIdsOn();
    surfaceFilter->Update();

    select_mapper = sMap;
    vtkSmartPointer<vtkHardwareSelector> sel =
        vtkSmartPointer<vtkHardwareSelector>::New();
    //sel->SetRenderer(render_win->GetRenderers()->GetFirstRenderer());
    sel->SetRenderer(select_render);
    double x0 = select_render->GetPickX1();
    double y0 = select_render->GetPickY1();
    double x1 = select_render->GetPickX2();
    double y1 = select_render->GetPickY2();
    sel->SetArea(static_cast<int>(x0),static_cast<int>(y0),static_cast<int>(x1),
                 static_cast<int>(y1));
    sel->SetFieldAssociation(vtkDataObject::FIELD_ASSOCIATION_CELLS);
    vtkSmartPointer<vtkSelection> res =
        vtkSmartPointer<vtkSelection>::New();
    res = sel->Select();
    if (!res)
      {
        cerr << "Selection not supported." << endl;
        return;
      }
    cout << "The number of selected nodes: " << res->GetNumberOfNodes() << endl;
    CREATE_NEW(vtkSelectionNode, cellids);

    cellids = res->GetNode(0);
    CREATE_NEW(vtkExtractSelectedPolyDataIds, extr);
    if (cellids)
      {
        cellids->SetFieldType(vtkSelectionNode::CELL);
        cellids->SetContentType(vtkSelectionNode::INDICES);
        //extr->SetInput(0, select_poly);
        extr->SetInput(0, surfaceFilter->GetOutput());
        vtkSmartPointer<vtkSelection> temp=
            vtkSmartPointer<vtkSelection>::New();
        temp->AddNode(cellids);
        extr->SetInput(1, temp);
        extr->Update();
        sMap->SetInput(extr->GetOutput());
        sMap->ScalarVisibilityOff();
        vtkSmartPointer<vtkDataSetSurfaceFilter> surface_filter =
            vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
        surface_filter->SetInputConnection(extr->GetOutputPort());
        surface_filter->PassThroughCellIdsOn();
        surface_filter->Update();
        vtkPolyData* selected = surface_filter->GetOutput();
        std::cout << "Selected " << selected->GetNumberOfPoints() << " points." << std::endl;
        std::cout << "Selected " << selected->GetNumberOfCells() << " cells." << std::endl;

        vtkIdTypeArray* ids = vtkIdTypeArray::SafeDownCast(selected->GetCellData()->GetArray("OriginalIds"));
        cell_index.clear();
        for(vtkIdType i = 0; i < ids->GetNumberOfTuples(); i++)
          {
            // std::cout << "Id " << i << " : " << ids->GetValue(i) << std::endl;
            cell_index.push_back(ids->GetValue(i));
          }
        //    vtkSmartPointer<vtkCell> tmp_cell;
        //    triangle tmp_tri;
        //    //select_poly->GetPointCells(1, index_list);
        //    cell_index.clear();
        //    for(size_t i = 0; i < selected->GetNumberOfCells(); ++i)
        //    {
        //      tmp_cell = selected->GetCell(i);
        //      vtkSmartPointer<vtkIdList>  index_list =
        //          vtkSmartPointer<vtkIdList>::New();
        //      for(size_t j = 0; j < tmp_cell->GetNumberOfPoints(); ++j)
        //      {
        //        //cout << "point index " << tmp_cell->GetPointId(j) << endl;
        //        select_poly->GetPointCells(tmp_cell->GetPointId(j), index_list);
        //        //extr->GetOutput()->GetPointCells(tmp_cell->GetPointId(j), index_list);
        //        //search_point_cell(tmp_cell->GetPointId(j), index_list);
        //        for(size_t k = 0; k < index_list->GetNumberOfIds(); ++k)
        //          tmp_tri.vertex.push_back(index_list->GetId(k));
        //        cell_list.push_back(tmp_tri);
        //        tmp_tri.vertex.clear();
        //      }

        //      cell_index.push_back(get_cell_index(cell_list));
        //      cell_list.clear();
        //    }
        num++;
        cout << "pick_type : " << pick_type << endl;
        /*< multi pick function, it allow user pick area and the area picked prvious is also
     display */
        if(pick_type == 2)
          {
            fstream fout;
            fout.open("multi_pick_surface_lq.txt", fstream::out | fstream::app);
            for(size_t i = 0; i < cell_index.size(); ++i)
              fout << cell_index[i] << endl;
            select_mapper_list.push_back(sMap);
            CREATE_NEW(vtkActor, multi_actor);
            multi_actor->SetMapper(sMap);
            multi_actor->GetProperty()->EdgeVisibilityOn();
            multi_actor->GetProperty()->SetEdgeColor(0, 0, 1);
            multi_actor->GetProperty()->SetLineWidth(3);
            select_render->AddActor(multi_actor);
          }
        else if(pick_type == 1)/*> single pick function, it allow user pick area and
                            the area picked previous is not exit*/
          {
            select_mapper_list.clear();
            select_mapper_list.push_back(sMap);
            select_actor->SetMapper(sMap);
            select_actor->GetProperty()->EdgeVisibilityOn();
            select_actor->GetProperty()->SetEdgeColor(0, 0, 1);
            select_actor->GetProperty()->SetLineWidth(3);
            // if(num > 1)
            // select_render->RemoveActor(select_actor);
            // select_render->ResetCamera();
            select_render->AddActor(select_actor);
          }
      }
    else
      {
        cerr << "Empty color buffer selection!" << endl;
        sMap->SetInput(emptyPD);
      }

  }

  void area_select_callback::search_point_cell(size_t point_index,
                                               vtkSmartPointer<vtkIdList> index_list)
  {
//    vtkSmartPointer<vtkCell> tmp_cell(vtkSmartPointer<vtkCell>::New());
    for(size_t i = 0; i < select_poly->GetNumberOfCells(); ++i)
      {
        vtkSmartPointer<vtkCell> tmp_cell(select_poly->GetCell(i));
        for(size_t j = 0; j < tmp_cell->GetNumberOfPoints(); ++j)
          {
            if(tmp_cell->GetPointId(j) == point_index)
              index_list->InsertNextId(i);
          }
      }
  }

}
