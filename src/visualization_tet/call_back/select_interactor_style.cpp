#include "select_interactor_style.h"
#include <vtkPlanes.h>
#include <vtkExtractPolyDataGeometry.h>
#include <vtkIdTypeArray.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPointData.h>
#include <vtkRendererCollection.h>
#include <vtkProperty.h>
#include <vtkExtractGeometry.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetSurfaceFilter.h>

namespace lq {

void select_interactor_style::set_poly_data(const vtkSmartPointer<vtkPolyData> &poly)
{
  this->selected_poly = poly;
}

// select_interactor_style::select_interactor_style()
// {
//   this->selected_poly = vtkSmartPointer<vtkPolyData>::New();
//   this->selected_actor = vtkSmartPointer<vtkActor>::New();
//   this->selected_mapper = vtkSmartPointer<vtkDataSetMapper>::New();
//   this->selected_actor->SetMapper(this->selected_mapper);
// }

void select_interactor_style::OnLeftButtonUp()
{
  // Forward events
  vtkInteractorStyleRubberBandPick::OnLeftButtonUp();
 
  if(this->CurrentMode == VTKISRBP_SELECT)
  {
    vtkPlanes* frustum = static_cast<vtkAreaPicker*>(this->GetInteractor()->GetPicker())->GetFrustum();
    cout << "plane number " << frustum->GetNumberOfPlanes() << endl;
    // vtkSmartPointer<vtkPoints>  tmp = vtkSmartPointer<vtkPoints>::New();
    // tmp->Initialize();
    // //tmp = frustum->GetPoints();
    // frustum->SetPoints(tmp);
    // cout << "points number " << tmp->GetNumberOfPoints() << endl;
    // for(size_t i = 0; i < tmp->GetNumberOfPoints(); ++i)
    // {
    //   double p[3];
    //   tmp->GetPoint(i, p);
    //   cout << "point: " << p[0] << " " << p[1] << " " << p[2] << endl;
    // }
    // double bounding[6];
    // for(size_t i = 0; i < 6; ++i)
    //   bounding[i] =0;
    vtkDataSet* data_set = static_cast<vtkAreaPicker*>(this->GetInteractor()->GetPicker())->GetDataSet();
    cout << "selected cell number " << data_set->GetNumberOfCells() << endl;
    vtkSmartPointer<vtkExtractGeometry> extract_poly_geometry =  
        vtkSmartPointer<vtkExtractGeometry>::New();
    // frustum->SetBounds(bounding);
    extract_poly_geometry->SetImplicitFunction(frustum);
    //  vtkSmartPointer<vtkExtractPolyDataGeometry> extract_poly_geometry =
//         vtkSmartPointer<vtkExtractPolyDataGeometry>::New();
//     extract_poly_geometry->SetImplicitFunction(frustum);
// #if
   //  VTK_MAJOR_VERSION <= 5
    // extract_poly_geometry->SetInputConnection(0, this->selected_poly->GetProducerPort());
// #else
//     extract_poly_geometry->SetInputData(this->selected_poly);
// #endif
//     extract_poly_geometry->Update();
//     vtkSmartPointer<vtkVertexGlyphFilter> glyph_filter =
//         vtkSmartPointer<vtkVertexGlyphFilter>::New();
//     glyph_filter->SetInputConnection(extract_poly_geometry->GetOutputPort());
//     glyph_filter->Update();
 
//     vtkPolyData* selected = glyph_filter->GetOutput();
//     std::cout << "Selected " << selected->GetNumberOfPoints() << " points." << std::endl;
    //     std::cout << "Selected " << selected->GetNumberOfCells() << " cells." << std::endl;
#if VTK_MAJOR_VERSION <= 5  
    extract_poly_geometry->SetInput(this->selected_poly);  
#else  
    extract_poly_geometry->SetInputData(this->selected_poly);  
#endif
// #if VTK_MAJOR_VERSION <= 5  
//     extract_poly_geometry->SetInput(1, data_set);  
// #else  
//     extract_poly_geometry->SetInputData(1, data_set);  
// #endif
    extract_poly_geometry->Update();
    vtkSmartPointer<vtkDataSetSurfaceFilter> surface_filter =
        vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surface_filter->SetInputConnection(extract_poly_geometry->GetOutputPort());
    surface_filter->PassThroughCellIdsOff();
    surface_filter->Update(); 
    // vtkSmartPointer<vtkVertexGlyphFilter> glyphFilter =
    //     vtkSmartPointer<vtkVertexGlyphFilter>::New();
    // glyphFilter->SetInputConnection(extract_poly_geometry->GetOutputPort());
    // glyphFilter->Update();
 
    vtkPolyData* selected_1 = surface_filter->GetOutput();
    std::cout << "Selected " << selected_1->GetNumberOfPoints() << " points." << std::endl;
    std::cout << "Selected " << selected_1->GetNumberOfCells() << " cells." << std::endl;
    vtkSmartPointer<vtkUnstructuredGrid> selected =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    selected->ShallowCopy(extract_poly_geometry->GetOutput());
    std::cout << "There are " << selected->GetNumberOfPoints()
            << " points in the selection." << std::endl;
    std::cout << "There are " << selected->GetNumberOfCells()
              << " cells in the selection." << std::endl;
    this->selected_mapper->SetInputConnection(extract_poly_geometry->GetOutputPort());
    this->selected_mapper->ScalarVisibilityOff();
// #if VTK_MAJOR_VERSION <= 5
//     this->selected_mapper->SetInputConnection(
//         selected->GetProducerPort());
// #else
//     this->selected_mapper->SetInputData(selected);
// #endif
//     this->selected_mapper->ScalarVisibilityOff();
    // vtkIdTypeArray* ids = vtkIdTypeArray::SafeDownCast(selected->GetPointData()->GetArray("OriginalIds"));
    // for(vtkIdType i = 0; i < ids->GetNumberOfTuples(); i++)
    // {
    //   std::cout << "Id " << i << " : " << ids->GetValue(i) << std::endl;
    // }
    //this->selected_actor->GetProperty()->SetColor(1.0, 0.0, 0.0); //(R,G,B)
    this->selected_actor->GetProperty()->EdgeVisibilityOn();
    this->selected_actor->GetProperty()->SetEdgeColor(1,0,0);
    this->selected_actor->GetProperty()->SetLineWidth(3);
    //    this->selected_actor->GetProperty()->SetPointSize(5);
    this->GetInteractor()->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(selected_actor);
    this->GetInteractor()->GetRenderWindow()->Render();
    this->HighlightProp(NULL);
  }
}



}
