#ifndef SELECT_INTERACTOR_STYLE_H
#define SELECT_INTERACTOR_STYLE_H
#include <vtkUnstructuredGrid.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkIdFilter.h>
#include <vtkInteractorStyleRubberBandPick.h>
#include <vtkAreaPicker.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkActor.h>
#include <vtkDataSetMapper.h>
#include <vtkObjectFactory.h>
#define VTKISRBP_ORIENT 0
#define VTKISRBP_SELECT 1
namespace lq {

  class select_interactor_style : public vtkInteractorStyleRubberBandPick
  {
 public :

    static select_interactor_style* New()
    {
      return new select_interactor_style;
    }
    vtkTypeMacro(select_interactor_style, vtkInteractorStyleRubberBandPick);
 select_interactor_style() : vtkInteractorStyleRubberBandPick()
    {
      this->selected_poly = vtkSmartPointer<vtkPolyData>::New();
      this->selected_actor = vtkSmartPointer<vtkActor>::New();
      this->selected_mapper = vtkSmartPointer<vtkDataSetMapper>::New();
      this->selected_actor->SetMapper(this->selected_mapper);
    }
    virtual void OnLeftButtonUp();
    void set_poly_data(const vtkSmartPointer<vtkPolyData> &poly);
 private:
    vtkSmartPointer<vtkPolyData> selected_poly;
    vtkSmartPointer<vtkActor> selected_actor;
    vtkSmartPointer<vtkDataSetMapper> selected_mapper;
  };
  //vtkStandardNewMacro(select_interactor_style);
}

#endif
