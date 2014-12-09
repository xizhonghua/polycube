#ifndef PLANE_WIDGET_H
#define PLANE_WIDGET_H

#include <vtkPlaneWidget.h>
namespace lq {

  class plane_widget : public vtkPlaneWidget
      
  {

 public:
    static plane_widget *New()
    {
      return new plane_widget;
    }
    plane_widget(){}
    bool is_selected();
    void set_selected(bool b);
 private:

    bool selected;
  };
}

#endif
