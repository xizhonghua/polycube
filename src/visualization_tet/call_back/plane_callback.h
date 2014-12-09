#ifndef PLANE_CALLBACK_H
#define PLANE_CALLBACK_H

#include "../widget/plane_widget.h"
#include "../widget/box_widget.h"

#include <vtkSmartPointer.h>
#include <vtkCallbackCommand.h>
#include <vector>

namespace lq {

  class plane_callback : public vtkCallbackCommand
  {

 public:
    void set_plane_list(std::vector<vtkSmartPointer<plane_widget> > *plane);
    void set_box_list(std::vector<vtkSmartPointer<box_widget> > *box);
    static plane_callback *New()
    {
      return new plane_callback;
    }
    virtual void Execute(vtkObject *caller, unsigned long eventId, void*);
    plane_callback(){}
    ~plane_callback();
    
 private:

    std::vector<vtkSmartPointer<box_widget> > *box_list;
    std::vector<vtkSmartPointer<plane_widget> > *plane_list;
    
    
  };
}


#endif
