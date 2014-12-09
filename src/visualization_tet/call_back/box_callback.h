/**
 *@class box_callback
 *@brief deal the callback event of box wodget 

 *@detail
 *-Execute reload the virtual function 
 *
 *@author Li Quan
 *@data 2012-11-20
 *@version 1.0
 */

#ifndef BOX_CALLBACK_H
#define BOX_CALLBACK_H
#include "../widget/box_widget.h"
#include "../widget/plane_widget.h"
#include <vtkSmartPointer.h>
#include <vtkCallbackCommand.h>
#include <vtkPoints.h>
#include <vector>

//class box_widget;

namespace lq {

  class box_callback : public vtkCallbackCommand 
  {
 public:
    int get_index();
    void set_index(const int num);
    void set_box_list(std::vector<vtkSmartPointer<box_widget> > *box);
    void set_plane_list(std::vector<vtkSmartPointer<plane_widget> > *plane);
    static box_callback *New()
    {
      return new box_callback;
    }
    virtual void Execute(vtkObject *caller, unsigned long eventId, void*);
    box_callback(){}
    ~box_callback();

 private:
    std::vector<vtkSmartPointer<box_widget> > *box_list;
    std::vector<vtkSmartPointer<plane_widget> > *plane_list;
    int index;

  };
  

}


#endif
