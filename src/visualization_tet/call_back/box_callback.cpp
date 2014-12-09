/**
 * @file  box_callback.cpp
 * @author Li Quan
 * @version 1.0
 * @date 2012-11-20
 * 
 * @ingroup visualization_tet
 * @brief deal the callback event of box wodget
 * @details
 *-Execute reload the virtual function 
 * 
 * @history
 *
 */

#include "../call_back/box_callback.h"
#include <vtkPolyData.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkBoxRepresentation.h>
#include <vtkTransform.h>
#include <cassert>
//#define TEST_BOX_TRANSFORM
namespace lq {

vtkSmartPointer<vtkSelectEnclosedPoints> select =
    vtkSmartPointer<vtkSelectEnclosedPoints>::New();

int box_callback::get_index()
{
  return this->index;
}

void box_callback::set_index(const int num)
{
  this->index = num;
}

void box_callback::set_box_list(std::vector<vtkSmartPointer<box_widget> > *box)
{
  this->box_list = box;
}

void box_callback::set_plane_list(std::vector<vtkSmartPointer<plane_widget> > *plane)
{
  this->plane_list = plane;
}
box_callback::~box_callback()
{
  delete box_list;
}

void clear_plane_widget(std::vector<vtkSmartPointer<plane_widget> > *plane)
{
  std::vector<vtkSmartPointer<plane_widget> >::iterator it;
  for(it = plane->begin(); it != plane->end(); ++it)
    (*it)->set_selected(false);
}


void box_callback::Execute(vtkObject *caller, unsigned long eventId, void*)
{
  assert(box_list != NULL);
  //std::cout << "box_list.size " << box_list->size() << std::endl;
  vtkSmartPointer<box_widget> box =
      vtkSmartPointer<box_widget>::New();
  //get the real box widget, You can get or set the data in the class box widget
  box =  reinterpret_cast<box_widget*>(caller);
  
  std::vector<vtkSmartPointer<box_widget> >::iterator box_it;
  clear_plane_widget(plane_list);

  vtkSmartPointer<box_widget> tmp;
  
  for(box_it = box_list->begin(); box_it != box_list->end(); ++box_it)
  {
    tmp = *(box_it);
    if(tmp == box)
      tmp->set_selected(true);
    else
      tmp->set_selected(false);
  }
  std::cout << "box_list size " << box_list->size() << endl;
  // box->set_selected(true);
  // std::cout << "event " << GetStringFromEventId(eventId) <<std::endl;
}

}
