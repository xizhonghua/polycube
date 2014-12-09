#include "../call_back/plane_callback.h"
#include <cassert>
namespace lq {

void plane_callback::set_plane_list(std::vector<vtkSmartPointer<plane_widget> > *plane)
{
  this->plane_list = plane;
}

void plane_callback::set_box_list(std::vector<vtkSmartPointer<box_widget> > *box)
{
  this->box_list = box;
}

plane_callback::~plane_callback()
{
  delete plane_list;
}

void clear_box_widget(std::vector<vtkSmartPointer<box_widget> > *box)
{
  std::vector<vtkSmartPointer<box_widget> >::iterator it;
  for(it = box->begin(); it != box->begin(); ++it)
    (*it)->set_selected(false);
}

void plane_callback::Execute(vtkObject *caller, unsigned long eventId, void*)
{
  assert(plane_list != NULL);
  vtkSmartPointer<plane_widget> plane =
      vtkSmartPointer<plane_widget>::New();
  plane = reinterpret_cast<plane_widget *>(caller);
  vtkSmartPointer<plane_widget> tmp;
  std::vector<vtkSmartPointer<plane_widget> >::iterator plane_it;
  clear_box_widget(box_list);
  for(plane_it = plane_list->begin(); plane_it != plane_list->end(); ++plane_it)
  {
    tmp = (*plane_it);
    if(tmp == plane)
      (*plane_it)->set_selected(true);
    else
      (*plane_it)->set_selected(false);
    cout << "plane is_selected " << (*plane_it)->is_selected() << endl;
  }
  
}


}
