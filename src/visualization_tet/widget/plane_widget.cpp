#include "plane_widget.h"
namespace lq {


bool plane_widget::is_selected()
{
  return this->selected;
}

void plane_widget::set_selected(bool b)
{
  this->selected = b;
}


}
