#include "box_widget.h"

namespace lq {

std::vector<int>& box_widget::get_point_index()
{
  return this->point_index;
}

size_t box_widget::get_num()
{
  return this->num;
}

bool box_widget::is_selected()
{
  return this->selected;
}

void box_widget::set_point_index(const std::vector<int> &index)
{
  this->point_index = index;
}


void box_widget::set_num(const size_t n)
{
  this->num = n;
}

void box_widget::set_selected(const bool &b)
{
  this->selected = b;
}
}
