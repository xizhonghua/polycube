#ifndef BOX_WIDGET_H
#define BOX_WIDGET_H
#include <vtkBoxWidget2.h>
#include <vector>

namespace lq {

  class box_widget : public vtkBoxWidget2 {

 public:
    static box_widget *New()
    {
      return new box_widget;
    }
    box_widget(){}
    std::vector<int>& get_point_index();
    bool is_selected();
    size_t get_num();
    
    void set_point_index(const std::vector<int> &index);
    void set_num(const size_t n);
    void set_selected(const bool &b);

 private:
    
    std::vector<int> point_index;
    size_t num;
    bool selected;
    
  };
}

#endif
