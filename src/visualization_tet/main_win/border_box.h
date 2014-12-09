#ifndef BorderWidgetQt_H
#define BorderWidgetQt_H

#include "ui_BorderWidgetQt.h"
#include "../widget/box_widget.h"
#include "../widget/plane_widget.h"
#include "../configure/data_type.h"
#include "../side_bar/side_bar.h"
#include <vtkSmartPointer.h>
#include <QMainWindow>
#include <vtkBoxWidget2.h>
#include <vtkBoxRepresentation.h>

#include <vector>
class vtkBorderWidget;
class vtkBoxWidget2;

namespace lq {

  class border_box : public QMainWindow, public Ui::BorderWidgetQt
  {
    Q_OBJECT
        public:
 
    // Constructor/Destructor
    border_box();
    ~border_box() {}

 signals:
    
    public slots:
    void open_file();
    void save();
    void import();
    void add_box();
    void delete_one_box();
    void delete_all_box();
    void get_tet();
    void rotate_x();
    void rotate_y();
    void rotate_z();
    void add_plane();
    void add_line();
    void set_blue();
    void set_white();
    void set_red();
    void show_surface_edge();
    void export_cell_data();
    void left_win_click();
    void right_win_click();
    void display_dual_screen();
    void display_single_screen();
    void link_camera();
    void single_pick();
    void multi_pick();
    void del_model(const QString &model_name);
    //void delete_model();
 protected:
    bool eventFilter(QObject *obj, QEvent *event);
    void modified_handler();
    void modified_handler_left();
 private:

    side_bar* side_panel;    
    matrixd point;
    matrixst tet;
    matrixst tri;
  
   
    std::vector<vtkSmartPointer<box_widget> > boxwidget;
    std::vector<vtkSmartPointer<plane_widget> > planewidget;
    std::vector<vtkSmartPointer<box_widget> >::iterator box_it;
    std::vector<vtkSmartPointer<plane_widget> >::iterator plane_it;
    bool load_file(const QString &file_name);
    void init();
    void create_render(const string &file_name, const string &str,
                       QVTKWidget* qvtk, const size_t flag,
                       vtkSmartPointer<vtkRenderer> &render);
    void create_action();
    int save_file(const QString &file_name);
    int import_file(const QString &file_name);
    void add_side_panel();
    void update_color(const size_t flag);
    void close_right_link_camera();
    void close_left_link_camera();
    void link_left_to_right();
    void link_right_to_left();
  };
}
#endif
