/**
 * @file  border_box.cpp
 * @author Li Quan
 * @version 1.0
 * @date 2012-11-20
 * 
 * @ingroup visualization_tet
 * @brief management the all operation of qt windows
 * @details
 *-load tet model control event
 *-load vtk model control event
 *-add box and plane control event
 *-delete one box and plane control event
 *-delete all box and plane control event
 *-get tet index control event
 *-save box position control event
 *-load the box position file control event
 * 
 * @history
 */


#include "border_box.h"
#include "../call_back/box_callback.h" 
#include "../call_back/plane_callback.h"
#include "../vtk_operate/vtk_operation.h"

#include "../call_back/area_select_callback.h"
#include "../side_bar/configure.h"
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>


#include <QtGui>
#include <QMessageBox>
#include <QErrorMessage>
#include <QEvent>
#include <QObject>

#include <string>
#include <cstring>

#include <vtkIndent.h>
#define CREATE_NEW(class, variable)                              \
  vtkSmartPointer<class> variable = vtkSmartPointer<class>::New();

namespace lq {
using namespace std;
const int max_angle = 1;
const int min_angle = 2;
const int add_angle = 1;
const int minus_angle = 2;
size_t count = 0;
size_t win_cnt_left = 0;
size_t win_cnt_right = 0;
bool left_win = false;
bool right_win = false;
bool is_link_camera = false;
bool left_render_loaded = false;
bool right_render_loaded = false;
vector<string> model_list;
string left_model_name, right_model_name;
vtk_operation vtk_op;
double axis[3] = {0, 0, 0};
// VTK Renderer
vtkSmartPointer<vtkRenderer> renderer;
vtkSmartPointer<vtkRenderer> left_render;
vtkSmartPointer<vtkRenderer> right_render;
vtkSmartPointer<vtkActor> actor;
vtkSmartPointer<vtkDataSetMapper> mapper_test;
vtkSmartPointer<vtkRenderer> render_test;
vtkSmartPointer<box_callback> callback;
vtkSmartPointer<plane_callback> callback_plane;
vtkSmartPointer<vtkPoints> points;

QString stripped_name(const QString &file_name)
{
    return QFileInfo(file_name).fileName();
}


void get_model_name(const std::string &file_path,
                           std::string &file_name)
{
  size_t flag;
  flag = file_path.find_last_of("/\\");
  file_name = file_path.substr(flag + 1);
}

int get_file_type(const string &file_name)
{
  string str;
  size_t flag;
  flag = file_name.find_last_of('.');
  if(flag != string::npos)
    str = file_name.substr(flag + 1);
  //cout << str << endl;
  if(strcmp(str.c_str(),"tet") == 0)
    return 1;
  else if(strcmp(str.c_str(), "vtk") == 0)
    return 2;
  else if(strcmp(str.c_str(), "obj") == 0)
    return 3;
  else
    return 0;
}

/*!\brief This function is to listen the event of qvtk widget, such as, keyboard event,
 *  mouse event.
 *
 *\param obj The caller object pointer.
 *\param event Caller event type and include some event information.
 *
 *\return bool
 */


bool border_box::eventFilter(QObject *obj, QEvent *event)
{
  QKeyEvent *keyEvent = static_cast<QKeyEvent *>(event);
  // if(event->type() == QEvent::KeyPress)
  // {
   
  //   if(keyEvent->key() == Qt::Key_R)
  //   {
  //     cout << "'R' is pressed !" << endl;
  //     close_link_camera();
  //   }
  //   else
  //     cout << "invilad " << endl;
  //   this->qvtk_widget_left->update();
  //   //this->qvtk_widget_right->update();
  //   return QObject::eventFilter(obj, event);
  // }
  if(obj == qvtk_widget_left)
  {
    if (event->type() == QEvent::KeyPress)
    {
      if(keyEvent->key() == Qt::Key_R)
      {
        cout << "win_cnt_left " << win_cnt_left << endl;
        ///when user press the 'r' key(the count of windwos is even),
        ///it mean user want to use the function of area pick.
        ///so, the camera share should be disconnected. If user press 'r' key again(windows count
        ///is odd), progrem auto to link the selected windows to the other one.
        close_left_link_camera();
        close_right_link_camera();
        if(win_cnt_left & 1)
          link_camera();
        win_cnt_left++;
      }
      for(box_it = boxwidget.begin(); box_it != boxwidget.end(); ++box_it)
      {
        if((*box_it)->is_selected())
          break;
      }
      if(box_it == boxwidget.end())
        return QObject::eventFilter(obj, event);
      if(keyEvent->key() == Qt::Key_Left)
        vtk_op.rotate_box(add_angle, max_angle, axis, (*box_it));
      else if(keyEvent->key() == Qt::Key_Right)
        vtk_op.rotate_box(minus_angle, max_angle, axis, (*box_it));
      else if(keyEvent->key() == Qt::Key_Up)
        vtk_op.rotate_box(add_angle, min_angle, axis, (*box_it));
      else if(keyEvent->key() == Qt::Key_Down)
        vtk_op.rotate_box(minus_angle, min_angle, axis, (*box_it));
   
      this->qvtk_widget_left->update();
      //this->qvtk_widget_right->update();
      return QObject::eventFilter(obj, event);
    }
    else
    {
      return QObject::eventFilter(obj, event);
    }
  }
  else if(obj == this->qvtk_widget_right)
  {
    if (event->type() == QEvent::KeyPress)
    {
      cout << "win_cnt_right " << win_cnt_right << endl;
      if(keyEvent->key() == Qt::Key_R)
      {
        close_left_link_camera();
        close_right_link_camera();
        if(win_cnt_right & 1)
          link_camera();
        win_cnt_right++;
      }
    }
    this->qvtk_widget_left->update();
    return QObject::eventFilter(obj, event);
  }
 
  
}

/*! 
 * \brief The call back function of share camera right to leftt 
 */


void border_box::modified_handler()
{
  this->qvtk_widget_left->GetRenderWindow()->Render();
}

/*!
 * \brief The call back function of share camera left to right 
 */

void border_box::modified_handler_left()
{
  this->qvtk_widget_right->GetRenderWindow()->Render();
}

void border_box::create_render(const string &file_name, const string &str,
                               QVTKWidget* qvtk, const size_t flag,
                               vtkSmartPointer<vtkRenderer> &render)
{
  this->side_panel->add_tree_leaf(file_name);
  if(get_file_type(file_name) == 1)
  {
    // vtk_op.set_qvtk(this->qvtk_widget_left);
    vtk_op.set_points(points);
    vtk_op.init_tet_render(str.c_str(), point, tet, tri, render);
  }
  else if( get_file_type(file_name) == 2)
    vtk_op.init_vtk_render(str.c_str(),  flag, render);
  else if( get_file_type(file_name) == 3)
    vtk_op.init_obj_render(qvtk, str.c_str(), render);
  else if( get_file_type(file_name) == 0)
  {
    cerr << "Can not open this file" << endl;
    cerr << __FILE__ << endl;
    cerr << __LINE__ << endl;
    return;
  }
  qvtk->GetRenderWindow()->AddRenderer(render);

  vtk_op.add_area_picker(qvtk, flag, render);
  qvtk->GetRenderWindow()->Render();
  qvtk->update();
}

bool border_box::load_file(const QString &file_name)
{
  
  //set_current_file(file_name);
  std::string str(qPrintable(file_name));
  string model_name;
  get_model_name(str, model_name);
  model_list.push_back(model_name);
  if(left_win && right_win)
  {
    QMessageBox messageBox;
    int ret = QMessageBox::warning(this, tr(WARNNING_TITLE),
                                   tr(NO_WIN_DISPLAY),
                                   QMessageBox::Cancel);
    return false;
  }
  else if(!left_win && !right_win)
  {
    QMessageBox messageBox;
    int ret = QMessageBox::warning(this, tr(WARNNING_TITLE),
                                   tr(NO_SELECT_WIN),
                                   QMessageBox::Cancel);
    return false;
  }
  else if(left_win)
  {
    create_render(model_name, str, qvtk_widget_left, LEFT, left_render);
    left_model_name = model_name;
    left_render_loaded = true;
  }
  else if(right_win)
  {
    create_render(model_name, str, qvtk_widget_right, RIGHT, right_render);
    right_model_name = model_name;
    right_render_loaded = true;
  }
 //vtk_op.add_axes(renderer);
  return true;
}

void border_box::init()
{
  renderer = vtkSmartPointer<vtkRenderer>::New();
  left_render = vtkSmartPointer<vtkRenderer>::New();
  right_render = vtkSmartPointer<vtkRenderer>::New();
  render_test = vtkSmartPointer<vtkRenderer>::New();
  callback = vtkSmartPointer<box_callback>::New();
  callback_plane = vtkSmartPointer<plane_callback>::New();
  points = vtkSmartPointer<vtkPoints>::New();
  callback->set_box_list(&boxwidget);
  callback->set_plane_list(&planewidget);
  callback_plane->set_box_list(&boxwidget);
  callback_plane->set_plane_list(&planewidget);
  this->link_camera_action->setEnabled(false);
  side_panel = NULL;
  add_side_panel();
}

border_box::border_box()
{
  
  this->setupUi(this);
  this->qvtk_widget_right->hide();
  qvtk_widget_left->installEventFilter(this);
  //  qvtk_widget_right->installEventFilter(this);
  init();
  create_action();
 
}

int border_box::save_file(const QString &file_name)
{
  string name;
  ofstream fout;
  fout.open(qPrintable(file_name));
  if(fout.fail())
  {
    cerr << "fail to save box!" << endl;
    cerr << "file:" << __FILE__ << endl;
    cerr << "line:" << __LINE__ << endl;
    return 1;
  }
  vtk_op.output_box_position(boxwidget, fout);
  fout.close();
  return 0;
}

int border_box::import_file(const QString &file_name)
{
  ifstream fin;
  fin.open(qPrintable(file_name));
  if(fin.fail())
  {
    cerr << "fail to import box!" << endl;
    cerr << "file:" << __FILE__ << endl;
    cerr << "line:" << __LINE__ << endl;
    return 1;
  }
  vtk_op.import_box_position(fin, this->qvtk_widget_left, boxwidget, callback);
  this->qvtk_widget_left->update();
}

void border_box::add_side_panel()
{
  if( !side_panel )
  {
    Qt::WindowFlags flags = Qt::SubWindow;
    side_panel = new side_bar(this);
    side_panel->setWindowFlags(flags);
    side_panel->setFixedWidth(230);
    addDockWidget(Qt::LeftDockWidgetArea, side_panel);
    //side_panel->show();
    //  side_panel->hide();
    QWidget::update();
  }
  else
    side_panel->hide();

}

void border_box::update_color(const size_t flag)
{
  if(!left_win && !right_win)
  {
    QMessageBox messageBox;
    int ret = QMessageBox::warning(this, tr(WARNNING_TITLE),
                                   tr(NO_SELECT_WIN),
                                   QMessageBox::Cancel);
    return ;
  }
  if(left_win)
  {
    vtk_op.change_color(flag, left_render, LEFT);
    this->qvtk_widget_left->update();
  }
  else if(right_win)
  {
    vtk_op.change_color(flag, right_render, RIGHT);
    this->qvtk_widget_right->update();
  }
}

void border_box::close_right_link_camera()
{
  
  cout << "Close right link camera " << endl; 
  this->qvtk_widget_right->GetRenderWindow()->RemoveAllObservers();
  left_render->ResetCameraClippingRange();
  right_render->ResetCameraClippingRange();
}

void border_box::close_left_link_camera()
{
  cout << "Close left link camera..." << endl;
  this->qvtk_widget_left->GetRenderWindow()->RemoveAllObservers();
  left_render->ResetCameraClippingRange();
  right_render->ResetCameraClippingRange();
}

/*! \brief this function is close right camera link and build a new
 *  camera share link in left to right.
 *
 *  This member function add the link_left_to_right in order to solve the problem of
 *  share the camera between two render windows.
 */


void border_box::link_left_to_right()
{
  close_right_link_camera();
  right_render->SetActiveCamera(left_render->GetActiveCamera());
  this->qvtk_widget_left->GetRenderWindow()->AddObserver(vtkCommand::AnyEvent,
                                                         this, &border_box::modified_handler_left);
}

/*! \brief this function is close left camera link and build a new
 *  camera share link in right to left.
 */
void border_box::link_right_to_left()
{
  close_left_link_camera();
  // right_render->ResetCamera();
  // left_render->ResetCamera();
 
  left_render->SetActiveCamera(right_render->GetActiveCamera());
  this->qvtk_widget_right->GetRenderWindow()->AddObserver(vtkCommand::AnyEvent,
                                                          this, &border_box::modified_handler);
}

void border_box::open_file()
{    
  QString fileName =QFileDialog::getOpenFileName(this,
                                                 tr("Open FILE"), ".",
                                                 tr("files (*.mesh *.vtk *.obj)"));
  if (!fileName.isEmpty())
  {
    load_file(fileName);
    status_bar->showMessage(tr("File open..."), 2000);
  }
  else
    status_bar->showMessage(tr("..."), 2000);
}

void border_box::save()
{
  QString file_name = QFileDialog::getSaveFileName(this,
                                                   tr("Save Box Position"), ".",
                                                   tr("Position File(*.pf)"));
  file_name += ".pf";
  if (file_name.isEmpty())
    return;
  if(!save_file(file_name))
    status_bar->showMessage(tr("File saved..."), 2000);
}

void border_box::import()
{
  QString file_name = QFileDialog::getOpenFileName(this,
                                                 tr("Open FILE"), ".",
                                                 tr("Position Files(*.pf)"));
  if(!file_name.isEmpty())
    import_file(file_name);
  
}

void border_box::add_box()
{
  vtk_op.add_box_widget(point, this->qvtk_widget_left, boxwidget, callback);
  
  this->qvtk_widget_left->update();
}

void border_box::add_plane()
{
  vtk_op.add_plane_widget(point, this->qvtk_widget_left, planewidget, callback_plane);
  this->qvtk_widget_left->update();
}

void border_box::add_line()
{
  std::vector<vtkSmartPointer<vtkLineWidget2> > line_list;
  vtk_op.add_line_widget(this->qvtk_widget_left, renderer, line_list);
  this->qvtk_widget_left->update();
}

void border_box::delete_one_box()
{
  
  for(box_it = boxwidget.begin(); box_it != boxwidget.end(); ++box_it)
  {
    if((*box_it)->is_selected())
    {
      boxwidget.erase(box_it);
      count--;
      break;
    }
  }
  for(plane_it = planewidget.begin(); plane_it != planewidget.end(); ++plane_it)
  {
    if((*plane_it)->is_selected())
    {
      (*plane_it)->Off();
      planewidget.erase(plane_it);
      break;
    }
  }
  this->qvtk_widget_left->update();
}

void border_box::get_tet()
{
  vtk_op.get_tet_data(boxwidget, tet);
  vtk_op.get_plane_tet(planewidget, tet, point);
  status_bar->showMessage(tr("tet index outputed..."), 2000);
}

void border_box::delete_all_box()
{
  boxwidget.clear();
  count = 0;
  this->qvtk_widget_left->update();
}

void border_box::rotate_x()
{
  axis[0] = 1;
  axis[1] = 0;
  axis[2] = 0;
}

void border_box::rotate_y()
{
  axis[0] = 0;
  axis[1] = 1;
  axis[2] = 0;
}

void border_box::rotate_z()
{
  axis[0] = 0;
  axis[1] = 0;
  axis[2] = 1;
}

void border_box::set_blue()
{
  update_color(BLUE);
}

void border_box::set_white()
{
  update_color(WHITE);
}

void border_box::set_red()
{
  update_color(RED);
}

void border_box::show_surface_edge()
{
  if(left_win)
  {
    vtk_op.display_surface_edge_mesh(LEFT);
    this->qvtk_widget_left->update();
  }
  else if(right_win)
  {
    vtk_op.display_surface_edge_mesh(RIGHT);
    this->qvtk_widget_right->update();
  }
}

void border_box::export_cell_data()
{
  if(!left_win && !right_win)
  {
    QMessageBox messageBox;
    int ret = QMessageBox::warning(this, tr(WARNNING_TITLE),
                                   tr(NO_SELECT_WIN),
                                   QMessageBox::Cancel);
    return ;
  }
  else if(left_win)
    vtk_op.output_cell_data(LEFT);
  else if(right_win)
    vtk_op.output_cell_data(RIGHT);
}

void border_box::display_dual_screen()
{
  this->qvtk_widget_left->show();
  this->qvtk_widget_left->update();
  this->qvtk_widget_right->show();
  this->qvtk_widget_right->update();
  this->qvtk_widget_right->installEventFilter(this);
  this->link_camera_action->setEnabled(true);
  status_bar->showMessage(tr(DISPLAY_DUAL_WIN), 2000);
}

void border_box::display_single_screen()
{
  this->qvtk_widget_right->hide();
  this->link_camera_action->setEnabled(false);
  status_bar->showMessage(tr(DISPLAY_SINGLE_WIN), 2000);
}

/*! \brief When click the left render windows 
 */
void border_box::left_win_click()
{
  right_win = false;
  left_win = true;
  if(left_render_loaded && right_render_loaded)
  {
    if(!(win_cnt_left & 1) )
      link_left_to_right();
  }
  status_bar->showMessage(tr(WIN_SELECTED_STATE_LEFT), 2000);
}

void border_box::right_win_click()
{
  left_win = false;
  right_win = true;
  status_bar->showMessage(tr(WIN_SELECTED_STATE_RIGHT), 2000);
}

void border_box::link_camera()
{
 
  // right_render->ResetCamera();
  // left_render->ResetCamera();
 
  left_render->SetActiveCamera(right_render->GetActiveCamera());
  this->qvtk_widget_right->GetRenderWindow()->AddObserver(vtkCommand::AnyEvent,
                                                          this, &border_box::modified_handler);
  is_link_camera = true;
}

void border_box::single_pick()
{
  if(left_win)
    vtk_op.switch_pick_type(1, LEFT);
  else if(right_win)
    vtk_op.switch_pick_type(1, RIGHT);
}

void border_box::multi_pick()
{
  if(left_win)
    vtk_op.switch_pick_type(2, LEFT);
  else if(right_win)
    vtk_op.switch_pick_type(2, RIGHT);
}

void border_box::del_model(const QString &model_name)
{
  string tmp;
  tmp = qPrintable(model_name);
  if(strcmp(tmp.c_str(), left_model_name.c_str()) == 0)
  {
    vtk_op.reset_render(LEFT);
    qvtk_widget_left->update();
  }
  else if(strcmp(tmp.c_str(), right_model_name.c_str()) == 0)
  {
    vtk_op.reset_render(RIGHT);
    qvtk_widget_right->update();
  }
}

void border_box::create_action()
{
   connect(this->open_file_action, SIGNAL(triggered()),
           this, SLOT(open_file()));
   connect(this->save_action, SIGNAL(triggered()),
           this, SLOT(save()));
   connect(this->import_action, SIGNAL(triggered()),
           this, SLOT(import()));
   connect(this->add_box_action, SIGNAL(triggered()),
           this, SLOT(add_box()));
   connect(this->add_plane_action, SIGNAL(triggered()),
           this, SLOT(add_plane()));
   connect(this->add_line_action, SIGNAL(triggered()),
           this, SLOT(add_line()));
   connect(this->delete_one_action, SIGNAL(triggered()),
           this, SLOT(delete_one_box()));
   connect(this->delete_all_action, SIGNAL(triggered()),
           this, SLOT(delete_all_box()));
   
   connect(this->get_tet_action, SIGNAL(triggered()),
           this, SLOT(get_tet()));
   
   connect(this->rotate_x_action, SIGNAL(triggered()),
           this, SLOT(rotate_x()));
   connect(this->rotate_y_action, SIGNAL(triggered()),
           this, SLOT(rotate_y()));
   connect(this->rotate_z_action, SIGNAL(triggered()),
           this, SLOT(rotate_z()));

   connect(this->blue_action, SIGNAL(triggered()),
           this, SLOT(set_blue()));
   connect(this->white_action, SIGNAL(triggered()),
           this, SLOT(set_white()));
   connect(this->red_action, SIGNAL(triggered()),
           this, SLOT(set_red()));
   connect(this->surface_edge_action, SIGNAL(triggered()),
           this, SLOT(show_surface_edge()));
   connect(this->export_cell_action, SIGNAL(triggered()),
           this, SLOT(export_cell_data()));
   connect(this->dual_screen_action, SIGNAL(triggered()),
           this, SLOT(display_dual_screen()));
   connect(this->single_screen_action, SIGNAL(triggered()),
           this, SLOT(display_single_screen()));

   //set connect with qvtkwidget in order to get the customer selection windows
   
   connect(this->qvtk_widget_left, SIGNAL(mouseEvent(QMouseEvent*)),
           this, SLOT(left_win_click()));
   connect(this->qvtk_widget_right, SIGNAL(mouseEvent(QMouseEvent*)),
           this, SLOT(right_win_click()));

   connect(this->link_camera_action, SIGNAL(triggered()),
           this, SLOT(link_camera()));

   connect(this->single_pick_action, SIGNAL(triggered()),
           this, SLOT(single_pick()));
   connect(this->multi_pick_action, SIGNAL(triggered()),
           this, SLOT(multi_pick()));

   connect(side_panel, SIGNAL(delete_model(const QString &)),
           this, SLOT(del_model(const QString &)));
   
}
}
