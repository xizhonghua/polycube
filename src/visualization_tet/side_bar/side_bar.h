#ifndef MAINPANEL_H
#define MAINPANEL_H

#include"ui_side_bar.h"
#include<string>
#include<QMainWindow>

namespace lq {
  
  class side_bar : public QDockWidget, public Ui::side_bar
  {
    Q_OBJECT

 public:
    side_bar(QWidget* parent = 0);
    ~side_bar();
    void set_tree(const std::vector<std::string> &model_list);
    void add_tree_leaf(const std::string &file_name);
    QTreeWidgetItem *root;

 signals:
    void delete_model(const QString &);
    void send_adjust_color(const int &);
    void send_color_legend();

 public slots:
    void select_item(QTreeWidgetItem *item, int column);
    void delete_btn_down();
    void adjust_color();
    void show_color_legend();

 private:
    // pointer to the main window (Optional)
    // QMainWindow* m_pMainWindow;
    void init();
    void create_action();
   
  };
}
#endif
