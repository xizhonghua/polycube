#include <QApplication>
#include "../main_win/border_box.h"

int main(int argc, char* argv[])
{
  QApplication app( argc, argv );
  lq::border_box borderWidgetQt;
  borderWidgetQt.show();
  
  return app.exec();
}
