#ifndef CONFIGURE_H
#define CONFIGURE_H

namespace lq {
  //use in main_panel.cpp
  const char TREE_TITLE[] = "Model List";
  const char MODEL[] = "model";
  const int STD_MAX_RANGE = 100;
  const int STD_ZERO = 0;
  const int STD_DEFAULT_VALUE = 50;
  //use in model_win.cpp
  const char ABOUT[] = "About";
  const char ABOUT_CONTENT[] =
      "This is surface quality inspection program\nversion 1.0\n";
  const char  WARNNING_TITLE[] = "warnning";
  const char  NO_WIN_DISPLAY[] = "This is no window to display file.";
  const	char  NO_SELECT_WIN[] = "There is no selected windows!";
  const char  NO_TWO_MODEL[] =
      "there is no two model to compare!\nplease load two model first!";
  const char  NO_COMPARE_TWO_MODEL[] = "You should compute error first!";
  const char  NO_COLOR_FILE[] = "Can't find color data!";
  //using in border_box.cpp
  const char FILE_OPEN[] = "File have been opened...";
  const char WIN_SELECTED_STATE_LEFT[] = "Left render windows is selected...";
  const char WIN_SELECTED_STATE_RIGHT[] = "Right render windows is selected...";
  const char DISPLAY_DUAL_WIN[] = "Start the dual render windows...";
  const char DISPLAY_SINGLE_WIN[] = "Close one render windows...";
  
}
#endif
