/*=============================================================================*\
 *                                                                             *
 *                                 CATM                                        *
 *         Copyright (C) 2011-2021 by State Key Lab of CAD&CG 410              *
 *			Zhejiang University, shizeyun                          *
 *                        http://www.cad.zju.edu.cn/                           *
 *                                                                             *
 *-----------------------------------------------------------------------------* 
 *  This file is part of QuadRangular                                          *
 *  Created by shizeyun                                                        *
 *                                                                             *
\*=============================================================================*/

#ifndef SZY_INLINE_H
#define SZY_INLINE_H

//== INCLUDES ===================================================================

//STANDARD
#include<iostream>
#include<string>
#include<vector>

//== NAMESPACES =================================================================

namespace szy {

//=== IMPLEMENTATION ============================================================

//! @brief get suffix of file name 
inline const std::string get_filename_suffix(const std::string &file_name) {
  const size_t pos = file_name.find_last_of(".");
  return (std::string::npos == pos)? "" : file_name.substr(pos+1);
}

//-------------------------------------------------------------------------------

//! @brief get prefix of file name
inline const std::string get_filename_prefix(const std::string & file_name) {
  const size_t pos = file_name.find_last_of(".");
  return (std::string::npos == pos)? "" : file_name.substr(0, pos);
}

//-------------------------------------------------------------------------------

//! @brief if this file line carries no information
inline bool is_line_invalid(const std::string &line) {
  return (line.empty() || 13 == line[0]);
}

//-------------------------------------------------------------------------------

inline void split(const std::string &s, char delim,
                  std::vector<std::string> &elems) {
  std::stringstream ss(s);
  std::string item;
  while(std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
}

//===============================================================================

} //namespace szy

//===============================================================================

#endif //SZY_INLINE_H

//===============================================================================
