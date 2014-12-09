#ifndef VOL_PARAM_COMMON_IO_H
#define VOL_PARAM_COMMON_IO_H

#include <vector>
#include <zjucad/matrix/matrix.h>

/// \brief load_equation_file load transition equation
/// \param file
/// \param variable_idx
/// \param coefficient
/// \return 0 if ok, or return 1
///
int load_equation_file(
    const char * file,
    std::vector<std::vector<size_t> > &variable_idx,
    std::vector<std::vector<double> > &coefficient,
    const zjucad::matrix::matrix<size_t> * node_mapping = 0);

///
/// \brief dump_equation_file
/// \param file
/// \param variable_idx
/// \param variable_coeff
/// \return
///
int dump_equation_file(
    const char * file,
    const std::vector<std::vector<size_t> > & variable_idx,
    const std::vector<std::vector<double> > & variable_coeff);

///
/// \brief load_variable_group_file
/// \param file
/// \param groups
/// \return
///
int load_variable_group_file(
    const char * file,
    std::vector<std::vector<size_t> > &groups);


///
/// \brief load_integer_constraint
/// \param filename
/// \param integer_variables
/// \return
///
int load_integer_constraint(
    const char * filename,
    std::vector<size_t> & integer_variables);

///
/// \brief load_restricted_path
/// \param filename
/// \param restricted_path
/// \return
///
int load_restricted_path(
    const char * filename,
    std::vector<std::pair<size_t,size_t> > & restricted_path);
#endif //

