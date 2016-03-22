/*    Copyright (C) 2016 University of Southern California and
 *                       Meng Zhou
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MATRIXUTIL_HPP
#define MATRIXUTIL_HPP

#include <initializer_list>
#include <vector>

#include "smithlab_utils.hpp"

typedef std::vector<std::vector<double> > matrix;

extern const double LOG_ZERO;

void
normalize_vec_inplace(std::vector<double> &v, const bool logged);

double
log_sum_log(const double p, const double q);

double
log_sum_log_list(const std::initializer_list<double> &list);

size_t
argmax_list(const std::initializer_list<double> &list);

size_t
argmax_vec(const std::vector<double> &v);

#endif
