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

#include "matrix_utils.hpp"

#include <cmath>
#include <algorithm>

using std::vector;

const double LOG_ZERO = -1e20;

void
normalize_vec_inplace(vector<double> &v, const bool logged = true) {
  double sum = logged ?
    smithlab::log_sum_log_vec(v, v.size())
    : accumulate(v.begin(), v.end(), 0.0);
  for (vector<double>::iterator i = v.begin();
      i < v.end(); ++i)
    if (logged)
      *i = *i - sum;
    else
      *i = sum < 1e-6 ? 0.0 : *i / sum;
}

double
log_sum_log(const double p, const double q) {
  if (p == 0) {return q;}
  else if (q == 0) {return p;}
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
}

double
log_sum_log_list(const std::initializer_list<double> &list) {
  const vector<double> vals(list);
  const vector<double>::const_iterator x = 
    std::max_element(vals.begin(), vals.end());
  const double max_val = *x;
  const size_t max_idx = x - vals.begin();
  double sum = 1.0;
  for (size_t i = 0; i < vals.size(); ++i) {
    if (i != max_idx) {
      sum += exp(vals[i] - max_val);
    }
  }
  return max_val + log(sum);
}

size_t
argmax_list(const std::initializer_list<double> &list) {
  const vector<double> v(list);
  const vector<double>::const_iterator x = 
     std::max_element(v.begin(), v.end());
  return x - v.begin();
}

size_t
argmax_vec(const vector<double> &v) {
  const vector<double>::const_iterator x = 
     std::max_element(v.begin(), v.end());
  return x - v.begin();
}
