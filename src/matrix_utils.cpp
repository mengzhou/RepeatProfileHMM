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
#include <numeric>
#include <iomanip>

using std::vector;
using std::pair;
using std::endl;
using std::accumulate;

const double LOG_ZERO = -1e20;

void
normalize_vec_inplace(vector<double> &v, const bool logged) {
  // potentially needed to add another parameter for the scaling factor,
  // since not in all cases is the vector required to be normalized to
  // the sum of 1.
  const double sum = logged ?
    smithlab::log_sum_log_vec(v, v.size())
    : std::accumulate(v.begin(), v.end(), 0.0);
  for (vector<double>::iterator i = v.begin();
       i < v.end(); ++i)
    if (logged)
      *i = *i - sum;
    else
      *i = sum < 1e-6 ? 0.0 : *i / sum;
}

double
log_sum_log(const double p, const double q) {
  //if (p == 0) {return q;}
  //else if (q == 0) {return p;}
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
}

double
log_minus_log(const double p, const double q) {
  if (p == 0) {return -q;}
  else if (q == 0) {return p;}
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  const double diff = larger + log(1.0 - exp(smaller - larger));
  return diff;
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

void
log_transform_matrix(matrix &m) {
  for (matrix::iterator i = m.begin();
       i < m.end(); ++i) {
    for (vector<double>:: iterator j = (*i).begin();
         j < (*i).end(); ++j) {
      if (*j < 1e-6)
        *j = LOG_ZERO;
      else
        *j = log(*j);
    }
  }
}

void
log_transform_vec(vector<double> &v) {
  for (vector<double>:: iterator i = v.begin();
       i < v.end(); ++i) {
    if (*i < 1e-6)
      *i = LOG_ZERO;
    else
      *i = log(*i);
  }
}

vector<double>
normalize_vec(const vector<double> &v, const bool logged) {
  static const double MAGIC_CUTOFF = 1e-6;
  vector<double> n;
  double sum = logged ?
    smithlab::log_sum_log_vec(v, v.size())
    : accumulate(v.begin(), v.end(), 0.0);
  for (vector<double>::const_iterator i = v.begin();
       i < v.end(); ++i)
    if (logged)
      n.push_back(*i - sum);
    else
      n.push_back(sum < MAGIC_CUTOFF ? 0.0 : *i / sum);
  return n;
}

vector<double>
combine_normalize(const vector<double> &v1, const vector<double> &v2,
                  const bool logged,
                  const pair<double, double> &weight) {
  assert(v1.size() == v2.size());
  vector<double> combined, n1, n2;
  n1 = normalize_vec(v1, logged);
  n2 = normalize_vec(v2, logged);
  for (size_t i = 0; i < n1.size(); ++i) {
    if (logged)
      combined.push_back(log_sum_log(n1[i] + log(weight.first),
                                     n2[i] + log(weight.second)));
    else
      combined.push_back(n1[i]*weight.first + n2[i]*weight.second);
  }
  combined = normalize_vec(combined, logged);
  return combined;
}

std::ostream&
operator<<(std::ostream &s, const matrix &m) {
  s << std::setprecision(4) << std::fixed;
  for (size_t k = 0; k < m.front().size(); ++k)
    s << "\t" << k;
  s << "\tRowSum" << endl;
  for (vector<vector<double> >::const_iterator i = m.begin();
      i < m.end(); ++i) {
    vector<double> list;
    s << i - m.begin();
    for (vector<double>::const_iterator j = i->begin();
        j < i->end(); ++j) {
      //s << "\t" << exp(*j);
      s << "\t" << -*j;
      list.push_back(*j);
    }
    //s << "\t" << exp(smithlab::log_sum_log_vec(list, list.size())) << endl;
    s << "\t" << -smithlab::log_sum_log_vec(list, list.size()) << endl;
  }
  return s;
}

std::ostream&
operator<<(std::ostream &s, const vector<double> &v) {
  //s << std::setprecision(4) << std::fixed;
  for (size_t k = 0; k < v.size(); ++k)
    s << k << "\t";
  s << endl;
  for (vector<double>::const_iterator i = v.begin();
      i < v.end(); ++i)
    s << *i << "\t";
  s << endl;
  return s;
}
