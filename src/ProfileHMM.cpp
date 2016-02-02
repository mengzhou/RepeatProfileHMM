/*    methcounts: a program for counting the methylated and
 *    unmethylated reads mapping over each CpG or C
 *
 *    Copyright (C) 2014 University of Southern California and
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

#include "ProfileHMM.hpp"

#include <cmath>

using std::vector;

inline double
ProfileHMM::log_sum_log(const double p, const double q) const {
  if (p == 0) {return q;}
  else if (q == 0) {return p;}
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
}

double
ProfileHMM::log_sum_log_vec(const vector<double> &vals, size_t limit) const {
  const vector<double>::const_iterator x = 
    std::max_element(vals.begin(), vals.begin() + limit);
  const double max_val = *x;
  const size_t max_idx = x - vals.begin();
  double sum = 1.0;
  for (size_t i = 0; i < limit; ++i) {
    if (i != max_idx) {
      sum += exp(vals[i] - max_val);
    }
  }
  return max_val + log(sum);
}

double
ProfileHMM::ViterbiDecoding(const vector<vector<double> > &transition,
    const vector<vector<double> > &emission, const vector<double> &initial,
    const vector<vector<double> > &bg,
    const size_t state, const size_t position,
    const vector<size_t> &observation, vector<size_t> &path) {

  // use a very small value as -Inf
  const double LOG_ZERO = -1000.0;
  vector<vector<double> > vm(num_model, vector<double>(num_seq, LOG_ZERO));
  vector<vector<double> > vi(num_model, vector<double>(num_seq, LOG_ZERO));
  vector<vector<double> > vd(num_model, vector<double>(num_seq, LOG_ZERO));

  // Initialization for vm, vi, vd
  // vm[j][i]: prob of given state M_j and observing sequence X_i
  // similarly for vi and vd.
  for (size_t j = 0; j < vm.size(); ++j) {
    vm[j][0] = LOG_ZERO;
    vi[j][0] = LOG_ZERO;
  }
  for (size_t i = 0; i < vm[0].size(); ++i) {
    // M_0 is the begin state with no emission
    vm[0][i] = LOG_ZERO;
    // D_0 does not exist
    vd[0][i] = LOG_ZERO;
  }
  vm[0][0] = 0;

  // I_0 is the background emission state
  vi[0][1] = transition[index_m[0]][index_i[1]];
  for (size_t i = 2; i < vi[0].size(); ++i) {
    vi[0][i] = vi[0][i-1] + transition[index_i[1]][index_i[1]];
  }
  // D_1 is the leading deletion state for local alignment
  vd[1][0] = transition[index_m[0]][index_d[1]];
  for (size_t i = 1; i < vd[1].size() - 1; ++i) {
    vd[1][i] = vi[0][i-1] + transition[index_i[1]][index_d[1]];
  }

  for (size_t j = 0; j < vm.size(); ++j) {
    for (size_t i =0; i < vm[0].size(); ++i) {
      vm;
    }
  }
}
