/*    methcounts: a program for counting the methylated and
 *    unmethylated reads mapping over each CpG or C
 *
 *    Copyright (C) 2016 University of Southern California and
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
#include <initializer_list>
#include <algorithm>
#include <vector>

using std::vector;
using std::max;
using std::pair;

inline double
ProfileHMM::log_sum_log(const double p, const double q) const {
  if (p == 0) {return q;}
  else if (q == 0) {return p;}
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
}

double
ProfileHMM::max_item(const std::initializer_list<double> list) {
  const vector<double> v(list);
  const vector<double>::const_iterator x = 
     std::max_element(v.begin(), v.end());
  return x - v.begin();
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

size_t
ProfileHMM::index_m(const size_t idx) const {
  // M_0 ~ M_L
  return idx;
}

size_t
ProfileHMM::index_i(const size_t idx) const {
  // I_0 ~ I_L-1 + (L+1)*M
  return idx + model_len + 1;
}

size_t
ProfileHMM::index_d(const size_t idx) const {
  // D_1 ~ D_L + (L+1)*M + L*I
  return idx + model_len * 2;
}

double
ProfileHMM::ViterbiDecoding(const bool VERBOSE,
    const vector<vector<double> > &transition,
    const vector<vector<double> > &emission,
    const vector<double> &initial, const vector<double> &ending,
    const vector<vector<double> > &bg,
    const size_t state, const size_t position,
    const vector<size_t> &observation,
    vector<pair<char, size_t> > &trace) {

  // use a very small value as -Inf
  const double LOG_ZERO = -1000.0;
  vector<vector<double> > vm(model_len, vector<double>(seq_len, LOG_ZERO));
  vector<vector<double> > vi(model_len, vector<double>(seq_len, LOG_ZERO));
  vector<vector<double> > vd(model_len, vector<double>(seq_len, LOG_ZERO));

  // Initialization for vm, vi, vd
  // vm[j][i]: prob of given state M_j and observing sequence X_i
  // similarly for vi and vd.
  for (size_t j = 0; j < model_len; ++j) {
    vm[j][0] = LOG_ZERO;
    vi[j][0] = LOG_ZERO;
  }
  for (size_t i = 0; i < seq_len; ++i) {
    // M_0 is the begin state with no emission
    vm[0][i] = LOG_ZERO;
    // D_0 does not exist
    vd[0][i] = LOG_ZERO;
  }
  vm[0][0] = 0;

  // I_0 is the background emission state
  vi[0][1] = transition[index_m(0)][index_i(1)];
  for (size_t i = 2; i < seq_len; ++i) {
    vi[0][i] = vi[0][i-1] + transition[index_i(1)][index_i(1)];
  }
  // D_1 is the leading deletion state for local alignment
  vd[1][0] = transition[index_m(0)][index_d(1)];
  for (size_t i = 1; i < seq_len - 1; ++i) {
    vd[1][i] = vi[0][i-1] + transition[index_i(1)][index_d(1)];
  }

  // some special cases before going to the main loop
  for (size_t i = 1; i < seq_len; ++i) {
    // D_2 only has two incoming transition from I_1 and M_1
    vd[2][i] = max( {vm[1][i] + transition[index_m(1)][index_d(2)],
      vi[1][i] + transition[index_i(1)][index_d(2)]});
  }

  // main loop
  for (size_t i = 1; i < seq_len; ++i) {
    // M_1 only has one incoming transition from D_1
    vm[1][i] = vd[1][i-1] + transition[index_d(1)][index_m(j)]
      + emission[index_m(1)][observation[i]];
    // I_1 only has two incoming transition including one from M_1
      ;
    vi[1][i] = max(
      {vm[1][i-1] + transition[index_m(1)][index_i(1)],
      vi[1][i-1] + transition[index_i(1)][index_i(1)]})
      + emission[index_i(1)][observation[i]]
      - bg[observation[i]];
    for (size_t j = 2; j < model_len - 1; ++j) {
      // M_i
      vm[j][i] = max(
        {vm[j-1][i-1] + transition[index_m[j-1]][index_m(j)],
        vi[j-1][i-1] + transition[index_i[j-1]][index_m(j)],
        vd[j-1][i-1] + transition[index_d[j-1]][index_m(j)],
        vd[1][i-1] + transition[index_d(1)][index_m(j)]})
        + emission[index_m(j)][observation[i]]
        - bg[observation[i]];
      // I_i
      vi[j][i] = max(
        {vm[j][i-1] + transition[index_m(j)][index_i(j)],
        vi[j][i-1] + transition[index_i(j)][index_i(j)],
        vd[j][i-1] + transition[index_d(j)][index_i(j)]})
        + emission[index_i(j)][observation[i]]
        - bg[observation[i]];
      // D_i
      vd[j][i] = max( 
        {vm[j-1][i] + transition[index_m[j-1]][index_d(j)],
        vi[j-1][i] + transition[index_i[j-1]][index_d(j)],
        vd[j-1][i] + transition[index_d[j-1]][index_d(j)]});
    }
    // M_model_len
    vm[model_len][i] = max( 
      {vm[model_len-1][i-1]
        + transition[index_m[model_len-1]][index_m[model_len]],
      vi[model_len-1][i-1]
        + transition[index_i[model_len-1]][index_m[model_len]],
      vd[model_len-1][i-1]
        + transition[index_d[model_len-1]][index_m[model_len]],
      vd[1][i-1] + transition[index_d(1)][index_m(L)]})
      + emission[index_m[model_len]][observation[i]]
      - bg[observation[i]];
    // I_model_len does not exist
    // D_model_len
    vector<double> list;
    for (size_t k = 1; k < model_len; ++k) {
      list.push_back(vm[k][i]);
    }
    const vector<double>::const_iterator x = 
      std::max_element(list.begin(), list.end());
    size_t max_idx = x - list.begin();
    vd[model_len][i] = *x + transition[index_m[max_idx+1]][index_d[model_len]];
    // I_0
    vi[0][i] = max( 
      {vi[1][i-1] + transition[index_i(1)][index_i(1)],
      vd[model_len][i] + transition[index_d[model_len]][index_i(1)]});
    // D_1
    vd[1][i] = vi[0][i] + transition[index_i(0)][index_d(1)];
  }

  //traceback
  trace.push_back('E', 0);
  size_t state_idx;
  char state;
  if (vd[model_len][seq_len] > vi[0][seq_len]) {
    state = 'D';
    state_idx = model_len;
  }
  else {
    state = 'I';
    state_idx = 0;
  }

  size_t seq_idx = seq_len;
  while (!(state == 'M' && state_idx == 0) && seq_idx > 0) {
    trace.push_back(state, state_idx);
    switch (state) {
      case 'M':
        // assuming state_idx = 1~L
        max_idx = max_item({vm[state_idx-1][seq_idx-1],
          vi[state_idx-1][seq_idx-1],
          vd[state_idx-1][seq_idx-1],
          vd[1][seq_idx-1]});
        switch (max_idx) {
          case 0:
            ;
          case 1:
            state = 'I';
          case 2:
            state = 'D';
          case 3:
            state = 'D';
        }
        --state_idx;
        --seq_idx;
      case 'I':
        if (state_idx == 0) {
          max_idx = max_item({vm[0][seq_idx-1],
            vi[0][seq_idx-1]},
            vd[model_len][seq_idx-1]);
          switch (max_idx) {
            case 0:
              state = 'M';
              state_idx = 0;
            case 1:
              ;
            case 2:
              state = 'D';
              state_idx = model_len;
          }
        }
        else {
          max_idx = max_item({vm[state_idx][seq_idx-1],
            vi[state_idx][seq_idx-1]},
            vd[model_len][state_idx-1]);
          switch (max_idx) {
            case 0:
              state = 'M';
            case 1:
              ;
            case 2:
              state = 'D';
          }
        --seq_idx;
        }
      case 'D':
        if (state_idx == model_len) {
          vector<double> list;
          for (size_t k = 1; k < model_len; ++k) {
            list.push_back(vm[k][seq_idx]);
          }
          const vector<double>::const_iterator x = 
            std::max_element(list.begin(), list.end());
          state_idx = x - list.begin();
        }
        else if (state_idx == 1) {
          max_idx = max_item({vm[0][seq_idx],
            vi[0][seq_idx]});
          switch (max_idx) {
            case 0:
              state = 'M';
            case 1:
              state = 'I';
          }
          state_idx = 0;
        }
        else {
          max_idx = max_item({vm[state_idx-1][seq_idx],
            vi[state_idx-1][seq_idx],
            vd[state_idx-1][seq_idx]});
          switch (max_idx) {
            case 0:
              state = 'M';
            case 1:
              state = 'D';
            case 2:
              ;
          }
          --state_idx;
        }
    }
  }
  trace.push_back('M', 0);
  reverse(trace.begin(), trace.end());

  return max(vd[model_len][seq_len] + ending[index_d[model_len]],
      vi[0][seq_len] + ending[index_i(0)];
}
