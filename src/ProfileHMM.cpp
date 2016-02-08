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
#include <iostream>

using std::vector;
using std::max;
using std::pair;
using std::endl;
using std::cout;

inline double
ProfileHMM::log_sum_log(const double p, const double q) const {
  if (p == 0) {return q;}
  else if (q == 0) {return p;}
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
}

double
ProfileHMM::max_item(const std::initializer_list<double> list) const {
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
ProfileHMM::ViterbiDecoding(const vector<vector<double> > &transition,
    const vector<vector<double> > &emission,
    const vector<double> &initial,
    const vector<int> &observation,
    vector<pair<char, size_t> > &trace) const {

  // use a very small value as -Inf
  const double LOG_ZERO = -1000.0;
  const size_t total_size = model_len * 3 + 2;
  // Note here: vm[i][j] corresponds to M_i on o_j, where i=0~L
  // vd[i][j]: i = 0~L-1 corresponds to D_1 ~ D_L
  // vi[i][j]: i = 0~L-1 corresponds to I_0 ~ I_L-1
  // crazy subscriptions!
  vector<vector<double> > vm(model_len+1, vector<double>(seq_len+1, LOG_ZERO));
  vector<vector<double> > vi(model_len, vector<double>(seq_len+1, LOG_ZERO));
  vector<vector<double> > vd(model_len, vector<double>(seq_len+1, LOG_ZERO));

  // Initialization for vm, vi, vd
  // vm[j][i]: prob of given state M_j and observing sequence X_i
  // similarly for vi and vd.
  vm[0][0] = 0;

  // I_0 is the background emission state
  //vi[0][1] = transition[index_m(0)][index_i(0)];
  //for (size_t i = 2; i < seq_len; ++i) {
  //  vi[0][i] = vi[0][i-1] + transition[index_i(0)][index_i(0)];
  //}
  // D_1 (vd[0]) is the leading deletion state for local alignment
  vd[0][0] = transition[index_m(0)][index_d(1)];
  //for (size_t i = 1; i < seq_len - 1; ++i) {
  //  vd[1-1][i] = vi[0][i-1] + transition[index_i(1)][index_d(1)];
  //}

  // some special cases before going to the main loop
  //for (size_t i = 1; i < seq_len; ++i) {
  //  // D_2 only has two incoming transition from I_1 and M_1
  //  vd[2][i] = max( {vm[1][i] + transition[index_m(1)][index_d(2)],
  //    vi[1][i] + transition[index_i(1)][index_d(2)]});
  //}

  // main loop
  for (size_t i = 1; i < seq_len + 1; ++i) {
    // M_1 only has one incoming transition from D_1
    vm[1][i] = vd[0][i-1] + transition[index_d(1)][index_m(1)]
      + emission[index_m(1)][observation[i-1]]
      - emission[index_i(0)][observation[i-1]];
    // I_1 only has two incoming transition including one from M_1
    vi[1][i] = max(
      {vm[1][i-1] + transition[index_m(1)][index_i(1)],
      vi[1][i-1] + transition[index_i(1)][index_i(1)]})
      + emission[index_i(1)][observation[i-1]]
      - emission[index_i(0)][observation[i-1]];
    for (size_t j = 2; j < model_len; ++j) {
      // M_i
      vm[j][i] = max(
        {vm[j-1][i-1] + transition[index_m(j-1)][index_m(j)],
        vi[j-1][i-1] + transition[index_i(j-1)][index_m(j)],
        vd[j-2][i-1] + transition[index_d(j-1)][index_m(j)],
        vd[0][i-1] + transition[index_d(1)][index_m(j)]})
        + emission[index_m(j)][observation[i-1]]
        - emission[index_i(0)][observation[i-1]];
      // I_i
      vi[j][i] = max(
        {vm[j][i-1] + transition[index_m(j)][index_i(j)],
        vi[j][i-1] + transition[index_i(j)][index_i(j)],
        vd[j-1][i-1] + transition[index_d(j)][index_i(j)]})
        + emission[index_i(j)][observation[i-1]]
        - emission[index_i(0)][observation[i-1]];
      // D_i
      vd[j-1][i] = max( 
        {vm[j-1][i] + transition[index_m(j-1)][index_d(j)],
        vi[j-1][i] + transition[index_i(j-1)][index_d(j)],
        vd[j-2][i] + transition[index_d(j-1)][index_d(j)]});
    }
    // M_model_len
    vm[model_len][i] = max( 
      {vm[model_len-1][i-1]
        + transition[index_m(model_len-1)][index_m(model_len)],
      vi[model_len-1][i-1]
        + transition[index_i(model_len-1)][index_m(model_len)],
      vd[model_len-2][i-1]
        + transition[index_d(model_len-1)][index_m(model_len)],
      vd[0][i-1] + transition[index_d(1)][index_m(model_len)]})
      + emission[index_m(model_len)][observation[i-1]]
      - emission[index_i(0)][observation[i-1]];
    // I_model_len does not exist
    // D_model_len
    vector<double> list;
    for (size_t k = 1; k < model_len; ++k) {
      list.push_back(vm[k][i]
        + transition[index_m(k)][index_d(model_len)]);
    }
    const vector<double>::const_iterator x = 
      std::max_element(list.begin(), list.end());
    vd[model_len-1][i] = *x;
    // I_0
    vi[0][i] = max( 
      {vi[0][i-1] + transition[index_i(0)][index_i(0)],
      vd[model_len-1][i-1] + transition[index_d(model_len)][index_i(0)]});
    // D_1
    vd[0][i] = vi[0][i] + transition[index_i(0)][index_d(1)];
  }

  if (VERBOSE) {
    cout << "V_M" << endl;
    for (size_t j = 0; j < seq_len + 1; ++j)
      cout << "\t" << j;
    cout << endl;
    for (vector<vector<double> >::const_iterator i = vm.begin(); i < vm.end(); ++i) {
      cout << i - vm.begin();
      for (size_t j = 0; j < seq_len +1; ++j)
        printf("\t%.4f", exp((*i)[j]));
        //printf("\t%.4f", (*i)[j]);
      cout << endl;
    }
    cout << endl << "V_I" << std::endl;
    for (size_t j = 0; j < seq_len + 1; ++j)
      cout << "\t" << j;
    cout << endl;
    for (vector<vector<double> >::const_iterator i = vi.begin(); i < vi.end(); ++i) {
      cout << i - vi.begin();
      for (size_t j = 0; j < seq_len +1; ++j)
        printf("\t%.4f", exp((*i)[j]));
        //printf("\t%.4f", (*i)[j]);
      cout << endl;
    }
    cout << endl << "V_D" << std::endl;
    for (size_t j = 0; j < seq_len + 1; ++j)
      cout << "\t" << j;
    cout << endl;
    for (vector<vector<double> >::const_iterator i = vd.begin(); i < vd.end(); ++i) {
      cout << i - vd.begin();
      for (size_t j = 0; j < seq_len +1; ++j)
        printf("\t%.4f", exp((*i)[j]));
        //printf("\t%.4f", (*i)[j]);
      cout << endl;
    }
  }

  //traceback
  std::pair<char, int> step('E', 0);
  trace.push_back(step);
  size_t state_idx;
  char state;
  size_t seq_idx = seq_len;
  if (vd[model_len-1][seq_len] + transition[index_d(model_len)][total_size-1]
      > vi[0][seq_len] + transition[index_i(0)][total_size-1]) {
    state = 'D';
    state_idx = model_len;
  }
  else {
    state = 'I';
    state_idx = 0;
  }

  while (!(state == 'M' && state_idx == 0) && seq_idx <= seq_len) {
    std::pair<char, int> step(state, state_idx);
    trace.push_back(step);
    if (VERBOSE)
      cout << state << state_idx << ", " << seq_idx << endl;
    switch (state) {
      case 'M': {
        // assuming state_idx = 1~L
        size_t max_idx = max_item({
          vm[state_idx-1][seq_idx-1] 
            + transition[index_m(state_idx-1)][index_m(state_idx)],
          vi[state_idx-1][seq_idx-1]
            + transition[index_i(state_idx-1)][index_m(state_idx)],
          vd[state_idx-2][seq_idx-1]
            + transition[index_d(state_idx-1)][index_m(state_idx)],
          vd[0][seq_idx-1]
            + transition[index_d(1)][index_m(state_idx)]});
        switch (max_idx) {
          case 0:
            --state_idx;
          case 1:
            state = 'I';
            --state_idx;
          case 2:
            state = 'D';
            --state_idx;
          case 3:
            state = 'D';
            state_idx = 1;
        }
        --seq_idx;
        break;
      }
      case 'I': {
        if (state_idx == 0) {
          size_t max_idx = max_item({
            vm[0][seq_idx-1]
              + transition[index_m(0)][index_i(0)],
            vi[0][seq_idx-1]
              + transition[index_i(0)][index_i(0)],
            vd[model_len-1][seq_idx-1]
              + transition[index_d(model_len)][index_i(0)],
            });
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
          size_t max_idx = max_item({
            vm[state_idx][seq_idx-1]
              + transition[index_m(state_idx)][index_i(state_idx)],
            vi[state_idx][seq_idx-1]
              + transition[index_i(state_idx)][index_i(state_idx)],
            vd[state_idx-1][seq_idx-1]
              + transition[index_d(state_idx)][index_i(state_idx)]});
          switch (max_idx) {
            case 0:
              state = 'M';
            case 1:
              ;
            case 2:
              state = 'D';
          }
        }
        --seq_idx;
        break;
      }
      case 'D': {
        if (state_idx == model_len) {
          vector<double> list;
          for (size_t k = 1; k < model_len + 1; ++k) {
            list.push_back(vm[k][seq_idx]
                + transition[index_m(k)][index_d(model_len)]);
          }
          const vector<double>::const_iterator x = 
            std::max_element(list.begin(), list.end());
          state_idx = x - list.begin() + 1;
          state = 'M';
        }
        else if (state_idx == 1) {
          size_t max_idx = max_item({
            vm[0][seq_idx]
              + transition[index_m(0)][index_d(1)],
            vi[0][seq_idx]
              + transition[index_i(0)][index_d(1)]});
          switch (max_idx) {
            case 0:
              state = 'M';
            case 1:
              state = 'I';
          }
          state_idx = 0;
        }
        else {
          size_t max_idx = max_item({
            vm[state_idx-1][seq_idx]
              + transition[index_m(state_idx-1)][index_d(state_idx)],
            vi[state_idx-1][seq_idx],
              + transition[index_i(state_idx-1)][index_d(state_idx)],
            vd[state_idx-2][seq_idx]
              + transition[index_d(state_idx-1)][index_d(state_idx)]});
          switch (max_idx) {
            case 0:
              state = 'M';
            case 1:
              state = 'I';
            case 2:
              ;
          }
          --state_idx;
        }
        break;
      }
    }
  }
  //step = std::make_pair ('M', 0);
  //trace.push_back(step);
  reverse(trace.begin(), trace.end());

  return max(vd[model_len-1][seq_len]
      + transition[index_d(model_len)][total_size-1],
      vi[0][seq_len] + transition[index_i(0)][total_size-1]);
}
