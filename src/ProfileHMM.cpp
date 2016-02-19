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

inline double
ProfileHMM::max_item(const std::initializer_list<double> list) const {
  const vector<double> v(list);
  const vector<double>::const_iterator x = 
     std::max_element(v.begin(), v.end());
  return x - v.begin();
}

inline double
ProfileHMM::log_sum_log_list(const std::initializer_list<double> list) const {
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

inline double
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

inline size_t
ProfileHMM::index_m(const size_t idx) const {
  // M_0 ~ M_L
  assert(idx >= 0 && idx <= model_len);
  return idx;
}

inline size_t
ProfileHMM::index_i(const size_t idx) const {
  // I_0 ~ I_L-1 + (L+1)*M
  assert(idx >= 0 && idx <= model_len - 1);
  return idx + model_len + 1;
}

inline size_t
ProfileHMM::index_d(const size_t idx) const {
  // D_1 ~ D_L + (L+1)*M + L*I
  assert(idx >= 1 && idx <= model_len);
  return idx + model_len * 2;
}

double
ProfileHMM::ViterbiDecoding(const bool VERBOSE,
    const vector<vector<double> > &transition,
    const vector<vector<double> > &emission,
    const vector<int> &observation,
    vector<pair<char, size_t> > &trace) const {
  const size_t seq_len = observation.size();
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

  // D_1 (vd[0]) is the leading deletion state for local alignment
  vd[0][0] = transition[index_m(0)][index_d(1)];

  // main loop
  for (size_t i = 1; i < seq_len + 1; ++i) {
    // M_1 only has one incoming transition from D_1
    vm[1][i] = vd[0][i-1] + transition[index_d(1)][index_m(1)]
      + emission[index_m(1)][observation[i-1]]
      - emission[index_i(0)][observation[i-1]];
    // I_1 only has two incoming transition including one from M_1
    vi[1][i] = max({
      vm[1][i-1] + transition[index_m(1)][index_i(1)],
      vi[1][i-1] + transition[index_i(1)][index_i(1)]})
      + emission[index_i(1)][observation[i-1]]
      - emission[index_i(0)][observation[i-1]];
    for (size_t j = 2; j < model_len; ++j) {
      // M_i
      vm[j][i] = max({
        vm[j-1][i-1] + transition[index_m(j-1)][index_m(j)],
        vi[j-1][i-1] + transition[index_i(j-1)][index_m(j)],
        vd[j-2][i-1] + transition[index_d(j-1)][index_m(j)],
        vd[0][i-1] + transition[index_d(1)][index_m(j)]})
        + emission[index_m(j)][observation[i-1]]
        - emission[index_i(0)][observation[i-1]];
      // I_i
      vi[j][i] = max({
        vm[j][i-1] + transition[index_m(j)][index_i(j)],
        vi[j][i-1] + transition[index_i(j)][index_i(j)],
        vd[j-1][i-1] + transition[index_d(j)][index_i(j)]})
        + emission[index_i(j)][observation[i-1]]
        - emission[index_i(0)][observation[i-1]];
      // D_i
      if (j > 2)
        vd[j-1][i] = max({
          vm[j-1][i] + transition[index_m(j-1)][index_d(j)],
          vi[j-1][i] + transition[index_i(j-1)][index_d(j)],
          vd[j-2][i] + transition[index_d(j-1)][index_d(j)]});
      else
        vd[j-1][i] = max({
          vm[j-1][i] + transition[index_m(j-1)][index_d(j)],
          vi[j-1][i] + transition[index_i(j-1)][index_d(j)]});
    }
    // M_model_len
    vm[model_len][i] = max({
      vm[model_len-1][i-1]
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
    vi[0][i] = max({
      vi[0][i-1] + transition[index_i(0)][index_i(0)],
      vd[model_len-1][i-1] + transition[index_d(model_len)][index_i(0)],
      vm[0][i-1] + transition[index_m(0)][index_i(0)]});
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
    cout << endl << "V_I" << endl;
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
    cout << endl << "V_D" << endl;
    for (size_t j = 0; j < seq_len + 1; ++j)
      cout << "\t" << j;
    cout << endl;
    for (vector<vector<double> >::const_iterator i = vd.begin(); i < vd.end(); ++i) {
      cout << i - vd.begin() + 1;
      for (size_t j = 0; j < seq_len +1; ++j)
        printf("\t%.4f", exp((*i)[j]));
        //printf("\t%.4f", (*i)[j]);
      cout << endl;
    }
  }

  //traceback
  std::pair<char, size_t> step('E', 0);
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
    std::pair<char, size_t> step(state, state_idx);
    trace.push_back(step);
    if (VERBOSE)
      cout << state << state_idx << ", " << seq_idx << endl;
    if (state == 'M') {
      // state_idx = 2~L
      if (state_idx > 1) {
        size_t max_idx = max_item({
          vm[state_idx-1][seq_idx-1] 
            + transition[index_m(state_idx-1)][index_m(state_idx)],
          vi[state_idx-1][seq_idx-1]
            + transition[index_i(state_idx-1)][index_m(state_idx)],
          vd[state_idx-2][seq_idx-1]
            + transition[index_d(state_idx-1)][index_m(state_idx)],
          vd[0][seq_idx-1]
            + transition[index_d(1)][index_m(state_idx)]});
        if (max_idx == 0) {
          --state_idx;
        } else if (max_idx == 1) {
          state = 'I';
          --state_idx;
        } else if (max_idx == 2) {
          state = 'D';
          --state_idx;
        } else if (max_idx == 3) {
          state = 'D';
          state_idx = 1;
        }
      } else if (state_idx == 1) {
        state = 'D';
        state_idx = 1;
      }
      --seq_idx;
    }
    else if (state == 'I') {
      if (state_idx == 0) {
        size_t max_idx = max_item({
          vm[0][seq_idx-1]
            + transition[index_m(0)][index_i(0)],
          vi[0][seq_idx-1]
            + transition[index_i(0)][index_i(0)],
          vd[model_len-1][seq_idx-1]
            + transition[index_d(model_len)][index_i(0)],
          });
        if (max_idx == 0) {
          state = 'M';
          state_idx = 0;
        } else if (max_idx == 1) {
          ;
        } else if (max_idx == 2) {
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
        if (max_idx == 0) {
          state = 'M';
        } else if (max_idx == 1) {
          ;
        } else if (max_idx == 2) {
          state = 'D';
        }
      }
      --seq_idx;
    }
    else if (state == 'D') {
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
        if (max_idx == 0) {
          state = 'M';
        } else if (max_idx == 1) {
          state = 'I';
        }
        state_idx = 0;
      }
      else {
        size_t max_idx = max_item({
          vm[state_idx-1][seq_idx]
            + transition[index_m(state_idx-1)][index_d(state_idx)],
          vi[state_idx-1][seq_idx]
            + transition[index_i(state_idx-1)][index_d(state_idx)],
          vd[state_idx-2][seq_idx]
            + transition[index_d(state_idx-1)][index_d(state_idx)]});
        if (max_idx == 0)
          state = 'M';
        else if (max_idx == 1)
          state = 'I';
        --state_idx;
      }
    }
  }
  step = std::make_pair (state, state_idx);
  trace.push_back(step);
  reverse(trace.begin(), trace.end());

  return max(vd[model_len-1][seq_len]
      + transition[index_d(model_len)][total_size-1],
      vi[0][seq_len] + transition[index_i(0)][total_size-1]);
}

void
ProfileHMM::forward_algorithm(const bool VERBOSE,
    const vector<vector<double> > &transition,
    const vector<vector<double> > &emission,
    const vector<int> &observation) {
  const size_t seq_len = observation.size();
  // Note here: fm[i][j] corresponds to M_i on o_j, where i=0~L
  // fd[i][j]: i = 0~L-1 corresponds to D_1 ~ D_L
  // fi[i][j]: i = 0~L-1 corresponds to I_0 ~ I_L-1
  // crazy subscriptions!
  fm.resize(model_len+1, vector<double>(seq_len+1, LOG_ZERO));
  fi.resize(model_len, vector<double>(seq_len+1, LOG_ZERO));
  fd.resize(model_len, vector<double>(seq_len+1, LOG_ZERO));

  // Initialization for fm, fi, fd
  // fm[j][i]: prob of given state M_j and observing sequence up to X_i
  // similarly for vi and vd.
  fm[0][0] = 0;

  // D_1 (fd[0]) is the leading deletion state for local alignment
  fd[0][0] = transition[index_m(0)][index_d(1)];

  // main loop
  for (size_t i = 1; i < seq_len + 1; ++i) {
    // M_1 only has one incoming transition from D_1
    fm[1][i] = fd[0][i-1] + transition[index_d(1)][index_m(1)]
      + emission[index_m(1)][observation[i-1]]
      - emission[index_i(0)][observation[i-1]];
    // I_1 only has two incoming transition including one from M_1
    fi[1][i] = log_sum_log_list(
      {fm[1][i-1] + transition[index_m(1)][index_i(1)],
      fi[1][i-1] + transition[index_i(1)][index_i(1)]})
      + emission[index_i(1)][observation[i-1]]
      - emission[index_i(0)][observation[i-1]];
    for (size_t j = 2; j < model_len; ++j) {
      // M_i
      fm[j][i] = log_sum_log_list({
        fm[j-1][i-1] + transition[index_m(j-1)][index_m(j)],
        fi[j-1][i-1] + transition[index_i(j-1)][index_m(j)],
        fd[j-2][i-1] + transition[index_d(j-1)][index_m(j)],
        fd[0][i-1] + transition[index_d(1)][index_m(j)]})
        + emission[index_m(j)][observation[i-1]]
        - emission[index_i(0)][observation[i-1]];
      // I_i
      fi[j][i] = log_sum_log_list({
        fm[j][i-1] + transition[index_m(j)][index_i(j)],
        fi[j][i-1] + transition[index_i(j)][index_i(j)],
        fd[j-1][i-1] + transition[index_d(j)][index_i(j)]})
        + emission[index_i(j)][observation[i-1]]
        - emission[index_i(0)][observation[i-1]];
      // D_i
      if (j > 2)
        fd[j-1][i] = log_sum_log_list({
          fm[j-1][i] + transition[index_m(j-1)][index_d(j)],
          fi[j-1][i] + transition[index_i(j-1)][index_d(j)],
          fd[j-2][i] + transition[index_d(j-1)][index_d(j)]});
      else
        fd[j-1][i] = log_sum_log_list({
          fm[j-1][i] + transition[index_m(j-1)][index_d(j)],
          fi[j-1][i] + transition[index_i(j-1)][index_d(j)]});
    }
    // M_model_len
    fm[model_len][i] = log_sum_log_list({
      fm[model_len-1][i-1]
        + transition[index_m(model_len-1)][index_m(model_len)],
      fi[model_len-1][i-1]
        + transition[index_i(model_len-1)][index_m(model_len)],
      fd[model_len-2][i-1]
        + transition[index_d(model_len-1)][index_m(model_len)],
      fd[0][i-1]
        + transition[index_d(1)][index_m(model_len)]})
      + emission[index_m(model_len)][observation[i-1]]
      - emission[index_i(0)][observation[i-1]];
    // I_model_len does not exist
    // D_model_len
    vector<double> list;
    for (size_t k = 1; k < model_len; ++k) {
      list.push_back(fm[k][i]
        + transition[index_m(k)][index_d(model_len)]);
    }
    fd[model_len-1][i] = log_sum_log_vec(list, list.size());
    // I_0
    fi[0][i] = log_sum_log_list({
      fi[0][i-1] + transition[index_i(0)][index_i(0)],
      fd[model_len-1][i-1] + transition[index_d(model_len)][index_i(0)],
      fm[0][i-1] + transition[index_m(0)][index_i(0)]});
    // D_1
    fd[0][i] = fi[0][i] + transition[index_i(0)][index_d(1)];
  }

  if (VERBOSE) {
    cout << "F_M" << endl;
    for (size_t j = 0; j < seq_len + 1; ++j)
      cout << "\t" << j;
    cout << endl;
    for (vector<vector<double> >::const_iterator i = fm.begin(); i < fm.end(); ++i) {
      cout << i - fm.begin();
      for (size_t j = 0; j < seq_len +1; ++j)
        printf("\t%.4f", exp((*i)[j]));
        //printf("\t%.4f", (*i)[j]);
      cout << endl;
    }
    cout << endl << "F_I" << endl;
    for (size_t j = 0; j < seq_len + 1; ++j)
      cout << "\t" << j;
    cout << endl;
    for (vector<vector<double> >::const_iterator i = fi.begin(); i < fi.end(); ++i) {
      cout << i - fi.begin();
      for (size_t j = 0; j < seq_len +1; ++j)
        printf("\t%.4f", exp((*i)[j]));
        //printf("\t%.4f", (*i)[j]);
      cout << endl;
    }
    cout << endl << "F_D" << endl;
    for (size_t j = 0; j < seq_len + 1; ++j)
      cout << "\t" << j;
    cout << endl;
    for (vector<vector<double> >::const_iterator i = fd.begin(); i < fd.end(); ++i) {
      cout << i - fd.begin() + 1;
      for (size_t j = 0; j < seq_len +1; ++j)
        printf("\t%.4f", exp((*i)[j]));
        //printf("\t%.4f", (*i)[j]);
      cout << endl;
    }
  }

}

double
ProfileHMM::forward_prob(const char state, const size_t state_idx,
    const size_t obs_pos) const {
  if (state == 'M')
    return fm[state_idx][obs_pos];
  else if (state == 'I')
    return fi[state_idx][obs_pos];
  else if (state == 'D')
    return fd[state_idx-1][obs_pos];
  else
    return log_sum_log_list({
      fi[0][fi[1].size()-1],
      fd[model_len-1][fd[1].size()-1]});
  return -1;
}

void
ProfileHMM::backward_algorithm(const bool VERBOSE,
    const vector<vector<double> > &transition,
    const vector<vector<double> > &emission,
    const vector<int> &observation) {
  const size_t seq_len = observation.size();
  bm.resize(model_len+1, vector<double>(seq_len+1, LOG_ZERO));
  bi.resize(model_len, vector<double>(seq_len+1, LOG_ZERO));
  bd.resize(model_len, vector<double>(seq_len+1, LOG_ZERO));

  // Initialization
  for (size_t i = 0; i < bm.size(); ++i)
    bm[i][seq_len] = 0.0;
  for (size_t i = 0; i < bi.size(); ++i)
    bi[i][seq_len] = 0.0;
  for (size_t i = 0; i < bd.size(); ++i)
    bd[i][seq_len] = 0.0;

  for (signed long i = seq_len - 1; i >= 0; --i) {
    // D_L
    bd[model_len-1][i] = 
      bi[0][i+1] + transition[index_d(model_len)][index_i(0)];
    // D_1
    vector<double> list;
    for (size_t k = 1; k <= model_len; ++k)
      list.push_back(bm[k][i+1] + transition[index_d(1)][index_m(k)]
        + emission[index_m(k)][observation[i]]
        - emission[index_i(0)][observation[i]]);
    bd[0][i] = log_sum_log_vec(list, list.size());
    // D_L-1
    bd[model_len-2][i] = log_sum_log(
      bm[model_len][i+1]
        + transition[index_d(model_len-1)][index_m(model_len)]
        + emission[index_m(model_len)][observation[i]]
        - emission[index_i(0)][observation[i]],
      bi[model_len-1][i+1]
        + transition[index_d(model_len-1)][index_i(model_len-1)]
        + emission[index_i(model_len-1)][observation[i]]
        - emission[index_i(0)][observation[i]]);
    // I_L-1
    bi[model_len-1][i] = log_sum_log(
      bi[model_len-1][i+1]
        + transition[index_i(model_len-1)][index_i(model_len-1)]
        + emission[index_i(model_len-1)][observation[i]]
        - emission[index_i(0)][observation[i]],
      bm[model_len][i+1]
        + transition[index_i(model_len-1)][index_m(model_len)]
        + emission[index_m(model_len)][observation[i]]
        - emission[index_i(0)][observation[i]]);
    // M_L
    bm[model_len][i] = bd[model_len-1][i]
      + transition[index_m(model_len)][index_d(model_len)];
    // M_L-1
    bm[model_len-1][i] = log_sum_log_list({
      bm[model_len][i+1]
        + transition[index_m(model_len-1)][index_m(model_len)]
        + emission[index_m(model_len)][observation[i]]
        - emission[index_i(0)][observation[i]],
      bi[model_len-1][i+1]
        + transition[index_m(model_len-1)][index_i(model_len-1)]
        + emission[index_i(model_len-1)][observation[i]]
        - emission[index_i(0)][observation[i]],
      bd[model_len-1][i]
        + transition[index_m(model_len-1)][index_d(model_len)]});

    // general loop
    for (size_t j = model_len - 2; j > 0; --j) {
      // D_i, i=2~L-2
      if (j > 1) {
        bd[j-1][i] = log_sum_log_list({
          bm[j+1][i+1]
            + transition[index_d(j)][index_m(j+1)]
            + emission[index_m(j+1)][observation[i]]
            - emission[index_i(0)][observation[i]],
          bi[j][i+1]
            + transition[index_d(j)][index_i(j)]
            + emission[index_i(j)][observation[i]]
            - emission[index_i(0)][observation[i]],
          bd[j][i]
            + transition[index_d(j)][index_d(j+1)]});
      }

      // I_i, i=1~L-2
      bi[j][i] = log_sum_log_list({
        bi[j][i+1]
          + transition[index_i(j)][index_i(j)]
          + emission[index_i(j)][observation[i]]
          - emission[index_i(0)][observation[i]],
        bm[j+1][i+1]
          + transition[index_i(j)][index_m(j+1)]
          + emission[index_m(j+1)][observation[i]]
          - emission[index_i(0)][observation[i]],
        bd[j+1][i]
          + transition[index_i(j)][index_d(j+1)]});

      // M_i, i=1~L-2
      bm[j][i] = log_sum_log_list({
        bi[j][i+1]
          + transition[index_m(j)][index_i(j)]
          + emission[index_i(j)][observation[i]]
          - emission[index_i(0)][observation[i]],
        bm[j+1][i+1]
          + transition[index_m(j)][index_m(j+1)]
          + emission[index_m(j+1)][observation[i]]
          - emission[index_i(0)][observation[i]],
        bd[j+1][i]
          + transition[index_m(j)][index_d(j+1)],
        bd[model_len-1][i]
          + transition[index_m(j)][index_d(model_len)]});
    }

    // I_0
    bi[0][i] = log_sum_log(
      bd[0][i] + transition[index_i(0)][index_d(1)],
      bi[0][i+1] + transition[index_i(0)][index_i(0)]);

    // M_0
    bm[0][i] = log_sum_log(
      bd[0][i] + transition[index_m(0)][index_d(1)],
      bi[0][i+1] + transition[index_m(0)][index_i(0)]);
  }

  if (VERBOSE) {
    cout << "B_M" << endl;
    for (size_t j = 0; j < seq_len+1; ++j)
      cout << "\t" << j;
    cout << endl;
    for (vector<vector<double> >::const_iterator i = bm.begin(); i < bm.end(); ++i) {
      cout << i - bm.begin();
      for (size_t j = 0; j < seq_len+1; ++j)
        printf("\t%.4f", exp((*i)[j]));
        //printf("\t%.4f", (*i)[j]);
      cout << endl;
    }
    cout << endl << "B_I" << endl;
    for (size_t j = 0; j < seq_len+1; ++j)
      cout << "\t" << j;
    cout << endl;
    for (vector<vector<double> >::const_iterator i = bi.begin(); i < bi.end(); ++i) {
      cout << i - bi.begin();
      for (size_t j = 0; j < seq_len+1; ++j)
        printf("\t%.4f", exp((*i)[j]));
        //printf("\t%.4f", (*i)[j]);
      cout << endl;
    }
    cout << endl << "B_D" << endl;
    for (size_t j = 0; j < seq_len+1; ++j)
      cout << "\t" << j;
    cout << endl;
    for (vector<vector<double> >::const_iterator i = bd.begin(); i < bd.end(); ++i) {
      cout << i - bd.begin() + 1;
      for (size_t j = 0; j < seq_len+1; ++j)
        printf("\t%.4f", exp((*i)[j]));
        //printf("\t%.4f", (*i)[j]);
      cout << endl;
    }
  }
}

double
ProfileHMM::backward_prob(const char state, const size_t state_idx,
    const size_t obs_pos) const {
  if (state == 'M')
    return bm[state_idx][obs_pos];
  else if (state == 'I')
    return bi[state_idx][obs_pos];
  else if (state == 'D')
    return bd[state_idx-1][obs_pos];
  else
    return -1;
}

void
ProfileHMM::BW_training(const bool VERBOSE,
    vector<vector<double> > &transition,
    vector<vector<double> > &emission,
    const vector<int> &observation) {

  const size_t seq_len = observation.size();
  double ll_new = 1e10;
  double ll = forward_prob('E', 1, 1);
  vector<vector<double> > e_trans(total_size,
      vector<double>(total_size, LOG_ZERO));
  vector<vector<double> > e_emiss(total_size,
      vector<double>(emission[0].size(), LOG_ZERO));
  for (size_t itr = 0; ll_new - ll > tolerance && itr < max_iterations;
      ++itr) {
    // get expectations
    // Transitions
    // M_0 to D_1 and I_0
    e_trans[index_m(0)][index_d(1)] = 
        fm[0][0] + transition[index_m(0)][index_d(1)]
        + bd[0][0] - ll;
    e_trans[index_m(0)][index_i(0)] = 
        fm[0][0] + transition[index_m(0)][index_i(0)]
        + bi[0][0] - ll;
    // I_0 to D_1, I_0, and E
    vector<double> list;
    for (size_t i = 1; i < seq_len; ++i) {
      list.push_back(fi[0][i] + transition[index_i(0)][index_d(1)]
        + bd[0][i]);
    }
    e_trans[index_i(0)][index_d(1)] = log_sum_log_vec(list, list.size())
      - ll;
    list.clear();
    for (size_t i = 1; i < seq_len; ++i) {
      // still use log-odds ratio here (no emission)?
      list.push_back(fi[0][i] + transition[index_i(0)][index_i(0)]
        + bi[0][i+1]);
    }
    e_trans[index_i(0)][index_i(0)] = log_sum_log_vec(list, list.size())
      - ll;
    list.clear();
    e_trans[index_i(0)][total_size-1] =
      fi[0][seq_len] + transition[index_i(0)][total_size-1] - ll;
    // D_1 to M_i
    for (size_t j = 1; j <= model_len; ++j) {
      for (size_t i = 0; i < seq_len; ++i) {
        list.push_back(fd[0][i] + transition[index_d(1)][index_m(j)]
          + emission[index_m(j)][observation[i]]
          - emission[index_i(0)][observation[i]]
          + bm[j][i+1]);
      }
      e_trans[index_d(1)][index_m(j)] = log_sum_log_vec(list, list.size()) - ll;
      list.clear();
    }
    // D_L to I_0 and E
    for (size_t i = 0; i < seq_len; ++i) {
      list.push_back(fd[model_len-1][i]
        + transition[index_d(model_len)][index_i(0)]
        + bi[0][i+1]);
    }
    e_trans[index_d(model_len)][index_i(0)] =
      log_sum_log_vec(list, list.size()) - ll;
    list.clear();
    e_trans[index_d(model_len)][total_size-1] =
      fd[model_len-1][seq_len] + transition[index_d(model_len)][total_size-1]
      - ll;
    // D_L-1 to M_L and I_L-1
    for (size_t i = 0; i < seq_len; ++i) {
      list.push_back(fd[model_len-2][i]
        + transition[index_d(model_len-1)][index_m(model_len)]
        + emission[index_m(model_len)][observation[i]]
        - emission[index_i(0)][observation[i]]
        + bm[model_len][i+1]);
    }
    e_trans[index_d(model_len-1)][index_m(model_len)] =
      log_sum_log_vec(list, list.size()) - ll;
    list.clear();
    for (size_t i = 0; i < seq_len; ++i) {
      list.push_back(fd[model_len-2][i]
        + transition[index_d(model_len-1)][index_i(model_len-1)]
        + emission[index_i(model_len-1)][observation[i]]
        - emission[index_i(0)][observation[i]]
        + bi[model_len-1][i+1]);
    }
    e_trans[index_d(model_len-1)][index_i(model_len-1)] =
      log_sum_log_vec(list, list.size()) - ll;
    list.clear();
    // I_L-1 to I_L-1 and M_L-1
    for (size_t i = 1; i < seq_len; ++i) {
      list.push_back(fi[model_len-1][i]
        + transition[index_i(model_len-1)][index_i(model_len-1)]
        + emission[index_i(model_len-1)][observation[i]]
        - emission[index_i(0)][observation[i]]
        + bi[model_len-1][i+1]);
    }
    e_trans[index_i(model_len-1)][index_i(model_len-1)] = 
      log_sum_log_vec(list, list.size()) - ll;
    list.clear();
    for (size_t i = 1; i < seq_len; ++i) {
      list.push_back(fi[model_len-1][i]
        + transition[index_i(model_len-1)][index_m(model_len)]
        + emission[index_m(model_len)][observation[i]]
        - emission[index_i(0)][observation[i]]
        + bm[model_len][i+1]);
    }
    e_trans[index_i(model_len-1)][index_m(model_len)] =
      log_sum_log_vec(list, list.size()) - ll;
    list.clear();
    // M_L-1 to I_L-1, M_L and D_L
    for (size_t i = 1; i < seq_len; ++i) {
      list.push_back(fm[model_len-1][i]
        + transition[index_m(model_len-1)][index_i(model_len-1)]
        + emission[index_i(model_len-1)][observation[i]]
        - emission[index_i(0)][observation[i]]
        + bi[model_len-1][i+1]);
    }
    e_trans[index_m(model_len-1)][index_i(model_len-1)] = 
      log_sum_log_vec(list, list.size()) - ll;
    list.clear();
    for (size_t i = 1; i < seq_len; ++i) {
      list.push_back(fm[model_len-1][i]
        + transition[index_m(model_len-1)][index_m(model_len)]
        + emission[index_m(model_len)][observation[i]]
        - emission[index_i(0)][observation[i]]
        + bm[model_len][i+1]);
    }
    e_trans[index_m(model_len-1)][index_m(model_len)] = 
      log_sum_log_vec(list, list.size()) - ll;
    list.clear();
    for (size_t i = 1; i <= seq_len; ++i) {
      list.push_back(fm[model_len-1][i]
        + transition[index_m(model_len-1)][index_d(model_len)]
        + bd[model_len-1][i]);
    }
    e_trans[index_m(model_len-1)][index_d(model_len)] = 
      log_sum_log_vec(list, list.size()) - ll;
    list.clear();
    // M_L to D_L is always 1, so no need to learn
    // general
    for (size_t j = 1; j < model_len - 1; ++j) {
      // I_1~L-2: to I_j, M_j+1, and D_j+1
      vector <double> list_m, list_i, list_d, list_l;
      for (size_t i = 1; i < seq_len; ++i) {
        list_i.push_back(fi[j][i]
          + transition[index_i(j)][index_i(j)]
          + emission[index_i(j)][observation[i]]
          - emission[index_i(0)][observation[i]]
          + bi[j][i+1]);
        list_m.push_back(fi[j][i]
          + transition[index_i(j)][index_m(j+1)]
          + emission[index_m(j+1)][observation[i]]
          - emission[index_i(0)][observation[i]]
          + bm[j+1][i+1]);
        list_d.push_back(fi[j][i]
          + transition[index_i(j)][index_d(j+1)]
          + bd[j][i]);
      }
      e_trans[index_i(j)][index_i(j)] =
        log_sum_log_vec(list_i, list_i.size()) - ll;
      e_trans[index_i(j)][index_m(j+1)] =
        log_sum_log_vec(list_m, list_m.size()) - ll;
      e_trans[index_i(j)][index_d(j+1)] =
        log_sum_log_vec(list_d, list_d.size()) - ll;
      list_i.clear();
      list_m.clear();
      list_d.clear();
      // M_1~L-2: to M_j+1, I_j, D_j+1 and D_L
      for (size_t i = 1; i < seq_len; ++i) {
        list_i.push_back(fm[j][i]
          + transition[index_m(j)][index_i(j)]
          + emission[index_i(j)][observation[i]]
          - emission[index_i(0)][observation[i]]
          + bi[j][i+1]);
        list_m.push_back(fm[j][i]
          + transition[index_m(j)][index_m(j+1)]
          + emission[index_m(j+1)][observation[i]]
          - emission[index_i(0)][observation[i]]
          + bm[j+1][i+1]);
        list_d.push_back(fm[j][i]
          + transition[index_m(j)][index_d(j+1)]
          + bd[j][i]);
        list_l.push_back(fm[j][i]
          + transition[index_m(j)][index_d(model_len)]
          + bd[model_len-1][i]);
      }
      e_trans[index_m(j)][index_i(j)] =
        log_sum_log_vec(list_i, list_i.size()) - ll;
      e_trans[index_m(j)][index_m(j+1)] =
        log_sum_log_vec(list_m, list_m.size()) - ll;
      e_trans[index_m(j)][index_d(j+1)] =
        log_sum_log_vec(list_d, list_d.size()) - ll;
      e_trans[index_m(j)][index_d(model_len)] =
        log_sum_log_vec(list_l, list_l.size()) - ll;
      list_i.clear();
      list_m.clear();
      list_d.clear();
      list_l.clear();
      // D_2~L-2: to D_j+1, M_j+1, I_j
      // since D_j are inner deletion states, fd[j][0] are 0,
      // so the loop bounds (1~seq_len-1) are same as M_j and I_j
      if (j > 1) {
        for (size_t i = 1; i < seq_len; ++i) {
          list_m.push_back(fd[j-1][i]
            + transition[index_d(j)][index_m(j+1)]
            + emission[index_m(j+1)][observation[i]]
            - emission[index_i(0)][observation[i]]
            + bm[j+1][i+1]);
          list_i.push_back(fd[j-1][i]
            + transition[index_d(j)][index_i(j)]
            + emission[index_i(j)][observation[i]]
            - emission[index_i(0)][observation[i]]
            + bi[j][i+1]);
          list_d.push_back(fd[j-1][i]
            + transition[index_d(j)][index_d(j+1)]
            + bd[j][i]);
        }
        e_trans[index_d(j)][index_m(j+1)] =
          log_sum_log_vec(list_m, list_m.size()) - ll;
        e_trans[index_d(j)][index_i(j)] = 
          log_sum_log_vec(list_i, list_i.size()) - ll;
        e_trans[index_d(j)][index_d(j+1)] =
          log_sum_log_vec(list_d, list_d.size()) - ll;
        list_m.clear();
        list_i.clear();
        list_d.clear();
      }
    }

    // Emissions
    // not sure how to train bg emission, since log-odds ratio is used
    // for all I_i:1~L-1
    // how to determine pseudo count?
    vector<vector<double> > list_emis(4, vector<double>(1, log(0.1)));
    for (size_t j = 1; j < model_len; ++j) {
      for (size_t i = 0; i < seq_len; ++i) {
        list_emis[observation[i]].push_back(fi[j][i+1] + bi[j][i+1]);
      }
      for (size_t k = 0; k < 4; ++k) {
        e_emiss[index_i(j)][k] =
          log_sum_log_vec(list_emis[k], list_emis[k].size()) - ll;
      }
    }
    list_emis.resize(4, vector<double>(1, log(0.01)));
    // for all M_i:1~L
    for (size_t j = 1; j <= model_len; ++j) {
      for (size_t i = 0; i < seq_len; ++i) {
        list_emis[observation[i]].push_back(fm[j][i+1] + bm[j][i+1]);
      }
      for (size_t k = 0; k < 4; ++k) {
        e_emiss[index_m(j)][k] =
          log_sum_log_vec(list_emis[k], list_emis[k].size()) - ll;
      }
    }

    // maximize parameters
    // Transitions
    // M_0
    double sum = log_sum_log(e_trans[index_m(0)][index_d(1)],
      e_trans[index_m(0)][index_i(0)]);
    transition[index_m(0)][index_d(1)] = e_trans[index_m(0)][index_d(1)]
      - sum;
    transition[index_m(0)][index_i(0)] = e_trans[index_m(0)][index_i(0)]
      - sum;
    // I_0
    sum = log_sum_log_list({e_trans[index_i(0)][index_d(1)],
      e_trans[index_i(0)][index_i(0)], e_trans[index_i(0)][total_size-1]});
    transition[index_i(0)][index_d(1)] = e_trans[index_i(0)][index_d(1)] - sum;
    transition[index_i(0)][index_i(0)] = e_trans[index_i(0)][index_i(0)] - sum;
    transition[index_i(0)][total_size-1] =
      e_trans[index_i(0)][total_size-1] - sum;
    // D1 to M_i
    vector<double>::const_iterator first = e_trans[index_d(1)].begin()
      + index_m(1);
    vector<double>::const_iterator last = e_trans[index_d(1)].begin()
      + index_m(model_len);
    vector<double> sum_list(first, last);
    sum = log_sum_log_vec(sum_list, sum_list.size());
    for (size_t i = 1; i <= model_len; ++i) {
      transition[index_d(1)][index_m(i)] = e_trans[index_d(1)][index_m(i)]
        - sum;
    }
    // D_L to I_0 and E
    sum = log_sum_log(e_trans[index_d(model_len)][index_i(0)],
      e_trans[index_d(model_len)][total_size-1]);
    transition[index_d(model_len)][index_i(0)] =
      e_trans[index_d(model_len)][index_i(0)] - sum;
    transition[index_d(model_len)][total_size-1] =
      e_trans[index_d(model_len)][total_size-1] - sum;
    // D_L-1 to M_L and I_L-1
    sum = log_sum_log(e_trans[index_d(model_len-1)][index_m(model_len)],
      e_trans[index_d(model_len-1)][index_i(model_len-1)]);
    transition[index_d(model_len-1)][index_m(model_len)] = 
      e_trans[index_d(model_len-1)][index_m(model_len)] - sum;
    transition[index_d(model_len-1)][index_i(model_len-1)] = 
      e_trans[index_d(model_len-1)][index_i(model_len-1)] - sum;
    // I_L-1 to I_L-1 and M_L
    sum = log_sum_log(e_trans[index_i(model_len-1)][index_i(model_len-1)],
      e_trans[index_i(model_len-1)][index_m(model_len)]);
    transition[index_i(model_len-1)][index_i(model_len-1)] = 
      e_trans[index_i(model_len-1)][index_i(model_len-1)] - sum;
    transition[index_i(model_len-1)][index_m(model_len)] = 
      e_trans[index_i(model_len-1)][index_m(model_len)] - sum;
    // M_L-1 to I_L-1, M_L and D_L
    sum = log_sum_log_list({
      e_trans[index_m(model_len-1)][index_i(model_len-1)],
      e_trans[index_m(model_len-1)][index_m(model_len)],
      e_trans[index_m(model_len-1)][index_d(model_len)]});
    transition[index_m(model_len-1)][index_i(model_len-1)] = 
      e_trans[index_m(model_len-1)][index_i(model_len-1)] - sum;
    transition[index_m(model_len-1)][index_m(model_len)] = 
      e_trans[index_m(model_len-1)][index_m(model_len)] - sum;
    transition[index_m(model_len-1)][index_d(model_len)] = 
      e_trans[index_m(model_len-1)][index_d(model_len)] - sum;
    // general
    for (size_t j = 1; j < model_len - 1; ++j) {
      // I_j: 1~L-2
      sum = log_sum_log_list({
        e_trans[index_i(j)][index_i(j)],
        e_trans[index_i(j)][index_m(j+1)],
        e_trans[index_i(j)][index_d(j+1)]});
      transition[index_i(j)][index_i(j)] = 
        e_trans[index_i(j)][index_i(j)] - sum;
      transition[index_i(j)][index_m(j+1)] = 
        e_trans[index_i(j)][index_m(j+1)] - sum;
      transition[index_i(j)][index_d(j+1)] = 
        e_trans[index_i(j)][index_d(j+1)] - sum;
      // M_j: 1~L-2
      sum = log_sum_log_list({
        e_trans[index_m(j)][index_m(j+1)],
        e_trans[index_m(j)][index_i(j)],
        e_trans[index_m(j)][index_d(j+1)],
        e_trans[index_m(j)][index_d(model_len)]});
      transition[index_m(j)][index_m(j+1)] =
        e_trans[index_m(j)][index_m(j+1)] - sum;
      transition[index_m(j)][index_i(j)] = 
        e_trans[index_m(j)][index_i(j)] - sum;
      transition[index_m(j)][index_d(j+1)] = 
        e_trans[index_m(j)][index_d(j+1)] - sum;
      transition[index_m(j)][index_d(model_len)] =
        e_trans[index_m(j)][index_d(model_len)] - sum;
      // D_j: 2~L-2
      if (j > 1) {
        sum = log_sum_log_list({
          e_trans[index_d(j)][index_d(j+1)],
          e_trans[index_d(j)][index_i(j)],
          e_trans[index_d(j)][index_m(j+1)]});
        transition[index_d(j)][index_d(j+1)] = 
          e_trans[index_d(j)][index_d(j+1)] - sum;
        transition[index_d(j)][index_i(j)] = 
          e_trans[index_d(j)][index_i(j)] - sum;
        transition[index_d(j)][index_m(j+1)] = 
          e_trans[index_d(j)][index_m(j+1)] - sum;
      }
    }

    // Emissions
    // I_i:1~L-1
    // M_i:1~L
    for (size_t j = 1; j <= model_len; ++j) {
      vector<double> list_i, list_m;
      for (size_t k =0; k < 4; ++k) {
        if (j < model_len) {
          for (size_t l = 0; l < 4; ++l) {
            list_i.push_back(emission[index_i(j)][l]);
          }
          list_i.erase(list_i.begin() + k);
          emission[index_i(j)][k] = e_emiss[index_i(j)][k]
            - log_sum_log_vec(list_i, list_i.size());
        }
        list_i.clear();
        for (size_t l = 0; l < 4; ++l) {
          list_m.push_back(emission[index_m(j)][l]);
        }
        list_m.erase(list_m.begin() + k);
        emission[index_m(j)][k] = e_emiss[index_m(j)][k]
          - log_sum_log_vec(list_m, list_m.size());
        list_m.clear();
      }
    }

    // update forward and backward prob using updated parameters
    forward_algorithm(false, transition, emission, observation);
    backward_algorithm(false, transition, emission, observation);
    if (itr > 1) ll = ll_new;
    ll_new = forward_prob('E', 1, 1);
    if (VERBOSE) {
      cout << "ITER:" << itr << "/" << max_iterations << endl;
      cout << ll << "\t" << ll_new << "\t" << ll_new - ll << endl;
    }
  }
}

size_t
random_weighted_sample(const gsl_rng* rng, const vector<double> &prob) {
  double cumulative = 0.0;
  double rand = gsl_rng_uniform(rng);
  size_t i;
  for (i = 0; i < prob.size() && cumulative < rand; ++i) {
    cumulative += exp(prob[i]);
  }
  return i-1;
}

void
ProfileHMM::sample_sequence(const bool VERBOSE,
    const gsl_rng* rng,
    const vector<vector<double> > &transition,
    const vector<vector<double> > &emission,
    vector<int> &seq,
    vector<size_t> &states) const {
  size_t idx = 0;
  states.push_back(idx);

  while (idx < total_size - 1) {
    idx = random_weighted_sample(rng, transition[idx]);
    states.push_back(idx);
    if (exp(log_sum_log_vec(emission[idx], emission[idx].size())) - 1.0
      > -tolerance) {
      seq.push_back(static_cast<int>(
            random_weighted_sample(rng, emission[idx])));
    }
  }

  if (VERBOSE) {
    for (vector<size_t>::const_iterator i = states.begin();
        i < states.end(); ++i)
      cout << *i << " ";
    cout << endl;
    for (vector<int>::const_iterator i = seq.begin();
        i < seq.end(); ++i)
      cout << *i << " ";
    cout << endl;
  }
}
