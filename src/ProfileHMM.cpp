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
#include <iomanip>

using std::vector;
using std::max;
using std::pair;
using std::string;
using std::endl;
using std::cout;
using std::cerr;
using std::map;

const double LOG_ZERO = -1e20;

size_t
baseint2stateint(const size_t &baseint,
    const bool marked) {
  if (marked && baseint < 4)
    return 0;
  else if (marked && baseint == 4)
    return 2;
  else if (marked && baseint == 5)
    return 3; // mark at outside of alignment, truncation
  else if (!marked && baseint < 4)
    return 1;
  else
    return 4; // gap in not marked column; this is not supposed to happen
}

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

state::state() {};

state::state(const char b, const size_t i, const bool m) {
  baseint = base2int(b);
  idx = i;
  stateint = baseint2stateint(base2int(b), m);
}

state::state(const state &obj) {
  baseint = obj.baseint;
  idx = obj.idx;
  stateint = obj.stateint;
}

state::state(const size_t model_len, const size_t matrix_idx) {
  baseint = 0;

  if (matrix_idx <= model_len) {
    stateint = 0;
    idx = matrix_idx;
  }
  else if (matrix_idx <= model_len * 2) {
    stateint = 1;
    idx = matrix_idx - model_len - 1;
  }
  else {
    stateint = 2;
    idx = matrix_idx - model_len * 2;
  }
}

bool
state::isvalid(void) const {
  return (stateint < 2 && baseint < 4) 
    || (stateint == 2 && baseint == 4);
}

size_t
state::index(const size_t model_len) const {
  if (stateint == 0)
    return idx;
  else if (stateint == 1)
    return idx + model_len + 1;
  else if (stateint == 2)
    return idx + model_len * 2;
  else
    return 0;
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

ProfileHMM::ProfileHMM() {};

ProfileHMM::ProfileHMM(const size_t ml) {
  model_len = ml;
  total_size = model_len * 3 + 2;
  transitions_to = get_viable_transitions_to();
  transitions_from = get_viable_transitions_from();
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

string
ProfileHMM::state_idx_to_str(const size_t idx) const {
  if (idx <= model_len)
    return "M" + std::to_string(idx);
  else if (idx <= model_len * 2)
    return "I" + std::to_string(idx - model_len - 1);
  else if (idx <= model_len * 3)
    return "D" + std::to_string(idx - model_len * 2);
  else
    return "E";
}

double
ProfileHMM::ViterbiDecoding(const bool VERBOSE,
    const matrix &transition,
    const matrix &emission,
    const vector<int> &observation,
    vector<pair<char, size_t> > &trace) const {
  const size_t seq_len = observation.size();
  // Note here: vm[i][j] corresponds to M_i on o_j, where i=0~L
  // vd[i][j]: i = 0~L-1 corresponds to D_1 ~ D_L
  // vi[i][j]: i = 0~L-1 corresponds to I_0 ~ I_L-1
  // crazy subscriptions!
  matrix vm(model_len+1, vector<double>(seq_len+1, LOG_ZERO));
  matrix vi(model_len, vector<double>(seq_len+1, LOG_ZERO));
  matrix vd(model_len, vector<double>(seq_len+1, LOG_ZERO));

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
    for (matrix::const_iterator i = vm.begin(); i < vm.end(); ++i) {
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
    for (matrix::const_iterator i = vi.begin(); i < vi.end(); ++i) {
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
    for (matrix::const_iterator i = vd.begin(); i < vd.end(); ++i) {
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
        size_t max_idx = argmax_list({
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
        size_t max_idx = argmax_list({
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
        size_t max_idx = argmax_list({
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
        size_t max_idx = argmax_list({
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
        size_t max_idx = argmax_list({
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
    const matrix &transition,
    const matrix &emission,
    const vector<int> &observation,
    matrix &forward) const {
  const size_t seq_len = observation.size();
  forward.resize(total_size, vector<double>(seq_len+1, LOG_ZERO));

  forward[0][0] = 0;
  forward[index_d(1)][0] = transition[index_m(0)][index_d(1)];

  for (size_t pos = 1; pos < seq_len + 1; ++pos) {
    for (size_t state_idx = 0; state_idx < index_d(1); ++state_idx) {
      // states with emission: M/I
      const map<size_t, vector<size_t> >::const_iterator i = 
        transitions_from.find(state_idx);
      if (i != transitions_from.end()) {
        vector<double> list;
        for (vector<size_t>::const_iterator j = i->second.begin();
            j < i->second.end(); ++j) {
          list.push_back(forward[*j][pos-1] + transition[*j][i->first]);
        }
        forward[i->first][pos] =
          smithlab::log_sum_log_vec(list, list.size())
          + emission[i->first][observation[pos-1]];
      }
    }
    for (size_t state_idx = index_d(1); state_idx < total_size; ++state_idx) {
      // states without emission: D
      const map<size_t, vector<size_t> >::const_iterator i = 
        transitions_from.find(state_idx);
      if (i != transitions_from.end()) {
        vector<double> list;
        for (vector<size_t>::const_iterator j = i->second.begin();
            j < i->second.end(); ++j) {
          list.push_back(forward[*j][pos] + transition[*j][i->first]);
        }
        forward[i->first][pos] =
          smithlab::log_sum_log_vec(list, list.size());
      }
    }
  }

  if (VERBOSE) {
    cout << endl << "#Forward" << endl;
    cout << "#M";
    for (size_t i = 1; i <= model_len; ++i) {
      cout << "\tM" << i;
    }
    cout << endl;
    for (size_t pos = 1; pos < seq_len +1; ++pos) {
      for (size_t i = 1; i <= model_len; ++i) {
        cout << pos << "\t" << forward[i][pos];
      }
      cout << endl;
    }

    cout << endl << "#I";
    for (size_t i = 0; i < model_len; ++i) {
      cout << "\tI" << i;
    }
    cout << endl;
    for (size_t pos = 1; pos < seq_len +1; ++pos) {
      for (size_t i = 0; i < model_len; ++i) {
        cout << pos << "\t" << forward[index_i(i)][pos];
      }
      cout << endl;
    }

    cout << endl << "#D";
    for (size_t i = 1; i <= model_len; ++i) {
      cout << "\tD" << i;
    }
    cout << endl;
    for (size_t pos = 1; pos < seq_len +1; ++pos) {
      for (size_t i = 1; i <= model_len; ++i) {
        cout << pos << "\t" << forward[index_d(i)][pos];
      }
      cout << endl;
    }
  }
}

void
ProfileHMM::backward_algorithm(const bool VERBOSE,
    const matrix &transition,
    const matrix &emission,
    const vector<int> &observation,
    matrix &backward) const {
  const size_t seq_len = observation.size();
  backward.resize(total_size, vector<double>(seq_len+1, LOG_ZERO));

  // Initialization
  for (size_t i = 0; i < total_size; ++i)
    backward[i].back() = 0.0;

  // loop is going backwards
  for (signed long pos = seq_len - 1; pos >= 0; --pos) {
    for (signed long state_idx = total_size - 1;
        state_idx >= 0; --state_idx) {
      const map<size_t, vector<size_t> >::const_iterator i = 
        transitions_to.find(state_idx);
      if (i != transitions_to.end()) {
        vector<double> list;
        for (vector<size_t>::const_reverse_iterator j = i->second.rbegin();
            j < i->second.rend(); ++j) {
          if (*j > 0 && *j < index_d(1)) {
            // *j is an M/I state with viable emission
            list.push_back(backward[*j][pos+1]
                + transition[i->first][*j]
                + emission[*j][observation[pos]]);
          }
          else {
            // *j is a D state without viable emission
            // this assert makes sure that
            // backward[*j][pos] has been calculated
            assert(*j > i->first);
            list.push_back(backward[*j][pos]
                + transition[i->first][*j]);
          }
        }
        backward[i->first][pos] =
          smithlab::log_sum_log_vec(list, list.size());
      }
    }
  }

  if (VERBOSE) {
    cout << endl << "#Backward" << endl;
    cout << "#M";
    for (size_t i = 1; i <= model_len; ++i) {
      cout << "\tM" << i;
    }
    cout << endl;
    for (size_t pos = 1; pos < seq_len +1; ++pos) {
      for (size_t i = 1; i <= model_len; ++i) {
        cout << pos << "\t" << backward[i][pos];
      }
      cout << endl;
    }

    cout << endl << "#I";
    for (size_t i = 0; i < model_len; ++i) {
      cout << "\tI" << i;
    }
    cout << endl;
    for (size_t pos = 1; pos < seq_len +1; ++pos) {
      for (size_t i = 0; i < model_len; ++i) {
        cout << pos << "\t" << backward[index_i(i)][pos];
      }
      cout << endl;
    }

    cout << endl << "#D";
    for (size_t i = 1; i <= model_len; ++i) {
      cout << "\tD" << i;
    }
    cout << endl;
    for (size_t pos = 1; pos < seq_len +1; ++pos) {
      for (size_t i = 1; i <= model_len; ++i) {
        cout << pos << "\t" << backward[index_d(i)][pos];
      }
      cout << endl;
    }
  }
}

void
ProfileHMM::Train(const bool VERBOSE,
    matrix &transition,
    matrix &emission,
    const vector<int> &observation) {
  const size_t seq_len = observation.size();
  matrix forward, backward;
  forward_algorithm(false, transition, emission, observation, forward);
  backward_algorithm(false, transition, emission, observation, backward);
  double ll = -1e10;
  double ll_new = posterior_prob(forward);
  if (VERBOSE)
    cout << "ITR\tDELTA" << endl;
  for (size_t itr = 0;
      //itr < max_iterations;
      std::abs((ll_new - ll)/ll_new) > tolerance && itr < max_iterations;
      ++itr) {
    std::swap(ll, ll_new);
    matrix e_trans(total_size,
        vector<double>(total_size, LOG_ZERO));
    matrix e_emiss(total_size,
        vector<double>(emission.front().size(), LOG_ZERO));
    // Expectation
    // transition
    for (map<size_t, vector<size_t> >::const_iterator trans =
        transitions_to.begin(); trans != transitions_to.end(); ++trans) {
      const size_t i = trans->first;
      for (vector<size_t>::const_iterator j = trans->second.begin();
          j < trans->second.end(); ++j) {
        vector<double> list;
        for (size_t pos = 0; pos <= seq_len; ++pos) {
          if (*j > index_m(0) && *j < index_d(1)) {
            if (pos < seq_len)
              list.push_back(forward[i][pos] + transition[i][*j]
                  + emission[*j][observation[pos]]
                  + backward[*j][pos + 1]);
          }
          else {
            list.push_back(forward[i][pos] + transition[i][*j]
                + backward[*j][pos]);
          }
        }
        e_trans[i][*j] = smithlab::log_sum_log_vec(list, list.size()) - ll;
      }
    }
    // emission
    for (size_t i = index_m(1); i < index_d(1); ++i) {
      matrix list(4, vector<double>());
      for (size_t pos = 0; pos < seq_len; ++pos) {
        list[observation[pos]].push_back(forward[i][pos+1] + backward[i][pos+1]);
      }
      for (size_t k = 0; k < 4; ++k) {
        e_emiss[i][k] =
          smithlab::log_sum_log_vec(list[k], list[k].size()) - ll;
      }
    }

    // Maximization
    // transition
    for (map<size_t, vector<size_t> >::const_iterator from
        = transitions_to.begin();
        from != transitions_to.end(); ++from) {
      vector<double> list;
      for (vector<size_t>::const_iterator to = from->second.begin();
          to < from->second.end(); ++to) {
        list.push_back(e_trans[from->first][*to]);
      }
      const double sum = smithlab::log_sum_log_vec(list, list.size());
      for (vector<size_t>::const_iterator to = from->second.begin();
          to < from->second.end(); ++to) {
        transition[from->first][*to] = e_trans[from->first][*to] - sum;
      }
    }
    // emission
    for (size_t i = index_m(1); i < index_d(1); ++i) {
      const double sum =
        smithlab::log_sum_log_vec(e_emiss[i], e_emiss[i].size());
      for (size_t j = 0; j < 4; ++j) {
        emission[i][j] = e_emiss[i][j] - sum;
      }
    }
    // add pseudocount and normalize
    add_pseudocount_uniform(transition, emission);
    for (matrix::iterator row = transition.begin();
        row < transition.end() - 1; ++row)
      normalize_vec_inplace(*row);
    for (size_t i = index_m(1); i < index_d(1); ++i)
      normalize_vec_inplace(emission[i]);

    // update forward and backward prob using updated parameters
    forward_algorithm(false, transition, emission, observation, forward);
    backward_algorithm(false, transition, emission, observation, backward);
    // Need some function for posterior probability here!
    ll_new = posterior_prob(forward);
    if (VERBOSE) {
      cout << itr + 1 << "/" << max_iterations << "\t"
        << ll << "\t" << ll_new << "\t"
        << std::abs((ll_new - ll)/ll_new) << endl;
    }
  }
}

void
ProfileHMM::SampleSequence(const bool VERBOSE,
    const gsl_rng* rng,
    const matrix &transition,
    const matrix &emission,
    vector<int> &seq,
    vector<size_t> &states) const {
  size_t idx = 0;
  states.push_back(idx);

  while (idx < total_size - 1) {
    idx = random_weighted_sample(rng, transition[idx]);
    states.push_back(idx);
    if (exp(smithlab::log_sum_log_vec(emission[idx], emission[idx].size())) - 1.0
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

void
ProfileHMM::PosteriorDecoding(const bool DEBUG,
    const matrix &transition,
    const matrix &emission,
    const vector<int> &observation,
    vector<size_t> &states) const {
  const size_t seq_len = observation.size();
  matrix forward, backward;

  forward_algorithm(false, transition, emission, observation, forward);
  backward_algorithm(false, transition, emission, observation, backward);

  vector<double> xi(index_d(1), LOG_ZERO);
  for (size_t i = 1; i <= seq_len; ++i) {
    for (size_t j = index_m(1); j < index_d(1); ++j) {
      xi[j] = forward[j][i] + backward[j][i];
    }
    size_t idx = argmax_vec(xi);
    states.push_back(idx);
  }
  if (DEBUG) {
    cerr << "State sequence:" << endl;
    for (vector<size_t>::const_iterator i = states.begin();
        i < states.end(); ++i) {
      cerr << i - states.begin() + 1 << "\t" << state_idx_to_str(*i) << endl;
    }
  }
}

double
ProfileHMM::PosteriorProb(const matrix &transition,
    const matrix &emission,
    const vector<int> &observation) const {
  matrix forward;

  forward_algorithm(false, transition, emission, observation, forward);
  return posterior_prob(forward);
}

double
ProfileHMM::posterior_prob(const matrix &forward) const {
  vector<double> list;
  for (matrix::const_iterator i = forward.begin(); i < forward.end(); ++i)
    list.push_back((*i).back());
  return smithlab::log_sum_log_vec(list, list.size());
}

void
ProfileHMM::add_pseudocount_uniform(matrix &transition,
    matrix &emission) const {
  assert(transition.size() == total_size);
  assert(emission.size() == total_size);
  const double PSEUDOCOUNT = -5;
  
  // transition
  // M_0 to I_0 and D_1
  transition[0][index_d(1)] =
    log_sum_log(transition[0][index_d(1)], PSEUDOCOUNT);
  transition[0][index_i(0)] =
    log_sum_log(transition[0][index_i(0)], PSEUDOCOUNT);
  // I_0 to I_0, D_1 and end
  transition[index_i(0)][index_i(0)] =
    log_sum_log(transition[index_i(0)][index_i(0)], PSEUDOCOUNT);
  transition[index_i(0)][index_d(1)] =
    log_sum_log(transition[index_i(0)][index_d(1)], PSEUDOCOUNT);
  transition[index_i(0)].back() =
    log_sum_log(transition[index_i(0)].back(), PSEUDOCOUNT);
  // D_1 to M_i
  for (size_t i = 1; i <= model_len; ++i)
    transition[index_d(1)][index_m(i)] =
      log_sum_log(transition[index_d(1)][index_m(i)], PSEUDOCOUNT);
  // general
  for (size_t i = 1; i < model_len; ++i) {
    // M_i
    transition[index_m(i)][index_m(i+1)] =
      log_sum_log(transition[index_m(i)][index_m(i+1)], PSEUDOCOUNT);
    transition[index_m(i)][index_i(i)] =
      log_sum_log(transition[index_m(i)][index_i(i)], PSEUDOCOUNT);
    transition[index_m(i)][index_d(i+1)] =
      log_sum_log(transition[index_m(i)][index_d(i+1)], PSEUDOCOUNT);
    // I_i
    transition[index_i(i)][index_m(i+1)] =
      log_sum_log(transition[index_i(i)][index_m(i+1)], PSEUDOCOUNT);
    transition[index_i(i)][index_i(i)] =
      log_sum_log(transition[index_i(i)][index_i(i)], PSEUDOCOUNT);
    if (i < model_len - 1)
      transition[index_i(i)][index_d(i+1)] =
        log_sum_log(transition[index_i(i)][index_d(i+1)], PSEUDOCOUNT);
    // D_i
    if (i > 1) {
      transition[index_d(i)][index_m(i+1)] =
        log_sum_log(transition[index_d(i)][index_m(i+1)], PSEUDOCOUNT);
      transition[index_d(i)][index_i(i)] =
        log_sum_log(transition[index_d(i)][index_i(i)], PSEUDOCOUNT);
      if (i < model_len - 1)
        transition[index_d(i)][index_d(i+1)] =
          log_sum_log(transition[index_d(i)][index_d(i+1)], PSEUDOCOUNT);
    }
  }
  // D_L to I_0 and end
  transition[index_d(model_len)][index_i(0)] =
    log_sum_log(transition[index_d(model_len)][index_i(0)], PSEUDOCOUNT);
  transition[index_d(model_len)].back() =
    log_sum_log(transition[index_d(model_len)].back(), PSEUDOCOUNT);

  // emission
  for (size_t i = index_m(1); i <= index_i(model_len - 1); ++i) {
    for (size_t j = 0; j < 4; ++j)
      emission[i][j] = log_sum_log(emission[i][j], PSEUDOCOUNT);
  }
}

map<size_t, vector<size_t> >
ProfileHMM::get_viable_transitions_to(void) const {
  map<size_t, vector<size_t> > t;
  
  // M_0 to I_0, D_1
  t[0] = vector<size_t>({index_i(0), index_d(1)});
  // I_0 to I_0, D_1, E
  t[index_i(0)] = vector<size_t>({index_i(0), index_d(1), total_size - 1});
  // D_1 to M_1...M_L
  t[index_d(1)] = vector<size_t>();
  for (size_t i = 1; i <= model_len; ++i)
    t[index_d(1)].push_back(index_m(i));
  // general: i~L-1
  for (size_t i = 1; i < model_len; ++i) {
    t[index_m(i)] = vector<size_t>({index_m(i+1), index_i(i),
        index_d(i+1), index_d(model_len)});
    if (i < model_len - 1) {
      t[index_i(i)] = vector<size_t>({index_m(i+1), index_i(i), index_d(i+1)});
      if (i > 1)
        t[index_d(i)] =
          vector<size_t>({index_m(i+1), index_i(i), index_d(i+1)});
    }
  }
  // I_L-1
  t[index_i(model_len-1)] =
    vector<size_t>({index_i(model_len-1), index_m(model_len)});
  // D_L-1
  t[index_d(model_len-1)] =
    vector<size_t>({index_i(model_len-1), index_m(model_len)});
  // M_L
  t[index_m(model_len)] = vector<size_t>(1, index_d(model_len));
  // D_L
  t[index_d(model_len)] = vector<size_t>({index_i(0), total_size - 1});

  return t;
}

map<size_t, vector<size_t> >
ProfileHMM::get_viable_transitions_from(void) const {
  map<size_t, vector<size_t> > t;

  // to I_0: I_0, M_0, D_L
  t[index_i(0)] = vector<size_t>({index_i(0), index_m(0), index_d(model_len)});
  // to D_1: M_0, I_0
  t[index_d(1)] = vector<size_t>({index_m(0), index_i(0)});
  // to M_1: D_1
  t[index_m(1)] = vector<size_t>(1, index_d(1));
  // to I_1: M_1, I_1
  t[index_i(1)] = vector<size_t>({index_m(1), index_i(1)});
  // to M_2: M_1, I_1, D_1
  t[index_m(2)] = vector<size_t>({index_m(1), index_i(1), index_d(1)});
  // to I_2: M_1, I_2, D_2
  t[index_i(2)] = vector<size_t>({index_m(2), index_i(2), index_d(2)});
  // to D_2: M_1, I_1
  t[index_d(2)] = vector<size_t>({index_m(1), index_i(1)});
  // general: i=3~L-1
  for (size_t i = 3; i < model_len; ++i) {
    t[index_m(i)] = vector<size_t>({index_m(i-1), index_i(i-1),
        index_d(i-1), index_d(1)});
    t[index_i(i)] = vector<size_t>({index_m(i), index_i(i), index_d(i)});
    t[index_d(i)] = vector<size_t>({index_m(i-1), index_i(i-1), index_d(i-1)});
  }
  // to M_L: M, I, D_L-1, D_1
  t[index_m(model_len)] = vector<size_t>({index_m(model_len-1),
      index_i(model_len-1), index_d(model_len-1), index_d(1)});
  // to D_L: M_1 to M_L
  t[index_d(model_len)] = vector<size_t>();
  for (size_t i = 1; i <= model_len; ++i)
    t[index_d(model_len)].push_back(index_m(i));
  // to E: I_0, D_L
  t[total_size - 1] = vector<size_t>({index_i(0), index_d(model_len)});

  return t;
}

void
print_transition(const matrix &transition) {
  const size_t model_len = (transition.size() - 2) / 3;

  cout << "#M_0->I_0, M_0->D_1 = "
    << exp(transition[0][state(0ul, 0, 1).index(model_len)]) << ", "
    << exp(transition[0][state(0ul, 1, 2).index(model_len)]) << endl;
  cout << "#M_i:" 
    << "\t" << "M->M" << "\t" << "M->I" << "\t" << "M->D"
    << "\t" << "I->M" << "\t" << "I->I" << "\t" << "I->D"
    << "\t" << "D->M" << "\t" << "D->I" << "\t" << "D->D"
    << endl;
  for (size_t i = 1; i < model_len; ++i) {
    cout << i;
    // M
    cout << std::setprecision(4) << std::fixed
      << "\t" << exp(transition[state(0ul,i,0).index(model_len)]
        [state(0ul,i+1,0).index(model_len)])
      << "\t" << exp(transition[state(0ul,i,0).index(model_len)]
        [state(0ul,i,1).index(model_len)])
      << "\t" << exp(transition[state(0ul,i,0).index(model_len)]
        [state(0ul,i+1,2).index(model_len)]);
    // I
    cout
      << "\t" << exp(transition[state(0ul,i,1).index(model_len)]
        [state(0ul,i+1,0).index(model_len)])
      << "\t" << exp(transition[state(0ul,i,1).index(model_len)]
        [state(0ul,i,1).index(model_len)])
      << "\t" << exp(transition[state(0ul,i,1).index(model_len)]
        [state(0ul,i+1,2).index(model_len)]);
    // D
    if (i > 1) {
      cout
        << "\t" << exp(transition[state(0ul,i,2).index(model_len)]
          [state(0ul,i+1,0).index(model_len)])
        << "\t" << exp(transition[state(0ul,i,2).index(model_len)]
          [state(0ul,i,1).index(model_len)]);
      if (i < model_len - 1)
        cout << "\t" << exp(transition[state(0ul,i,2).index(model_len)]
            [state(0ul,i+1,2).index(model_len)]);
      else
        cout << "\tnan";
    }
    else
      cout << "\tnan\tnan\tnan";
    cout << endl;
  }
  cout << endl << "#LOG\tT->M\tM->T" << endl;
  for (size_t i = 1; i <= model_len; ++i) {
    cout << i << std::setprecision(2)
      << "\t" << transition[state(0ul,1,2).index(model_len)]
        [state(0ul, i, 0).index(model_len)]
      << "\t" << transition[state(0ul,i,0).index(model_len)]
        [state(0ul,model_len,2).index(model_len)]
      << endl;
  }
  cout << std::setprecision(4);
  cout << endl << "#I_0\tI_0\tD_1\tE" << endl
    << "\t" << exp(transition[state(0ul,0,1).index(model_len)]
      [state(0ul,0,1).index(model_len)])
    << "\t" << exp(transition[state(0ul,0,1).index(model_len)]
      [state(0ul,1,2).index(model_len)])
    << "\t" << exp(transition[state(0ul,0,1).index(model_len)].back())
    << endl;
  cout << endl << "#D_L\tI_0\tE" << endl
    << "\t" << exp(transition[state(0ul,model_len,2).index(model_len)]
      [state(0ul,0,1).index(model_len)])
    << "\t" << exp(transition[state(0ul,model_len,2).index(model_len)].back())
    << endl;
}

void
print_emission(const matrix &emission) {
  const size_t model_len = (emission.size() - 2) / 3;
  cout << endl << "#M_i:" 
    << "\t" << "A" << "\t" << "C" << "\t" << "G" << "\t" << "T" << endl;
  for (size_t i = 1; i <= model_len; ++i) {
    cout << i;
    cout << std::setprecision(4) << std::fixed
      << "\t" << exp(emission[state(0ul,i,0).index(model_len)][0])
      << "\t" << exp(emission[state(0ul,i,0).index(model_len)][1])
      << "\t" << exp(emission[state(0ul,i,0).index(model_len)][2])
      << "\t" << exp(emission[state(0ul,i,0).index(model_len)][3])
      << endl;
  }

  cout << endl << "#I_i:" 
    << "\t" << "A" << "\t" << "C" << "\t" << "G" << "\t" << "T" << endl;
  for (size_t i = 0; i < model_len; ++i) {
    cout << i;
    cout << std::setprecision(4)
      << "\t" << exp(emission[state(0ul,i,1).index(model_len)][0])
      << "\t" << exp(emission[state(0ul,i,1).index(model_len)][1])
      << "\t" << exp(emission[state(0ul,i,1).index(model_len)][2])
      << "\t" << exp(emission[state(0ul,i,1).index(model_len)][3])
      << endl;
  }
}

void
log_odds_transform(matrix &emission) {
  const size_t model_len = (emission.size() - 2) / 3;
  vector<double> bg = emission[state(0ul, 0, 1).index(model_len)];
  for (size_t i = 0; i < model_len; ++i) {
    for (size_t j = 0; j < emission.front().size(); ++j) {
      emission[state(0ul, i+1, 0).index(model_len)][j] =
        emission[state(0ul, i+1, 0).index(model_len)][j] - bg[j];
      emission[state(0ul, i, 1).index(model_len)][j] =
        emission[state(0ul, i, 1).index(model_len)][j] - bg[j];
    }
  }
}
