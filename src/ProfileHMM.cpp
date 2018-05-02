/*    ProfileHMM: implemenation of profile-HMM with internal loop
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

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

using std::vector;
using std::pair;
using std::max;
using std::string;
using std::to_string;
using std::endl;
using std::cout;
using std::cerr;
using std::ostream;
using std::cerr;
using std::unordered_map;

size_t
baseint2stateint(const size_t &baseint,
    const bool marked) {
  if (marked && baseint < 4)
    return 0; // match
  else if (marked && baseint == 4)
    return 2; // deletion
  else if (marked && baseint == 5)
    return 3; // mark at outside of alignment, truncation
  else if (!marked && baseint < 4)
    return 1; // insertion
  else
    return 4; // gap in not marked column; this is not supposed to happen
}

string
state_type_to_str(const size_t model_len, const size_t idx) {
  // M_0~M_L; I_0~I_L-1; D_1~D_L
  if (idx <= model_len)
    return "M";
  else if (idx == model_len+1)
    return "N";
  else if (idx <= model_len*2)
    return "I";
  else if (idx <= model_len*3)
    return "D";
  else
    return "*";
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
  else if (matrix_idx > model_len * 2 && matrix_idx <= model_len*3) {
    stateint = 2;
    idx = matrix_idx - model_len * 2;
  }
  else {
    // the E state
    stateint = 3;
    idx = 0;
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

ProfileHMM::ProfileHMM() {
  model_len = 0;
  total_size = 0;
}

ProfileHMM::ProfileHMM(const matrix &t,
    const matrix &e) {
  transition = t;
  emission = e;
  total_size = t.size();
  model_len = (total_size - 2) / 3;
  get_viable_transitions_to();
  get_viable_transitions_from();
  name = "NULL";
  //collapse_states();
}

ProfileHMM::ProfileHMM(const string inf) {
  load_from_file(inf);
  get_viable_transitions_to();
  get_viable_transitions_from();
  name = inf.substr(0, inf.find_last_of('.'));
  //collapse_states();
}

ProfileHMM::ProfileHMM(const ProfileHMM &other) {
  model_len = other.model_len;
  total_size = other.total_size;
  transition = other.transition;
  emission = other.emission;
  transition_c = other.transition_c;
  emission_c = other.emission_c;
  get_viable_transitions_to();
  get_viable_transitions_from();
  name = other.name;
}

ProfileHMM&
ProfileHMM::operator=(const ProfileHMM &rhs) {
  ProfileHMM tmp(rhs);
  std::swap(model_len, tmp.model_len);
  std::swap(total_size, tmp.total_size);
  std::swap(name, tmp.name);
  std::swap(transition, tmp.transition);
  std::swap(emission, tmp.emission);
  std::swap(transition_c, tmp.transition_c);
  std::swap(emission_c, tmp.emission_c);
  get_viable_transitions_to();
  get_viable_transitions_from();
  return *this;
}

void
ProfileHMM::ComplementBackground(void) {
  // When doing decoding for the reverse-complement strand,
  // the background emission also needs to be complemented.
  std::swap(emission[index_i(0)][0], emission[index_i(0)][3]);
  std::swap(emission[index_i(0)][1], emission[index_i(0)][2]);
}

size_t
ProfileHMM::index_m(const size_t idx) const {
  // M_0 ~ M_L
  assert(idx <= model_len);
  return idx;
}

size_t
ProfileHMM::index_i(const size_t idx) const {
  // I_0 ~ I_L-1 + (L+1)*M
  assert(idx <= model_len - 1);
  return idx + model_len + 1;
}

size_t
ProfileHMM::index_d(const size_t idx) const {
  // D_1 ~ D_L + (L+1)*M + L*I
  assert(idx >= 1 && idx <= model_len);
  return idx + model_len * 2;
}

string
ProfileHMM::state_idx_to_str(const size_t idx) const {
  if (idx == 0)
    return "B";
  else if (idx <= model_len)
    return "M" + to_string(idx);
  else if (idx <= model_len * 2)
    return "I" + to_string(idx - model_len - 1);
  else if (idx <= model_len * 3)
    return "D" + to_string(idx - model_len * 2);
  else
    return "E";
}

string
ProfileHMM::state_idx_to_str_c(const size_t idx) const {
  if (idx < model_len)
    return "M" + to_string(idx+1);
  else if (idx < model_len * 2)
    return "I" + to_string(idx - model_len);
  else
    return "N";
}

double
ProfileHMM::ViterbiDecoding(const bool VERBOSE,
    const string &observation,
    vector<pair<char, size_t> > &trace) const {
  const size_t seq_len = observation.length();
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
      + emission[index_m(1)][base2int(observation[i-1])]
      - emission[index_i(0)][base2int(observation[i-1])];
    // I_1 only has two incoming transition including one from M_1
    vi[1][i] = max({
      vm[1][i-1] + transition[index_m(1)][index_i(1)],
      vi[1][i-1] + transition[index_i(1)][index_i(1)]})
      + emission[index_i(1)][base2int(observation[i-1])]
      - emission[index_i(0)][base2int(observation[i-1])];
    for (size_t j = 2; j < model_len; ++j) {
      // M_i
      vm[j][i] = max({
        vm[j-1][i-1] + transition[index_m(j-1)][index_m(j)],
        vi[j-1][i-1] + transition[index_i(j-1)][index_m(j)],
        vd[j-2][i-1] + transition[index_d(j-1)][index_m(j)],
        vd[0][i-1] + transition[index_d(1)][index_m(j)]})
        + emission[index_m(j)][base2int(observation[i-1])]
        - emission[index_i(0)][base2int(observation[i-1])];
      // I_i
      vi[j][i] = max({
        vm[j][i-1] + transition[index_m(j)][index_i(j)],
        vi[j][i-1] + transition[index_i(j)][index_i(j)],
        vd[j-1][i-1] + transition[index_d(j)][index_i(j)]})
        + emission[index_i(j)][base2int(observation[i-1])]
        - emission[index_i(0)][base2int(observation[i-1])];
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
      + emission[index_m(model_len)][base2int(observation[i-1])]
      - emission[index_i(0)][base2int(observation[i-1])];
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
    const bool USE_LOG_ODDS,
    const string &observation,
    matrix &forward) const {
  if (VERBOSE)
    cerr << "FORWARD ALGORITHM" << endl;
  const size_t seq_len = observation.length();
  forward.resize(seq_len+1, vector<double>(total_size, LOG_ZERO));

  forward[0][0] = 0;
  forward[0][index_d(1)] = transition[0][index_d(1)];

  for (size_t pos = 1; pos < seq_len + 1; ++pos) {
    if (VERBOSE && pos % 10000 == 0)
      cerr << "\tPROCESSED " << 100 * pos / seq_len << "%" << endl;
    // iterations must follow a specific order, such that internal loops
    // can be processed properly.
    vector<size_t> itr_order;
    for (size_t idx = 1; idx < total_size - 1; ++idx)
      if (idx != index_d(1) && idx != index_d(model_len))
        itr_order.push_back(idx);
    itr_order.push_back(index_d(model_len));
    itr_order.push_back(index_d(1));
    for (vector<size_t>::const_iterator idx = itr_order.begin();
      idx < itr_order.end(); ++idx) {
      const size_t state_idx = *idx;
      const unordered_map<size_t, vector<size_t> >::const_iterator i = 
        transitions_from.find(state_idx);
      vector<double> list;
      // use observation[pos-1] because in this vector the position
      // is 0-based
      const int curr_baseint = base2int(observation[pos-1]);
      if (state_can_emit(state_idx)) {
        for (vector<size_t>::const_iterator j = i->second.begin();
            j < i->second.end(); ++j) {
          list.push_back(forward[pos-1][*j] + transition[*j][state_idx]);
        }
        double emission_score = LOG_ZERO;
        if (curr_baseint < 4) {
          // ACGT
          emission_score = emission[state_idx][curr_baseint];
          if (USE_LOG_ODDS)
            emission_score = emission_score - emission[index_i(0)][curr_baseint];
        }
        else {
          // an N is seen
          if (state_can_emit_n(state_idx))
            emission_score = 0.0;
          else
            emission_score = LOG_ZERO;
        }
        forward[pos][state_idx] = smithlab::log_sum_log_vec(list, list.size())
          + emission_score;
      }
      else {
        for (vector<size_t>::const_iterator j = i->second.begin();
            j < i->second.end(); ++j) {
          list.push_back(forward[pos][*j] + transition[*j][state_idx]);
        }
        forward[pos][state_idx] =
          smithlab::log_sum_log_vec(list, list.size());
      }
    }
  }
}

void
ProfileHMM::backward_algorithm(const bool VERBOSE,
    const bool USE_LOG_ODDS,
    const string &observation,
    matrix &backward) const {
  if (VERBOSE)
    cerr << "BACKWARD ALGORITHM" << endl;
  const size_t seq_len = observation.length();
  backward.resize(seq_len+1, vector<double>(total_size, LOG_ZERO));
  // Initialization conditions
  for (size_t state_idx = 0; state_idx < total_size; ++state_idx)
    backward.back()[state_idx] = 0.0;

  // loop is going backwards
  for (signed long pos = seq_len - 1; pos >= 0; --pos) {
    if (VERBOSE && (seq_len - pos) % 10000 == 0)
      cerr << "\tPROCESSED " << 100 - 100 * pos / seq_len << "%" << endl;
    vector<signed long> itr_order;
    itr_order.push_back(index_d(1));
    itr_order.push_back(index_d(model_len));
    for (signed long state_idx = total_size - 1; state_idx >= 0; --state_idx)
      if (state_idx != static_cast<signed long>(index_d(1))
        && state_idx != static_cast<signed long>(index_d(model_len)))
        itr_order.push_back(state_idx);
    for (vector<signed long>::const_iterator idx = itr_order.begin();
        idx < itr_order.end(); ++idx) {
      const signed long state_idx = *idx;
      const unordered_map<size_t, vector<size_t> >::const_iterator i = 
        transitions_to.find(static_cast<size_t>(state_idx));
      if (i != transitions_to.end()) {
        vector<double> list;
        const int next_baseint = base2int(observation[pos]);
        for (vector<size_t>::const_reverse_iterator j = i->second.rbegin();
            j < i->second.rend(); ++j) {
          if (state_can_emit(*j)) {
            // *j is an M/I state with viable emission
            const double transition_score = backward[pos+1][*j] 
              + transition[state_idx][*j];
            double emission_score = LOG_ZERO;
            if (next_baseint < 4) {
              // ACGT
              emission_score = emission[*j][next_baseint];
              if (USE_LOG_ODDS)
                emission_score = emission_score 
                  - emission[index_i(0)][next_baseint];
            }
            else {
              // an N is seen
              if (state_can_emit_n(*j))
                emission_score = 0.0;
              else
                emission_score = LOG_ZERO;
            }
            list.push_back(transition_score + emission_score);
          }
          else {
            // *j is a D state without viable emission
            list.push_back(backward[pos][*j]
                + transition[state_idx][*j]);
          }
        }
        backward[pos][state_idx] =
          smithlab::log_sum_log_vec(list, list.size());
      }
    }
  }
}

void
ProfileHMM::forward_algorithm_c(const bool VERBOSE,
    const bool USE_LOG_ODDS,
    const string &observation,
    matrix &forward) const {
  const size_t seq_len = observation.length();
  const size_t state_num = model_len*2;
  if (VERBOSE)
    cerr << "FORWARD ALGORITHM " << seq_len << "x" << state_num << endl;
  forward.resize(seq_len, vector<double>(state_num, LOG_ZERO));
  for (size_t state = 0; state < state_num; ++state) {
    forward.front()[state] = initial_prob[state]
      + emission_c[state][base2int(observation.front())];
  }

  for (size_t pos = 1; pos < seq_len; ++pos) {
    if (VERBOSE && pos % 10000 == 0)
      cerr << "\tPROCESSED " << 100 * pos / seq_len << "%" << endl;
    const int curr_baseint = base2int(observation[pos]);
    for (size_t state = 0; state < state_num; ++state) {
      vector<double> list;
      for (size_t prev_state = 0; prev_state < state_num; ++prev_state) {
        if (transition_c[prev_state][state] > LOG_ZERO
            && forward[pos-1][prev_state] > LOG_ZERO) {
          list.push_back(forward[pos-1][prev_state] +
            transition_c[prev_state][state]);
        }
      }
      if (USE_LOG_ODDS) {
        forward[pos][state] = emission_c[state][curr_baseint]
          - emission_c[model_len][curr_baseint] // supposed to be I_0
          + smithlab::log_sum_log_vec(list, list.size());
      }
      else {
        forward[pos][state] = emission_c[state][curr_baseint]
          + smithlab::log_sum_log_vec(list, list.size());
      }
    }
  }
}

void
ProfileHMM::backward_algorithm_c(const bool VERBOSE,
    const bool USE_LOG_ODDS,
    const string &observation,
    matrix &backward) const {
  const size_t seq_len = observation.length();
  const size_t state_num = model_len*2;
  backward.resize(seq_len, vector<double>(state_num, LOG_ZERO));
  if (VERBOSE)
    cerr << "BACKWARD ALGORITHM " << seq_len << "x" << state_num << endl;
  for (size_t state = 0; state < state_num; ++state)
    backward.back()[state] = 0.0;

  for (signed long pos = seq_len - 2; pos >= 0; --pos) {
    if (VERBOSE && (seq_len - pos) % 10000 == 0)
      cerr << "\tPROCESSED " << 100 - 100 * pos / seq_len << "%" << endl;
    const int next_baseint = base2int(observation[pos+1]);
    for (size_t state = 0; state < state_num; ++state) {
      vector<double> list;
      for (size_t next_state = 0; next_state < state_num; ++next_state) {
        const double transition_prob = transition_c[state][next_state];
        if (transition_prob > LOG_ZERO) {
          if (USE_LOG_ODDS)
            list.push_back(backward[pos+1][next_state]
              + emission_c[next_state][next_baseint]
              - emission_c[model_len][next_baseint] // supposed to be I_0
              + transition_prob);
          else
            list.push_back(backward[pos+1][next_state]
              + emission_c[next_state][next_baseint]
              + transition_prob);
        }
      }
      backward[pos][state] = smithlab::log_sum_log_vec(list, list.size());
    }
  }
}

void
ProfileHMM::Train(const bool VERBOSE,
    const double tolerance,
    const size_t max_iterations,
    const vector<string> &observations) {
  double ll = -1e10, ll_new = -1.0;
  if (VERBOSE)
    cerr << "ITR\tCURR_NLL\tDELTA" << endl;
  matrix e_trans, e_emiss;
  for (size_t itr = 0; itr < max_iterations; ++itr) {
    e_trans.clear();
    e_emiss.clear();
    e_trans.resize(total_size, vector<double>(total_size, LOG_ZERO));
    e_emiss.resize(total_size, vector<double>(4, LOG_ZERO));
    vector<double> ll_list;
    // Expectation
    for (vector<string>::const_iterator seq = observations.begin();
        seq < observations.end(); ++seq) {
      matrix forward, backward;
      forward_algorithm(false, false, *seq, forward);
      backward_algorithm(false, false, *seq, backward);
      ll_list.push_back(posterior_prob(forward));

      train_expectation(VERBOSE, *seq, forward, backward, e_trans, e_emiss);
    }
    ll_new = smithlab::log_sum_log_vec(ll_list, ll_list.size());
    if (VERBOSE)
      cerr << itr + 1 << "/" << max_iterations << "\t"
        << -1.0 * ll_new << "\t" << std::abs((ll_new - ll)/ll_new) << endl;
    if (std::abs((ll_new - ll)/ll_new) < tolerance) break;
    std::swap(ll, ll_new);

    // Maximization
    // transition
    for (unordered_map<size_t, vector<size_t> >::const_iterator from
        = transitions_to.begin();
        from != transitions_to.end(); ++from) {
      vector<double> list;
      for (vector<size_t>::const_iterator to = from->second.begin();
          to < from->second.end(); ++to)
        list.push_back(e_trans[from->first][*to]);
      const double sum = smithlab::log_sum_log_vec(list, list.size());
      for (vector<size_t>::const_iterator to = from->second.begin();
          to < from->second.end(); ++to)
        transition[from->first][*to] = e_trans[from->first][*to] - sum;
    }
    // emission
    for (size_t i = index_m(1); i < index_d(1); ++i) {
      const double sum =
        smithlab::log_sum_log_vec(e_emiss[i], e_emiss[i].size());
      for (size_t j = 0; j < 4; ++j)
        emission[i][j] = e_emiss[i][j] - sum;
    }
    // add pseudocount and normalize
    pseudo_count();
    for (matrix::iterator row = transition.begin();
        row < transition.end() - 1; ++row)
      normalize_vec_inplace(*row, true);
    for (size_t i = index_m(1); i < index_d(1); ++i)
      normalize_vec_inplace(emission[i], true);
  }
  //for (size_t i = index_m(1); i < index_d(1); ++i) {
  //  const double sum =
  //    smithlab::log_sum_log_vec(e_emiss[i], e_emiss[i].size());
  //  cout << i << "\t" << exp(sum) << endl;
  //}
  //collapse_states();
}

void
ProfileHMM::train_expectation(const bool VERBOSE,
    const string &observation,
    const matrix &forward, const matrix &backward,
    matrix &e_trans, matrix &e_emiss) const {
  const size_t seq_len = observation.length();
  const double ll = posterior_prob(forward);
  // transition
  // N's are excluded in training process because they
  // do not provide any information
  for (unordered_map<size_t, vector<size_t> >::const_iterator trans =
      transitions_to.begin(); trans != transitions_to.end(); ++trans) {
    const size_t i = trans->first;
    for (vector<size_t>::const_iterator j = trans->second.begin();
        j < trans->second.end(); ++j) {
      vector<double> list;
      for (size_t pos = 0; pos < seq_len - 1; ++pos) {
        if (base2int(observation[pos])<4 && base2int(observation[pos+1])<4) {
          if (state_can_emit(*j)) {
            list.push_back(forward[pos+1][i] + transition[i][*j]
                + emission[*j][base2int(observation[pos+1])]
                + backward[pos+2][*j]);
          }
          else {
            list.push_back(forward[pos+1][i] + transition[i][*j]
                + backward[pos+1][*j]);
          }
        }
      }
      e_trans[i][*j] = log_sum_log(e_trans[i][*j],
        smithlab::log_sum_log_vec(list, list.size()) - ll);
    }
  }
  // emission
  for (size_t i = index_m(1); i < index_d(1); ++i) {
    matrix list(4, vector<double>());
    for (size_t pos = 0; pos < seq_len - 1; ++pos) {
      const int curr_baseint = base2int(observation[pos]);
      if (curr_baseint<4)
        list[curr_baseint].push_back(
            forward[pos+1][i] + backward[pos+1][i]);
    }
    for (size_t k = 0; k < 4; ++k)
      e_emiss[i][k] = log_sum_log(e_emiss[i][k],
        smithlab::log_sum_log_vec(list[k], list[k].size()) - ll);
  }
}

void
ProfileHMM::SampleSequence(const bool VERBOSE,
    const gsl_rng* rng,
    vector<int> &seq,
    vector<size_t> &states) const {
  size_t idx = 0;
  states.push_back(idx);

  while (idx < total_size - 1) {
    idx = random_weighted_sample(rng, transition[idx]);
    states.push_back(idx);
    if (state_can_emit(idx)) {
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
ProfileHMM::PosteriorDecoding(const bool VERBOSE,
    const bool DEBUG,
    const bool USE_LOG_ODDS,
    const string &observation,
    vector<size_t> &states) const {
  const size_t seq_len = observation.length();
  states.clear();
  matrix forward, backward;

  forward_algorithm(VERBOSE, USE_LOG_ODDS, observation, forward);
  backward_algorithm(VERBOSE, USE_LOG_ODDS, observation, backward);

  // only compute posterior for states M_0 to I_L-1
  vector<double> xi(index_d(1), LOG_ZERO);
  for (size_t pos = 1; pos <= seq_len; ++pos) {
    for (size_t state = index_m(1); state < index_d(1); ++state) {
      xi[state] = forward[pos][state] + backward[pos][state];
    }
    size_t idx = argmax_vec(xi);
    states.push_back(idx);
  }
  if (DEBUG) {
    //cerr << "State sequence:" << endl;
    for (vector<size_t>::const_iterator i = states.begin();
        i < states.end(); ++i)
      cerr << i - states.begin()
        //<< "\t" << *i
        << "\t" << state_idx_to_str(*i) << "\t" << *i << endl;

    //cerr << std::setprecision(4);
    //cerr << endl << "Forward matrix:" << endl;
    //for (size_t state = 0; state < forward.front().size(); ++state)
    //  cerr << "\t" << state_idx_to_str(state);
    //cerr << endl;
    //for (size_t pos = 0; pos < forward.size(); ++pos) {
    //  cerr << pos;
    //  for (size_t state = 0; state < forward.front().size(); ++state)
    //    //cerr << "\t" << exp(forward[pos][state]);
    //    cerr << "\t" << forward[pos][state];
    //  cerr << endl;
    //}
    //cerr << endl << "Backward matrix:" << endl;
    //for (size_t state = 0; state < backward.front().size(); ++state)
    //  cerr << "\t" << state_idx_to_str(state);
    //cerr << endl;
    //for (size_t pos = 0; pos < backward.size(); ++pos) {
    //  cerr << pos;
    //  for (size_t state = 0; state < backward.front().size(); ++state)
    //    //cerr << "\t" << exp(backward[pos][state]);
    //    cerr << "\t" << backward[pos][state];
    //  cerr << endl;
    //}
  }
}

void
ProfileHMM::PosteriorDecoding_c(const bool VERBOSE,
    const bool DEBUG,
    const bool USE_LOG_ODDS,
    const string &observation,
    vector<size_t> &states) const {
  const size_t seq_len = observation.length();
  states.clear();
  matrix forward, backward;

  forward_algorithm_c(VERBOSE, USE_LOG_ODDS, observation, forward);
  backward_algorithm_c(VERBOSE, USE_LOG_ODDS, observation, backward);

  const size_t state_num = model_len*2;
  vector<double> xi(state_num, LOG_ZERO);
  //if (DEBUG) {
  //  cerr << std::setprecision(4);
  //  cerr << "Xi matrix:" << endl;
  //  for (size_t state = 0; state < state_num; ++state)
  //    cerr << "\t" << state;
  //  cerr << endl;
  //}
  for (size_t pos = 0; pos < seq_len; ++pos) {
    //if (DEBUG) cerr << pos;
    for (size_t state = 0; state < state_num; ++state) {
      xi[state] = forward[pos][state] + backward[pos][state];
      //if (DEBUG) cerr << "\t" << xi[state];
    }
    //if (DEBUG) cerr << endl;
    size_t idx = argmax_vec(xi);
    // idx+1 is because in collased matrix, state[0] is M_1, while in the
    // full matrix state[0] is M_0 which is not included in decoding.
    states.push_back(idx);
  }
  if (DEBUG) {
    cerr << endl << "State sequence:" << endl;
    for (vector<size_t>::const_iterator i = states.begin();
        i < states.end(); ++i)
      cerr << i - states.begin() + 1
        << "\t" << *i
        << "\t" << state_idx_to_str_c(*i) << endl;

    //cerr << endl << "Forward matrix:" << endl;
    //for (size_t state = 0; state < forward.front().size(); ++state)
    //  cerr << "\t" << state_idx_to_str_c(state);
    //cerr << endl;
    //for (size_t pos = 0; pos < forward.size(); ++pos) {
    //  cerr << pos;
    //  for (size_t state = 0; state < forward.front().size(); ++state)
    //    //cerr << "\t" << exp(forward[pos][state]);
    //    cerr << "\t" << forward[pos][state];
    //  cerr << endl;
    //}
    //cerr << endl << "Backward matrix:" << endl;
    //for (size_t state = 0; state < backward.front().size(); ++state)
    //  cerr << "\t" << state_idx_to_str_c(state);
    //cerr << endl;
    //for (size_t pos = 0; pos < backward.size(); ++pos) {
    //  cerr << pos;
    //  for (size_t state = 0; state < backward.front().size(); ++state)
    //    //cerr << "\t" << exp(backward[pos][state]);
    //    cerr << "\t" << backward[pos][state];
    //  cerr << endl;
    //}
  }
}

double
ProfileHMM::PosteriorProb(const bool USE_LOG_ODDS,
    const string &observation) const {
  matrix forward;

  forward_algorithm(false, USE_LOG_ODDS, observation, forward);
  return posterior_prob(forward);
}

double
ProfileHMM::PosteriorProb_c(const bool USE_LOG_ODDS,
    const string &observation) const {
  matrix forward;

  forward_algorithm_c(false, USE_LOG_ODDS, observation, forward);
  return posterior_prob(forward);
}

double
ProfileHMM::posterior_prob(const matrix &forward) const {
  vector<double> list;
  for (vector<double>::const_iterator i = forward.back().begin();
      i < forward.back().end(); ++i)
    list.push_back(*i);
  return smithlab::log_sum_log_vec(list, list.size());
}

void
ProfileHMM::pseudo_count(void) {
  const double PSEUDOCOUNT = -5;
  for (unordered_map<size_t, vector<size_t> >::const_iterator from =
      transitions_to.begin();
      from != transitions_to.end(); ++from) {
    for (vector<size_t>::const_iterator to = from->second.begin();
        to < from->second.end(); ++to) {
      transition[from->first][*to] =
        log_sum_log(transition[from->first][*to], PSEUDOCOUNT);
    }
  }
  for (size_t i = index_m(1); i < index_d(1); ++i) {
    for (size_t j = 0; j < 4; ++j) {
      emission[i][j] = log_sum_log(emission[i][j], PSEUDOCOUNT);
    }
  }
}

void
ProfileHMM::get_viable_transitions_to(void) {
  transitions_to.clear();
  // M_0 to I_0, D_1
  transitions_to[0] = vector<size_t>({index_i(0), index_d(1)});
  // I_0 to I_0, D_1, E
  transitions_to[index_i(0)] =
    vector<size_t>({index_i(0), index_d(1)});
  // D_1 to M_1...M_L
  transitions_to[index_d(1)] = vector<size_t>();
  for (size_t i = 1; i <= model_len; ++i)
    transitions_to[index_d(1)].push_back(index_m(i));
  // general: i=1~L-2
  for (size_t i = 1; i < model_len - 1; ++i) {
    transitions_to[index_m(i)] = vector<size_t>({index_m(i+1), index_i(i),
        index_d(i+1), index_d(model_len)});
    transitions_to[index_i(i)] =
      vector<size_t>({index_m(i+1), index_i(i), index_d(i+1)});
    if (i > 1)
      transitions_to[index_d(i)] =
        vector<size_t>({index_m(i+1), index_i(i), index_d(i+1)});
  }
  // M_L-1
  transitions_to[index_m(model_len-1)] =
    vector<size_t>({index_m(model_len), index_i(model_len-1),
        index_d(model_len)});
  // I_L-1
  transitions_to[index_i(model_len-1)] =
    vector<size_t>({index_i(model_len-1), index_m(model_len)});
  // D_L-1
  transitions_to[index_d(model_len-1)] =
    vector<size_t>({index_i(model_len-1), index_m(model_len)});
  // M_L
  transitions_to[index_m(model_len)] = vector<size_t>(1, index_d(model_len));
  // D_L
  transitions_to[index_d(model_len)] =
    vector<size_t>({index_i(0), index_d(1)});
}

void
ProfileHMM::get_viable_transitions_from(void) {
  transitions_from.clear();
  // to I_0: I_0, M_0, D_L
  transitions_from[index_i(0)] = vector<size_t>({index_m(0), index_i(0)
    , index_d(model_len)});
  // to D_L: M_1 to M_L
  transitions_from[index_d(model_len)] = vector<size_t>();
  for (size_t i = 1; i <= model_len; ++i)
    transitions_from[index_d(model_len)].push_back(index_m(i));
  // to D_1: M_0, I_0, D_L
  transitions_from[index_d(1)] = vector<size_t>({index_m(0), index_i(0)
    , index_d(model_len)});
  // to M_1: D_1
  transitions_from[index_m(1)] = vector<size_t>(1, index_d(1));
  // to I_1: M_1, I_1
  transitions_from[index_i(1)] = vector<size_t>({index_m(1), index_i(1)});
  // to M_2: M_1, I_1, D_1
  transitions_from[index_m(2)] =
    vector<size_t>({index_m(1), index_i(1), index_d(1)});
  // to I_2: M_1, I_2, D_2
  transitions_from[index_i(2)] =
    vector<size_t>({index_m(2), index_i(2), index_d(2)});
  // to D_2: M_1, I_1
  transitions_from[index_d(2)] = vector<size_t>({index_m(1), index_i(1)});
  // general: i=3~L-1
  for (size_t i = 3; i < model_len; ++i) {
    transitions_from[index_m(i)] = vector<size_t>({index_m(i-1), index_i(i-1),
        index_d(i-1), index_d(1)});
    transitions_from[index_i(i)] =
      vector<size_t>({index_m(i), index_i(i), index_d(i)});
    transitions_from[index_d(i)] =
      vector<size_t>({index_m(i-1), index_i(i-1), index_d(i-1)});
  }
  // to M_L: M, I, D_L-1, D_1
  transitions_from[index_m(model_len)] = vector<size_t>({index_m(model_len-1),
      index_i(model_len-1), index_d(model_len-1), index_d(1)});
  // to E: I_0, D_L
  transitions_from[total_size - 1] =
    vector<size_t>({index_i(0), index_d(model_len)});
}

bool
ProfileHMM::state_can_emit(const size_t idx) const {
  return idx >= index_m(1) && idx < index_d(1);
}

bool
ProfileHMM::state_can_emit_n(const size_t idx) const {
  return idx >= index_i(0) && idx < index_d(1);
}

void
ProfileHMM::Print(ostream& out, const bool HUM_READABLE) const {
  out << "#Model length: " << model_len << endl;
  out << "#Consensus: " << Consensus() << endl;
  out << "#Transition" << endl;
  print_transition(out, HUM_READABLE);
  out << endl << "#Emission" << endl;
  print_emission(out, HUM_READABLE);
  out << "//" << endl;
}

string
ProfileHMM::Consensus(void) const {
  vector<char> bases;
  for (size_t i = 1; i <= model_len; ++i) {
    const vector<double>::const_iterator max_itr =
      std::max_element(emission[i].begin(), emission[i].end());
    bases.push_back(int2base(max_itr - emission[i].begin()));
  }
  const string consensus(bases.begin(), bases.end());
  return consensus;
}

void
ProfileHMM::DebugOutput(void) const {
  cerr << "Trantision:" << endl;
  cerr << transition;
  cerr << "Trantision collapsed:" << endl;
  cerr << transition_c;
}

void
ProfileHMM::print_transition(ostream& out, const bool HUM_READABLE) const {
  const std::streamsize ss = out.precision();
  if (HUM_READABLE)
    out << std::setprecision(4) << std::fixed;
  out << "#Initial\tI_0\tD_1" << endl << "\t";
  if (HUM_READABLE)
    out
      << exp(transition[0][state(0ul, 0, 1).index(model_len)]) << "\t"
      << exp(transition[0][state(0ul, 1, 2).index(model_len)]) << endl;
  else
    out
      << transition[0][state(0ul, 0, 1).index(model_len)] << "\t"
      << transition[0][state(0ul, 1, 2).index(model_len)] << endl;
  out << "#pos_i:" 
    << "\t" << "M->M" << "\t" << "M->I" << "\t" << "M->D"
    << "\t" << "I->M" << "\t" << "I->I" << "\t" << "I->D"
    << "\t" << "D->M" << "\t" << "D->I" << "\t" << "D->D"
    << endl;
  for (size_t i = 1; i < model_len - 1; ++i) {
    out << i;
    // M
    if (HUM_READABLE)
      out
        << "\t" << exp(transition[state(0ul,i,0).index(model_len)]
          [state(0ul,i+1,0).index(model_len)])
        << "\t" << exp(transition[state(0ul,i,0).index(model_len)]
          [state(0ul,i,1).index(model_len)])
        << "\t" << exp(transition[state(0ul,i,0).index(model_len)]
          [state(0ul,i+1,2).index(model_len)]);
    else
      out
        << "\t" << transition[state(0ul,i,0).index(model_len)]
          [state(0ul,i+1,0).index(model_len)]
        << "\t" << transition[state(0ul,i,0).index(model_len)]
          [state(0ul,i,1).index(model_len)]
        << "\t" << transition[state(0ul,i,0).index(model_len)]
          [state(0ul,i+1,2).index(model_len)];
    // I
    if (HUM_READABLE)
      out
        << "\t" << exp(transition[state(0ul,i,1).index(model_len)]
          [state(0ul,i+1,0).index(model_len)])
        << "\t" << exp(transition[state(0ul,i,1).index(model_len)]
          [state(0ul,i,1).index(model_len)])
        << "\t" << exp(transition[state(0ul,i,1).index(model_len)]
          [state(0ul,i+1,2).index(model_len)]);
    else
      out
        << "\t" << transition[state(0ul,i,1).index(model_len)]
          [state(0ul,i+1,0).index(model_len)]
        << "\t" << transition[state(0ul,i,1).index(model_len)]
          [state(0ul,i,1).index(model_len)]
        << "\t" << transition[state(0ul,i,1).index(model_len)]
          [state(0ul,i+1,2).index(model_len)];
    // D
    if (i > 1) {
      if (HUM_READABLE)
        out
          << "\t" << exp(transition[state(0ul,i,2).index(model_len)]
            [state(0ul,i+1,0).index(model_len)])
          << "\t" << exp(transition[state(0ul,i,2).index(model_len)]
            [state(0ul,i,1).index(model_len)])
          << "\t" << exp(transition[state(0ul,i,2).index(model_len)]
            [state(0ul,i+1,2).index(model_len)]);
      else
        out
          << "\t" << transition[state(0ul,i,2).index(model_len)]
            [state(0ul,i+1,0).index(model_len)]
          << "\t" << transition[state(0ul,i,2).index(model_len)]
            [state(0ul,i,1).index(model_len)]
          << "\t" << transition[state(0ul,i,2).index(model_len)]
            [state(0ul,i+1,2).index(model_len)];
    }
    else
      out << "\tnan\tnan\tnan";
    out << endl;
  }
  out << model_len - 1;
  if (HUM_READABLE)
    out
      // M_L-1
      << "\t" << exp(transition[state(0ul,model_len-1,0).index(model_len)]
        [state(0ul,model_len,0).index(model_len)])
      << "\t" << exp(transition[state(0ul,model_len-1,0).index(model_len)]
        [state(0ul,model_len-1,1).index(model_len)])
      << "\t" << "nan"
      // I_L-1
      << "\t" << exp(transition[state(0ul,model_len-1,1).index(model_len)]
        [state(0ul,model_len,0).index(model_len)])
      << "\t" << exp(transition[state(0ul,model_len-1,1).index(model_len)]
        [state(0ul,model_len-1,1).index(model_len)])
      << "\t" << "nan"
      // D_L-1
      << "\t" << exp(transition[state(0ul,model_len-1,2).index(model_len)]
        [state(0ul,model_len,0).index(model_len)])
      << "\t" << exp(transition[state(0ul,model_len-1,2).index(model_len)]
        [state(0ul,model_len-1,1).index(model_len)])
      << "\t" << "nan";
  else
    out
      << "\t" << transition[state(0ul,model_len-1,0).index(model_len)]
        [state(0ul,model_len,0).index(model_len)]
      << "\t" << transition[state(0ul,model_len-1,0).index(model_len)]
        [state(0ul,model_len-1,1).index(model_len)]
      << "\t" << "nan"
      << "\t" << transition[state(0ul,model_len-1,1).index(model_len)]
        [state(0ul,model_len,0).index(model_len)]
      << "\t" << transition[state(0ul,model_len-1,1).index(model_len)]
        [state(0ul,model_len-1,1).index(model_len)]
      << "\t" << "nan"
      << "\t" << transition[state(0ul,model_len-1,2).index(model_len)]
        [state(0ul,model_len,0).index(model_len)]
      << "\t" << transition[state(0ul,model_len-1,2).index(model_len)]
        [state(0ul,model_len-1,1).index(model_len)]
      << "\t" << "nan";

  out << endl << "#LOG\t5'T->M\tM->3'T" << endl;
  for (size_t i = 1; i <= model_len; ++i) {
    out << i
      << "\t" << transition[state(0ul,1,2).index(model_len)]
        [state(0ul, i, 0).index(model_len)]
      << "\t" << transition[state(0ul,i,0).index(model_len)]
        [state(0ul,model_len,2).index(model_len)]
      << endl;
  }
  if (HUM_READABLE) {
    out << "#I_0\tI_0\tD_1" << endl
      << "\t" << exp(transition[state(0ul,0,1).index(model_len)]
        [state(0ul,0,1).index(model_len)])
      << "\t" << exp(transition[state(0ul,0,1).index(model_len)]
        [state(0ul,1,2).index(model_len)])
      << endl;
    out << "#D_L\tI_0\tD_1" << endl
      << "\t" << exp(transition[state(0ul,model_len,2).index(model_len)]
        [state(0ul,0,1).index(model_len)])
      << "\t" << exp(transition[state(0ul,model_len,2).index(model_len)]
        [state(0ul,1,2).index(model_len)])
      << endl;
  }
  else {
    out << "#I_0\tI_0\tD_1" << endl
      << "\t" << transition[state(0ul,0,1).index(model_len)]
        [state(0ul,0,1).index(model_len)]
      << "\t" << transition[state(0ul,0,1).index(model_len)]
        [state(0ul,1,2).index(model_len)]
      << endl;
    out << "#D_L\tI_0\tD_1" << endl
      << "\t" << transition[state(0ul,model_len,2).index(model_len)]
        [state(0ul,0,1).index(model_len)]
      << "\t" << transition[state(0ul,model_len,2).index(model_len)]
        [state(0ul,1,2).index(model_len)]
      << endl;
  }
  out << std::setprecision(ss);
}

void
ProfileHMM::print_emission(ostream& out, const bool HUM_READABLE) const {
  const std::streamsize ss = out.precision();
  if (HUM_READABLE)
    out << std::setprecision(4) << std::fixed;
  out << "#M_i:" 
    << "\t" << "A" << "\t" << "C" << "\t" << "G" << "\t" << "T" << endl;
  for (size_t i = 1; i <= model_len; ++i) {
    out << i;
    if (HUM_READABLE)
      out
        << "\t" << exp(emission[state(0ul,i,0).index(model_len)][0])
        << "\t" << exp(emission[state(0ul,i,0).index(model_len)][1])
        << "\t" << exp(emission[state(0ul,i,0).index(model_len)][2])
        << "\t" << exp(emission[state(0ul,i,0).index(model_len)][3])
        << endl;
    else
      out
        << "\t" << emission[state(0ul,i,0).index(model_len)][0]
        << "\t" << emission[state(0ul,i,0).index(model_len)][1]
        << "\t" << emission[state(0ul,i,0).index(model_len)][2]
        << "\t" << emission[state(0ul,i,0).index(model_len)][3]
        << endl;
  }

  out << endl << "#I_i:" 
    << "\t" << "A" << "\t" << "C" << "\t" << "G" << "\t" << "T" << endl;
  for (size_t i = 0; i < model_len; ++i) {
    out << i;
    if (HUM_READABLE)
      out
        << "\t" << exp(emission[state(0ul,i,1).index(model_len)][0])
        << "\t" << exp(emission[state(0ul,i,1).index(model_len)][1])
        << "\t" << exp(emission[state(0ul,i,1).index(model_len)][2])
        << "\t" << exp(emission[state(0ul,i,1).index(model_len)][3])
        << endl;
    else
      out
        << "\t" << emission[state(0ul,i,1).index(model_len)][0]
        << "\t" << emission[state(0ul,i,1).index(model_len)][1]
        << "\t" << emission[state(0ul,i,1).index(model_len)][2]
        << "\t" << emission[state(0ul,i,1).index(model_len)][3]
        << endl;
  }
  out << std::setprecision(ss);
}

vector<string>
split(const string s, const char delim = '\t') {
  vector<string> v;
  std::stringstream ss(s);
  string item;
  while(getline(ss, item, delim))
    v.push_back(item);
  return v;
}

void
ProfileHMM::load_from_file(const string filename) {
  std::ifstream in(filename.c_str());
  if (!in)
    throw SMITHLABException("cannot open input file " + filename);

  vector<string> f;
  string line;
  getline(in, line);
  f = split(line, ' ');
  model_len = stoi(f.back());
  total_size = model_len * 3 + 2;
  transition.resize(total_size, vector<double>(total_size, LOG_ZERO));
  emission.resize(total_size, vector<double>(4, LOG_ZERO));

  while(getline(in, line) && (line[0] == '#' || line.empty()));
  f = split(line);
  assert(stof(f[1]) < 0); // make sure the values are log transformed
  transition[0][state(0ul, 0, 1).index(model_len)] = stof(f[1]);
  transition[0][state(0ul, 1, 2).index(model_len)] = stof(f[2]);

  while(getline(in, line) && (line[0] == '#' || line.empty()));
  // general transition
  while(line[0] != '#' && !line.empty()) {
    f = split(line);
    assert(f.size() == 10);
    const size_t i = stoi(f[0]);
    // M_i
    transition[state(0ul,i,0).index(model_len)]
      [state(0ul,i+1,0).index(model_len)] = stof(f[1]);
    transition[state(0ul,i,0).index(model_len)]
      [state(0ul,i,1).index(model_len)] = stof(f[2]);
    transition[state(0ul,i,0).index(model_len)]
      [state(0ul,i+1,2).index(model_len)] =
        f[3] == "nan" ? LOG_ZERO : stof(f[3]);
    // I_i
    transition[state(0ul,i,1).index(model_len)]
      [state(0ul,i+1,0).index(model_len)] = stof(f[4]);
    transition[state(0ul,i,1).index(model_len)]
      [state(0ul,i,1).index(model_len)] = stof(f[5]);
    transition[state(0ul,i,1).index(model_len)]
      [state(0ul,i+1,2).index(model_len)] =
        f[6] == "nan" ? LOG_ZERO : stof(f[6]);
    // D_i
    transition[state(0ul,i,2).index(model_len)]
      [state(0ul,i+1,0).index(model_len)] =
        f[7] == "nan" ? LOG_ZERO : stof(f[7]);
    transition[state(0ul,i,2).index(model_len)]
      [state(0ul,i,1).index(model_len)] =
        f[8] == "nan" ? LOG_ZERO : stof(f[8]);
    transition[state(0ul,i,2).index(model_len)]
      [state(0ul,i+1,2).index(model_len)] =
        f[9] == "nan" ? LOG_ZERO : stof(f[9]);
    
    getline(in, line);
  }

  while(getline(in, line) && (line[0] == '#' || line.empty()));
  // truncation
  while(line[0] != '#' && !line.empty()) {
    f = split(line);
    assert(f.size() == 3);
    assert(stof(f[1]) < 0); // make sure the values are log transformed
    const size_t i = stoi(f[0]);
    transition[state(0ul,1,2).index(model_len)]
      [state(0ul, i, 0).index(model_len)] = stof(f[1]);
    transition[state(0ul,i,0).index(model_len)]
      [state(0ul,model_len,2).index(model_len)] = stof(f[2]);

    getline(in, line);
  }

  while(getline(in, line) && (line[0] == '#' || line.empty()));
  // I_0 transition probability
  f = split(line);
  transition[state(0ul,0,1).index(model_len)]
    [state(0ul,0,1).index(model_len)] = stof(f[1]);
  transition[state(0ul,0,1).index(model_len)]
    [state(0ul,1,2).index(model_len)] = stof(f[2]);
  //transition[state(0ul,0,1).index(model_len)].back() = stof(f[3]);

  while(getline(in, line) && (line[0] == '#' || line.empty()));
  // D_L transition probability
  f = split(line);
  transition[state(0ul,model_len,2).index(model_len)]
    [state(0ul,0,1).index(model_len)] = stof(f[1]);
  transition[state(0ul,model_len,2).index(model_len)]
    [state(0ul,1,2).index(model_len)] = stof(f[2]);

  while(getline(in, line) && (line[0] == '#' || line.empty()));
  // emission
  while(line[0] != '#' && !line.empty()) {
    f = split(line);
    assert(f.size() == 5);
    assert(stof(f[1]) < 0);
    const size_t i = stoi(f[0]);
    for (size_t j = 1; j < f.size(); j++) {
      emission[state(0ul,i,0).index(model_len)][j-1] = stof(f[j]);
    }

    getline(in, line);
  }
  while(getline(in, line) && (line[0] == '#' || line.empty()));
  while(line[0] != '#' && !line.empty() && line[0] != '/') {
    f = split(line);
    assert(f.size() == 5);
    assert(stof(f[1]) < 0);
    const size_t i = stoi(f[0]);
    for (size_t j = 1; j < f.size(); j++)
      emission[state(0ul,i,1).index(model_len)][j-1] = stof(f[j]);

    getline(in, line);
  }
}

void
ProfileHMM::FisherScoreVector(const string &sequence,
    vector<double> &score) const {
  // Organization of the output values (the score vector):
  // 4 states for emission for each state of M_1~M_L and I_0;
  // plus expected counts for each of those states (model_len+1).
  matrix forward, backward;
  forward_algorithm(false, false, sequence, forward);
  backward_algorithm(false, false, sequence, backward);
  // model_len*4: all M states; +4: I_0
  score.resize(model_len*4+4, 0.0);
  const double posterior_p = posterior_prob(forward);

  // scores for states M_1 to M_L, plus I_0
  for (size_t state_index = index_m(1); state_index < index_i(0)+1;
      ++state_index) {
    // posterior count of staying at state_idx
    vector<double> state_count;
    for (size_t i = 1; i < forward.size(); ++i) {
      // the loop is from position 1 to N of the string (1-based coordinates)
      // see forward_algorithm()
      state_count.push_back(forward[i][state_index] + backward[i][state_index]
        - posterior_p);
    }
    const double expected_state_count =
      exp(smithlab::log_sum_log_vec(state_count, state_count.size()));

    vector<vector<double> > nt_count(4, vector<double>(1, LOG_ZERO));
    //vector<double> m_to_next_m, m_to_curr_i, m_to_next_d;
    //const size_t idx_next_m = state_index + 1;
    //const size_t idx_curr_i = state_index + model_len;
    //const size_t idx_next_d = idx_next_m + 1;
    for (size_t i = 0; i < sequence.length()-1; ++i) {
      // posterior emission count per nt per state
      const int nt_idx = base2int(sequence[i]);
      if (nt_idx < 4)
        nt_count[nt_idx].push_back(state_count[i]);
      // posterior transition count per state of interest
      //if (state_index < model_len - 1) {
      //  // M_1 to M_L-1
      //  m_to_next_m.push_back(forward[i][state_index]
      //    + transition_c[state_index][idx_next_m]
      //    + emission_c[idx_next_m][base2int(sequence[i+1])]
      //    + backward[i+1][idx_next_m]);
      //  m_to_curr_i.push_back(forward[i][state_index]
      //    + transition_c[state_index][idx_curr_i]
      //    + emission_c[idx_curr_i][base2int(sequence[i+1])]
      //    + backward[i+1][idx_curr_i]);
      //  if (state_index < model_len - 2)
      //    // M_L-1 does not have next D
      //    m_to_next_d.push_back(forward[i][state_index]
      //      + transition_c[state_index][idx_next_d]
      //      + emission_c[idx_next_d][base2int(sequence[i+1])]
      //      + backward[i+1][idx_next_d]);
      //}
    }
    // scores for emission
    for (size_t nt_idx = 0; nt_idx < 4; ++nt_idx) {
      const double expected_emission_count =
        smithlab::log_sum_log_vec(nt_count[nt_idx], nt_count[nt_idx].size());
      // use state_index-1 is because forward[0] corresponds to M_0
      score[(state_index-1)*4+nt_idx] = exp(expected_emission_count
          - emission[state_index][nt_idx])
        - expected_state_count;
    }
    // scores for expected state count (normalization purpose)
    score.push_back(expected_state_count);
    // scores for transition, M states only
    //if (state_index < model_len - 1) {
    //  score.push_back(exp(smithlab::log_sum_log_vec(m_to_next_m,
    //      m_to_next_m.size()) - posterior_p
    //      - transition_c[state_index][idx_next_m])
    //    - expected_state_count);
    //  score.push_back(exp(smithlab::log_sum_log_vec(m_to_curr_i,
    //      m_to_curr_i.size()) - posterior_p
    //      - transition_c[state_index][idx_curr_i])
    //    - expected_state_count);
    //  if (state_index < model_len - 2)
    //    score.push_back(exp(smithlab::log_sum_log_vec(m_to_next_d,
    //        m_to_next_d.size()) - posterior_p
    //        - transition_c[state_index][idx_next_d])
    //      - expected_state_count);
    //}
  }
  // posterior transition count per state of interest for I_1 to I_L-1
  //for (size_t state_index = model_len+1;state_index < model_len*2;
  //    ++state_index) {
  //  vector<double> state_count;
  //  for (size_t i = 0; i < forward.size() - 1; ++i) {
  //    state_count.push_back(forward[i][state_index] + backward[i][state_index]
  //      - posterior_p);
  //  }
  //  assert(state_count.size() == sequence.length() - 1);
  //  const double expected_state_count =
  //    exp(smithlab::log_sum_log_vec(state_count, state_count.size()));
  //  vector<double> i_to_next_m, i_to_curr_i, i_to_next_d;
  //  const size_t idx_next_m = state_index - model_len;
  //  const size_t idx_curr_i = state_index;
  //  const size_t idx_next_d = idx_next_m + 1;
  //  for (size_t i = 0; i < sequence.length()-1; ++i) {
  //    i_to_next_m.push_back(forward[i][state_index]
  //      + transition_c[state_index][idx_next_m]
  //      + emission_c[idx_next_m][base2int(sequence[i+1])]
  //      + backward[i+1][idx_next_m]);
  //    i_to_curr_i.push_back(forward[i][state_index]
  //      + transition_c[state_index][idx_curr_i]
  //      + emission_c[idx_curr_i][base2int(sequence[i+1])]
  //      + backward[i+1][idx_curr_i]);
  //    if (state_index < model_len * 2 - 1)
  //      // I_L-1 does not have next D
  //      i_to_next_d.push_back(forward[i][state_index]
  //        + transition_c[state_index][idx_next_d]
  //        + emission_c[idx_next_d][base2int(sequence[i+1])]
  //        + backward[i+1][idx_next_d]);
  //  }
  //  score.push_back(exp(smithlab::log_sum_log_vec(i_to_next_m,
  //      i_to_next_m.size()) - posterior_p
  //      - transition_c[state_index][idx_next_m])
  //    - expected_state_count);
  //  score.push_back(exp(smithlab::log_sum_log_vec(i_to_curr_i,
  //      i_to_curr_i.size()) - posterior_p
  //      - transition_c[state_index][idx_curr_i])
  //    - expected_state_count);
  //  if (state_index < model_len * 2 - 1)
  //    score.push_back(exp(smithlab::log_sum_log_vec(i_to_next_d,
  //        i_to_next_d.size()) - posterior_p
  //        - transition_c[state_index][idx_next_d])
  //      - expected_state_count);
  //}
}

void
ProfileHMM::redistribute_prob(matrix &input,
    const size_t row_idx, const vector<size_t> &collapse_cols) {
  // Collapsing a row:
  // Given a matrix, a row to be collapsed and a list of columns to be
  // collapsed, going over the list, first set the value of the cell to 0,
  // then redistribute its outgoing transition probabilities to the 
  // corresponding cells in the same row.
  // E.g. to collapse col 8 of row 5, set row_5 = row_5 + row_8 * cell_5,8;
  // then set cell_5,8 = 0
  // Note: assume there is no self loop for any state to be collapsed.
  for (vector<size_t>::const_iterator col_itr = collapse_cols.begin();
      col_itr < collapse_cols.end(); ++col_itr) {
    if (input[row_idx][*col_itr] > LOG_ZERO) {
      for (size_t i = 0; i < input[row_idx].size(); ++i) {
        if (input[*col_itr][i] > LOG_ZERO)
          input[row_idx][i] = log_sum_log(input[row_idx][i],
              input[row_idx][*col_itr] + input[*col_itr][i]);
      }
      input[row_idx][*col_itr] = LOG_ZERO;
    }
  }
  for (vector<size_t>::const_iterator col_itr = collapse_cols.begin();
      col_itr < collapse_cols.end(); ++col_itr) {
    if (input[row_idx][*col_itr] > LOG_ZERO)
      redistribute_prob(input, row_idx, collapse_cols);
  }
  return;
}

void
ProfileHMM::collapse_states(void) {
  transition_c = transition;
  const size_t first_non_emis = index_d(1);

  vector<size_t> collapse_cols(1, 0);
  for (size_t i = first_non_emis; i < transition[0].size(); ++i)
    collapse_cols.push_back(i);

  for (size_t row = 0; row < first_non_emis; ++row) {
    redistribute_prob(transition_c, row, collapse_cols);
    // delete non-emitting columns (M_0, D_1 - D_L, E)
    transition_c[row].erase(transition_c[row].begin()+first_non_emis,
        transition_c[row].end());
    transition_c[row].erase(transition_c[row].begin());
  }
  // delete non-emitting rows (M_0, D_1 -> D_L), and use row M_0 as initial
  transition_c.erase(transition_c.begin()+first_non_emis,
      transition_c.end());
  initial_prob = transition_c.front();
  transition_c.erase(transition_c.begin());
  // now the new transition matrix should be 2*L x 2*L
  assert(transition_c.size() == model_len*2);
  assert(transition_c[0].size() == model_len*2);
  assert(initial_prob.size() == model_len*2);

  emission_c.erase(emission_c.begin(), emission_c.end());
  for (size_t i = 1; i < first_non_emis; ++i) {
    emission_c.push_back(emission[i]);
  }
  assert(emission_c.size() == model_len * 2);
}
