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

#include "ProfileHMM.hpp"

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

using std::vector;
using std::pair;
using std::max;
using std::string;
using std::endl;
using std::cout;
using std::cerr;
using std::ostream;
using std::cerr;
using std::map;

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
}

ProfileHMM::ProfileHMM(const string inf) {
  load_from_file(inf);
  get_viable_transitions_to();
  get_viable_transitions_from();
  name = inf.substr(0, inf.find_last_of('.'));
}

ProfileHMM::ProfileHMM(const ProfileHMM &other) {
  model_len = other.model_len;
  total_size = other.total_size;
  transition = other.transition;
  emission = other.emission;
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
  assert(idx >= 0 && idx <= model_len);
  return idx;
}

size_t
ProfileHMM::index_i(const size_t idx) const {
  // I_0 ~ I_L-1 + (L+1)*M
  assert(idx >= 0 && idx <= model_len - 1);
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
  // the forward matrix f[state_idx][pos] is an L x N matrix
  if (VERBOSE)
    cerr << "FORWARD ALGORITHM" << endl;
  const size_t seq_len = observation.length();
  forward.resize(seq_len+1, vector<double>(total_size, LOG_ZERO));

  forward[0][0] = 0;

  for (size_t pos = 1; pos < seq_len + 1; ++pos) {
    if (VERBOSE && pos % 10000 == 0)
      cerr << "\tPROCESSED " << 100 * pos / seq_len << "%" << endl;
    for (size_t state_idx = 0; state_idx < total_size; ++state_idx) {
      // states with emission: M/I
      const map<size_t, vector<size_t> >::const_iterator i = 
        transitions_from.find(state_idx);
      if (i != transitions_from.end()) {
        vector<double> list;
        // use observation[pos-1] because in this vector the position
        // is 0-based
        const int curr_baseint = base2int(observation[pos-1]);
        const size_t curr_state = i->first;
        if (state_can_emit(curr_state)) {
          // if pos == 0, then only the non-emitting states should
          // be calculated
          if (pos > 0) {
            for (vector<size_t>::const_iterator j = i->second.begin();
                j < i->second.end(); ++j) {
              list.push_back(forward[pos-1][*j] + transition[*j][curr_state]);
            }
            if (USE_LOG_ODDS)
              forward[pos][curr_state] =
                smithlab::log_sum_log_vec(list, list.size())
                + emission[curr_state][curr_baseint]
                - emission[index_i(0)][curr_baseint];
            else
              forward[pos][curr_state] =
                smithlab::log_sum_log_vec(list, list.size())
                + emission[curr_state][curr_baseint];
          }
        }
        else {
          for (vector<size_t>::const_iterator j = i->second.begin();
              j < i->second.end(); ++j) {
            list.push_back(forward[pos][*j] + transition[*j][curr_state]);
          }
          forward[pos][curr_state] =
            smithlab::log_sum_log_vec(list, list.size());
        }
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
  for (signed long pos = seq_len - 1; pos > 0; --pos) {
    if (VERBOSE && pos % 10000 == 0)
      cerr << "\tPROCESSED " << 100 - 100 * pos / seq_len << "%" << endl;
    for (signed long state_idx = total_size - 1;
        state_idx >= 0; --state_idx) {
      const map<size_t, vector<size_t> >::const_iterator i = 
        transitions_to.find(state_idx);
      if (i != transitions_to.end()) {
        vector<double> list;
        const int next_baseint = base2int(observation[pos]);
        const size_t curr_state = i->first;
        for (vector<size_t>::const_reverse_iterator j = i->second.rbegin();
            j < i->second.rend(); ++j) {
          if (state_can_emit(*j)) {
            // *j is an M/I state with viable emission
            if (USE_LOG_ODDS)
              list.push_back(backward[pos+1][*j]
                  + transition[curr_state][*j]
                  + emission[*j][next_baseint]
                  - emission[index_i(0)][next_baseint]);
            else
              list.push_back(backward[pos+1][*j]
                  + transition[curr_state][*j]
                  + emission[*j][next_baseint]);
          }
          else {
            // *j is a D state without viable emission
            // this assert makes sure that
            // backward[pos][*j] has been calculated
            assert(*j > curr_state);
            list.push_back(backward[pos][*j]
                + transition[curr_state][*j]);
          }
        }
        backward[pos][curr_state] =
          smithlab::log_sum_log_vec(list, list.size());
      }
    }
  }
}

void
ProfileHMM::Train(const bool VERBOSE,
    const double tolerance,
    const size_t max_iterations,
    const string &observation) {
  const size_t seq_len = observation.length();
  matrix forward, backward;
  forward_algorithm(VERBOSE, false, observation, forward);
  backward_algorithm(VERBOSE, false, observation, backward);
  double ll = -1e10;
  double ll_new = posterior_prob(forward);
  if (VERBOSE)
    cerr << "ITR\tPREV_LL\tCURR_LL\tDELTA" << endl;
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
          if (state_can_emit(*j)) {
            if (pos < seq_len)
              list.push_back(forward[pos][i] + transition[i][*j]
                  + emission[*j][base2int(observation[pos])]
                  + backward[pos + 1][*j]);
          }
          else {
            list.push_back(forward[pos][i] + transition[i][*j]
                + backward[pos][*j]);
          }
        }
        e_trans[i][*j] = smithlab::log_sum_log_vec(list, list.size()) - ll;
      }
    }
    // emission
    for (size_t i = index_m(1); i < index_d(1); ++i) {
      matrix list(4, vector<double>());
      for (size_t pos = 0; pos < seq_len; ++pos) {
        list[base2int(observation[pos])].push_back(
            forward[pos+1][i] + backward[pos+1][i]);
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
    pseudo_count();
    for (matrix::iterator row = transition.begin();
        row < transition.end() - 1; ++row)
      normalize_vec_inplace(*row, true);
    for (size_t i = index_m(1); i < index_d(1); ++i)
      normalize_vec_inplace(emission[i], true);

    // update forward and backward prob using updated parameters
    forward_algorithm(false, false, observation, forward);
    backward_algorithm(false, false, observation, backward);
    // Need some function for posterior probability here!
    ll_new = posterior_prob(forward);
    if (VERBOSE) {
      cerr << itr + 1 << "/" << max_iterations << "\t"
        << ll << "\t" << ll_new << "\t" 
        << std::abs((ll_new - ll)/ll_new) << endl;
    }
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

  vector<double> xi(index_d(1), LOG_ZERO);
  for (size_t pos = 1; pos <= seq_len; ++pos) {
    for (size_t state = index_m(1); state < index_d(1); ++state) {
      xi[state] = forward[pos][state] + backward[pos][state];
    }
    size_t idx = argmax_vec(xi);
    states.push_back(idx);
  }
  if (DEBUG) {
    cerr << "State sequence:" << endl;
    for (vector<size_t>::const_iterator i = states.begin();
        i < states.end(); ++i)
      cerr << i - states.begin() + 1 << "\t" << state_idx_to_str(*i) << endl;

    cerr << std::setprecision(4);
    cerr << endl << "Forward matrix:" << endl;
    for (size_t state = 0; state < forward.front().size(); ++state)
      cerr << "\t" << state_idx_to_str(state);
    cerr << endl;
    for (size_t pos = 0; pos < forward.size(); ++pos) {
      cerr << pos;
      for (size_t state = 0; state < forward.front().size(); ++state)
        //cerr << "\t" << exp(forward[pos][state]);
        cerr << "\t" << forward[pos][state];
      cerr << endl;
    }
    cerr << endl << "Backward matrix:" << endl;
    for (size_t state = 0; state < backward.front().size(); ++state)
      cerr << "\t" << state_idx_to_str(state);
    cerr << endl;
    for (size_t pos = 0; pos < backward.size(); ++pos) {
      cerr << pos;
      for (size_t state = 0; state < backward.front().size(); ++state)
        //cerr << "\t" << exp(backward[pos][state]);
        cerr << "\t" << backward[pos][state];
      cerr << endl;
    }
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
  for (map<size_t, vector<size_t> >::const_iterator from =
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
    vector<size_t>({index_i(0), index_d(1), total_size - 1});
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
    vector<size_t>({index_i(0), total_size - 1});
}

void
ProfileHMM::get_viable_transitions_from(void) {
  transitions_from.clear();
  // to I_0: I_0, M_0, D_L
  transitions_from[index_i(0)] =
    vector<size_t>({index_i(0), index_m(0), index_d(model_len)});
  // to D_1: M_0, I_0
  transitions_from[index_d(1)] = vector<size_t>({index_m(0), index_i(0)});
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
  // to D_L: M_1 to M_L
  transitions_from[index_d(model_len)] = vector<size_t>();
  for (size_t i = 1; i <= model_len; ++i)
    transitions_from[index_d(model_len)].push_back(index_m(i));
  // to E: I_0, D_L
  transitions_from[total_size - 1] =
    vector<size_t>({index_i(0), index_d(model_len)});
}

bool
ProfileHMM::state_can_emit(const size_t idx) const {
  return idx >= index_m(1) && idx < index_d(1);
}

void
ProfileHMM::Print(ostream& out, const bool HUM_READABLE) const {
  out << "#Model length: " << model_len << endl;
  out << "#Transition" << endl;
  print_transition(out, HUM_READABLE);
  out << endl << "#Emission" << endl;
  print_emission(out, HUM_READABLE);
  out << "//" << endl;
}

void
ProfileHMM::print_transition(ostream& out, const bool HUM_READABLE) const {
  const std::streamsize ss = out.precision();
  if (HUM_READABLE)
    out << std::setprecision(4) << std::fixed;
  out << "#M_0\tI_0\tD_1" << endl << "\t";
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
    out << "#I_0\tI_0\tD_1\tE" << endl
      << "\t" << exp(transition[state(0ul,0,1).index(model_len)]
        [state(0ul,0,1).index(model_len)])
      << "\t" << exp(transition[state(0ul,0,1).index(model_len)]
        [state(0ul,1,2).index(model_len)])
      << "\t" << exp(transition[state(0ul,0,1).index(model_len)].back())
      << endl;
    out << "#D_L\tI_0\tE" << endl
      << "\t" << exp(transition[state(0ul,model_len,2).index(model_len)]
        [state(0ul,0,1).index(model_len)])
      << "\t" << exp(transition[state(0ul,model_len,2).index(model_len)].back())
      << endl;
  }
  else {
    out << "#I_0\tI_0\tD_1\tE" << endl
      << "\t" << transition[state(0ul,0,1).index(model_len)]
        [state(0ul,0,1).index(model_len)]
      << "\t" << transition[state(0ul,0,1).index(model_len)]
        [state(0ul,1,2).index(model_len)]
      << "\t" << transition[state(0ul,0,1).index(model_len)].back()
      << endl;
    out << "#D_L\tI_0\tE" << endl
      << "\t" << transition[state(0ul,model_len,2).index(model_len)]
        [state(0ul,0,1).index(model_len)]
      << "\t" << transition[state(0ul,model_len,2).index(model_len)].back()
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
  transition[0][state(0ul, 0, 1).index(model_len)] = stof(f[1]);
  transition[0][state(0ul, 1, 2).index(model_len)] = stof(f[2]);

  while(getline(in, line) && (line[0] == '#' || line.empty()));
  // general transition
  while(line[0] != '#' && !line.empty()) {
    f = split(line);
    assert(f.size() == 10);
    assert(stof(f[1]) < 0); // make sure the values are log transformed
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
  // background transition probability
  f = split(line);
  transition[state(0ul,0,1).index(model_len)]
    [state(0ul,0,1).index(model_len)] = stof(f[1]);
  transition[state(0ul,0,1).index(model_len)]
    [state(0ul,1,2).index(model_len)] = stof(f[2]);
  transition[state(0ul,0,1).index(model_len)].back() = stof(f[3]);

  while(getline(in, line) && (line[0] == '#' || line.empty()));
  // ending transition probability
  f = split(line);
  transition[state(0ul,model_len,2).index(model_len)]
    [state(0ul,0,1).index(model_len)] = stof(f[1]);
  transition[state(0ul,model_len,2).index(model_len)].back() = stof(f[2]);

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
