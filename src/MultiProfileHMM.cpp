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

#include <map>
#include <vector>
#include <utility>
#include <string>

#include "MultiProfileHMM.hpp"

using std::vector;
using std::pair;
using std::make_pair;
using std::map;
using std::cout;
using std::cerr;
using std::endl;
using std::string;

MultiProfileHMM::MultiProfileHMM(vector<ProfileHMM> &v) {
  num_states = 2;
  for (vector<ProfileHMM>::iterator i = v.begin();
      i < v.end(); ++i) {
    models.push_back(&*i);
    num_states += (*i).total_size;
  }
  num_models = models.size();
  transition.resize(num_models + 2,
      vector<double>(num_models + 2, log(1.0/num_models)));
}

void
MultiProfileHMM::Print(void) {
  for (vector<ProfileHMM*>::const_iterator i = models.begin();
      i < models.end(); ++i)
    cout << (**i).total_size << endl;
}

void
MultiProfileHMM::get_viable_transitions_to(void) {
  transitions_to.clear();
  pair<ProfileHMM*, size_t> begin = make_pair(&dummy_this, 0);
  pair<ProfileHMM*, size_t> end = make_pair(&dummy_this, 1);
  transitions_to[begin] = vector<multihmm_state>();
  for (vector<ProfileHMM*>::const_iterator i = models.begin();
      i < models.end(); ++i) {
    // B_0 to B of each model
    transitions_to[begin].push_back(make_pair(*i, (**i).index_m(0)));
    // E of each model to E_0
    pair<ProfileHMM*, size_t> e_i = make_pair(*i, (**i).total_size);
    transitions_to[e_i] = vector<multihmm_state>(1, end);
  }
  // E_0 back to B_0
  transitions_to[end] = vector<multihmm_state>(1, begin);
}

void
MultiProfileHMM::get_viable_transitions_from(void) {
  transitions_from.clear();
  pair<ProfileHMM*, size_t> begin = make_pair(&dummy_this, 0);
  pair<ProfileHMM*, size_t> end = make_pair(&dummy_this, 1);
  // E_0 to B_0
  transitions_from[begin] = vector<multihmm_state>(1, end);
  transitions_to[end] = vector<multihmm_state>();
  for (vector<ProfileHMM*>::const_iterator i = models.begin();
      i < models.end(); ++i) {
    // E of each model to E_0
    transitions_from[end].push_back(make_pair(*i, (**i).total_size));
    // B_0 to B of each model
    pair<ProfileHMM*, size_t> b_i = make_pair(*i, (**i).index_m(0));
    transitions_from[b_i] = vector<multihmm_state>(1, begin);
  }
}

void
MultiProfileHMM::set_transition(void) {
  transition.resize(num_states, vector<double>(num_states, LOG_ZERO));
}

bool
MultiProfileHMM::state_can_emit(multihmm_state state) const {
  if (state.first == &dummy_this) {
    return false;
  }
  else
    return state.first->state_can_emit(state.second);
}

void
MultiProfileHMM::forward_algorithm(const bool VERBOSE,
    const bool USE_LOG_ODDS,
    const string &observation,
    matrix &forward) const {
  if (VERBOSE)
    cerr << "FORWARD ALGORITHM" << endl;
  const size_t seq_len = observation.length();
  forward.resize(seq_len+1, vector<double>(num_states, LOG_ZERO));
  forward[0][0] = 0.0;
  for (size_t pos = 0; pos < seq_len + 1; ++pos) {
    if (VERBOSE && pos % 10000 == 0)
      cerr << "\tPROCESSED " << 100 * pos / seq_len << "%" << endl;
    size_t offset = 1;
    // E_0 to B_0
    if (pos > 0)
      forward[pos][0] = forward[pos].back() + transition.back()[0];
    vector<double> ending_list;
    for (vector<ProfileHMM*>::const_iterator model = models.begin();
        model < models.end(); ++model) {
      // B_0 to B_i
      forward[pos][offset] = forward[pos][0]
        + transition[0][model - models.begin() + 1];
      for (size_t state_idx = 1;
          state_idx < (**model).total_size - 1; ++state_idx) {
        // states with emission: M/I
        const map<size_t, vector<size_t> >::const_iterator i = 
          (**model).transitions_from.find(state_idx);
        if (i != (**model).transitions_from.end()) {
          vector<double> list;
          // use observation[pos-1] because in this vector the position
          // is 0-based
          const int curr_baseint = base2int(observation[pos-1]);
          const size_t curr_state = i->first;
          if ((**model).state_can_emit(curr_state)) {
            // if pos == 0, then only the non-emitting states should
            // be calculated
            if (pos > 0) {
              for (vector<size_t>::const_iterator j = i->second.begin();
                  j < i->second.end(); ++j) {
                list.push_back(forward[pos-1][*j + offset]
                    + (**model).transition[*j][curr_state]);
              }
              if (USE_LOG_ODDS)
                forward[pos][curr_state + offset] =
                  smithlab::log_sum_log_vec(list, list.size())
                  + (**model).emission[curr_state][curr_baseint]
                  - (**model).emission[(**model).index_i(0)][curr_baseint];
              else
                forward[pos][curr_state + offset] =
                  smithlab::log_sum_log_vec(list, list.size())
                  + (**model).emission[curr_state][curr_baseint];
            }
          }
          else {
            for (vector<size_t>::const_iterator j = i->second.begin();
                j < i->second.end(); ++j) {
              list.push_back(forward[pos][*j + offset]
                  + (**model).transition[*j][curr_state]);
            }
            forward[pos][curr_state + offset] =
              smithlab::log_sum_log_vec(list, list.size());
          }
        }
      }
      offset += (**model).total_size;
      // E_i to E_0
      ending_list.push_back(forward[pos][(**model).total_size - 1 + offset]
          + transition[model - models.begin() + 1].back());
    }
    // E_i to E_0
    forward[pos].back() =
      smithlab::log_sum_log_vec(ending_list, ending_list.size());
  }
}

void
MultiProfileHMM::backward_algorithm(const bool VERBOSE,
    const bool USE_LOG_ODDS,
    const string &observation,
    matrix &backward) const {
  if (VERBOSE)
    cerr << "BACKWARD ALGORITHM" << endl;
  const size_t seq_len = observation.length();
  backward.resize(seq_len+1, vector<double>(num_states, LOG_ZERO));
  // Initialization conditions
  for (size_t state_idx = 0; state_idx < num_states; ++state_idx)
    backward.back()[state_idx] = 0.0;

  // loop of backward algorithm starts from D_1^i for all model i,
  // which are the states that have no outgoing transition to any
  // non-emitting states. i.e.
  // backward[pos][D_1^i] = sum_j(backward[pos-1][M_j^i] * ...)
  // otherwise, the backward[pos][state] for an non-emitting state will
  // depend on some other backward[pos][state'] which possibly has not
  // been calculated.
  for (signed long pos = seq_len - 1; pos >= 0; --pos) {
    if (VERBOSE && pos % 10000 == 0)
      cerr << "\tPROCESSED " << 100 - 100 * pos / seq_len << "%" << endl;
    const int next_baseint = base2int(observation[pos]);
    size_t offset = 1;
    vector<double> beginning_list;
    for (vector<ProfileHMM*>::const_iterator model = models.begin();
        model < models.end(); ++model) {
      const size_t idx_d_1 = (**model).index_d(1);
      const size_t idx_i_0 = (**model).index_i(0);
      const map<size_t, vector<size_t> >::const_iterator i = 
        (**model).transitions_to.find(idx_d_1);
      vector<double> list;
      // D_1^i to all M_j^i
      for (vector<size_t>::const_reverse_iterator j = i->second.rbegin();
          j < i->second.rend(); ++j) {
        if (USE_LOG_ODDS)
          list.push_back(backward[pos+1][offset + *j]
              + (**model).transition[idx_d_1][*j]
              + (**model).emission[*j][next_baseint]
              - (**model).emission[idx_i_0][next_baseint]);
        else
          list.push_back(backward[pos+1][offset + *j]
              + (**model).transition[idx_d_1][*j]
              + (**model).emission[*j][next_baseint]);
      }
      backward[pos][idx_d_1 + offset] =
        smithlab::log_sum_log_vec(list, list.size());
      // B^i to D_1^i and I_0^i
      if (USE_LOG_ODDS)
        backward[pos][offset] = log_sum_log(
          backward[pos][idx_d_1 + offset] + (**model).transition[0][idx_d_1],
          backward[pos+1][idx_i_0 + offset]
            + (**model).transition[0][idx_i_0]);
      else
        backward[pos][offset] = log_sum_log(
          backward[pos][idx_d_1 + offset] + (**model).transition[0][idx_d_1],
          backward[pos+1][idx_i_0 + offset] + (**model).transition[0][idx_i_0]
            + (**model).emission[idx_i_0][next_baseint]);

      beginning_list.push_back(backward[pos][offset]
        + transition[0][model - models.begin() + 1]);
      offset += (**model).total_size;
    }
    // B_0 to all B^i
    backward[pos][0] =
      smithlab::log_sum_log_vec(beginning_list, beginning_list.size());
    // E_0 to B_0
    backward[pos].back() = backward[pos][0] + transition.back()[0];
    // E^i to E_0 and D_L^i to E^i
    offset = 1;
    for (vector<ProfileHMM*>::const_iterator model = models.begin();
        model < models.end(); ++model) {
      const size_t idx_d_L = (**model).index_d((**model).model_len);
      const size_t idx_e = (**model).total_size - 1;
      backward[pos][offset + idx_e] =
        backward[pos].back() + transition[model - models.begin() + 1].back();
      backward[pos][offset + idx_d_L] = backward[pos][offset + idx_e]
        + (**model).transition[idx_d_L].back();

      offset += (**model).total_size;
    }
    // other internal states of individual models
    offset = 1;
    for (vector<ProfileHMM*>::const_iterator model = models.begin();
        model < models.end(); ++model) {
      for (signed long state_idx = (**model).total_size - 1;
          state_idx >= 0; --state_idx) {
        const map<size_t, vector<size_t> >::const_iterator i = 
          (**model).transitions_to.find(state_idx);
        // skip states that have been calculated
        if (i != (**model).transitions_to.end()
            && i->first != (**model).index_d(1)
            && i->first != (**model).index_d((**model).model_len)
            && i->first != 0
           ) {
          vector<double> list;
          const size_t curr_state = i->first;
          for (vector<size_t>::const_reverse_iterator j = i->second.rbegin();
              j < i->second.rend(); ++j) {
            if ((**model).state_can_emit(*j)) {
              // *j is an M/I state with viable emission
              if (USE_LOG_ODDS)
                list.push_back(backward[pos+1][offset + *j]
                    + (**model).transition[curr_state][*j]
                    + (**model).emission[*j][next_baseint]
                    - (**model).emission[(**model).index_i(0)][next_baseint]);
              else
                list.push_back(backward[pos+1][offset + *j]
                    + (**model).transition[curr_state][*j]
                    + (**model).emission[*j][next_baseint]);
            }
            else {
              // *j is a D state without viable emission
              list.push_back(backward[pos][offset + *j]
                  + (**model).transition[curr_state][*j]);
            }
          }
          backward[pos][offset + curr_state] =
            smithlab::log_sum_log_vec(list, list.size());
        }
      }
      offset += (**model).total_size;
    }
  }
}
