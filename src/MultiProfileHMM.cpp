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

#include <vector>
#include <utility>
#include <string>

#include "MultiProfileHMM.hpp"

using std::vector;
using std::pair;
using std::make_pair;
using std::unordered_map;
using std::cout;
using std::cerr;
using std::endl;
using std::string;

MultiProfileHMM::MultiProfileHMM(vector<ProfileHMM> &v) {
  num_states = 2;
  num_columns = 0;
  size_t offset = 0;
  for (vector<ProfileHMM>::iterator i = v.begin();
      i < v.end(); ++i) {
    models.push_back(&*i);
    num_states += (*i).total_size;
    num_columns += (*i).model_len;
    offset += (*i).model_len*2;
    xi_idx_upper.push_back(offset);
  }
  num_models = models.size();
  set_transition();
}

void
MultiProfileHMM::Print(void) const {
  for (vector<ProfileHMM*>::const_iterator i = models.begin();
      i < models.end(); ++i)
    cout << (**i).total_size << endl;
}

void
MultiProfileHMM::get_viable_transitions_to(void) {
  transitions_to.clear();
  const multihmm_state begin = make_pair(&dummy_this, 0);
  const multihmm_state end = make_pair(&dummy_this, 1);
  transitions_to[begin] = vector<multihmm_state>();
  for (vector<ProfileHMM*>::const_iterator i = models.begin();
      i < models.end(); ++i) {
    // B_0 to B of each model
    transitions_to[begin].push_back(make_pair(*i, (**i).index_m(0)));
    // E of each model to E_0
    multihmm_state e_i = make_pair(*i, (**i).total_size);
    transitions_to[e_i] = vector<multihmm_state>(1, end);
  }
  // E_0 back to B_0
  transitions_to[end] = vector<multihmm_state>(1, begin);
}

void
MultiProfileHMM::get_viable_transitions_from(void) {
  transitions_from.clear();
  const multihmm_state begin = make_pair(&dummy_this, 0);
  const multihmm_state end = make_pair(&dummy_this, 1);
  // E_0 to B_0
  transitions_from[begin] = vector<multihmm_state>(1, end);
  transitions_to[end] = vector<multihmm_state>();
  for (vector<ProfileHMM*>::const_iterator i = models.begin();
      i < models.end(); ++i) {
    // E of each model to E_0
    transitions_from[end].push_back(make_pair(*i, (**i).total_size));
    // B_0 to B of each model
    multihmm_state b_i = make_pair(*i, (**i).index_m(0));
    transitions_from[b_i] = vector<multihmm_state>(1, begin);
  }
}

void
MultiProfileHMM::set_transition(void) {
  transition.resize(num_models+2, vector<double>(num_models+2, LOG_ZERO));
  // Set transitions of B_0 to B^i
  for (vector<double>::iterator i = transition.front().begin()+1;
      i < transition.front().end()-1; ++i)
    *i = log(1.0/num_models);
  // Set transitions of E^i to E_0
  for (matrix::iterator i = transition.begin()+1;
      i < transition.end()-1; ++i)
    (*i).back() = 0.0;
  // Set transition of E_0 to B_0
  transition.back().front() = 0.0;
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
  /* Notations
   * 1. B_0, E_0: the structural states that wraps all states of individual
   *    family models;
   * 2. B^i, E^i: the B and E states of individual model i;
   * 3. M_j^i, I_j^i, D_j^i: the M, I, D states of individual model i.
   *
   * Matrix organization
   * The forward matrix is a NxM matrix, where N is the length of
   * the observation sequence, and M is calculated by:
   * M = 2 + \sum_i L_i, i.e.
   * B_0, B^1, M_1^1, ..., E^1, B^2, ..., E^k, E_0
   */
  if (VERBOSE)
    cerr << "FORWARD ALGORITHM" << endl;
  const size_t seq_len = observation.length();
  forward.resize(seq_len+1, vector<double>(num_states, LOG_ZERO));
  forward[0][0] = 0.0;
  // calculate B^i, D_1^i as initialization
  size_t offset = 1;
  for (vector<ProfileHMM*>::const_iterator model = models.begin();
      model < models.end(); ++model) {
    forward[0][offset] = forward[0][0] + transition[0][model-models.begin()+1];
    const size_t idx_d_1 = (**model).index_d(1);
    forward[0][offset+idx_d_1] = forward[0][offset]
      + (**model).transition[0][idx_d_1];
    offset += (**model).total_size;
  }

  // loop of forward algorithm starts from M_1^i for all model i,
  // going through all emitting states, then non-emitting states:
  // forward[pos][M_1^i] = sum_j(forward[pos-1][B^i] * transition...)
  // otherwise, the forward[pos][state] for a state might depend on
  // some other non-emitting state forward[pos][state'], which possibly
  // has not been calculated.
  for (size_t pos = 0; pos < seq_len + 1; ++pos) {
    if (VERBOSE && pos % 10000 == 0)
      cerr << "\tPROCESSED " << 100 * pos / seq_len << "%" << endl;
    offset = 1;
    vector<double> ending_list;
    for (vector<ProfileHMM*>::const_iterator model = models.begin();
        model < models.end(); ++model) {
      for (size_t state_idx = 1;
          state_idx < (**model).total_size - 1; ++state_idx) {
        // states with emission: M/I
        const unordered_map<size_t, vector<size_t> >::const_iterator i = 
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
      // D_L and I_0 to E^i
      const size_t idx_d_l = (**model).index_d((**model).model_len);
      const size_t idx_i_0 = (**model).index_i(0);
      const size_t idx_e_i = (**model).total_size-1;
      forward[pos][idx_e_i+offset] = log_sum_log(
          forward[pos][idx_d_l+offset] + (**model).transition[idx_d_l][idx_e_i],
          forward[pos][idx_i_0+offset] + (**model).transition[idx_i_0][idx_e_i]
          );
      // E^i to E_0
      ending_list.push_back(forward[pos][idx_e_i+offset]
          + transition[model-models.begin()+1].back());
      offset += (**model).total_size;
    }
    // E^i to E_0
    forward[pos].back() =
      smithlab::log_sum_log_vec(ending_list, ending_list.size());
    // E_0 to B_0
    if (pos > 0)
      forward[pos][0] = forward[pos].back() + transition.back()[0];
    // B_0 to B^i and D_1^i
    offset = 1;
    for (vector<ProfileHMM*>::const_iterator model = models.begin();
        model < models.end(); ++model) {
      const size_t idx_d_1 = (**model).index_d(1);
      const size_t idx_i_0 = (**model).index_i(0);
      forward[pos][offset] = forward[pos][0]
        + transition[0][model-models.begin()+1];
      forward[pos][idx_d_1+offset] = log_sum_log(
          forward[pos][offset] + (**model).transition[0][idx_d_1],
          forward[pos][idx_i_0+offset] + (**model).transition[idx_i_0][idx_d_1]
          );
      offset += (**model).total_size;
    }
  }
}

void
MultiProfileHMM::backward_algorithm(const bool VERBOSE,
    const bool USE_LOG_ODDS,
    const string &observation,
    matrix &backward) const {
  /* Notations:
   * 1. B_0, E_0: the structural states that wraps all states of individual
   *    family models;
   * 2. B^i, E^i: the B and E states of individual model i;
   * 3. M_j^i, I_j^i, D_j^i: the M, I, D states of individual model i.
   */
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
      const unordered_map<size_t, vector<size_t> >::const_iterator i = 
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
        const unordered_map<size_t, vector<size_t> >::const_iterator i = 
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

void
MultiProfileHMM::PosteriorDecoding(const bool VERBOSE,
    const bool DEBUG,
    const bool USE_LOG_ODDS,
    const string &observation,
    vector<size_t> &states) const {
  const size_t seq_len = observation.length();
  states.clear();
  matrix forward, backward;

  forward_algorithm(VERBOSE, USE_LOG_ODDS, observation, forward);
  backward_algorithm(VERBOSE, USE_LOG_ODDS, observation, backward);

  vector<double> xi(num_columns*2, LOG_ZERO);
  for (size_t pos = 1; pos <= seq_len; ++pos) {
    size_t state_offset = 2;
    size_t xi_offset = 0;
    for (vector<ProfileHMM*>::const_iterator model = models.begin();
        model < models.end(); ++model) {
      for (size_t state = (**model).index_m(1);
          state < (**model).index_d(1); ++state) {
        xi[state+xi_offset] = forward[pos][state+state_offset]
          + backward[pos][state+state_offset];
      }
      state_offset += (**model).total_size;
      xi_offset += (**model).model_len*2;
    }
    size_t idx = argmax_vec(xi);
    states.push_back(idx);
  }
  if (DEBUG) {
    cerr << "State sequence:" << endl;
    for (vector<size_t>::const_iterator i = states.begin();
        i < states.end(); ++i)
      //cerr << i - states.begin() + 1 << "\t" << *i << endl;
      cerr << i - states.begin() + 1 << "\t" << xi_idx_to_str(*i) << endl;

    //cerr << endl << "Forward matrix:" << endl;
    //for (size_t state = 0; state < forward.front().size(); ++state)
    //  cerr << "\t" << state;
    //cerr << endl;
    //for (size_t pos = 0; pos < forward.size(); ++pos) {
    //  cerr << pos;
    //  for (size_t state = 0; state < forward.front().size(); ++state)
    //    cerr << "\t" << forward[pos][state];
    //  cerr << endl;
    //}
    //cerr << endl << "Backward matrix:" << endl;
    //for (size_t state = 0; state < backward.front().size(); ++state)
    //  cerr << "\t" << state;
    //cerr << endl;
    //for (size_t pos = 0; pos < backward.size(); ++pos) {
    //  cerr << pos;
    //  for (size_t state = 0; state < backward.front().size(); ++state)
    //    cerr << "\t" << backward[pos][state];
    //  cerr << endl;
    //}
  }
}

string
MultiProfileHMM::xi_idx_to_str(const size_t idx) const {
  size_t offset = 0;
  vector<ProfileHMM*>::const_iterator model = models.begin();
  while (model < models.end() && offset+(**model).model_len*2 < idx) {
    offset += (**model).model_len*2;
    ++model;
  }
  // the idx here is the xi matrix idx. The xi matrix is composed of all
  // emitting states from individual models, with each model having 2*model_len
  // number of states (M_1~M_L+I_0~I_L-1). Therefore the offset is added up
  // by 2*model_len per model.
  // To obtain the state name of individual models from this idx, 1 needs to
  // be added because in each model there is the B state in the transition
  // matrix before any other states.
  return (**model).name + "." + (**model).state_idx_to_str(idx-offset+1);
}

size_t
MultiProfileHMM::which_model(const size_t idx) const {
  // the idx is the index of the xi matric, not the grand forward matrix.
  vector<size_t>::const_iterator which =
    std::upper_bound(xi_idx_upper.begin(), xi_idx_upper.end(), idx);
  return which - xi_idx_upper.begin();
}

void
MultiProfileHMM::ComplementBackground(void) {
  for (vector<ProfileHMM*>::iterator i = models.begin();
      i < models.end(); ++i) {
    (**i).ComplementBackground();
  }
}
