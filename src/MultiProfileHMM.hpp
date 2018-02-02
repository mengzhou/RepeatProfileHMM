/*    MultiProfileHMM: class implementation of combining multiple
 *    profile-HMMs.
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

#ifndef MULTIPHMM_HPP
#define MULTIPHMM_HPP

#include <unordered_map>

#include "ProfileHMM.hpp"

typedef std::pair<ProfileHMM*, size_t> multihmm_state;

namespace std {
  template <>
  struct hash<multihmm_state> : public unary_function<multihmm_state, size_t> {
    size_t operator()(const multihmm_state& multi_state) const {
      const size_t h1(std::hash<std::string>{}(multi_state.first->Name()));
      const size_t h2(std::hash<size_t>{}(multi_state.first->Length()));
      // just some simple combination of the hash values
      return h1 ^ (h2 << 2) ^ (multi_state.second << 1);
    }
  };
}

class MultiProfileHMM {
public:
  MultiProfileHMM();
  MultiProfileHMM(std::vector<ProfileHMM> &v);

  void
  PosteriorDecoding(const bool VERBOSE,
      const bool DEBUG,
      const bool USE_LOG_ODDS,
      const std::string &observation,
      std::vector<size_t> &states) const;

  void
  ComplementBackground(void);

  void
  Print() const;

  size_t
  which_model(const size_t xi_idx) const;

  std::vector<ProfileHMM*>::const_iterator
  begin() const { return models.cbegin(); }
  std::vector<ProfileHMM*>::const_iterator
  end() const { return models.cend(); }

private:
  void
  forward_algorithm(const bool VERBOSE,
      const bool USE_LOG_ODDS,
      const std::string &observation,
      matrix &forward) const;

  void
  backward_algorithm(const bool VERBOSE,
      const bool USE_LOG_ODDS,
      const std::string &observation,
      matrix &backward) const;

  void
  get_viable_transitions_to(void);

  void
  get_viable_transitions_from(void);

  void
  set_transition(void);

  bool
  state_can_emit(const multihmm_state state) const;

  std::string
  xi_idx_to_str(const size_t idx) const;

  std::vector<ProfileHMM*> models;
  // this vector stores *xi* indices of the last emitting state for individual
  // models, e.g.
  // (model[0].model_len*2, model[0].model_len*2+model[1].model_len*2, ...)
  // used for posterior decoding
  std::vector<size_t> xi_idx_upper;
  matrix transition;
  size_t num_models;
  size_t num_states;
  size_t num_columns; // sum of individual model lengths

  std::unordered_map<multihmm_state, std::vector<multihmm_state> > transitions_to;
  std::unordered_map<multihmm_state, std::vector<multihmm_state> > transitions_from;
  ProfileHMM dummy_this;
};
#endif
