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

#ifndef MULTIPHMM_HPP
#define MULTIPHMM_HPP

#include "ProfileHMM.hpp"

typedef std::pair<ProfileHMM*, size_t> multihmm_state;

class MultiProfileHMM {
public:
  MultiProfileHMM();
  MultiProfileHMM(std::vector<ProfileHMM> &v);

  void
  Print(void);

  //operator+;
private:
  void
  get_viable_transitions_to(void);

  void
  get_viable_transitions_from(void);

  bool
  state_can_emit(const multihmm_state state) const;

  std::vector<ProfileHMM*> models;
  matrix transition;
  size_t num_models;
  size_t num_states;

  std::map<multihmm_state, std::vector<multihmm_state> > transitions_to;
  std::map<multihmm_state, std::vector<multihmm_state> > transitions_from;
  ProfileHMM dummy_this;
};
#endif
