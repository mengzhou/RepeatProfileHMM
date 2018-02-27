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

#ifndef PHMM_HPP
#define PHMM_HPP

#include <unordered_map>
#include <gsl/gsl_randist.h>

#include "smithlab_utils.hpp"
#include "matrix_utils.hpp"

size_t
baseint2stateint(const size_t &baseint, const bool marked);

std::string
state_type_to_str(const size_t model_len, const size_t idx);

struct state {
  state();
  state(const size_t b, const size_t i, const size_t s) :
    baseint(b), idx(i), stateint(s) {}
  state(const char b, const size_t i, const bool m);
  state(const state &obj);
  state(const size_t model_len, const size_t matrix_idx);

  size_t baseint;
  size_t idx;
  size_t stateint;

  bool isvalid(void) const;
  size_t index(const size_t model_len) const;
};

class ProfileHMM {
  // profile-HMM allowing local alignment and re-occurrence of M_i;
  // see Durbin book p.114 (bottom figure)
friend class MultiProfileHMM;
public:
  ProfileHMM();
  ProfileHMM(const matrix &t, const matrix &e);
  ProfileHMM(const std::string inf);
  ProfileHMM(const ProfileHMM &other);
  ProfileHMM &operator=(const ProfileHMM &rhs);

  double
  ViterbiDecoding(const bool VERBOSE,
      const std::string &observation,
      std::vector<std::pair<char, size_t> > &trace) const;

  void
  Train(const bool VERBOSE,
      const double tolerance,
      const size_t max_iterations,
      const std::vector<std::string> &observations);

  void
  SampleSequence(const bool VERBOSE,
      const gsl_rng* rng,
      std::vector<int> &seq,
      std::vector<size_t> &states) const;

  void
  PosteriorDecoding(const bool VERBOSE,
      const bool DEBUG,
      const bool USE_LOG_ODDS,
      const std::string &observation,
      std::vector<size_t> &states) const;

  void
  PosteriorDecoding_c(const bool VERBOSE,
      const bool DEBUG,
      const bool USE_LOG_ODDS,
      const std::string &observation,
      std::vector<size_t> &states) const;

  double
  PosteriorProb(const bool USE_LOG_ODDS,
      const std::string &observation) const;

  double
  PosteriorProb_c(const bool USE_LOG_ODDS,
      const std::string &observation) const;

  void
  FisherScoreVector(const std::string &sequence,
      std::vector<double> &score) const;

  void
  Print(std::ostream& out, const bool HUM_READABLE) const;

  size_t
  Length(void) const {return model_len;}

  std::string
  Name(void) const {return name;}

  void
  ComplementBackground(void);
  
  std::string
  Consensus(void) const;

  void
  DebugOutput() const;

private:
  size_t
  index_m(const size_t idx) const;

  size_t
  index_i(const size_t idx) const;

  size_t
  index_d(const size_t idx) const;

  std::string
  state_idx_to_str(const size_t idx) const;

  std::string
  state_idx_to_str_c(const size_t idx) const;

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
  forward_algorithm_c(const bool VERBOSE,
      const bool USE_LOG_ODDS,
      const std::string &observation,
      matrix &forward) const;

  void
  backward_algorithm_c(const bool VERBOSE,
      const bool USE_LOG_ODDS,
      const std::string &observation,
      matrix &backward) const;

  double
  posterior_prob(const matrix &forward) const;

  void
  pseudo_count(void);

  void
  get_viable_transitions_to(void);

  void
  get_viable_transitions_from(void);

  bool
  state_can_emit(const size_t idx) const;

  bool
  state_can_emit_n(const size_t idx) const;

  void
  print_transition(std::ostream& out, const bool HUM_READABLE) const;
  
  void
  print_emission(std::ostream& out, const bool HUM_READABLE) const;

  void
  load_from_file(const std::string inf);

  void
  redistribute_prob(matrix &input,
    const size_t row_idx, const std::vector<size_t> &collapse_cols);

  void
  collapse_states(void);

  void
  train_expectation(const bool VERBOSE, const std::string &seq,
    const matrix &forward, const matrix &backward,
    matrix &e_trans, matrix &e_emiss) const;

  size_t model_len;
  size_t total_size;
  std::string name;
  std::unordered_map<size_t, std::vector<size_t> > transitions_to;
  std::unordered_map<size_t, std::vector<size_t> > transitions_from;
  matrix transition, transition_c;
  matrix emission, emission_c;
  std::vector<double> initial_prob;
};
#endif
