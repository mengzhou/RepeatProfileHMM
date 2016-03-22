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

#ifndef PHMM_HPP
#define PHMM_HPP

#include <string>
#include <initializer_list>
#include <gsl/gsl_randist.h>
#include <map>

#include "smithlab_utils.hpp"

typedef std::vector<std::vector<double> > matrix;

extern const double LOG_ZERO;

size_t
baseint2stateint(const size_t &baseint, const bool marked);

void
print_transition(const matrix &transition);

void
print_emission(const matrix &emission);

void
log_odds_transform(matrix &emission);

void
normalize_vec_inplace(std::vector<double> &v, const bool logged);

double
log_sum_log_list(const std::initializer_list<double> &list);

double
log_sum_log(const double p, const double q);

size_t
argmax_list(const std::initializer_list<double> &list);

size_t
argmax_vec(const std::vector<double> &v);

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
public:
  ProfileHMM();
  ProfileHMM(const size_t ml);

  double
  ViterbiDecoding(const bool VERBOSE,
      const matrix &transition,
      const matrix &emission,
      const std::vector<int> &observation,
      std::vector<std::pair<char, size_t> > &trace) const;

  void
  Train(const bool VERBOSE,
      matrix &transition,
      matrix &emission,
      const std::vector<int> &observation);

  void
  SampleSequence(const bool VERBOSE,
      const gsl_rng* rng,
      const matrix &transition,
      const matrix &emission,
      std::vector<int> &seq,
      std::vector<size_t> &states) const;

  void
  PosteriorDecoding(const bool VERBOSE,
      const matrix &transition,
      const matrix &emission,
      const std::vector<int> &observation,
      std::vector<size_t> &states) const;

  double
  PosteriorProb(const matrix &transition,
      const matrix &emission,
      const std::vector<int> &observation) const;

private:
  size_t
  index_m(const size_t idx) const;

  size_t
  index_i(const size_t idx) const;

  size_t
  index_d(const size_t idx) const;

  std::string
  state_idx_to_str(const size_t idx) const;

  void
  forward_algorithm(const bool VERBOSE,
      const matrix &transition,
      const matrix &emission,
      const std::vector<int> &observation,
      matrix &forward) const;

  void
  backward_algorithm(const bool VERBOSE,
      const matrix &transition,
      const matrix &emission,
      const std::vector<int> &observation,
      matrix &backward) const;

  double
  posterior_prob(const matrix &forward) const;

  void
  add_pseudocount_uniform(matrix &transition,
      matrix &emission) const;

  std::map<size_t, std::vector<size_t> >
  get_viable_transitions_to(void) const;

  std::map<size_t, std::vector<size_t> >
  get_viable_transitions_from(void) const;

  size_t model_len;
  size_t total_size;
  std::map<size_t, std::vector<size_t> > transitions_to;
  std::map<size_t, std::vector<size_t> > transitions_from;
  const double tolerance = 1e-9;
  const size_t max_iterations = 50;
};
#endif
