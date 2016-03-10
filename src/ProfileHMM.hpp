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

#include <initializer_list>
#include <gsl/gsl_randist.h>
#include "smithlab_utils.hpp"

class ProfileHMM {
  // profile-HMM allowing local alignment and re-occurrence of M_i;
  // see Durbin book p.114 (bottom figure)
public:
  ProfileHMM(const size_t ml) : model_len(ml) {}

  static double
  log_sum_log_list(const std::initializer_list<double> &list);

  static double
  log_sum_log(const double p, const double q);

  static size_t
  argmax_list(const std::initializer_list<double> &list);

  static size_t
  argmax_vec(const std::vector<double> &v);

  double
  ViterbiDecoding(const bool VERBOSE,
      const std::vector<std::vector<double> > &transition,
      const std::vector<std::vector<double> > &emission,
      const std::vector<int> &observation,
      std::vector<std::pair<char, size_t> > &trace) const;

  void
  BW_training(const bool VERBOSE,
      std::vector<std::vector<double> > &transition,
      std::vector<std::vector<double> > &emission,
      const std::vector<int> &observation);

  void
  SampleSequence(const bool VERBOSE,
      const gsl_rng* rng,
      const std::vector<std::vector<double> > &transition,
      const std::vector<std::vector<double> > &emission,
      std::vector<int> &seq,
      std::vector<size_t> &states) const;

  void
  PosteriorDecoding(const bool VERBOSE,
      const std::vector<std::vector<double> > &transition,
      const std::vector<std::vector<double> > &emission,
      const std::vector<int> &observation,
      std::vector<size_t> &states);

  double
  PosteriorProb(const std::vector<std::vector<double> > &transition,
      const std::vector<std::vector<double> > &emission,
      const std::vector<int> &observation);

private:
  size_t
  index_m(const size_t idx) const;

  size_t
  index_i(const size_t idx) const;

  size_t
  index_d(const size_t idx) const;

  void
  forward_algorithm(const bool VERBOSE,
      const std::vector<std::vector<double> > &transition,
      const std::vector<std::vector<double> > &emission,
      const std::vector<int> &observation,
      std::vector<std::vector<double> > &fm,
      std::vector<std::vector<double> > &fi,
      std::vector<std::vector<double> > &fd);

  void
  backward_algorithm(const bool VERBOSE,
      const std::vector<std::vector<double> > &transition,
      const std::vector<std::vector<double> > &emission,
      const std::vector<int> &observation,
      std::vector<std::vector<double> > &bm,
      std::vector<std::vector<double> > &bi,
      std::vector<std::vector<double> > &bd);

  size_t model_len;
  const double LOG_ZERO = -1e7;
  const size_t total_size = model_len * 3 + 2;
  const double tolerance = 1e-10;
  const size_t max_iterations = 100;
  std::vector<std::vector<double> > all_fm, all_fi, all_fd;
  std::vector<std::vector<double> > all_bm, all_bi, all_bd;
};

#endif
