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
#include "smithlab_utils.hpp"

class ProfileHMM {
  // profile-HMM allowing local alignment and re-occurrence of M_i;
  // see Durbin book p.114 (bottom figure)
public:
  ProfileHMM (const bool v, const size_t ml, const size_t sl) :
    VERBOSE(v), model_len(ml),
    seq_len(sl) {}

  double
  ViterbiDecoding(const std::vector<std::vector<double> > &transition,
      const std::vector<std::vector<double> > &emission,
      const std::vector<double> &initial,
      const std::vector<int> &observation,
      std::vector<std::pair<char, size_t> > &trace) const;

  /*
  double
  BaumWelchTraining() const;

  double
  PosteriorDecoding() const;
  */

private:
  double
  log_sum_log_vec(const std::vector<double> &vals, size_t limit) const;

  double
  log_sum_log(const double p, const double q) const;

  double
  max_item(const std::initializer_list<double> list) const;

  size_t
  index_m(const size_t idx) const;

  size_t
  index_i(const size_t idx) const;

  size_t
  index_d(const size_t idx) const;

  bool VERBOSE;
  size_t model_len, seq_len;
};

#endif
