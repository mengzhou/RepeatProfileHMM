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

#include <iostream>
#include <string>
#include <cmath>
#include <stdio.h>
#include <random>
#include <gsl/gsl_randist.h>
#include <sys/types.h>
#include <unistd.h>

#include "ProfileHMM.hpp"
#include "OptionParser.hpp"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::pair;

size_t
index_i(const size_t model_len, const size_t idx) {
  // I_0 ~ I_L-1 + (L+1)*M
  return idx + model_len + 1;
}

size_t
index_d(const size_t model_len, const size_t idx) {
  // D_1 ~ D_L + (L+1)*M + L*I
  return idx + model_len * 2;
}

void seq_to_int(const string &seq, vector<int> &observation) {
  for (string::const_iterator i = seq.begin(); i < seq.end(); ++i) {
    switch (*i) {
      case 'A':
        observation.push_back(0);
      case 'C':
        observation.push_back(1);
      case 'G':
        observation.push_back(2);
      case 'T':
        observation.push_back(3);
    }
  }
}

void make_hmm_parameter(const bool VERBOSE, vector<vector<double> > &transition,
    vector<vector<double> > &emission, vector<double> &initial,
    const size_t model_len) {
  // total size 
  // M_0(B) ~ M_L + I_0 ~ I_L-1 + D_1 ~ D_L + E
  const size_t total_size = model_len * 3 + 2;
  const int alphabet_size = 4;
  const double LOG_ZERO = -1000.0;

  // setting transition
  transition.resize(total_size);
  vector<vector<double> >::iterator i;
  for (i = transition.begin(); i < transition.end(); ++i) {
    (*i).resize(total_size);
    for (vector<double>::iterator j = (*i).begin(); j < (*i).end(); ++j) {
      *j = LOG_ZERO;
    }
  }
  for (size_t x = 1; x < model_len; ++x) {
    if (x < model_len - 1) {
      // M_1 ~ M_L-2
      // next D
      transition[x][index_d(model_len, x+1)] = log(0.1);
      // next M
      transition[x][x+1] = log(0.6);
      // I
      transition[x][index_i(model_len, x)] = log(0.1);
      // D_L
      transition[x][index_d(model_len, model_len)] = log(0.2);
      // I_1 ~ I_L-2
      // self
      transition[index_i(model_len, x)][index_i(model_len, x)] = log(0.3);
      // next M
      transition[index_i(model_len, x)][x+1] = log(0.6);
      // next D
      transition[index_i(model_len, x)][index_d(model_len, x+1)] = log(0.1);
    }
    else {
      // M_L-1
      // next M
      transition[x][x+1] = log(0.7);
      // I
      transition[x][index_i(model_len, x)] = log(0.1);
      // D_L
      transition[x][index_d(model_len, model_len)] = log(0.2);
      // I_L-1
      // self
      transition[index_i(model_len, x)][index_i(model_len, x)] = log(0.35);
      // next M
      transition[index_i(model_len, x)][x+1] = log(0.65);
    }

    if (x > 1 && x < model_len - 1) {
      // D_2 ~ D_L-2
      // next D
      transition[index_d(model_len, x)][index_d(model_len, x+1)] = log(0.4);
      // I
      transition[index_d(model_len, x)][index_i(model_len, x)] = log(0.1);
      // next M
      transition[index_d(model_len, x)][x+1] = log(0.5);
    }
    else if (x > 1) {
      // D_L-1
      // I
      transition[index_d(model_len, x)][index_i(model_len, x)] = log(0.2);
      // next M
      transition[index_d(model_len, x)][x+1] = log(0.8);
    }
  }

  // B
  transition[0][index_d(model_len, 1)] = log(0.4);
  transition[0][index_i(model_len, 0)] = log(0.6);
  // D_1
  double first_base_weight = 0.5;
  transition[index_d(model_len, 1)][1] = log(first_base_weight);
  for (size_t x = 2; x <= model_len; ++x) {
    transition[index_d(model_len, 1)][x] =
      log(1.0 - first_base_weight) - log(model_len - 1);
  }

  // I_0
  // D_1
  transition[index_i(model_len, 0)][index_d(model_len, 1)] = log(0.4);
  // self
  transition[index_i(model_len, 0)][index_i(model_len, 0)] = log(0.4);
  // end
  transition[index_i(model_len, 0)][total_size - 1] = log(0.2);

  // D_L
  // I_0
  transition[index_d(model_len, model_len)][index_i(model_len, 0)] = log(0.6);
  // end
  transition[index_d(model_len, model_len)][total_size - 1] = log(0.4);

  // M_L
  transition[model_len][index_d(model_len, model_len)] = 0.0;

  if (VERBOSE) {
    cout << "Transition" << endl;
    for (size_t x = 0; x < total_size; ++x)
      cout << "\t" << x;
    cout << endl;
    for (size_t x = 0; x < total_size; ++x) {
      cout << x;
      double sum = 0.0;
      for (size_t y = 0; y < total_size; ++y) {
        printf("\t%.3f", exp(transition[x][y]));
        sum += exp(transition[x][y]);
      }
      cout << "\t" << sum << endl;
    }
  }

  // setting emission
  gsl_rng *rng;
  rng = gsl_rng_alloc(gsl_rng_default);
  size_t seed = time(0) * getpid();
  gsl_rng_set(rng, seed);
  emission.resize(total_size);
  for (i = emission.begin(); i < emission.end(); ++i) {
    (*i).resize(alphabet_size);
    for (vector<double>::iterator j = (*i).begin(); j < (*i).end(); ++j) {
      *j = LOG_ZERO;
    }
  }
  // bg A, C, G, T
  double insert[4] = {0.49, 0.01, 0.01, 0.49};
  double match[4] = {0.10, 0.70, 0.10, 0.10};
  emission[index_i(model_len, 0)][0] = log(insert[0]);
  emission[index_i(model_len, 0)][1] = log(insert[1]);
  emission[index_i(model_len, 0)][2] = log(insert[2]);
  emission[index_i(model_len, 0)][3] = log(insert[3]);
  
  const bool RANDOMIZE = true;
  for (size_t x = 1; x <= model_len; ++x) {
    double prob[4];
    if (RANDOMIZE) {
      gsl_ran_dirichlet(rng, 4, match, prob);
      emission[x][0] = log(prob[0]);
      emission[x][1] = log(prob[1]);
      emission[x][2] = log(prob[2]);
      emission[x][3] = log(prob[3]);
    } else {
      emission[x][0] = log(match[0]);
      emission[x][1] = log(match[1]);
      emission[x][2] = log(match[2]);
      emission[x][3] = log(match[3]);
    }

    if (x < model_len) {
      if (RANDOMIZE) {
        gsl_ran_dirichlet(rng, 4, insert, prob);
        emission[index_i(model_len, x)][0] = log(prob[0]);
        emission[index_i(model_len, x)][1] = log(prob[1]);
        emission[index_i(model_len, x)][2] = log(prob[2]);
        emission[index_i(model_len, x)][3] = log(prob[3]);
      } else {
        emission[index_i(model_len, x)][0] = log(insert[0]);
        emission[index_i(model_len, x)][1] = log(insert[1]);
        emission[index_i(model_len, x)][2] = log(insert[2]);
        emission[index_i(model_len, x)][3] = log(insert[3]);
      }
    }
  }

  if (VERBOSE) {
    cout << "Emission" << endl;
    cout << "\tA\tC\tG\tT" << endl;
    for (size_t x = 0; x < total_size; ++x) {
      cout << x;
      double sum = 0.0;
      for (size_t y = 0; y < 4; ++y) {
        printf("\t%.3f", exp(emission[x][y]));
        sum += exp(emission[x][y]);
      }
      cout << "\t" << sum << endl;
    }
  }

  // setting initial prob
  initial.resize(total_size);
  for (vector<double>::iterator j = initial.begin(); j < initial.end(); ++j) {
    *j = LOG_ZERO;
  }
  initial[0] = 0.0;
}

int main (int argc, const char **argv) {
  bool VERBOSE = false;
  string input;
  size_t model_len = 10;
  vector<vector<double> > transition, emission;
  vector<double> initial;

  OptionParser opt_parse(argv[0], "Program for profile-HMM.");

  opt_parse.add_opt("in", 'i', "Input sequence.", true, input);
  opt_parse.add_opt("len", 'l', "Model length.", false, model_len);
  opt_parse.add_opt("verbose", 'v', "Verbose mode.", false, VERBOSE);

  vector<string> leftover_args;
  opt_parse.parse(argc, argv, leftover_args);

  if (input.empty()) {
    cerr << opt_parse.help_message() << endl
      << opt_parse.about_message() << endl;
    return EXIT_SUCCESS;
  }

  cout << "Input: " << input << endl;
  cout << "Model length: " << model_len << endl;
  make_hmm_parameter(VERBOSE, transition, emission, initial, model_len);
  vector<int> observation;
  vector<pair<char, size_t> > trace;
  seq_to_int(input, observation);
  ProfileHMM hmm(VERBOSE, model_len, input.size());
  const double lh = 
  hmm.ViterbiDecoding(transition, emission, initial, observation, trace);

  cout << endl << "Result: " << lh << endl;
  for (vector<pair<char, size_t> >::const_iterator i = trace.begin();
      i < trace.end(); ++i) {
    cout << (*i).first << (*i).second << " ";
  }
  cout << endl;
}
