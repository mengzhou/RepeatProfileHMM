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
#include <fstream>
#include <sstream>
#include <random>
#include <gsl/gsl_randist.h>
#include <sys/types.h>
#include <unistd.h>

#include "ProfileHMM.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "smithlab_os.hpp"

using namespace std;

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

template <typename T>
void
print_matrix(const vector<vector<T> > &matrix, const size_t offset = 0,
    const bool print_row_sum = true) {
  typename vector<vector<T> >::const_iterator first = matrix.begin();
  for (typename vector<T>::const_iterator j = (*first).begin();
    j < (*first).end(); ++j)
    cout << "\t" << j - (*first).begin() + offset;
  cout << endl;
  for (typename vector<vector<T> >::const_iterator i = matrix.begin();
    i < matrix.end(); ++i) {
    cout << i - matrix.begin();
    double sum = 0.0;
    for(typename vector<T>::const_iterator j = (*i).begin();
        j < (*i).end(); ++j) {
      printf("\t%.3f", exp(*j));
      sum += exp(*j);
    }
    if (print_row_sum)
      cout << "\t" << sum << endl;
    else
      cout << endl;
  }
}

void
seq_to_int(const string &seq, vector<int> &observation) {
  for (string::const_iterator i = seq.begin(); i < seq.end(); ++i) {
    if (*i == 'A' || *i == 'a')
      observation.push_back(0);
    else if (*i == 'C' || *i == 'c')
      observation.push_back(1);
    else if (*i == 'G' || *i == 'g')
      observation.push_back(2);
    else if (*i == 'T' || *i == 't')
      observation.push_back(3);
    else {
      // no rng version
      observation.push_back('A');
    }
  }
}

void
seq_to_int(const gsl_rng* rng, const string &seq, vector<int> &observation) {
  for (string::const_iterator i = seq.begin(); i < seq.end(); ++i) {
    if (*i == 'A' || *i == 'a')
      observation.push_back(0);
    else if (*i == 'C' || *i == 'c')
      observation.push_back(1);
    else if (*i == 'G' || *i == 'g')
      observation.push_back(2);
    else if (*i == 'T' || *i == 't')
      observation.push_back(3);
    else {
      // usually this case is N; then randomly assign one base
      observation.push_back(gsl_rng_uniform_int(rng, 4));
    }
  }
}

void
int_to_seq(const vector<int> &observation, string &seq) {
  for (vector<int>::const_iterator i = observation.begin();
      i < observation.end(); ++i) {
    if (*i == 0)
      seq.append("A");
    else if (*i == 1)
      seq.append("C");
    else if (*i == 2)
      seq.append("G");
    else if (*i == 3)
      seq.append("T");
  }
}

void
print_trace(const vector<pair<char, size_t> > &trace) {
  size_t counter = 0;
  for (vector<pair<char, size_t> >::const_iterator i = trace.begin();
      i < trace.end(); ++i) {
    cout << counter << "\t" << (*i).first << (*i).second << endl;
    if (((*i).first == 'M' && (*i).second > 0) or (*i).first == 'I')
      ++counter;
  }
  cout << endl;
}

void
state_to_trace(const vector<size_t> &states,
    const size_t model_len,
    vector<pair<char, size_t> > &trace) {
  char state;
  size_t idx;
  trace.clear();
  for (vector<size_t>::const_iterator i = states.begin();
      i < states.end(); ++i) {
    if (*i <= model_len) {
      state = 'M';
      idx = *i;
    }
    else if (*i <= model_len * 2) {
      state = 'I';
      idx = *i - model_len - 1;
    }
    else if (*i <= model_len * 3) {
      state = 'D';
      idx = *i - model_len * 2;
    }
    else {
      state = 'E';
      idx = 0;
    }
  pair<char, size_t> step(state, idx);
  trace.push_back(step);
  }
}
    
void
load_hmm_parameter(const bool VERBOSE,
    const string &input_file,
    vector<vector<double> > &transition,
    vector<vector<double> > &emission,
    size_t &model_len) {
  /* Format of input file:
   transition_value1 transition_value2 ...
   //
   emission_value1 ...
   //
  */
  ifstream in(input_file.c_str());
  // skip header with #
  string line;
  getline(in, line);
  while((line.substr(0,1)) == "#"){
    getline(in, line);
  }

  assert(transition.empty());
  double val;
  while (!line.empty() && line.compare("//")) {
    transition.push_back(vector<double>());
    istringstream iss(line);
    while (iss >> val) {
      (*(transition.end()-1)).push_back(val);
    }
    getline(in, line);
  }
  model_len = (transition.size() - 2) / 3;
  assert(emission.empty());
  getline(in, line);
  while (!line.empty() && line.compare("//")) {
    emission.push_back(vector<double>());
    istringstream iss(line);
    while (iss >> val) {
      (*(emission.end()-1)).push_back(val);
    }
    getline(in, line);
  }

  const size_t total_size = model_len * 3 + 2;
  assert(transition.size() == total_size);
  assert((*transition.begin()).size() == total_size);
  assert(emission.size() == total_size);
}

void
write_hmm_parameter(const bool VERBOSE,
    const string &output_file,
    const vector<vector<double> > &transition,
    const vector<vector<double> > &emission) {
  ofstream out(output_file.c_str());
  for (vector<vector<double> >::const_iterator i = transition.begin();
      i < transition.end(); ++i) {
    for (vector<double>::const_iterator j = (*i).begin();
        j < (*i).end(); ++j) {
      out << *j << "\t";
    }
    out << endl;
  }
  out << "//" << endl;
  for (vector<vector<double> >::const_iterator i = emission.begin();
      i < emission.end(); ++i) {
    for (vector<double>::const_iterator j = (*i).begin();
        j < (*i).end(); ++j) {
      out << *j << "\t";
    }
    out << endl;
  }
  out << "//" << endl;
}

void
make_hmm_parameter(const bool VERBOSE, const gsl_rng* rng,
    vector<vector<double> > &transition,
    vector<vector<double> > &emission,
    const size_t model_len) {
  // total size 
  // M_0(B) ~ M_L + I_0 ~ I_L-1 + D_1 ~ D_L + E
  const size_t total_size = model_len * 3 + 2;
  const int alphabet_size = 4;

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
    print_matrix(transition);
  }

  // setting emission
  emission.resize(total_size);
  for (i = emission.begin(); i < emission.end(); ++i) {
    (*i).resize(alphabet_size);
    for (vector<double>::iterator j = (*i).begin(); j < (*i).end(); ++j) {
      *j = LOG_ZERO;
    }
  }
  // bg A, C, G, T
  double insert[4] = {0.40, 0.10, 0.10, 0.40};
  double match[4] = {0.05, 0.40, 0.50, 0.05};
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
    print_matrix(emission);
  }
}

void identify_repeats(const vector<vector<double> > &transition,
    const vector<vector<double> > &emission,
    const string &observation,
    const vector<size_t> &states,
    ProfileHMM &hmm,
    const string chr_name,
    vector<GenomicRegion> &coordinates) {
  const size_t model_len = (transition.size() - 1)/3;
  const size_t bg_state = index_i(model_len, 0);
  vector<size_t>::const_iterator i = states.begin();
  size_t start = 0, end = 0;
  for (;i < states.end() - 1; ++i) {
    vector<size_t>::const_iterator j = next(i);
    if (*i == bg_state && *j != bg_state)
      start = j - states.begin();
    else if (*i != bg_state && *j == bg_state) {
      end = j - states.begin();
      const string obs = observation.substr(start, end - start - 1);
      double score =
        hmm.PosteriorProb(true, obs) - log(end - start);
      string name = score > 50 ? "COPY" : "X";
      GenomicRegion new_copy(chr_name, start,
        end, name, score, '+');
      coordinates.push_back(new_copy);
    }
  }
}

int
main (int argc, const char **argv) {
  bool VERBOSE = false;
  string input, in_par, out_par;
  size_t model_len = 10;
  size_t seed = time(0) * getpid();
  vector<vector<double> > transition, emission;

  OptionParser opt_parse(argv[0], "Program for profile-HMM.");

  opt_parse.add_opt("in", 'i', "Input sequence (string) or a fasta file.",
    false, input);
  opt_parse.add_opt("len", 'l', "Model length.", false, model_len);
  opt_parse.add_opt("seed", 's', "Random seed.", false, seed);
  opt_parse.add_opt("in-params", 'p', "Input parameters file.", false, in_par);
  opt_parse.add_opt("out-params", 'q', "Output parameters file.", false, out_par);
  opt_parse.add_opt("verbose", 'v', "Verbose mode.", false, VERBOSE);

  vector<string> leftover_args;
  opt_parse.parse(argc, argv, leftover_args);

  if (in_par.compare(out_par) == 0) {
    cout << opt_parse.help_message() << endl
      << opt_parse.about_message() << endl;
    return EXIT_SUCCESS;
  }

  gsl_rng *rng;
  rng = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(rng, seed);

  if (in_par.empty())
    make_hmm_parameter(VERBOSE, rng, transition, emission, model_len);
  else
    load_hmm_parameter(VERBOSE, in_par, transition, emission, model_len);

  if (!out_par.empty())
    write_hmm_parameter(VERBOSE, out_par, transition, emission);

  ProfileHMM hmm(transition, emission);
  vector<string> chrom_names, ref_chroms;
  // decoding test
  if (!input.empty()) {
    if (VERBOSE) {
      cout << "Input: " << input << endl;
      cout << "Model length: " << model_len << endl;
    }
    const string suffix(input.substr(input.find_last_of(".")+1));
    string input_seq;
    if (suffix.compare("fa") == 0) {
      read_fasta_file(input, chrom_names, ref_chroms);
      input_seq = ref_chroms[0];
    }
    else {
      input_seq = input;
    }
    // Decoding test
    vector<pair<char, size_t> > trace;
    vector<size_t> states;
    //log_odds_transform(emission);
    //const double lh = 
    //hmm.ViterbiDecoding(VERBOSE, transition, emission, observation, trace);
    //print_trace(trace);
    hmm.PosteriorDecoding(false, true, input_seq, states);
    //state_to_trace(states, model_len, trace);
    //print_trace(trace);
    //vector<GenomicRegion> coordinates;
    //identify_repeats(transition, emission, observation, states,
    //  hmm, chrom_names[0], coordinates);
    //for (vector<GenomicRegion>::const_iterator i = coordinates.begin();
    //  i < coordinates.end(); ++i)
    //{
    //  cout << (*i) << endl;
    //}

    // Learning test
    //if (!leftover_args.empty()) {
    //  ProfileHMM h(leftover_args.front());
    //  h.Print(cout, true);
    //}
    hmm.Print(cout, true);
    hmm.Train(VERBOSE, 1e-4, 20, input_seq);
    hmm.Print(cout, true);
    //hmm.PosteriorDecoding(false, true, observation, states);
  }
  else {
    for (size_t i = 1; i <= 3; ++i) {
      // sampling test
      vector<pair<char, size_t> > trace;
      vector<size_t> states;
      vector<int> seq;
      string output;
      hmm.SampleSequence(VERBOSE, rng, seq, states);
      cout << "#Trial " << i << endl;
      int_to_seq(seq, output);
      cout << "  Sampled seq:\t" << output << endl;
      state_to_trace(states, model_len, trace);
      cout << "  Trace:\t";
      print_trace(trace);
      trace.clear();
      //hmm.ViterbiDecoding(VERBOSE, seq, trace);
      cout << "  Decoded:\t";
      print_trace(trace);
      cout << endl;
    }
  }

  // learning test
  //hmm.forward_algorithm(VERBOSE, transition, emission, observation);
  //cout << "P=" << exp(hmm.forward_prob('E', 2, 2)) << endl;
  //hmm.backward_algorithm(VERBOSE, transition, emission, observation);
  //hmm.BW_training(VERBOSE, transition, emission, observation);
  //print_matrix(transition,0,true);
  //print_matrix(emission,0,true);
}
