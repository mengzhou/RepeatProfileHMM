/*    scan: a program for finding repeats in the genome using profile-HMM
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
#include <tr1/unordered_map>
#include <sys/types.h>
#include <unistd.h>

#include "smithlab_os.hpp"
#include "OptionParser.hpp"
#include "GenomicRegion.hpp"
#include "ProfileHMM.hpp"

using std::vector;
using std::string;
using std::tr1::unordered_map;
using std::ifstream;
using std::istringstream;
using std::ostream;
using std::ofstream;
using std::cout;
using std::cerr;
using std::endl;

typedef unordered_map<string, string> chrom_file_map;

void
seq_to_int(const string &seq, vector<int> &observation) {
  // no rng version
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
load_hmm_parameter(const string &input_file,
    matrix &transition,
    matrix &emission,
    size_t &model_len) {
  /* Format of input file:
   # some header
   transition_value1 transition_value2 ...
   ...
   //
   emission_value1 ...
   ...
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

void identify_repeats(const matrix &transition,
    const matrix &emission,
    const vector<int> &observation,
    const vector<size_t> &states,
    ProfileHMM &hmm,
    const string chr_name,
    vector<GenomicRegion> &coordinates) {
  const size_t model_len = (transition.size() - 1)/3;
  const size_t bg_state = state(0ul, 0, 1).index(model_len);
  size_t start = 0, end = 0;
  for (vector<size_t>::const_iterator i = states.begin()
      ;i < states.end() - 1; ++i) {
    vector<size_t>::const_iterator j = next(i);
    if (*i == bg_state && *j != bg_state)
      start = j - states.begin();
    else if (*i != bg_state && *j == bg_state) {
      end = j - states.begin();
      const vector<int> obs(observation.begin() + start,
        observation.begin() + end);
      double score =
        hmm.PosteriorProb(transition, emission, obs) - log(end - start);
      string name = score > 50 ? "COPY" : "X";
      GenomicRegion new_copy(chr_name, start,
        end, name, score, '+');
      coordinates.push_back(new_copy);
    }
  }
}

int
main (int argc, const char **argv) {
  try {
    bool VERBOSE = false;
    string chrom_file, in_par, out_file;
    string fasta_suffix = "fa";
    size_t seed = time(0) * getpid();

    OptionParser opt_parse(strip_path(argv[0]), "Program for finding repeats.",
        "-c <chroms> <profile-HMM params file>");

    opt_parse.add_opt("chrom", 'c',
      "File or directory of chroms (FASTA format; .fa suffix)",
      true, chrom_file);
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)",
                      false, out_file);
    opt_parse.add_opt("verbose", 'v', "Verbose mode.", false, VERBOSE);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);

    if (leftover_args.empty()) {
      cout << opt_parse.help_message() << endl
        << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cout << opt_parse.option_missing_message() << endl
        << opt_parse.help_message() << endl;
      return EXIT_FAILURE;
    }
    in_par = leftover_args.front();

    std::ofstream of;
    if (!out_file.empty()) of.open(out_file.c_str());
    std::ostream out(out_file.empty() ? std::cout.rdbuf() : of.rdbuf());

    gsl_rng *rng;
    rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, seed);

    matrix transition, emission;
    size_t model_len;
    if (VERBOSE)
      cerr << "[LOADING HMM]" << endl;
    load_hmm_parameter(in_par, transition, emission, model_len);
    // use log-odds for emission
    log_odds_transform(emission);
    if (VERBOSE)
      cerr << "\tMODEL LENGTH=" << model_len << endl;

    if (VERBOSE)
      cerr << "[LOADING GENOME]" << endl;
    chrom_file_map chrom_files;
    identify_chromosomes(chrom_file, fasta_suffix, chrom_files);
    if (VERBOSE)
      cerr << "\tCHROMS_FOUND=" << chrom_files.size() << endl;


    for (chrom_file_map::const_iterator chrom = chrom_files.begin();
        chrom != chrom_files.end(); ++chrom) {
      string chr_seq;
      read_fasta_file(chrom->second, chrom->first, chr_seq);
      if (chr_seq.empty())
        throw SMITHLABException("could not find chrom: " + chrom->first);

      vector<int> observation;
      vector<size_t> states;
      seq_to_int(rng, chr_seq, observation);
      matrix fm(model_len+1, vector<double>(observation.size()+1, LOG_ZERO));
      matrix fi(model_len, vector<double>(observation.size()+1, LOG_ZERO));
      matrix fd(model_len, vector<double>(observation.size()+1, LOG_ZERO));
      matrix bm(model_len+1, vector<double>(observation.size()+1, LOG_ZERO));
      matrix bi(model_len, vector<double>(observation.size()+1, LOG_ZERO));
      matrix bd(model_len, vector<double>(observation.size()+1, LOG_ZERO));
      ProfileHMM hmm(model_len);
      hmm.PosteriorDecoding(false, transition, emission, observation, states,
          fm, fi, fd, bm, bi, bd);

      vector<GenomicRegion> coordinates;
      identify_repeats(transition, emission, observation, states,
        hmm, chrom->first, coordinates);
      for (vector<GenomicRegion>::const_iterator i = coordinates.begin();
        i < coordinates.end(); ++i)
      {
        out << (*i) << endl;
      }
    }
  }
  catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
