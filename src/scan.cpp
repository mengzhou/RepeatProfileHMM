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
//#include <gsl/gsl_randist.h>
#include <unordered_map>
#include <sys/types.h>
#include <unistd.h>

#include "smithlab_os.hpp"
#include "OptionParser.hpp"
#include "GenomicRegion.hpp"
#include "ProfileHMM.hpp"

using std::vector;
using std::string;
using std::unordered_map;
using std::ifstream;
using std::istringstream;
using std::ostream;
using std::ofstream;
using std::cout;
using std::cerr;
using std::endl;

typedef unordered_map<string, string> chrom_file_map;

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

string
get_state_bits(const vector<size_t> &copy_states,
    const size_t width) {
  // states are indexed as 0=M_1~M_L; L=I_0~I_L-1
  string bits(width, '0');
  for (vector<size_t>::const_iterator i = copy_states.begin();
      i < copy_states.end(); ++i)
    if (*i >= 0 && *i < width)
      bits[*i] = '1';
  return bits;
}

void
identify_repeats(const ProfileHMM &hmm,
    const string &chr_seq,
    const vector<size_t> &states,
    const string chr_name,
    const bool SENSE_STRAND,
    vector<GenomicRegion> &coordinates,
    vector<string> &state_bits) {
  const size_t model_len = hmm.Length();
  //const size_t bg_state = state(0ul, 0, 1).index(model_len);
  const size_t bg_state = model_len;
  const size_t chr_len = states.size();
  size_t start = 0, end = 0;
  for (vector<size_t>::const_iterator i = states.begin();
      i < states.end(); ++i) {
    vector<size_t>::const_iterator j = next(i);
    // find start of region
    if (*i != bg_state && i == states.begin())
      start = 0;
    else if (*i == bg_state && *j != bg_state)
      start = j - states.begin();
    // find end of region
    else if ((*i != bg_state && *j == bg_state)
        || (*i != bg_state && j == states.end())) {
      end = j - states.begin();
      const string obs = chr_seq.substr(start, end - start);
      // this is the score using log-odds
      double score =
        hmm.PosteriorProb_c(true, obs) / obs.length();
      // this is the score using reverse sequenc as normalization
      //string obs_rev(obs.rbegin(), obs.rend());
      //double score =
      //  hmm.PosteriorProb(true, obs) - hmm.PosteriorProb(true, obs_rev);
      string name = "X";
      if (SENSE_STRAND) {
        GenomicRegion new_copy(chr_name, start,
          end, name, score, '+');
        coordinates.push_back(new_copy);
      }
      else {
        GenomicRegion new_copy(chr_name, chr_len - end + 1,
          chr_len - start + 1, name, score, '-');
        coordinates.push_back(new_copy);
      }

      // get bit vector for matching states
      vector<size_t>::const_iterator first = states.begin() + start;
      vector<size_t>::const_iterator last = states.begin() + end;
      vector<size_t> copy_states(first, last);
      string bits = get_state_bits(copy_states, model_len);
      state_bits.push_back(bits);
    }
  }
}

void
zscore_filter(vector<GenomicRegion> &coordinates,
    const double CUT_OFF) {
  // apply z-score filter to find most likely true copies
  vector<double> scores;
  for (vector<GenomicRegion>::const_iterator i = coordinates.begin();
      i < coordinates.end(); ++i) {
    if ((*i).get_name() == "X")
      scores.push_back((*i).get_score());
  }
  const double mean = std::accumulate(scores.begin(),
      scores.end(), 0.0) / scores.size();
  const double stdev = std::inner_product(scores.begin(),
      scores.end(), scores.begin(), 0.0) / scores.size()
    - mean * mean;
  for (vector<GenomicRegion>::iterator i = coordinates.begin();
      i < coordinates.end(); ++i) {
    (*i).set_score(std::abs((*i).get_score() - mean) / stdev);
    if ((*i).get_score() > CUT_OFF)
      (*i).set_name("COPY");
    else
      (*i).set_name("X");
  }
}

int
main (int argc, const char **argv) {
  try {
    bool VERBOSE = false;
    bool DEBUG = false;
    string chrom_file, in_par, out_file;
    string fasta_suffix = "fa";
    //size_t seed = time(0) * getpid();

    OptionParser opt_parse(strip_path(argv[0]), "Program for finding repeats.",
        "-c <chroms> <profile-HMM params file>");

    opt_parse.add_opt("chrom", 'c',
      "File or directory of chroms (FASTA format; .fa suffix)",
      true, chrom_file);
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)",
                      false, out_file);
    opt_parse.add_opt("verbose", 'v', "Verbose mode.", false, VERBOSE);
    opt_parse.add_opt("debug", 'd', "Print debug information.", false, DEBUG);

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

    // currently not used. in the furture needed for handling N in the genome
    //gsl_rng *rng;
    //rng = gsl_rng_alloc(gsl_rng_default);
    //gsl_rng_set(rng, seed);

    if (VERBOSE)
      cerr << "[LOADING HMM]" << endl;
    ProfileHMM hmm(in_par);
    if (VERBOSE)
      cerr << "\tMODEL LENGTH=" << hmm.Length() << endl;

    if (VERBOSE)
      cerr << "[LOADING GENOME]" << endl;
    chrom_file_map chrom_files;
    identify_chromosomes(chrom_file, fasta_suffix, chrom_files);
    if (VERBOSE)
      cerr << "\tCHROMS_FOUND=" << chrom_files.size() << endl;

    for (chrom_file_map::const_iterator chrom = chrom_files.begin();
        chrom != chrom_files.end(); ++chrom) {
      vector<string> chr_name, chr_seq;
      read_fasta_file(chrom->second, chr_name, chr_seq);
      if (chr_name.size() < 1)
        throw SMITHLABException("could not find any sequence in: "
            + chrom->second);

      vector<size_t> states;
      vector<string> state_bits;

      if (VERBOSE)
        cerr << "[SCANNING " << chrom->second << "]" << endl;
      vector<GenomicRegion> coordinates;
      for (size_t i = 0; i < chr_seq.size(); ++i) {
        // the forward strand
        hmm.PosteriorDecoding_c(false, DEBUG, false, chr_seq[i], states);
        identify_repeats(hmm, chr_seq[i], states,
          chr_name[i], true, coordinates, state_bits);
        // the reverse strand
        revcomp_inplace(chr_seq[i]);
        hmm.ComplementBackground();
        hmm.PosteriorDecoding_c(false, DEBUG, false, chr_seq[i], states);
        identify_repeats(hmm, chr_seq[i], states,
          chr_name[i], false, coordinates, state_bits);
      }
      zscore_filter(coordinates, 3.0);
      for (vector<GenomicRegion>::const_iterator i = coordinates.begin();
        i < coordinates.end(); ++i)
      {
        out << (*i) << "\t" << state_bits[i-coordinates.begin()] << endl;
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
