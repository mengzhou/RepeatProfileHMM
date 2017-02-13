/*    multiscan: a program for finding repeats in the genome using profile-HMM
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
#include <unordered_map>
#include <sys/types.h>
#include <unistd.h>

#include "smithlab_os.hpp"
#include "OptionParser.hpp"
#include "GenomicRegion.hpp"
#include "ProfileHMM.hpp"
#include "MultiProfileHMM.hpp"

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

bool
is_bg_state(const vector<size_t> bg_states, const size_t xi_idx) {
  vector<size_t>::const_iterator which =
    std::find(bg_states.begin(), bg_states.end(), xi_idx);
  return !(which == bg_states.end());
}

void
identify_repeats_multifamily(const MultiProfileHMM &multihmm,
    const string &chr_seq,
    const vector<size_t> &states,
    const string chr_name,
    const bool SENSE_STRAND,
    vector<GenomicRegion> &coordinates) {
  // get indeices of background states for individual models
  vector<size_t> bg_states;
  size_t offset = 0;
  for (vector<ProfileHMM*>::const_iterator i = multihmm.begin();
      i < multihmm.end(); ++i) {
    const size_t bg_state_xi_idx = (**i).Length() + offset;
    bg_states.push_back(bg_state_xi_idx);
    offset += (**i).Length() * 2;
  }
  const size_t chr_len = states.size();
  size_t start = 0, end = 0;
  for (vector<size_t>::const_iterator i = states.begin();
      i < states.end() - 1; ++i) {
    vector<size_t>::const_iterator j = next(i);
    const size_t first_model_idx = multihmm.which_model(*i);
    const size_t second_model_idx = multihmm.which_model(*j);
    const bool is_bg_state_i = is_bg_state(bg_states, *i);
    const bool is_bg_state_j = is_bg_state(bg_states, *j);
    if (is_bg_state_i && !is_bg_state_j)
      start = j - states.begin();
    else if ((!is_bg_state_i && is_bg_state_j)
        || (first_model_idx != second_model_idx
          && !is_bg_state_i && !is_bg_state_j)) {
      end = j - states.begin();
      const string obs = chr_seq.substr(start, end - start);
      vector<ProfileHMM*>::const_iterator model =
        multihmm.begin() + first_model_idx;
      double score =
        (**model).PosteriorProb(true, obs) / obs.length();
      string name = (**model).Name();
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
      if (!is_bg_state_i && !is_bg_state_j)
        start = j - states.begin();
    }
  }
}

void
zscore_filter(vector<GenomicRegion> &coordinates,
    const double CUT_OFF) {
  // apply z-score filter to find most likely true copies
  vector<double> scores;
  for (vector<GenomicRegion>::const_iterator i = coordinates.begin();
      i < coordinates.end(); ++i)
    scores.push_back((*i).get_score());
  const double mean = std::accumulate(scores.begin(),
      scores.end(), 0.0) / scores.size();
  const double stdev = std::inner_product(scores.begin(),
      scores.end(), scores.begin(), 0.0) / scores.size()
    - mean * mean;
  for (vector<GenomicRegion>::iterator i = coordinates.begin();
      i < coordinates.end(); ++i) {
    (*i).set_score(std::abs((*i).get_score() - mean) / stdev);
    //if ((*i).get_score() > CUT_OFF)
    //  (*i).set_name("COPY");
    //else
    //  (*i).set_name("X");
  }
}

int
main (int argc, const char **argv) {
  try {
    bool VERBOSE = false;
    bool DEBUG = false;
    string chrom_file, in_par, out_file;
    string fasta_suffix = "fa";

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
    if (leftover_args.size() < 2) {
      cout << "Must provide more than 1 model parameter files!" << endl;
      return EXIT_FAILURE;
    }

    std::ofstream of;
    if (!out_file.empty()) of.open(out_file.c_str());
    std::ostream out(out_file.empty() ? std::cout.rdbuf() : of.rdbuf());

    if (VERBOSE)
      cerr << "[LOADING HMM]" << endl;
    vector<ProfileHMM> v;
    for (vector<string>::const_iterator file = leftover_args.begin();
        file < leftover_args.end(); ++file) {
      v.push_back(ProfileHMM(*file));
    }
    MultiProfileHMM multihmm(v);
    if (VERBOSE)
      cerr << "\tMODEL NUMBER=" << v.size() << endl;

    if (VERBOSE)
      cerr << "[LOADING GENOME]" << endl;
    chrom_file_map chrom_files;
    identify_chromosomes(chrom_file, fasta_suffix, chrom_files);
    if (VERBOSE)
      cerr << "\tCHROMS_FOUND=" << chrom_files.size() << endl;

    size_t file_counter = 1;
    for (chrom_file_map::const_iterator chrom = chrom_files.begin();
        chrom != chrom_files.end(); ++chrom) {
      vector<string> chr_name, chr_seq;
      read_fasta_file(chrom->second, chr_name, chr_seq);
      if (chr_name.size() < 1)
        throw SMITHLABException("could not find any sequence in: "
            + chrom->second);

      vector<size_t> states;

      if (VERBOSE)
        cerr << "[SCANNING " << chrom->second
          << " " << file_counter++ << " OF "
          << chrom_files.size() << "]" << endl;
      vector<GenomicRegion> coordinates;
      for (size_t i = 0; i < chr_seq.size(); ++i) {
        // the forward strand
        if (VERBOSE)
          cerr << "\t" << i+1 << "/" << chr_seq.size()
            << "\t" << chr_name[i] << endl;
        if (DEBUG)
          cerr << chr_name[i] << endl;
        multihmm.PosteriorDecoding(VERBOSE, DEBUG, true, chr_seq[i], states);
        identify_repeats_multifamily(multihmm, chr_seq[i], states,
          chr_name[i], true, coordinates);
        // the reverse strand
        revcomp_inplace(chr_seq[i]);
        multihmm.ComplementBackground();
        multihmm.PosteriorDecoding(VERBOSE, DEBUG, true, chr_seq[i], states);
        identify_repeats_multifamily(multihmm, chr_seq[i], states,
          chr_name[i], false, coordinates);
      }
      zscore_filter(coordinates, 1.0);
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
