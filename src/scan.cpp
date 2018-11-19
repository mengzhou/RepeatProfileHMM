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
#include <unordered_map>
#include <sys/types.h>
#include <unistd.h>
#ifdef OPENMP
#include <omp.h>
#endif

#include "smithlab_os.hpp"
#include "OptionParser.hpp"
#include "GenomicRegion.hpp"
#include "ProfileHMM.hpp"

using std::vector;
using std::string;
using std::to_string;
using std::unordered_map;
using std::ifstream;
using std::istringstream;
using std::ostream;
using std::ofstream;
using std::cout;
using std::cerr;
using std::endl;

typedef unordered_map<string, string> chrom_file_map;

string
get_state_bits(const vector<size_t> &copy_states,
    const size_t offset,
    const size_t width) {
  // states are indexed as 0=M_1~M_L; L=I_0~I_L-1
  // plus offset
  string bits(width, '0');
  for (vector<size_t>::const_iterator i = copy_states.begin();
      i < copy_states.end(); ++i) {
    assert(*i >= offset);
    const size_t idx = *i - offset;
    if (idx < width)
      bits[idx] = '1';
  }
  return bits;
}

string
get_state_CIGAR(const vector<size_t> &copy_states, const size_t width) {
  // This function works only for complete matrix indices, therefore no offset
  // needed.
  string cigar = "";
  size_t query_counter = 1;
  string prev_type = state_type_to_str(width, copy_states.front());
  state prev_state = state(width, copy_states.front());
  // 5' truncation
  if (prev_state.idx > 1)
    cigar = cigar + to_string(prev_state.idx-1) + "N";
  for (vector<size_t>::const_iterator i = copy_states.begin()+1;
      i < copy_states.end(); ++i) {
    state curr_state = state(width, *i);
    if (curr_state.stateint == prev_state.stateint) {
      if (curr_state.idx - prev_state.idx > 1 && curr_state.idx > prev_state.idx) {
        // deletion
        cigar = cigar + to_string(query_counter) + prev_type;
        cigar = cigar + to_string(curr_state.idx-prev_state.idx-1) + "D";
        query_counter = 0;
      }
    }
    else {
      // change of type
      cigar = cigar + to_string(query_counter) + prev_type;
      query_counter = 0;
    }
    ++query_counter;
    prev_type = state_type_to_str(width, *i);
    prev_state = curr_state;
  }
  cigar = cigar + to_string(query_counter) + prev_type;
  // 3' truncation
  if (prev_state.idx < width)
    cigar = cigar + to_string(width-prev_state.idx) + "N";
  return cigar;
}

void
identify_repeats(const ProfileHMM &hmm,
    const size_t bg_state,
    const size_t offset,
    const string &chr_seq,
    const vector<size_t> &states,
    const string chr_name,
    const bool SENSE_STRAND,
    const bool NO_LOG_ODDS,
    const bool USE_BITS,
    vector<GenomicRegion> &coordinates,
    vector<string> &state_bits) {
  // use bg_state = model_len, offset = 0 for collapsed matrix decoding;
  // use bg_state = model_len+1, offset = 1 for complete matrix decoding.
  const size_t INTERNAL_LOOP_LEN = 10;
  const size_t model_len = hmm.Length();
  const size_t chr_len = states.size();
  size_t start = 0, end = 0;
  for (vector<size_t>::const_iterator i = states.begin();
      i < states.end()-1; ++i) {
    vector<size_t>::const_iterator j = next(i);
    // find start and end of region
    if (*i != bg_state && i == states.begin())
      start = 0;
    else if (*i == bg_state && *j != bg_state)
      start = j - states.begin();
    else if ((*i != bg_state && *j == bg_state)
        || (*i != bg_state && j == states.end()))
      end = j - states.begin();
    else if (*i != bg_state && *j != bg_state) {
      // check if internal loop exists
      state prev(model_len, *i);
      state curr(model_len, *j);
      if (prev.idx > curr.idx + INTERNAL_LOOP_LEN)
        end = j - states.begin();
    }
    if (end > start) {
      const string obs = chr_seq.substr(start, end - start);
      // this is the score using log-odds
      //double score =
      //  hmm.PosteriorProb_c(true, obs) / obs.length();
      // this is the score using reverse sequenc as normalization
      string obs_rev(obs.rbegin(), obs.rend());
      double score =
        hmm.PosteriorProb(!NO_LOG_ODDS, obs)
          - hmm.PosteriorProb(!NO_LOG_ODDS, obs_rev);
      string name = "X";
      if (SENSE_STRAND) {
        GenomicRegion new_copy(chr_name, start,
          end, name, score, '+');
        coordinates.push_back(new_copy);
      }
      else {
        GenomicRegion new_copy(chr_name, chr_len - end,
          chr_len - start, name, score, '-');
        coordinates.push_back(new_copy);
      }

      // get bit vector for matching states
      vector<size_t>::const_iterator first = states.begin() + start;
      vector<size_t>::const_iterator last = states.begin() + end;
      vector<size_t> copy_states(first, last);
      string bits;
      if (USE_BITS)
        bits = get_state_bits(copy_states, offset, model_len);
      else
        bits = get_state_CIGAR(copy_states, model_len);
      state_bits.push_back(bits);
      start = end;
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
  const double stdev = std::sqrt(std::inner_product(scores.begin(),
    scores.end(), scores.begin(), 0.0) / (scores.size()-1) - mean * mean);
  for (vector<GenomicRegion>::iterator i = coordinates.begin();
      i < coordinates.end(); ++i) {
    (*i).set_score(((*i).get_score() - mean) / stdev);
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
    bool NO_LOG_ODDS = false;
    bool USE_BITS = false;
    double Z_CUTOFF = 1.0;
    string chrom_file, in_par, out_file;
    string fasta_suffix = "fa";
    size_t NUM_THREAD = 1;

    OptionParser opt_parse(strip_path(argv[0]), "Program for finding repeats.",
        "-c <chroms> <profile-HMM params file>");

    opt_parse.add_opt("chrom", 'c',
      "File or directory of chroms (FASTA format; .fa suffix)",
      true, chrom_file);
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)",
      false, out_file);
    opt_parse.add_opt("no-log", 'n', "Do not use log odds ratio as posterior.",
      false, NO_LOG_ODDS);
    opt_parse.add_opt("z-cutoff", 'z',
      "Z-score cutoff for occurrence identifiation. Default: 1.0;\
        setting to 0 will disable this functionality.",
      false, Z_CUTOFF);
#ifdef OPENMP
    opt_parse.add_opt("process", 'p', "Set the number of processes for parallelization. \
        Default: 1.", false, NUM_THREAD);
#endif
    opt_parse.add_opt("bits", 'b', "Use state bits instead of CIGAR.",
      false, USE_BITS);
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
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_FAILURE;
    }
    if (Z_CUTOFF < 0) {
      cerr << "Z-score cutoff must be positive." << endl;
      return EXIT_FAILURE;
    }
    in_par = leftover_args.front();
#ifdef OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(NUM_THREAD);
#endif

    std::ofstream of;
    if (!out_file.empty()) of.open(out_file.c_str());
    std::ostream out(out_file.empty() ? std::cout.rdbuf() : of.rdbuf());

    if (VERBOSE)
      cerr << "[LOADING HMM]" << endl;
    ProfileHMM hmm(in_par);
    if (VERBOSE)
      cerr << "\tMODEL LENGTH=" << hmm.Length() << endl;

    if (VERBOSE)
      cerr << "[LOADING GENOME]" << endl;
    chrom_file_map chrom_files;
    identify_chromosomes(chrom_file, fasta_suffix, chrom_files);
    vector<string> chr_name, chr_seq;
    for (chrom_file_map::const_iterator chrom = chrom_files.begin();
        chrom != chrom_files.end(); ++chrom) {
      read_fasta_file(chrom->second, chr_name, chr_seq);
      if (chr_name.size() < 1)
        throw SMITHLABException("could not find any sequence in: "
            + chrom->second);
    }
    if (VERBOSE)
      cerr << "[SCANNING " << chr_seq.size() << " SEQUENCES LOADED FROM "
        << chrom_files.size() << " FILES]" << endl;

    vector<GenomicRegion> coordinates;
    vector<string> state_bits;
    size_t progress = 0;
    const size_t TICK = NUM_THREAD*4+1>50 ? NUM_THREAD*4+1 : 50;
#ifdef OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < chr_seq.size(); ++i) {
      vector<size_t> states;
      if (DEBUG) {
#ifdef OPENMP
#pragma omp critical (debug_info)
#endif
        cerr << chr_name[i] << endl;
      }
      hmm.PosteriorDecoding(false, DEBUG, !NO_LOG_ODDS, chr_seq[i], states);
#ifdef OPENMP
#pragma omp atomic
#endif
      ++progress;
#ifdef OPENMP
#pragma omp critical (update_results1)
#endif
      {
      identify_repeats(hmm, hmm.Length()+1, 1, chr_seq[i], states,
        chr_name[i], true, NO_LOG_ODDS, USE_BITS, coordinates, state_bits);
      }
      if (VERBOSE && progress % TICK == 0) {
#ifdef OPENMP
#pragma omp critical (progress1)
#endif
        cerr << "\r\tPROCESSED " << 100*progress/chr_seq.size()/2 << "%";
      }
    }
#ifdef OPENMP
#pragma omp barrier
#endif
    hmm.ComplementBackground();
#ifdef OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < chr_seq.size(); ++i) {
      revcomp_inplace(chr_seq[i]);
      vector<size_t> states;
      hmm.PosteriorDecoding(false, false, !NO_LOG_ODDS, chr_seq[i], states);
#ifdef OPENMP
#pragma omp atomic
#endif
      ++progress;
#ifdef OPENMP
#pragma omp critical (update_results2)
#endif
      {
      identify_repeats(hmm, hmm.Length()+1, 1, chr_seq[i], states,
        chr_name[i], false, NO_LOG_ODDS, USE_BITS, coordinates, state_bits);
      }
      if (VERBOSE && progress % TICK == 0) {
#ifdef OPENMP
#pragma omp critical (progress2)
#endif
        cerr << "\r\tPROCESSED " << 100*progress/chr_seq.size()/2 << "%";
      }
    }
#ifdef OPENMP
#pragma omp barrier
#endif
    if (Z_CUTOFF - 0.0 > 1e-10) {
      if (VERBOSE)
        cerr << endl << "[CALCULATING Z-SCORES]" << endl;
      zscore_filter(coordinates, Z_CUTOFF);
    }
    for (vector<GenomicRegion>::const_iterator i = coordinates.begin();
      i < coordinates.end(); ++i)
    {
      if (Z_CUTOFF - 0.0 < 1e-10 || i->get_score() > Z_CUTOFF)
        out << *i << "\t" << state_bits[i-coordinates.begin()] << endl;
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
