/*    posterior: a program for computing posteriors using profile-HMM
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

int
main (int argc, const char **argv) {
  try {
    bool VERBOSE = false;
    bool DEBUG = false;
    bool NO_LOG_ODDS = false;
    double z_cutoff = 1.0;
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
    opt_parse.add_opt("no-log", 'n', "Do not use log odds ratio as posterior.",
      false, NO_LOG_ODDS);
    opt_parse.add_opt("z-cutoff", 'z',
      "Z-score cutoff for occurrence identifiation. Default: 1.0;\
        setting to 0 will disable this functionality.",
      false, z_cutoff);
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
    if (z_cutoff < 0) {
      cerr << "Z-score cutoff must be positive." << endl;
      return EXIT_FAILURE;
    }
    in_par = leftover_args.front();

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
      vector<string> state_bits;

      if (VERBOSE)
        cerr << "[SCANNING " << chrom->second
          << " " << file_counter++ << " OF "
          << chrom_files.size() << "]" << endl;
      vector<GenomicRegion> coordinates;
      for (size_t i = 0; i < chr_seq.size(); ++i) {
        if (VERBOSE)
          cerr << "\t" << i+1 << "/" << chr_seq.size()
            << "\t" << chr_name[i] << endl;
        if (DEBUG)
          cerr << chr_name[i] << endl;
        hmm.PosteriorDecoding(false, DEBUG, !NO_LOG_ODDS, chr_seq[i], states);
        identify_repeats(hmm, hmm.Length()+1, 1, chr_seq[i], states,
          chr_name[i], true, coordinates, state_bits);
      }
      for (vector<GenomicRegion>::const_iterator i = coordinates.begin();
        i < coordinates.end(); ++i)
      {
        out << *i << "\t" << state_bits[i-coordinates.begin()] << endl;
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
