/*    fisher: a program for finding repeats in the genome using profile-HMM
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
    size_t NUM_THREAD = 1;
    bool VERBOSE = false;
    bool DEBUG = false;
    string chrom_file, in_par, out_file;
    string fasta_suffix = "fa";

    OptionParser opt_parse(strip_path(argv[0]), "Program for finding repeats.",
        "-c <fasta> <profile-HMM params file>");

    opt_parse.add_opt("chrom", 'c',
      "File or directory of sequences (FASTA format; .fa suffix)",
      true, chrom_file);
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)",
                      false, out_file);
#ifdef OPENMP
    opt_parse.add_opt("process", 'p', "Set the number of processes for parallelization. \
        Default: 1.", false, NUM_THREAD);
#endif
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
      cerr << "[LOADING SEQUENCES]" << endl;
    vector<string> copies, copy_seq;
    read_fasta_file(chrom_file, copies, copy_seq);
    if (VERBOSE)
      cerr << "\tSEQUENCES FOUND=" << copies.size() << endl;
    vector<vector<double> > scores(copies.size(), vector<double>());
    size_t progress = 0;
    const size_t TICK = NUM_THREAD*4+1>50 ? NUM_THREAD*4+1 : 50;

#ifdef OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0;i < copies.size(); ++i) {
      vector<double> score;
      hmm.FisherScoreVector(copy_seq[i], score);
      scores[i] = score;
#ifdef OPENMP
#pragma omp atomic
#endif
      ++progress;
      if (VERBOSE && progress%TICK==0) {
#ifdef OPENMP
#pragma omp critical (progress)
#endif
        cerr << "\tProcessed "
          << 100*progress/copies.size() << "%" << endl;
      }
    }
#ifdef OPENMP
#pragma omp barrier
#endif
    for (size_t i = 0; i < scores.size(); ++i) {
      out << copies[i] << "\t" << copy_seq[i].size();
      for (vector<double>::const_iterator j = scores[i].begin();
          j < scores[i].end(); ++j)
        out << "\t" << *j;
      out << endl;
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
