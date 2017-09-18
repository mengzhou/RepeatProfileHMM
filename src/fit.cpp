/*    fit: a program for fitting a profile-HMM based on input data using
 *    Baum-Welch training algorithm
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
#include <fstream>

#include "ProfileHMM.hpp"
#include "OptionParser.hpp"
#include "smithlab_os.hpp"

using std::cout;
using std::endl;

using namespace std;
typedef unordered_map<string, string> chrom_file_map;

int
main (int argc, const char **argv) {
  try {
    bool VERBOSE = false;
    string chrom_file, input_file, in_par, out_par;
    size_t max_itr = 50;
    double tolerance = 1e-4;
    string fasta_suffix = "fa";

    OptionParser opt_parse(argv[0], "Program for profile-HMM.");

    opt_parse.add_opt("seq", 's', "FASTA file of input sequences.",
      true, input_file);
    opt_parse.add_opt("itr", 'r', "Max iteration.", false, max_itr);
    opt_parse.add_opt("tolerance", 't', "EM tolerance.", false, tolerance);
    opt_parse.add_opt("in-params", 'i', "Input parameters file.", true, in_par);
    opt_parse.add_opt("out-params", 'o', "Output parameters file.", false, out_par);
    opt_parse.add_opt("verbose", 'v', "Verbose mode.", false, VERBOSE);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);

    if (in_par.compare(out_par) == 0) {
      cout << opt_parse.help_message() << endl
        << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    ofstream outfs;
    if (!out_par.empty()) outfs.open(out_par.c_str());
    ostream out(out_par.empty() ? cout.rdbuf() : outfs.rdbuf());

    if (VERBOSE)
      cerr << "[LOADING HMM]" << endl;
    ProfileHMM hmm = ProfileHMM(in_par);
    if (VERBOSE)
      cerr << "\tMODEL LENGTH=" << hmm.Length() << endl;

    if (VERBOSE)
      cerr << "[LOADING SEQUENCES]" << endl;
    chrom_file_map chrom_files;
    identify_chromosomes(input_file, fasta_suffix, chrom_files);
    vector<string> chr_name, chr_seq;
    for (chrom_file_map::const_iterator chrom = chrom_files.begin();
        chrom != chrom_files.end(); ++chrom) {
      read_fasta_file(chrom->second, chr_name, chr_seq);
    }
    if (VERBOSE)
      cerr << "\tSEQUENCES_LOADED=" << chr_seq.size() << endl;

    if (VERBOSE)
      cerr << "[FITTING MODEL PARAMETERS]" << endl;
    hmm.Train(VERBOSE, tolerance, max_itr, chr_seq);

    hmm.Print(out, false);
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
