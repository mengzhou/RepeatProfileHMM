/*    align: a program for pair-wise global alignment
 *    Copyright (C) 2017 University of Southern California and
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
#include <fstream>
#include <iomanip>

#include "smithlab_os.hpp"
#include "OptionParser.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::istringstream;
using std::ostream;
using std::ofstream;

size_t
g_align_hamming(const string &seq1, const string &seq2) {
  const size_t L1 = seq1.size();
  const size_t L2 = seq2.size();
  vector<size_t> prev(L1+1, 0);
  vector<size_t> curr(L1+1, 0);
  //gap = 1;

  for (size_t i=1; i<L1+1; ++i) {
    prev[i] = i;
  }
  for (size_t j=1; j<L2+1; ++j) {
    curr[0] = j;
    for (size_t i=1; i<L1+1; ++i) {
      const size_t left = curr[i-1] + 1;
      const size_t up = prev[i] + 1;
      const size_t w_gap = left < up ? left : up;
      const size_t wo_gap = prev[i-1] + (seq1[i-1] != seq2[j-1]);
      curr[i] = wo_gap < w_gap ? wo_gap : w_gap;
    }
    prev = curr;
  }
  return prev.back();
}

int
main(int argc, const char **argv) {
  try {
    bool VERBOSE = false;
    string in_file, out_file;

    OptionParser opt_parse(strip_path(argv[0]),
        "Program pair-wise sequence global alignment.",
        "[Options] -o <Output.txt> <Input.fa | stdin>");

    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)",
      false, out_file);
    opt_parse.add_opt("verbose", 'v',
      "Verbose mode. Only use this for interactive enviroment.", false, VERBOSE);

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
    in_file = leftover_args.front();
    std::ofstream of;
    if (!out_file.empty()) of.open(out_file.c_str());
    std::ostream out(out_file.empty() ? std::cout.rdbuf() : of.rdbuf());

    if (VERBOSE)
      cerr << "[LOADING SEQUENCES]" << endl;
    vector<string> seq_name, seq_seq;
    read_fasta_file(in_file, seq_name, seq_seq);
    if (VERBOSE)
      cerr << "\tSEQ_FOUND=" << seq_name.size() << endl;
    if (seq_name.size() < 1)
      throw SMITHLABException("could not find any sequence in: "
          + in_file);

    const size_t total = seq_name.size()*(seq_name.size()-1)/2;
    vector<vector<size_t> > sc(seq_name.size(), vector<size_t>(seq_name.size(), 0));
    size_t counter = 0;
    for (size_t i = 0; i < seq_name.size(); ++i)
      out << "\t" << seq_name[i];
    out << endl;
    for (size_t i = 0; i < seq_seq.size(); ++i) {
      out << seq_name[i];
      for (size_t j = 0; j < seq_seq.size(); ++j) {
        if (j>i) {
          sc[i][j] = g_align_hamming(seq_seq[i], seq_seq[j]);
          ++counter;
          if (VERBOSE && counter % 100 == 1) {
            cerr << "\r[";
            for (int t=0; t < 20; ++t) {
              const int progress = 20.0*counter/total;
              if (t < progress)
                cerr << "=";
              else if (t == progress)
                cerr << ">";
              else
                cerr << " ";
            }
            const double perc = 100.0*counter/total;
            cerr << "] " << std::fixed << std::setprecision(2) << perc << "%";
            cerr.flush();
          }
          out << "\t" << sc[i][j];
        }
        else if (j==i)
          out << "\tNA";
        else
          out << "\t" << sc[j][i];
      }
      out << endl;
    }
    if (VERBOSE) cerr << endl;
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
