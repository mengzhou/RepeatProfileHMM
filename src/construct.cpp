/*    construct: a program for constructing profile-HMM from multiple
 *    sequence alignment
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
#include <sstream>
#include <string>
#include <initializer_list>
#include <numeric>
#include <cmath>
#include <iomanip>
#include <tr1/unordered_map>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "ProfileHMM.hpp"
#include "matrix_utils.hpp"

using std::vector;
using std::pair;
using std::string;
using std::tr1::unordered_map;
using std::ifstream;
using std::istringstream;
using std::ostream;
using std::ofstream;
using std::cout;
using std::cerr;
using std::endl;

void
load_alignment(bool VERBOSE, const string infile,
    unordered_map<string, string> &msa) {
  ifstream infs(infile.c_str());
  string line;

  while (getline(infs, line)) {
    istringstream iss(line);
    string buffer, name, seq;
    vector<string> tokens;

    while (iss >> buffer)
      tokens.push_back(buffer);
    if (tokens.size() == 2) {
      name = tokens.front();
      seq = tokens.back();
      const unordered_map<string, string>::const_iterator repeat(msa.find(name));
      if (repeat == msa.end()) {
        msa[name] = seq;
      }
      else {
        msa[name] = msa[name] + seq;
      }
    }
  }
  if (VERBOSE)
    cerr << "\t" << msa.size() << " SEQUENCES LOADED." << endl;
}

void
set_transition_prior(matrix &tran_prior) {
  // dimensionality of transition prior:
  // 4x4 - M, I, D_i (internal deletion), D_1 (truncation)
  // M to M, I, D_i, D_1
  tran_prior.push_back(vector<double>({0.8, 0.1, 0.1, 0.0}));
  // I to M, I, D_i, D_1
  tran_prior.push_back(vector<double>({0.4, 0.5, 0.1, 0.0}));
  // D_i to M, I, D_i, D_1
  tran_prior.push_back(vector<double>({0.4, 0.1, 0.5, 0.0}));
  // D_1 to M, I, D_i, D_1
  tran_prior.push_back(vector<double>({1.0, 0.0, 0.0, 0.0}));
}

void
set_emission_prior(matrix &emis_prior) {
  // dimensionality of emission prior:
  // 2x4 - M, I for each nucleotide A,C,G,T
  // M
  // human L1PA6 ORF2
  emis_prior.push_back(
      vector<double>({0.4138828,0.2082860,0.1713843,0.2064469}));
  // I
  // human chrM
  emis_prior.push_back(
      vector<double>({0.308551,0.313318,0.131555,0.246575}));
}

double
get_5prime_truncation_prior(const size_t model_len,
    const size_t idx) {
  /*
  const double first_col_weight = 0.5;
  if (idx == 1)
    return first_col_weight;
  else
    return (1.0 - first_col_weight) / (model_len - 1);
    */
  return -log(model_len);
}

double
get_3prime_truncation_prior(const size_t model_len,
    const size_t idx) {
  /*
  const double last_col_weight = 0.5;
  if (idx == model_len)
    return last_col_weight;
  else
    return (1.0 - last_col_weight) / (model_len - 1);
    */
  return -log(model_len);
}

double
posterior_transition_score(const unordered_map<string, string> &msa,
    const size_t start, const size_t end,
    const matrix &tran_prior) {
  // there are 3 types of states: M, I, D, so the number of possible
  // transitions between position start and end is 9=3x3.
  // index for states: M - 0, I - 1, D - 2 corresponding to the
  // transition prior matrix
  vector<vector<size_t> > count(4, vector<size_t>(4, 0));
  for (unordered_map<string, string>::const_iterator i = msa.begin();
      i != msa.end(); ++i) {
    // convert base to int with - converted to 4, skipping - in insertion
    // columns (start < nt < end)
    vector<size_t> seq_int;
    string::const_iterator left_bound;
    if (start == 0) {
      seq_int.push_back(5);
      left_bound = (i->second).begin();
    }
    else
      left_bound = (i->second).begin() + start - 1;
    for (string::const_iterator nt = left_bound;nt < (i->second).begin() + end
        && nt < (i->second).end(); ++nt) {
      if ((start > 0 && nt == left_bound)
          || nt == (i->second).begin() + end - 1
          || base2int(*nt) < 4)
        seq_int.push_back(base2int(*nt));
    }
    if (end > (i->second).size())
      seq_int.push_back(5);
    // now get transition counts
    for (vector<size_t>::const_iterator s = seq_int.begin();
        s < seq_int.end() - 1; ++s) {
      const size_t state_first = baseint2stateint(*s, s == seq_int.begin());
      const size_t state_second = baseint2stateint(*(s+1),
          s+1 == seq_int.end() - 1);
      assert(state_first < 4 && state_second < 4);
      ++count[state_first][state_second];
    }
  }

  // get posterior probabilities
  double total = 0.0;
  for (size_t i = 0; i < 4; ++i) {
    vector<double> row;
    for (vector<size_t>::const_iterator k = count[i].begin();
        k < count[i].end(); ++k)
      row.push_back(static_cast<double>(*k));
    row = combine_normalize(row, tran_prior[i], false);
    log_transform_vec(row);
    for (size_t j = 0; j < 4; ++j)
      total = total + count[i][j] * row[j];
  }
  return total;
}

double
posterior_emission_score(const unordered_map<string, string> &msa,
    const size_t start, const size_t end,
    const vector<double> &emis_prior) {
  if (start == ((msa.begin())->second).size())
    return 0.0;
  vector<size_t> count(4, 0);
  for (unordered_map<string, string>::const_iterator i = msa.begin();
      i != msa.end(); ++i) {
    for (size_t pos = start; pos < end; ++pos) {
      const size_t baseint = base2int((i->second)[pos]);
      if (baseint < 4) ++count[baseint];
    }
  }
  double sum = 0.0;
  for (size_t i = 0; i < 4; ++i)
    sum = sum + count[i] + emis_prior[i];
  double total = 0.0;
  for (size_t i = 0; i < 4; ++i)
    total = total + count[i] * (log(count[i] + emis_prior[i]) - log(sum));

  return total;
}

void
find_marked_cols_map(const bool VERBOSE,
    vector<bool> &marked,
    const unordered_map<string, string> &msa,
    const matrix &tran_prior,
    const matrix &emis_prior,
    const double lambda) {
  const size_t seq_len = marked.size();
  vector<double> s(seq_len+2, LOG_ZERO);
  vector<size_t> trace_back(seq_len+2, 0);
  s[0] = 0.0;

  for (size_t j = 1; j <= seq_len + 1; ++j) {
    vector<double> v;
    if (VERBOSE && j % 10 == 0)
      cerr << "\tPROCESSED " << 100 * j / seq_len << "%" << endl;
    for (size_t i = 0; i < j; ++i) {
      v.push_back(s[i] + posterior_transition_score(msa, i, j, tran_prior)
          + posterior_emission_score(msa, j, j+1, emis_prior.front())
          + posterior_emission_score(msa, i+1, j, emis_prior.back())
          + lambda);
    }
    s[j] = *max_element(v.begin(), v.end());
    trace_back[j] = argmax_vec(v);
  }

  size_t i = trace_back.back();
  marked[i-1] = true;
  while (i > 0) {
    i = trace_back[i];
    if (i > 0)
      marked[i-1] = true;
  }
}

void
find_marked_cols_heu(vector<bool> &marked,
    const unordered_map<string, string> &msa) {
  const size_t seq_len = marked.size();
  for (size_t col = 0; col < seq_len; ++col) {
    size_t match_count = 0;
    for (unordered_map<string, string>::const_iterator i = msa.begin();
        i != msa.end(); ++i) {
      if (base2int(i->second[col]) < 4) ++match_count;
    }
    if (match_count > msa.size() / 2)
      marked[col] = true;
  }
}

void
write_hmm_parameter(const bool VERBOSE,
    ostream &out,
    const matrix &transition,
    const matrix &emission) {
  const size_t model_len = (transition.size() - 2) / 3;
  out << "# Model length: " << model_len << endl;
  for (matrix::const_iterator i = transition.begin();
      i < transition.end(); ++i) {
    for (vector<double>::const_iterator j = (*i).begin();
        j < (*i).end(); ++j) {
      out << *j << "\t";
    }
    out << endl;
  }
  out << "//" << endl;
  for (matrix::const_iterator i = emission.begin();
      i < emission.end(); ++i) {
    for (vector<double>::const_iterator j = (*i).begin();
        j < (*i).end(); ++j) {
      out << *j << "\t";
    }
    out << endl;
  }
  out << "//" << endl;
}

int
main (int argc, const char **argv) {
  try {
    bool VERBOSE = false;
    bool DEBUG = false;
    bool USE_HEU = false;
    double lambda = 0.0;
    string outf;

    OptionParser opt_parse(strip_path(argv[0]),
        "Construct profile-HMM from multiple alignment.");
    opt_parse.add_opt("out", 'o', "Output parameters file.", false, outf);
    opt_parse.add_opt("lambda", 'l',
        "Model length adjusting parameter. Default: 0", false, lambda);
    opt_parse.add_opt("heuristic", 'u', "Heuristic mode.", false, USE_HEU);
    opt_parse.add_opt("verbose", 'v', "Verbose mode.", false, VERBOSE);
    opt_parse.add_opt("deubg", 'd', "Print debug information.", false, DEBUG);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cout << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cout << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested() || leftover_args.size() != 1) {
      cout << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    const string inf = leftover_args.front();

    // 1. load alignment
    if (VERBOSE)
      cerr << "[LOADING ALIGNMENT]" << endl;
    unordered_map<string, string> msa;
    load_alignment(VERBOSE, inf, msa);

    // 2. compute recursion to find marked (match) columns
    if (VERBOSE)
      cerr << "[FINDING OPTIMAL MATCH COLUMNS]" << endl;
    vector<bool> marked_cols((msa.begin()->second).length(), false);
    matrix tran_prior, emis_prior;
    set_transition_prior(tran_prior);
    set_emission_prior(emis_prior);
    if (USE_HEU)
      find_marked_cols_heu(marked_cols, msa);
    else
      find_marked_cols_map(VERBOSE,
          marked_cols, msa, tran_prior, emis_prior, lambda);
    
    // 3. estimate probabilities for individual states
    if (VERBOSE)
      cerr << "[ESTIMATING PARAMETERS]" << endl;
    matrix transition, emission;
    const size_t model_len =
      count(marked_cols.begin(), marked_cols.end(), true);
    transition.resize(model_len*3+2, vector<double>(model_len*3+2, 0.0));
    emission.resize(model_len*3+2, vector<double>(4, 0.0));

    // 3.1 transitions
    // Seuqnecs out of the first and last marked columns is not considered.
    // They belong to the background
    const vector<bool>::const_iterator first_marked =
      find(marked_cols.begin(), marked_cols.end(), true);
    const vector<bool>::const_iterator last_marked =
      find(marked_cols.rbegin(), marked_cols.rend(), true).base() - 1;
    state left_end(0ul, 1, 2);
    state right_end(0ul, model_len, 2);
    for (unordered_map<string, string>::const_iterator i = msa.begin();
        i != msa.end(); ++i) {
      // get truncation by jumping to the first and the last non-gap character
      // for the marked columns
      vector<bool>::const_iterator curr_itr = first_marked;
      while(curr_itr < marked_cols.begin() && 
          (!*curr_itr || i->second[curr_itr - marked_cols.begin()] == '-'))
        ++curr_itr;
      const vector<bool>::const_iterator first_nontruncated = curr_itr;
      curr_itr = last_marked; // last_marked is a reverse iterator
      while(curr_itr > first_nontruncated &&
          (!*curr_itr || i->second[curr_itr - marked_cols.begin()] == '-'))
        --curr_itr;
      const vector<bool>::const_iterator last_nontruncated = curr_itr;
      // get counts for the first marked column associating with truncation
      state curr_state(i->second[first_nontruncated - marked_cols.begin()],
          1, true);
      transition[left_end.index(model_len)][curr_state.index(model_len)] =
        transition[left_end.index(model_len)][curr_state.index(model_len)] + 1;
      emission[curr_state.index(model_len)][curr_state.baseint] =
        emission[curr_state.index(model_len)][curr_state.baseint] + 1;

      vector<state> seq_state(1, curr_state);
      size_t col_idx = 1;
      for (vector<bool>::const_iterator j = first_nontruncated + 1;
          j <= last_nontruncated; ++j) {
        if (*j) ++col_idx;
        curr_state = state(i->second[j - marked_cols.begin()],
            col_idx, *j);
        if (curr_state.isvalid())
          seq_state.push_back(curr_state);
      }
      for (vector<state>::const_iterator j = seq_state.begin(), k = j+1;
          k < seq_state.end(); ++j,++k) {
        transition[(*j).index(model_len)][(*k).index(model_len)] = 
          transition[(*j).index(model_len)][(*k).index(model_len)] + 1;
        if ((*k).stateint < 2)
          emission[(*k).index(model_len)][(*k).baseint] =
            emission[(*k).index(model_len)][(*k).baseint] + 1;
      }

      // 3' truncation
      transition[seq_state.back().index(model_len)][right_end.index(model_len)]
        = transition[seq_state.back().index(model_len)]
          [right_end.index(model_len)]
        + 1;
    }

    // combine with priors and normalize
    set_transition_prior(tran_prior);
    set_emission_prior(emis_prior);
    log_transform_matrix(transition);
    log_transform_matrix(emission);
    state curr_m, curr_i, curr_d, next_m, next_i, next_d;
    vector<double> row, prior;
    size_t idx = 1;
    for (idx = 1; idx < model_len - 1; ++idx) {
      curr_m = state(0ul, idx, 0);
      curr_i = state(0ul, idx, 1);
      curr_d = state(0ul, idx, 2);
      next_m = state(0ul, idx+1, 0);
      next_d = state(0ul, idx+1, 2);
      // M_1 ~ M_L-2
      row = vector<double>({
          transition[curr_m.index(model_len)][next_m.index(model_len)],
          transition[curr_m.index(model_len)][curr_i.index(model_len)],
          transition[curr_m.index(model_len)][next_d.index(model_len)],
          transition[curr_m.index(model_len)][right_end.index(model_len)]
          });
      prior = tran_prior[0];
      log_transform_vec(prior);
      prior[3] = get_3prime_truncation_prior(model_len, idx);
      row = combine_normalize(row, prior, true);
      transition[curr_m.index(model_len)][next_m.index(model_len)] = row[0];
      transition[curr_m.index(model_len)][curr_i.index(model_len)] = row[1];
      transition[curr_m.index(model_len)][next_d.index(model_len)] = row[2];
      transition[curr_m.index(model_len)][right_end.index(model_len)] = row[3];
      // I_1 ~ I_L-2
      row = vector<double>({
          transition[curr_i.index(model_len)][next_m.index(model_len)],
          transition[curr_i.index(model_len)][curr_i.index(model_len)],
          transition[curr_i.index(model_len)][next_d.index(model_len)]
          });
      prior = tran_prior[1];
      prior.erase(prior.begin() + 3);
      log_transform_vec(prior);
      row = combine_normalize(row, prior, true);
      transition[curr_i.index(model_len)][next_m.index(model_len)] = row[0];
      transition[curr_i.index(model_len)][curr_i.index(model_len)] = row[1];
      transition[curr_i.index(model_len)][next_d.index(model_len)] = row[2];
      // D_2 ~ D_L-2
      if (idx > 1) {
        row = vector<double>({
            transition[curr_d.index(model_len)][next_m.index(model_len)],
            transition[curr_d.index(model_len)][curr_i.index(model_len)],
            transition[curr_d.index(model_len)][next_d.index(model_len)]
            });
        prior = tran_prior[2];
        prior.erase(prior.begin() + 3);
        log_transform_vec(prior);
        row = combine_normalize(row, prior, true);
        transition[curr_d.index(model_len)][next_m.index(model_len)] = row[0];
        transition[curr_d.index(model_len)][curr_i.index(model_len)] = row[1];
        transition[curr_d.index(model_len)][next_d.index(model_len)] = row[2];
      }
    }
    // M_L-1
    curr_m = state(0ul, model_len - 1, 0);
    curr_i = state(0ul, model_len - 1, 1);
    next_m = state(0ul, model_len, 0);
    row = vector<double>({
        transition[curr_m.index(model_len)][next_m.index(model_len)],
        transition[curr_m.index(model_len)][curr_i.index(model_len)],
        transition[curr_m.index(model_len)][right_end.index(model_len)]
        });
    prior = tran_prior[0];
    prior.erase(prior.begin() + 2); // no next_d
    log_transform_vec(prior);
    prior.back() = get_3prime_truncation_prior(model_len, idx);
    row = combine_normalize(row, prior, true);
    transition[curr_m.index(model_len)][next_m.index(model_len)] = row[0];
    transition[curr_m.index(model_len)][curr_i.index(model_len)] = row[1];
    transition[curr_m.index(model_len)][right_end.index(model_len)] = row[2];
    // I_L-1
    curr_i = state(0ul, model_len - 1, 1);
    next_m = state(0ul, model_len, 0);
    row = vector<double>({
        transition[curr_i.index(model_len)][next_m.index(model_len)],
        transition[curr_i.index(model_len)][curr_i.index(model_len)]
        });
    prior = tran_prior[1];
    prior.erase(prior.begin() + 2, prior.end()); // no next_d
    log_transform_vec(prior);
    row = combine_normalize(row, prior, true);
    transition[curr_i.index(model_len)][next_m.index(model_len)] = row[0];
    transition[curr_i.index(model_len)][curr_i.index(model_len)] = row[1];
    // D_L-1
    curr_i = state(0ul, model_len - 1, 1);
    curr_d = state(0ul, model_len - 1, 2);
    next_m = state(0ul, model_len, 0);
    row = vector<double>({
        transition[curr_d.index(model_len)][next_m.index(model_len)],
        transition[curr_d.index(model_len)][curr_i.index(model_len)]
        });
    prior = tran_prior[2];
    prior.erase(prior.begin() + 2, prior.end()); // no next_d
    log_transform_vec(prior);
    row = combine_normalize(row, prior, true);
    transition[curr_d.index(model_len)][next_m.index(model_len)] = row[0];
    transition[curr_d.index(model_len)][curr_i.index(model_len)] = row[1];
    // M_L
    transition[state(0ul, model_len, 0).index(model_len)]
      [right_end.index(model_len)] = 0.0;
    // M_0
    transition[0][left_end.index(model_len)] = log(0.05);
    transition[0][state(0ul, 0, 1).index(model_len)] = log(0.95);
    // D_1
    transition[left_end.index(model_len)] =
      normalize_vec(transition[left_end.index(model_len)], true);
    for (size_t i = 1; i <= model_len; ++i) {
      transition[left_end.index(model_len)][i] = 
        log_sum_log(transition[left_end.index(model_len)][i],
            get_5prime_truncation_prior(model_len, i));
    }
    transition[left_end.index(model_len)] =
      normalize_vec(transition[left_end.index(model_len)], true);
    // I_0
    transition[state(0ul, 0, 1).index(model_len)]
      [state(0ul, 0, 1).index(model_len)] = log(0.8);
    transition[state(0ul, 0, 1).index(model_len)]
      [left_end.index(model_len)] = log(0.18);
    transition[state(0ul, 0, 1).index(model_len)].back() = log(0.02);
    // D_L
    transition[right_end.index(model_len)]
      [state(0ul, 0, 1).index(model_len)] = log(0.98);
    transition[right_end.index(model_len)].back() = log(0.02);

    // 3.2 emissions
    prior = emis_prior[1];
    log_transform_vec(prior);
    for (size_t i = 1; i <= model_len; ++i) {
      emission[state(0ul, i, 0).index(model_len)] =
        combine_normalize(emission[state(0ul, i, 0).index(model_len)],
            prior, true, std::make_pair(0.9, 0.1));
      if (i < model_len)
      emission[state(0ul, i, 1).index(model_len)] =
        combine_normalize(emission[state(0ul, i, 1).index(model_len)],
            prior, true, std::make_pair(0.9, 0.1));
    }
    emission[state(0ul, 0, 1).index(model_len)] = prior;

    ofstream outfs;
    if (!outf.empty()) outfs.open(outf.c_str());
    ostream out(outf.empty() ? cout.rdbuf() : outfs.rdbuf());

    write_hmm_parameter(VERBOSE, out, transition, emission);
    if (VERBOSE)
      cerr << "[DONE]" << endl;

    // display alignment and marked columns for debug purposes
    if (DEBUG) {
      for (unordered_map<string, string>::const_iterator i = msa.begin();
          i != msa.end(); ++i)
        cout << i->first << "\t" << i->second << endl;
      cout << std::setw(25);
      for (size_t i = 0; i < marked_cols.size(); ++i) {
        if (marked_cols[i])
          cout << "*";
        else
          cout << " ";
      }
      cout << endl;
      vector<bool> marked_heu((msa.begin()->second).length(), false);
      find_marked_cols_heu(marked_heu, msa);
      cout << std::setw(25);
      for (size_t i = 0; i < marked_heu.size(); ++i) {
        if (marked_heu[i])
          cout << "#";
        else
          cout << " ";
      }
      cout << endl;

      ProfileHMM hmm(transition, emission);
      hmm.Print();
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
