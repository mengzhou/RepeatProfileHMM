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

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "ProfileHMM.hpp"
#include "matrix_utils.hpp"

using std::vector;
using std::pair;
using std::string;
using std::ifstream;
using std::istringstream;
using std::ostream;
using std::ofstream;
using std::cout;
using std::cerr;
using std::endl;

void
load_alignment(bool VERBOSE, const string infile,
    vector<string> &names, vector<string> &seqs) {
  // sequnce names are truncated in the output of clustalw and muscle,
  // therefore duplicated names are possible and they cannot be used
  // as lookup keys
  ifstream infs(infile.c_str());
  string line;
  size_t index = 0;

  while (getline(infs, line)) {
    istringstream iss(line);
    string buffer, name, seq;
    vector<string> tokens;
    if (line.length() == 0)
      index = 0;

    while (iss >> buffer)
      tokens.push_back(buffer);
    if (tokens.size() == 2 && buffer.find('*') == std::string::npos) {
      name = tokens.front();
      seq = tokens.back();
      ++index;
      if (index > names.size()) {
        names.push_back(name);
        seqs.push_back(seq);
      }
      else {
        seqs[index-1] = seqs[index-1] + seq;
      }
    }
  }
  assert(names.size() == seqs.size());
  if (VERBOSE)
    cerr << "\t" << names.size() << " SEQUENCES LOADED." << endl;
}

void
set_transition_prior(matrix &tran_prior) {
  // dimensionality of transition prior:
  // 4x4 - M, I, D_i (internal deletion), D_1 (truncation)
  // M to M, I, D_i, D_1
  tran_prior.push_back(vector<double>({0.8, 0.1, 0.1, 0.0}));
  // I to M, I, D_i, D_1
  tran_prior.push_back(vector<double>({0.5, 0.4, 0.1, 0.0}));
  // D_i to M, I, D_i, D_1
  tran_prior.push_back(vector<double>({0.5, 0.1, 0.4, 0.0}));
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
  const double first_col_weight = 0.5;
  if (idx == 1)
    return log(first_col_weight);
  else
    return log(1.0 - first_col_weight) - log(model_len - 1);
  //return -log(model_len);
}

double
get_3prime_truncation_prior(const size_t model_len,
    const size_t idx) {
  //const double last_col_weight = 0.5;
  //if (idx == model_len)
  //  return last_col_weight;
  //else
  //  return log(1.0 - last_col_weight) - log(model_len - 1);
  return -log(model_len-1);
}

double
posterior_transition_score(const vector<string> &seqs,
    const size_t start, const size_t end,
    const matrix &tran_prior) {
  // there are 3 types of states: M, I, D, so the number of possible
  // transitions between position start and end is 9=3x3.
  // index for states: M - 0, I - 1, D - 2 corresponding to the
  // transition prior matrix
  vector<vector<size_t> > count(4, vector<size_t>(4, 0));
  for (vector<string>::const_iterator i = seqs.begin();
      i < seqs.end(); ++i) {
    // convert base to int with - converted to 4, skipping - in insertion
    // columns (start < nt < end)
    vector<size_t> seq_int;
    string::const_iterator left_bound;
    if (start == 0) {
      seq_int.push_back(5);
      left_bound = (*i).begin();
    }
    else
      left_bound = (*i).begin() + start - 1;
    for (string::const_iterator nt = left_bound;nt < (*i).begin() + end
        && nt < (*i).end(); ++nt) {
      if ((start > 0 && nt == left_bound)
          || nt == (*i).begin() + end - 1
          || base2int(*nt) < 4)
        seq_int.push_back(base2int(*nt));
    }
    if (end > (*i).length())
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
posterior_emission_score(const vector<string> &seqs,
    const size_t start, const size_t end,
    const vector<double> &emis_prior) {
  if (start == seqs.front().size())
    return 0.0;
  vector<size_t> count(4, 0);
  for (vector<string>::const_iterator i = seqs.begin();
      i < seqs.end(); ++i) {
    for (size_t pos = start; pos < end; ++pos) {
      const size_t baseint = base2int((*i)[pos]);
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
    const vector<string> &seqs,
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
      v.push_back(s[i] + posterior_transition_score(seqs, i, j, tran_prior)
          + posterior_emission_score(seqs, j, j+1, emis_prior.front())
          + posterior_emission_score(seqs, i+1, j, emis_prior.back())
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
    const vector<string> &seqs) {
  const size_t seq_len = marked.size();
  for (size_t col = 0; col < seq_len; ++col) {
    size_t match_count = 0;
    for (vector<string>::const_iterator i = seqs.begin();
        i < seqs.end(); ++i) {
      if (base2int((*i)[col]) < 4) ++match_count;
    }
    if (match_count > seqs.size() / 2)
      marked[col] = true;
  }
}

void
find_marked_cols_manual(vector<bool> &marked,
    const vector<string> &names, const vector<string> &seqs) {
  vector<string>::const_iterator ref = std::find(names.begin(),
      names.end(), "consensus");
  if (ref == names.end()) {
    throw SMITHLABException("The consensus is not found in input alignment!");
  }
  const size_t ref_idx = ref - names.begin();
  const size_t seq_len = marked.size();
  for (size_t col = 0; col < seq_len; ++col) {
    if (base2int(seqs[ref_idx][col]) < 4) marked[col] = true;
  }
}

string
get_consensus(const vector<bool> &marked,
    const vector<string> &seqs) {
  const size_t seq_len = marked.size();
  string consensus = "";
  for (size_t col = 0; col < seq_len; ++col) {
    if (marked[col]) {
      vector<size_t> nt_count(4, 0);
      for (vector<string>::const_iterator i = seqs.begin();
          i < seqs.end(); ++i) {
        const int nt = base2int((*i)[col]);
        if (nt < 4) ++nt_count[nt];
      }
      vector<size_t>::iterator majority_nt = std::max_element(nt_count.begin(),
          nt_count.end());
      consensus += int2base(majority_nt - nt_count.begin());
    }
  }
  return consensus;
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
    bool MANUAL = false;
    bool NOT_TRAIN_TRUNC = false;
    double lambda = 0.0;
    string outf;

    OptionParser opt_parse(strip_path(argv[0]),
        "Construct profile-HMM from multiple alignment.",
        "<multiple alignment ALN file>");
    opt_parse.add_opt("out", 'o', "Output parameters file.", false, outf);
    opt_parse.add_opt("lambda", 'l',
        "Model length adjusting parameter. Default: 0", false, lambda);
    opt_parse.add_opt("heuristic", 'u', "Heuristic mode.", false, USE_HEU);
    opt_parse.add_opt("manual", 'm'
        , "Manual mode. The input must have a sequence named as consensus."
        , false, MANUAL);
    opt_parse.add_opt("no-trunc", 'n', "Do not train truncation probabilities."
        , false, NOT_TRAIN_TRUNC);
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
    if (USE_HEU && MANUAL) {
      cout << "Heuristic mode (-u) is not compatible with manual mode (-m)!.";
      return EXIT_SUCCESS;
    }
    const string inf = leftover_args.front();

    // 1. load alignment
    if (VERBOSE)
      cerr << "[LOADING ALIGNMENT]" << endl;
    vector<string> names, seqs;
    load_alignment(VERBOSE, inf, names, seqs);
    if (seqs.size() == 0) {
      cerr << "NO SEQUENCE WAS LOADED. PLEASE CHECK INPUT." << endl;
      return EXIT_FAILURE;
    }

    // 2. compute recursion to find marked (match) columns
    if (VERBOSE)
      cerr << "[FINDING OPTIMAL MATCH COLUMNS]" << endl;
    vector<bool> marked_cols(seqs.front().length(), false);
    matrix tran_prior, emis_prior;
    set_transition_prior(tran_prior);
    set_emission_prior(emis_prior);
    if (USE_HEU)
      find_marked_cols_heu(marked_cols, seqs);
    else if (MANUAL)
      find_marked_cols_manual(marked_cols, names, seqs);
    else
      find_marked_cols_map(VERBOSE,
          marked_cols, seqs, tran_prior, emis_prior, lambda);
    
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
    for (vector<string>::const_iterator i = seqs.begin();
        i < seqs.end(); ++i) {
      // get truncation by jumping to the first and the last non-gap character
      // for the marked columns
      vector<bool>::const_iterator curr_itr = first_marked;
      while(curr_itr < marked_cols.begin() && 
          (!*curr_itr || (*i)[curr_itr - marked_cols.begin()] == '-'))
        ++curr_itr;
      const vector<bool>::const_iterator first_nontruncated = curr_itr;
      curr_itr = last_marked; // last_marked is a reverse iterator
      while(curr_itr > first_nontruncated &&
          (!*curr_itr || (*i)[curr_itr - marked_cols.begin()] == '-'))
        --curr_itr;
      const vector<bool>::const_iterator last_nontruncated = curr_itr;
      // get counts for the first marked column associating with truncation
      state curr_state((*i)[first_nontruncated - marked_cols.begin()],
          1, true);
      transition[left_end.index(model_len)][curr_state.index(model_len)] =
        transition[left_end.index(model_len)][curr_state.index(model_len)] + 1;
      if (curr_state.baseint < 4)
        emission[curr_state.index(model_len)][curr_state.baseint] =
          emission[curr_state.index(model_len)][curr_state.baseint] + 1;

      vector<state> seq_state(1, curr_state);
      size_t col_idx = 1;
      for (vector<bool>::const_iterator j = first_nontruncated + 1;
          j <= last_nontruncated; ++j) {
        if (*j) ++col_idx;
        curr_state = state((*i)[j - marked_cols.begin()],
            col_idx, *j);
        if (curr_state.isvalid())
          seq_state.push_back(curr_state);
      }
      for (vector<state>::const_iterator j = seq_state.begin(), k = j+1;
          k < seq_state.end(); ++j,++k) {
        transition[(*j).index(model_len)][(*k).index(model_len)] = 
          transition[(*j).index(model_len)][(*k).index(model_len)] + 1;
        if ((*k).stateint < 2)
          if ((*k).baseint < 4)
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
      if (NOT_TRAIN_TRUNC) {
        row[3] = prior[3];
        normalize_vec_inplace(row);
      }
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
      if (NOT_TRAIN_TRUNC) {
        row.back() = prior.back();
        normalize_vec_inplace(row);
      }
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
    transition[0][left_end.index(model_len)] = log(0.30);
    transition[0][state(0ul, 0, 1).index(model_len)] = log(0.70);
    // D_1
    transition[left_end.index(model_len)] =
      normalize_vec(transition[left_end.index(model_len)], true);
    for (size_t i = 1; i <= model_len; ++i) {
      if (NOT_TRAIN_TRUNC)
        transition[left_end.index(model_len)][i] =
          get_5prime_truncation_prior(model_len, i);
      else
        transition[left_end.index(model_len)][i] = 
          log_sum_log(transition[left_end.index(model_len)][i],
              get_5prime_truncation_prior(model_len, i));
    }
    transition[left_end.index(model_len)] =
      normalize_vec(transition[left_end.index(model_len)], true);
    // I_0
    transition[state(0ul, 0, 1).index(model_len)]
      [state(0ul, 0, 1).index(model_len)] = log(0.80);
    transition[state(0ul, 0, 1).index(model_len)]
      [left_end.index(model_len)] = log(0.20);
    //transition[state(0ul, 0, 1).index(model_len)].back() = log(0.02);
    // D_L
    // the two lines below allow transition only to I_0 and E
    //transition[right_end.index(model_len)]
    //  [state(0ul, 0, 1).index(model_len)] = log(0.98);
    //transition[right_end.index(model_len)].back() = log(0.02);
    // the two lines below allow transitions only to I_0 and D_1
    transition[right_end.index(model_len)]
      [state(0ul, 0, 1).index(model_len)] = log(0.80);
    transition[right_end.index(model_len)]
      [left_end.index(model_len)] = log(0.20);

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

    //write_hmm_parameter(VERBOSE, out, transition, emission);
    ProfileHMM hmm(transition, emission);
    hmm.Print(out, false);
    if (VERBOSE)
      cerr << "[DONE]" << endl;

    // display alignment and marked columns for debug purposes
    if (DEBUG) {
      // clustalW keeps first 30 letters of the seq name, but
      // muscle keeps first 32. This is just to make sure the
      // debug output is compatible with both.
      const int name_width = (names.front().length() + 8) / 8 * 8;
      for (size_t i = 0; i < seqs.size(); ++i)
        cerr << names[i] << "\t" << seqs[i] << endl;
      //cerr << "MAP columns\t\t\t\t";
      size_t col_count = 0;
      for (vector<bool>::const_iterator i = marked_cols.begin();
          i < marked_cols.end(); ++i)
        if (*i) ++col_count;
      cerr << std::left << std::setw(name_width) << "MAP columns ("
        + std::to_string(col_count) + ")";
      for (size_t i = 0; i < marked_cols.size(); ++i) {
        if (marked_cols[i])
          cerr << "*";
        else
          cerr << " ";
      }
      cerr << endl;
      vector<bool> marked_heu(seqs.front().length(), false);
      find_marked_cols_heu(marked_heu, seqs);
      col_count = 0;
      for (vector<bool>::const_iterator i = marked_heu.begin();
          i < marked_heu.end(); ++i)
        if (*i) ++col_count;
      cerr << std::left << std::setw(name_width) << "Heuristic columns ("
        + std::to_string(col_count) + ")";
      for (size_t i = 0; i < marked_heu.size(); ++i) {
        if (marked_heu[i])
          cerr << "#";
        else
          cerr << " ";
      }
      cerr << endl;
      cerr << "Consensus: " << endl
        << get_consensus(marked_cols, seqs) << endl;
      //hmm.DebugOutput();
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
