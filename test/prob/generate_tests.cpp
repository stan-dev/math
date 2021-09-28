#include <ostream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <utility>
#include <vector>
#include <iomanip>
#include <boost/algorithm/string.hpp>
#include <cmath>
#include <cstdlib>

using std::endl;
using std::pair;
using std::string;
using std::stringstream;
using std::vector;

#ifdef STAN_TEST_ROW_VECTORS
int ROW_VECTORS = 1;
#else
int ROW_VECTORS = 0;
#endif

void push_args(vector<string>& args, const string& type) {
  if (type.compare("varmat") == 0) {
    args.push_back("var");
    args.push_back("std::vector<var>");
    args.push_back(
        "stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>");
    if (ROW_VECTORS == 1)
      args.push_back(
          "stan::math::var_value<Eigen::Matrix<double, 1, Eigen::Dynamic>>");
  } else {
    args.push_back(type);
    args.push_back("std::vector<" + type + ">");
    args.push_back("Eigen::Matrix<" + type + ", Eigen::Dynamic, 1>");
    if (ROW_VECTORS == 1)
      args.push_back("Eigen::Matrix<" + type + ", 1, Eigen::Dynamic>");
  }
}

vector<string> lookup_argument(const string& argument, const int& ind) {
  using boost::iequals;
  vector<string> args;
  if (iequals(argument, "int")) {
    args.push_back("int");
  } else if (iequals(argument, "ints")) {
    push_args(args, "int");
  } else if (iequals(argument, "double")) {
    args.push_back("double");
    args.push_back("var");
  } else if (iequals(argument, "doubles")) {
    push_args(args, "double");
    if (ind == 1) {
      push_args(args, "var");
    } else if (ind == 2) {
      push_args(args, "fvar<double>");
    } else if (ind == 3) {
      push_args(args, "fvar<var>");
    } else if (ind == 4) {
      push_args(args, "fvar<fvar<double> >");
    } else if (ind == 5) {
      push_args(args, "fvar<fvar<var> >");
    } else if (ind == 6) {
      push_args(args, "varmat");
    }
  }
  return args;
}

std::ostream& operator<<(std::ostream& o, pair<string, string>& p) {
  o << "<" << p.first << ", " << p.second << ">" << endl;
  return o;
}

template <class T>
std::ostream& operator<<(std::ostream& o, vector<T>& vec) {
  o << "vector size: " << vec.size() << endl;
  for (size_t n = 0; n < vec.size(); n++) {
    o << "  \'" << vec[n] << "\'" << endl;
  }
  return o;
}

void write_includes(vector<std::ostream*>& outs, const string& include) {
  for (size_t n = 0; n < outs.size(); n++) {
    std::ostream* out = outs[n];
    *out << "#include <gtest/gtest.h>" << endl;
    *out << "#include <tuple>" << endl;
    *out << "#include <test/prob/test_fixture_distr.hpp>" << endl;
    *out << "#include <test/prob/test_fixture_cdf.hpp>" << endl;
    *out << "#include <test/prob/test_fixture_cdf_log.hpp>" << endl;
    *out << "#include <test/prob/test_fixture_ccdf_log.hpp>" << endl;
    *out << "#include <" << include << ">" << endl;
    *out << endl;
  }
}

vector<string> tokenize_arguments(const string& arguments) {
  vector<string> tokens;
  string delimiters = ", ";
  string args_only_string = arguments.substr(arguments.find(":") + 1);
  boost::algorithm::trim(args_only_string);
  boost::algorithm::split(tokens, args_only_string,
                          boost::is_any_of(delimiters),
                          boost::token_compress_on);
  return tokens;
}

size_t size(const vector<vector<string> >& sequences) {
  if (sequences.size() == 0)
    return 0;
  size_t N = 1;
  for (size_t n = 0; n < sequences.size(); n++)
    N *= sequences[n].size();
  return N;
}

bool is_argument_list(const string& line) {
  size_t comment = line.find("// ");
  if (comment == string::npos)
    return false;
  size_t keyword = line.find("Arguments:", comment + 1);
  if (keyword == string::npos)
    return false;
  return true;
}

string read_file(const string& in_name) {
  std::ifstream in(in_name.c_str());

  string file;
  in.seekg(0, std::ios::end);
  file.resize(in.tellg());
  in.seekg(0, std::ios::beg);
  in.read(&file[0], file.size());
  in.close();

  return file;
}

string read_arguments_from_file(const string& file) {
  string arguments;
  std::istringstream in(file);
  if (!in.eof()) {
    std::getline(in, arguments);
    while (in.good() && !is_argument_list(arguments)) {
      std::getline(in, arguments);
    }
  }
  if (!is_argument_list(arguments))
    arguments = "";
  return arguments;
}

pair<string, string> read_test_name_from_file(const string& file) {
  pair<string, string> name;

  size_t pos = 0;
  string class_keyword = "class ";
  string public_keyword = "public ";
  pos = file.find(class_keyword, pos);
  if (pos < file.size()) {
    pos += class_keyword.size();
    size_t pos2 = file.find(":", pos);
    string test_name = file.substr(pos, pos2 - pos);
    pos = file.find(public_keyword, pos) + public_keyword.size();
    pos2 = file.find("{", pos);
    string fixture_name = file.substr(pos, pos2 - pos);
    pos = file.find("};", pos) + 2;
    boost::algorithm::trim(test_name);
    boost::algorithm::trim(fixture_name);

    if (fixture_name.find("Test") != string::npos) {
      fixture_name += "Fixture";
      name = pair<string, string>(test_name, fixture_name);
    }
  }
  return name;
}

vector<vector<string> > build_argument_sequence(const string& arguments,
                                                const int& ind) {
  vector<string> argument_list = tokenize_arguments(arguments);
  vector<vector<string> > argument_sequence;
  for (size_t n = 0; n < argument_list.size(); n++)
    argument_sequence.push_back(lookup_argument(argument_list[n], ind));
  return argument_sequence;
}

bool check_all_double(string base, string arg) {
  string arguments = base + arg;
  bool result = true;
  bool temp;

  vector<string> tokens;
  string delimiters = ", ";
  boost::algorithm::trim(arguments);
  boost::algorithm::split(tokens, arguments, boost::is_any_of(delimiters),
                          boost::token_compress_on);

  for (size_t i = 0; i < tokens.size(); i++) {
    if (tokens[i] == "1" || tokens[i] == "Eigen::Dynamic>" || tokens[i] == "1>"
        || tokens[i] == "Eigen::Dynamic")
      result = result && true;
    else {
      temp = (tokens[i] == "double") || (tokens[i] == "std::vector<double>")
             || (tokens[i] == "Eigen::Matrix<double") || (tokens[i] == "int")
             || (tokens[i] == "std::vector<int>")
             || (tokens[i] == "Eigen::Matrix<int");
      result = result && temp;
    }
  }
  return result;
}

int num_doubles(string arguments) {
  vector<string> tokens;
  string delimiters = ", ";
  boost::algorithm::trim(arguments);
  boost::algorithm::split(tokens, arguments, boost::is_any_of(delimiters),
                          boost::token_compress_on);

  int num = 0;
  for (size_t i = 0; i < tokens.size(); i++) {
    if (tokens[i] == "Doubles")
      ++num;
  }
  return num;
}
int num_ints(string arguments) {
  vector<string> tokens;
  string delimiters = ", ";
  boost::algorithm::trim(arguments);
  boost::algorithm::split(tokens, arguments, boost::is_any_of(delimiters),
                          boost::token_compress_on);

  int num = 0;
  for (size_t i = 0; i < tokens.size(); i++) {
    if (tokens[i] == "Ints")
      ++num;
  }
  return num;
}

void write_types_typedef(vector<std::ostream*>& outs, string base, size_t& N,
                         vector<vector<string> > argument_sequence,
                         const size_t depth, const int& index,
                         const int& N_TESTS) {
  vector<string> args = argument_sequence.front();
  argument_sequence.erase(argument_sequence.begin());
  if (argument_sequence.size() > 0) {
    for (size_t n = 0; n < args.size(); n++)
      write_types_typedef(outs, base + args[n] + ", ", N, argument_sequence,
                          depth, index, N_TESTS);
  } else {
    string extra_args;
    for (size_t n = depth; n < 6; n++) {
      extra_args += ", empty";
    }
    for (size_t n = 0; n < args.size(); n++) {
      std::ostream* out = outs[int(N / N_TESTS)];
      if (index == 1) {
        *out << "typedef std::tuple<" << base << args[n] << extra_args;
        if (extra_args.size() == 0)
          *out << " ";
        *out << "> type_v_" << N << ";" << endl;
        N++;
      } else {
        if (check_all_double(base, args[n]) == false) {
          *out << "typedef std::tuple<" << base << args[n] << extra_args;
          if (extra_args.size() == 0)
            *out << " ";
          else if (index == 2)
            *out << "> type_fd_" << N << ";" << endl;
          else if (index == 3)
            *out << "> type_fv_" << N << ";" << endl;
          else if (index == 4)
            *out << "> type_ffd_" << N << ";" << endl;
          else if (index == 5)
            *out << "> type_ffv_" << N << ";" << endl;
          else if (index == 6)
            *out << "> type_vv_" << N << ";" << endl;
          N++;
        }
      }
    }
  }
}

size_t write_types(vector<std::ostream*>& outs,
                   const vector<vector<string> >& argument_sequence,
                   const int& index, const int& N_TESTS) {
  size_t N = 0;
  write_types_typedef(outs, "", N, argument_sequence, argument_sequence.size(),
                      index, N_TESTS);
  for (size_t n = 0; n < outs.size(); n++)
    *outs[n] << endl;
  return N;
}

void write_test(vector<std::ostream*>& outs, const string& test_name,
                const string& fixture_name, const size_t N, const int& index,
                const int& N_TESTS) {
  for (size_t n = 0; n < N; n++) {
    std::ostream* out = outs[int(n / N_TESTS)];
    if (index == 1)
      *out << "typedef std::tuple<" << test_name << ", type_v_" << n << "> "
           << test_name << "_v_" << n << ";" << endl;
    else if (index == 2)
      *out << "typedef std::tuple<" << test_name << ", type_fd_" << n << "> "
           << test_name << "_fd_" << n << ";" << endl;
    else if (index == 3)
      *out << "typedef std::tuple<" << test_name << ", type_fv_" << n << "> "
           << test_name << "_fv_" << n << ";" << endl;
    else if (index == 4)
      *out << "typedef std::tuple<" << test_name << ", type_ffd_" << n << "> "
           << test_name << "_ffd_" << n << ";" << endl;
    else if (index == 5)
      *out << "typedef std::tuple<" << test_name << ", type_ffv_" << n << "> "
           << test_name << "_ffv_" << n << ";" << endl;
    else if (index == 6)
      *out << "typedef std::tuple<" << test_name << ", type_vv_" << n << "> "
           << test_name << "_vv_" << n << ";" << endl;
  }
  for (size_t i = 0; i < outs.size(); i++) {
    *outs[i] << endl;
  }
  for (size_t n = 0; n < N; n++) {
    std::ostream* out = outs[int(n / N_TESTS)];
    if (index == 1)
      *out << "INSTANTIATE_TYPED_TEST_SUITE_P(" << test_name << "_v_" << n
           << ", " << fixture_name << ", " << test_name << "_v_" << n << ");"
           << endl;
    else if (index == 2)
      *out << "INSTANTIATE_TYPED_TEST_SUITE_P(" << test_name << "_fd_" << n
           << ", " << fixture_name << ", " << test_name << "_fd_" << n << ");"
           << endl;
    else if (index == 3)
      *out << "INSTANTIATE_TYPED_TEST_SUITE_P(" << test_name << "_fv_" << n
           << ", " << fixture_name << ", " << test_name << "_fv_" << n << ");"
           << endl;
    else if (index == 4)
      *out << "INSTANTIATE_TYPED_TEST_SUITE_P(" << test_name << "_ffd_" << n
           << ", " << fixture_name << ", " << test_name << "_ffd_" << n << ");"
           << endl;
    else if (index == 5)
      *out << "INSTANTIATE_TYPED_TEST_SUITE_P(" << test_name << "_ffv_" << n
           << ", " << fixture_name << ", " << test_name << "_ffv_" << n << ");"
           << endl;
    else if (index == 6)
      *out << "INSTANTIATE_TYPED_TEST_SUITE_P(" << test_name << "_vv_" << n
           << ", " << fixture_name << ", " << test_name << "_vv_" << n << ");"
           << endl;
  }
  for (size_t i = 0; i < outs.size(); i++) {
    *outs[i] << endl;
  }
}

void write_test_cases(vector<std::ostream*>& outs, const string& file,
                      const vector<vector<string> >& argument_sequence,
                      const int& index, const int& N_TESTS) {
  pair<string, string> name = read_test_name_from_file(file);
  string test_name = name.first;
  string fixture_name = name.second;

  size_t num_tests = write_types(outs, argument_sequence, index, N_TESTS);
  write_test(outs, test_name, fixture_name, num_tests, index, N_TESTS);
}

int create_files(const int& argc, const char* argv[], const int& index,
                 const int& start, const int& N_TESTS) {
  if (argc != 3)
    return -1;
  string in_suffix = "_test.hpp";

  string in_name = argv[1];

  size_t last_in_suffix
      = in_name.find_last_of(in_suffix) + 1 - in_suffix.length();
  string out_name_base = in_name.substr(0, last_in_suffix);

  string file = read_file(in_name);

  string arguments = read_arguments_from_file(file);
  vector<vector<string> > argument_sequence
      = build_argument_sequence(arguments, index);

  int num_tests;
  if (index == 1)
    num_tests = size(argument_sequence);
  else
    num_tests = size(argument_sequence)
                - std::pow(3 + ROW_VECTORS, num_ints(arguments))
                      * std::pow(3 + ROW_VECTORS, num_doubles(arguments));

  vector<std::ostream*> outs;
  const double BATCHES = N_TESTS > 0 ? num_tests / N_TESTS : -N_TESTS;
  for (int n = start + 1; n < start + 1 + BATCHES + 1; n++) {
    stringstream out_name;
    out_name << out_name_base;
    out_name << "_" << std::setw(5) << std::setfill('0') << n;
    if (index == 1)
      out_name << "_generated_v_test.cpp";
    else if (index == 2)
      out_name << "_generated_fd_test.cpp";
    else if (index == 3)
      out_name << "_generated_fv_test.cpp";
    else if (index == 4)
      out_name << "_generated_ffd_test.cpp";
    else if (index == 5)
      out_name << "_generated_ffv_test.cpp";
    else if (index == 6)
      out_name << "_generated_vv_test.cpp";
    std::string tmp(out_name.str());
    outs.push_back(new std::ofstream(tmp.c_str()));
  }

  write_includes(outs, in_name);
  if (N_TESTS > 0)
    write_test_cases(outs, file, argument_sequence, index, N_TESTS);
  else if (num_tests > 0)
    write_test_cases(outs, file, argument_sequence, index,
                     ceil(num_tests / BATCHES));

  for (size_t n = 0; n < outs.size(); n++) {
    static_cast<std::ofstream*>(outs[n])->close();
    delete (outs[n]);
  }
  outs.clear();
  return start + BATCHES;
}

/**
 * Generate test cases.
 *
 * @param argc Number of arguments
 * @param argv Arguments. Should contain one argument with a filename and
 * the number of tests per file or if non-positive, the number of files - 1
 *
 * @return 0 for success, negative number otherwise.
 */
int main(int argc, const char* argv[]) {
  int N_TESTS = atoi(argv[2]);

  create_files(argc, argv, 1, -1, N_TESTS);  // create var tests
  create_files(argc, argv, 5, -1, N_TESTS);  // create ffv tests
  create_files(argc, argv, 6, -1, N_TESTS);  // create varmat tests
#ifdef STAN_PROB_TEST_ALL
  create_files(argc, argv, 2, -1, N_TESTS);  // create fd tests
  create_files(argc, argv, 3, -1, N_TESTS);  // create fv tests
  create_files(argc, argv, 4, -1, N_TESTS);  // create ffd tests
#endif

  return 0;
}
