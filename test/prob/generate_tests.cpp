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
#include <filesystem>
namespace fs = std::filesystem;

using std::endl;
using std::pair;
using std::string;
using std::stringstream;
using std::vector;

enum class ad_test_type { V, FD, FV, FFD, FFV, VV };

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

vector<string> lookup_argument(const string& argument, const ad_test_type& test_type) {
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
    if (test_type == ad_test_type::V) {
      push_args(args, "var");
    } else if (test_type == ad_test_type::FD) {
      push_args(args, "fvar<double>");
    } else if (test_type == ad_test_type::FV) {
      push_args(args, "fvar<var>");
    } else if (test_type == ad_test_type::FFD) {
      push_args(args, "fvar<fvar<double> >");
    } else if (test_type == ad_test_type::FFV) {
      push_args(args, "fvar<fvar<var> >");
    } else if (test_type == ad_test_type::VV) {
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

inline std::string trim_front_path(const std::string& path) {
  size_t pos = path.find("test/prob/");
  if (pos == std::string::npos) {
    return path;
  }
  return std::string(path.substr(pos));
}

void write_header_includes(std::ostream& out, const string& include) {
    out << "#include <gtest/gtest.h>" << endl;
    out << "#include <tuple>" << endl;
    out << "#include <test/prob/test_fixture_distr.hpp>" << endl;
    out << "#include <test/prob/test_fixture_cdf.hpp>" << endl;
    out << "#include <test/prob/test_fixture_cdf_log.hpp>" << endl;
    out << "#include <test/prob/test_fixture_ccdf_log.hpp>" << endl;
    out << endl;
}

void write_includes(std::ostream& out, const std::string& type_header, const string& include) {
  // TODO: Add arg header needed
    out << "#include <" << trim_front_path(type_header) << ">" << endl;
    out << "#include <" << trim_front_path(include) << ">" << endl;
    out << endl;
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
                                                const ad_test_type& test_type) {
  vector<string> argument_list = tokenize_arguments(arguments);
  vector<vector<string> > argument_sequence;
  for (size_t n = 0; n < argument_list.size(); n++)
    argument_sequence.push_back(lookup_argument(argument_list[n], test_type));
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

void write_types_typedef(std::ostream& out, string base, size_t& arg_seq_size,
                         vector<vector<string> > argument_sequence,
                         const size_t depth, const ad_test_type& test_type,
                         const std::string& scalar_type) {
  vector<string> args = argument_sequence.front();
  argument_sequence.erase(argument_sequence.begin());
  if (argument_sequence.size() > 0) {
    for (size_t n = 0; n < args.size(); n++)
      write_types_typedef(out, base + args[n] + ", ", arg_seq_size, argument_sequence,
                          depth, test_type, scalar_type);
  } else {
    string extra_args;
    for (size_t n = depth; n < 6; n++) {
      extra_args += ", empty";
    }
    for (size_t n = 0; n < args.size(); n++) {
      if (test_type == ad_test_type::V) {
        out << "typedef std::tuple<" << base << args[n] << extra_args;
        if (extra_args.size() == 0)
          out << " ";
        out << "> type_v_" << scalar_type << arg_seq_size << ";" << endl;
        arg_seq_size++;
      } else {
        if (check_all_double(base, args[n]) == false) {
          out << "typedef std::tuple<" << base << args[n] << extra_args;
          if (extra_args.size() == 0)
            out << " ";
          else if (test_type == ad_test_type::FD)
            out << "> type_fd_" << scalar_type << arg_seq_size << ";" << endl;
          else if (test_type == ad_test_type::FV)
            out << "> type_fv_" << scalar_type << arg_seq_size << ";" << endl;
          else if (test_type == ad_test_type::FFD)
            out << "> type_ffd_" << scalar_type << arg_seq_size << ";" << endl;
          else if (test_type == ad_test_type::FFV)
            out << "> type_ffv_" << scalar_type << arg_seq_size << ";" << endl;
          else if (test_type == ad_test_type::VV)
            out << "> type_vv_" << scalar_type << arg_seq_size << ";" << endl;
          arg_seq_size++;
        }
      }
    }
  }
}

size_t write_types(std::ostream& out,
                   const vector<vector<string> >& argument_sequence,
                   const ad_test_type& test_type, const std::string& scalar_types) {
  size_t arg_seq_size = 0;
  write_types_typedef(out, "", arg_seq_size, argument_sequence, argument_sequence.size(),
                      test_type, scalar_types);
    out << endl;
  return arg_seq_size;
}

void write_test(std::ostream& out, const string& test_name,
                const string& fixture_name, const ad_test_type& test_type,
                const std::string& scalar_type,
                const size_t test_start,
                const int& test_end) {
    std::string test_type_str;
    if (test_type == ad_test_type::V) {
      test_type_str = "_v_";
    } else if (test_type == ad_test_type::FD) {
      test_type_str = "_fd_";
    } else if (test_type == ad_test_type::FV) {
      test_type_str = "_fv_";
    } else if (test_type == ad_test_type::FFD) {
      test_type_str = "_ffd_";
    } else if (test_type == ad_test_type::FFV) {
      test_type_str = "_ffv_";
    } else if (test_type == ad_test_type::VV) {
      test_type_str = "_vv_";
    }
    for (int n = test_start; n < test_end; ++n) {
      out << "typedef std::tuple<" << test_name << ", type" << test_type_str << scalar_type << n << "> "
            << test_name << test_type_str << scalar_type << n << ";" << endl;
      out << endl;
      out << "INSTANTIATE_TYPED_TEST_SUITE_P(" << test_name << test_type_str << scalar_type << n
          << ", " << fixture_name << ", " << test_name << test_type_str << scalar_type << n << ");"
          << endl;
      out << endl;
    }
}

void write_test_cases(std::ostream& out, const string& file,
                      const vector<vector<string> >& argument_sequence,
                      const ad_test_type& test_type, const std::string& scalar_type, const int start_test_num, const int end_test_num) {
  pair<string, string> name = read_test_name_from_file(file);
  string test_name = name.first;
  string fixture_name = name.second;
  write_test(out, test_name, fixture_name, test_type, scalar_type, start_test_num, end_test_num);
}

int count_lines(const string& file) {
  std::ifstream in(file);
  if (in.is_open()) {
    int count = 0;
    std::string line;
    while (std::getline(in, line)) {  // Loop through each line in the file
      count++;  // Incrementing line count for each line read
    }
    in.close();
    return count;
  } else {
    std::cout << "BAD!!!\n";
    return -1;
  }
}

int create_files(const int& argc, const std::filesystem::path& path,
const std::filesystem::path& parent_path,
                 const ad_test_type& test_type,
                 const int& start, const int& N_TESTS) {
  if (argc != 3)
    return -1;
  string in_suffix = "_test.hpp";
  string in_name = path;
  std::cout << in_name << std::endl;
  // ???
  size_t last_in_suffix
      = in_name.find_last_of(in_suffix) + 1 - in_suffix.length();
  string out_name_base = in_name.substr(0, last_in_suffix);
  std::cout << "cleaned name: " << in_name << std::endl;
  std::cout << "out name: " << out_name_base << std::endl;
  string file = read_file(in_name);

  string arguments = read_arguments_from_file(file);
  vector<vector<string> > argument_sequence
      = build_argument_sequence(arguments, test_type);
    std::cout << "Arg Seq: \n";
    for (auto&& args : argument_sequence) {
      for (auto&& arg : args) {
        std::cout << arg << ", ";
      }
      std::cout << "\n";
    }
  // We always have 8 lines in the header for includes then an EOL
  //const double BATCHES = N_TESTS > 0 ? num_tests / N_TESTS : -N_TESTS;
  stringstream arg_header_tmp;
  arg_header_tmp << std::string(parent_path) << "/args/arg_generated_";
  if (test_type == ad_test_type::V) {
    arg_header_tmp << "v_";
  } else if (test_type == ad_test_type::FD) {
    arg_header_tmp << "fd_";
  } else if (test_type == ad_test_type::FV) {
    arg_header_tmp << "fv_";
  } else if (test_type == ad_test_type::FFD) {
    arg_header_tmp << "ffd_";
  } else if (test_type == ad_test_type::FFV) {
    arg_header_tmp << "ffv_";
  } else if (test_type == ad_test_type::VV) {
    arg_header_tmp << "vv_";
  }
  stringstream arg_scalar_types_tmp;
  for (auto&& args : argument_sequence) {
    if (args[0] == "int") {
      arg_scalar_types_tmp << "int_";
    } else if (args[0] == "double") {
      arg_scalar_types_tmp << "real_";
    }
  }
  std::string string_arg_scalar_types_tmp = arg_scalar_types_tmp.str();
  arg_header_tmp << string_arg_scalar_types_tmp;
  arg_header_tmp << "pch.hpp";
  std::string arg_header(arg_header_tmp.str());
  std::cout << "arg header name: " << arg_header << std::endl;
  if (false) {
//    if (fs::exists(arg_header)) {
    std::cout << "header file exists: " << arg_header << std::endl;
  } else {
    std::cout << "header file does not exist: " << arg_header << std::endl;
    std::ofstream arg_header_stream(arg_header.c_str());
    write_header_includes(arg_header_stream, in_name);
    write_types(arg_header_stream, argument_sequence, test_type, string_arg_scalar_types_tmp);
    //write_header_typedefs(test_stream, file, argument_sequence, test_type);
  }

  int num_tests = count_lines(arg_header) - 9;
  std::cout << "num tests: " << num_tests << "\n";
  int n = 0;
  int file_count = 0;
  // NOTE: DOES NOT WORK FOR num_tests < N_TESTS
  for (; n < num_tests; n += N_TESTS, file_count++) {
    std::cout << "==========\n";
    stringstream out_name;
    out_name << out_name_base;
    out_name << "_" << std::setw(5) << std::setfill('0') << file_count;
    if (test_type == ad_test_type::V) {
      out_name << "_generated_v_test.cpp";
    } else if (test_type == ad_test_type::FD) {
      out_name << "_generated_fd_test.cpp";
    } else if (test_type == ad_test_type::FV) {
      out_name << "_generated_fv_test.cpp";
    } else if (test_type == ad_test_type::FFD) {
      out_name << "_generated_ffd_test.cpp";
    } else if (test_type == ad_test_type::FFV) {
      out_name << "_generated_ffv_test.cpp";
    } else if (test_type == ad_test_type::VV) {
      out_name << "_generated_vv_test.cpp";
    }
    std::string tmp(out_name.str());
    std::cout << "subtest name: " << tmp << std::endl;
    std::ofstream test_stream(tmp.c_str());
    write_includes(test_stream, arg_header, in_name);
    write_test_cases(test_stream, file, argument_sequence, test_type, string_arg_scalar_types_tmp, n, std::min(n + N_TESTS, num_tests));

    /*
    if (N_TESTS > 0) {
    } else if (num_tests > 0) {
      write_test_cases(test_stream, file, argument_sequence, test_type,
                      ceil(num_tests / BATCHES));
    }
    std::cout << "==========\n";
  */
  }
  return 0;
}

template <typename T1, typename T2, typename T3>
void recurse_directories(T1 argc, T2&& path, T3&& parent_path, int N_TESTS = 0) {
  for (const auto & entry : fs::directory_iterator(path)) {
    auto inner_path = entry.path();
//      std::cout << "entry1: " << inner_path << std::endl;
    if (fs::is_directory(inner_path)) {
      recurse_directories(argc, inner_path, parent_path, N_TESTS);
    } else {
      if (inner_path.extension() != ".hpp") {
        if (std::string(inner_path.filename()).find("test") != std::string::npos) {
          continue;
        }
      }
      if (std::string(inner_path.filename()).find("generated") != std::string::npos) {
        continue;
      }
      create_files(argc, inner_path, parent_path, ad_test_type::V, -1, N_TESTS);  // create var tests
      create_files(argc, inner_path, parent_path, ad_test_type::FFV, -1, N_TESTS);  // create ffv tests
      create_files(argc, inner_path, parent_path, ad_test_type::VV, -1, N_TESTS);  // create varmat tests
      #ifdef STAN_PROB_TEST_ALL
      create_files(argc, inner_path, parent_path, ad_test_type::FD, -1, N_TESTS);  // create fd tests
      create_files(argc, inner_path, parent_path, ad_test_type::FV, -1, N_TESTS);  // create fv tests
      create_files(argc, inner_path, parent_path, ad_test_type::FFD, -1, N_TESTS);  // create ffd tests
      #endif
    }
  }
}

void make_precompiled_header(const std::string& path_name) {
  std::vector<std::string> arg_headers;
  for (auto&& arg_file : fs::directory_iterator(path_name + "/args")) {
    arg_headers.push_back(arg_file.path());
  }
  std::ofstream pch_stream(path_name + "/generated_pch.hpp");
  pch_stream << R"(#include <stan/math/mix.hpp>
#include <stan/math/rev.hpp>
#include <test/prob/test_fixture_ccdf_log.hpp>
#include <test/prob/test_fixture_cdf_log.hpp>
#include <test/prob/test_fixture_cdf.hpp>
#include <test/prob/test_fixture_distr.hpp>
#include <test/prob/utility.hpp>
#include <gtest/gtest.h>
#include <stdexcept>
#include <tuple>
#include <type_traits>
)";

  for (auto&& arg_header : arg_headers) {
    pch_stream << "#include <" << trim_front_path(arg_header) << ">" << endl;
  }
  pch_stream << endl;



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
  std::cout << "RUNNING PROB BUILD STEP" << std::endl;
  int N_TESTS = atoi(argv[2]);
  for (const auto & entry : fs::directory_iterator(argv[1])) {
          std::cout << "entry0: " << entry.path() << std::endl;
    if (fs::is_directory(entry.path())) {
      recurse_directories(argc, entry.path(), argv[1], N_TESTS);
    }
  }
  make_precompiled_header(argv[1]);
  return 0;
}
