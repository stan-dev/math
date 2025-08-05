#pragma once
#include <stan/math/prim/fun/Eigen.hpp>   // or <Eigen/Dense>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace stan::math::test::laplace {

/**
 * Reads a CSV of doubles into a dynamically sized Eigen::Matrix.
 *
 * @param file_path  Path to a CSV file with R-style numeric matrix
 *                   (no header, no row names, comma-separated).
 * @return           Matrix(rows, cols) filled from the file.
 * @throws std::runtime_error on I/O error or inconsistent column counts.
 */
inline Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> read_matrix_csv(
    const std::string& file_path) {
  std::ifstream in(file_path);
  if (!in.is_open()) {
    throw std::runtime_error("Could not open CSV file: " + file_path);
  }

  std::vector<std::vector<double>> buffer;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;  // skip blank lines
    }
    std::vector<double> row;
    std::stringstream ss(line);
    std::string cell;
    while (std::getline(ss, cell, ',')) {
      try {
        row.push_back(std::stod(cell));
      } catch (const std::invalid_argument&) {
        throw std::runtime_error(
            "Non-numeric value in CSV at line: " + line);
      }
    }
    if (!buffer.empty() && row.size() != buffer[0].size()) {
      throw std::runtime_error(
          "Inconsistent column count in CSV: expected " +
          std::to_string(buffer[0].size()) +
          " but got " + std::to_string(row.size()));
    }
    buffer.push_back(std::move(row));
  }

  // If empty file, return a 0×0 matrix
  if (buffer.empty()) {
    return Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(0, 0);
  }

  const int rows = static_cast<int>(buffer.size());
  const int cols = static_cast<int>(buffer[0].size());
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mat(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      mat(i, j) = buffer[i][j];
    }
  }
  return mat;
}

}  // namespace stan::math::test::roaches
