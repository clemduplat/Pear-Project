#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "eigen-3.4.0/Eigen/Eigen"
#include "eigen-3.4.0/Eigen/Dense"
using namespace std;
using namespace Eigen;

MatrixXd read_file(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Error: '" + filename + "' is not open");
    }
    std::vector<std::vector<double>> data;
    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::vector<double> row;
        std::string cell;
        while (std::getline(ss, cell, ',')) {
            row.push_back(std::stod(cell));
        }
        data.push_back(row);
    }
    MatrixXd matrix_output(data.size(), data[0].size());
    for (size_t i = 0; i < data.size(); ++i) {
        for (size_t j = 0; j < data[i].size(); ++j) {
            matrix_output(i, j) = data[i][j];
        }
    }
    return matrix_output;
}