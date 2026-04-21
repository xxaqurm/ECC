#include <vector>
#include <fstream>
#include <filesystem>
#include <ctime>
#include <random>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace std;

using Matrix = vector<vector<double>>;

const double MATRIX_MIN_VALUE = -1000.0;
const double MATRIX_MAX_VALUE = 1000.0;

void createMatrix(Matrix &matrix, int matrixSize) {
    if (matrixSize <= 0) {
        throw invalid_argument("Matrix size must be positive");
    }
    
    random_device seed;
    mt19937 gen(seed());
    uniform_real_distribution<double> dist(MATRIX_MIN_VALUE, MATRIX_MAX_VALUE);

    for (int i = 0; i < matrixSize; i++) {
        for (int j = 0; j < matrixSize; j++) {
            matrix[i][j] = dist(gen);
        }
    }
}

void displayMatrix(const Matrix &matrix) {
    if (matrix.empty()) {
        cout << "[Empty matrix]" << endl;
        return;
    }
    
    for (size_t i = 0; i < matrix.size(); i++) {
        for (size_t j = 0; j < matrix[i].size(); j++) {
            cout << setw(12) << fixed << setprecision(2) << matrix[i][j];
        }
        cout << endl;
    }
}

Matrix dgemmBlas(const Matrix &m1, const Matrix &m2) {
    size_t n = m1.size();
    Matrix result(n, vector<double>(n, 0.0));
    
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            double sum = 0.0;
            for (size_t k = 0; k < n; k++) {
                sum += m1[i][k] * m2[k][j];
            }
            result[i][j] = sum;
        }
    }
    return result;
}

Matrix dgemmOpt1(const Matrix &m1, const Matrix &m2) {
    size_t n = m1.size();
    Matrix result(n, vector<double>(n, 0.0));

    for (size_t i = 0; i < n; i++) {
        for (size_t r = 0; r < n; r++) {
            double cachedValue = m1[i][r];
            for (size_t j = 0; j < n; j++) {
                result[i][j] += cachedValue * m2[r][j];
            }
        }
    }

    return result;
}

Matrix dgemmOpt2(const Matrix &m1, const Matrix &m2, int blockSize) {
    size_t n = m1.size();
    Matrix result(n, vector<double>(n, 0.0));
    
    // TODO: реализовать блочное умножение матриц
    // Пока возвращаем результат базового алгоритма
    return dgemmBlas(m1, m2);
}

void checkDgemmxPerformance(int maxValue, int step, const string &name) {
    try {
        filesystem::create_directories("./data/plots");
        filesystem::create_directories("./data/stats");
        
        ofstream dataFile("./data/stats/" + name + ".dat");
        if (!dataFile.is_open()) {
            throw runtime_error("Failed to open file for writing: ./data/stats/" + name + ".dat");
        }
        
        for (int matrixSize = 1; matrixSize < maxValue; matrixSize += step) {
            Matrix firstMatrix(matrixSize, vector<double>(matrixSize));
            Matrix secondMatrix(matrixSize, vector<double>(matrixSize));

            createMatrix(firstMatrix, matrixSize);
            createMatrix(secondMatrix, matrixSize);

            clock_t start = clock();
            Matrix result;
            if (name == "dgemmBlas") {
                result = dgemmBlas(firstMatrix, secondMatrix);
            } else if (name == "dgemmOpt1") {
                result = dgemmOpt1(firstMatrix, secondMatrix);
            } else if (name == "dgemmOpt2") {
                result = dgemmOpt2(firstMatrix, secondMatrix, 32);
            }
            clock_t end = clock();
            double duration = double(end - start) / CLOCKS_PER_SEC;

            dataFile << matrixSize << " " << duration << "\n";
        }
        dataFile.close();

        // Создание графика с помощью gnuplot
        string command = "gnuplot -e \"set terminal png; "
                          "set output './data/plots/" + name + ".png'; "
                          "set xlabel 'Matrix Size (N)'; "
                          "set ylabel 'Execution Time (s)'; "
                          "set grid; "
                          "plot './data/stats/" + name + ".dat' with linespoints title '" + name + "'\"";
        
        int ret = system(command.c_str());
        if (ret != 0) {
            cerr << "Warning: gnuplot command failed with code " << ret << endl;
        }
    } catch (const exception &e) {
        cerr << "Error in checkDgemmxPerformance: " << e.what() << endl;
    }
}