#include <vector>
#include <fstream>
#include <filesystem>
#include <ctime>
#include <random>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <immintrin.h>
#include <cmath>

using namespace std;

using Matrix = vector<vector<double>>;
using LinearMatrix = vector<double>;

const double MATRIX_MIN_VALUE = -1000.0;
const double MATRIX_MAX_VALUE = 1000.0;

void createMatrix(Matrix &matrix, int matrixSize) {
    if (matrixSize < 0) {
        throw invalid_argument("Matrix size must be positive");
    }
    
    random_device seed;
    mt19937 gen(seed());
    uniform_real_distribution<double> dist(MATRIX_MIN_VALUE, MATRIX_MAX_VALUE);

    for (size_t i = 0; i < matrixSize; i++) {
        for (size_t j = 0; j < matrixSize; j++) {
            matrix[i][j] = dist(gen);
        }
    }
}

void createLinearMatrix(LinearMatrix &matrix, int matrixSize) {
    if (matrixSize < 0) {
        throw invalid_argument("Matrix size must be positive");
    }

    random_device seed;
    mt19937 gen(seed());
    uniform_real_distribution<double> dist(MATRIX_MIN_VALUE, MATRIX_MAX_VALUE);

    for (size_t i = 0; i < matrixSize * matrixSize; i++) {
        matrix[i] = dist(gen);
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

void displayLinearMatrix(const LinearMatrix &matrix, int matrixSize) {
    if (matrix.empty()) {
        cout << "[Empty matrix]" << endl;
        return;
    }

    for (size_t i = 0; i < matrixSize; i++) {
        for (size_t j = i * matrixSize; j < (i + 1) * matrixSize; j++) {
            cout << setw(12) << fixed << setprecision(2) << matrix[j];
        }
        cout << endl;
    }
}

Matrix dgemmBlas(const Matrix &m1, const Matrix &m2) {
    // Последовательное умножение двух квадратных матриц
    size_t n = m1.size();
    Matrix result(n, vector<double>(n, 0.0));
    
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            for (size_t k =  0; k < n; k++) {
                result[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }
    return result;
}

Matrix dgemmOpt1(const Matrix &m1, const Matrix &m2) {
    // Оптимизация доступа к памяти засчет построчного перебора элементов
    size_t n = m1.size();
    Matrix result(n, vector<double>(n, 0.0));

    for (size_t i = 0; i < n; i++) {
        for (size_t k = 0; k < n; k++) {
            for (size_t j = 0; j < n; j++) {
                result[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }
    return result;
}

Matrix dgemmOpt2(const Matrix &m1, const Matrix &m2, int blockSize) {
    // Оптимизация доступа к памяти за счет блочного перебора
    size_t n = m1.size();
    Matrix result(n, vector<double>(n, 0.0));

    for (size_t bi = 0; bi < n; bi += blockSize) {
        for (size_t bk = 0; bk < n; bk += blockSize) {
            for (size_t bj = 0; bj < n; bj += blockSize) {
                for (size_t i = bi; i < bi + blockSize && i < n; i++) {
                    for (size_t k = bk; k < bk + blockSize && k < n; k++) {
                        for (size_t j = bj; j < bj + blockSize && j < n; j++) {
                            result[i][j] += m1[i][k] * m2[k][j];
                        }
                    }
                }
            }
        }
    }
    return result;
}

LinearMatrix dgemmOpt3(const LinearMatrix &m1, const LinearMatrix &m2, int n, int blockSize) {
    LinearMatrix result(n * n, 0.0);

    LinearMatrix m2_T(n * n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            m2_T[j * n + i] = m2[i * n + j];
        }
    }

    for (int bi = 0; bi < n; bi += blockSize) {
        for (int bj = 0; bj < n; bj += blockSize) {
            for (int bk = 0; bk < n; bk += blockSize) {

                int i_end = std::min(bi + blockSize, n);
                int j_end = std::min(bj + blockSize, n);
                int k_end = std::min(bk + blockSize, n);

                for (int i = bi; i < i_end; i++) {

                    int j = bj;

                    for (; j + 3 < j_end; j += 4) {

                        __m256d c = _mm256_loadu_pd(&result[i * n + j]);

                        for (int k = bk; k < k_end; k++) {

                            __m256d a = _mm256_set1_pd(m1[i * n + k]);
                            __m256d b = _mm256_loadu_pd(&m2_T[j * n + k]);

                            c = _mm256_fmadd_pd(a, b, c);
                        }

                        _mm256_storeu_pd(&result[i * n + j], c);
                    }

                    for (; j < j_end; j++) {
                        for (int k = bk; k < k_end; k++) {
                            result[i * n + j] += m1[i * n + k] * m2[k * n + j];
                        }
                    }
                }
            }
        }
    }

    return result;
}

void checkDgemmxPerformance(int maxValue, int step, const string &name) {
    try {
        filesystem::create_directories("./data/plots");
        filesystem::create_directories("./data/stats");
        
        ofstream dataFile("./data/stats/" + name + ".dat");
        if (!dataFile.is_open()) {
            throw runtime_error("Failed to open file for writing: ./data/stats/" + name + ".dat");
        }
        
        for (int matrixSize = 0; matrixSize <= maxValue; matrixSize += step) {
            double duration;
            if (name == "dgemmBlas" || name == "dgemmOpt1" || name == "dgemmOpt2") {
                Matrix firstMatrix(matrixSize, vector<double>(matrixSize));
                Matrix secondMatrix(matrixSize, vector<double>(matrixSize));

                createMatrix(firstMatrix, matrixSize);
                createMatrix(secondMatrix, matrixSize);

                Matrix result;
                clock_t start = clock();
                if (name == "dgemmBlas") {
                result = dgemmBlas(firstMatrix, secondMatrix);
                } else if (name == "dgemmOpt1") {
                    result = dgemmOpt1(firstMatrix, secondMatrix);
                } else if (name == "dgemmOpt2") {
                    result = dgemmOpt2(firstMatrix, secondMatrix, 32);
                }
                clock_t end = clock();
                duration = double(end - start) / CLOCKS_PER_SEC;
            } else if (name == "dgemmOpt3") {
                LinearMatrix firstLinearMatrix(matrixSize * matrixSize);
                LinearMatrix secondLinearMatrix(matrixSize * matrixSize);
                
                createLinearMatrix(firstLinearMatrix, matrixSize);
                createLinearMatrix(secondLinearMatrix, matrixSize);
                
                LinearMatrix linearResult;
                clock_t start = clock();
                linearResult = dgemmOpt3(firstLinearMatrix, secondLinearMatrix, matrixSize, 32);
                clock_t end = clock();
                duration = double(end - start) / CLOCKS_PER_SEC;
            }
            dataFile << matrixSize << " " << duration << "\n";
        }
        dataFile.close();

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