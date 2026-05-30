#include <vector>
#include <fstream>
#include <filesystem>
#include <ctime>
#include <random>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <cmath>

#include <immintrin.h>
#include <pthread.h>
#include <omp.h>
#include <tbb/tbb.h>

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
    // Оптимизация за счет векторизации кода
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

                int i_end = min(bi + blockSize, n);
                int j_end = min(bj + blockSize, n);
                int k_end = min(bk + blockSize, n);

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

struct ThreadData {
    const LinearMatrix* m1;
    const LinearMatrix* m2_T;
    LinearMatrix* result;

    int n;
    int blockSize;

    int rowStart;
    int rowEnd;
};

void* worker(void* arg) {
    ThreadData* data = (ThreadData*)arg;

    const LinearMatrix& m1 = *data->m1;
    const LinearMatrix& m2_T = *data->m2_T;
    LinearMatrix& result = *data->result;

    int n = data->n;
    int blockSize = data->blockSize;

    for (int bi = data->rowStart; bi < data->rowEnd; bi += blockSize) {
        for (int bj = 0; bj < n; bj += blockSize) {
            for (int bk = 0; bk < n; bk += blockSize) {

                int i_end = min(bi + blockSize, data->rowEnd);
                int j_end = min(bj + blockSize, n);
                int k_end = min(bk + blockSize, n);

                for (int i = bi; i < i_end; i++) {

                    int j = bj;

                    for (; j + 3 < j_end; j += 4) {

                        __m256d c =
                            _mm256_loadu_pd(&result[i * n + j]);

                        for (int k = bk; k < k_end; k++) {

                            __m256d a =
                                _mm256_set1_pd(m1[i * n + k]);

                            __m256d b =
                                _mm256_loadu_pd(&m2_T[j * n + k]);

                            c = _mm256_fmadd_pd(a, b, c);
                        }

                        _mm256_storeu_pd(
                            &result[i * n + j], c);
                    }

                    for (; j < j_end; j++) {
                        for (int k = bk; k < k_end; k++) {
                            result[i * n + j] +=
                                m1[i * n + k] *
                                m2_T[j * n + k];
                        }
                    }
                }
            }
        }
    }

    return nullptr;
}

LinearMatrix dgemmOpt4(const LinearMatrix& m1, const LinearMatrix& m2, int n, int blockSize, int numThreads) {
    // Оптимизация за счет распараллеливания вычислений при помощи POSIX Threads
    LinearMatrix result(n * n, 0.0);

    LinearMatrix m2_T(n * n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            m2_T[j * n + i] = m2[i * n + j];
        }
    }

    vector<pthread_t> threads(numThreads);
    vector<ThreadData> threadData(numThreads);

    int rowsPerThread = n / numThreads;

    for (int t = 0; t < numThreads; t++) {

        int start = t * rowsPerThread;

        int end = (t == numThreads - 1)
            ? n
            : start + rowsPerThread;

        threadData[t] = {
            &m1,
            &m2_T,
            &result,
            n,
            blockSize,
            start,
            end
        };

        pthread_create(
            &threads[t],
            nullptr,
            worker,
            &threadData[t]
        );
    }

    for (int t = 0; t < numThreads; t++) {
        pthread_join(threads[t], nullptr);
    }

    return result;
}

LinearMatrix dgemmOpt5(const LinearMatrix &m1,const LinearMatrix &m2,int n,int blockSize) {
    // Оптимизация за счет распараллеливания вычислений при помощи OpenMP
    LinearMatrix result(n * n, 0.0);

    LinearMatrix m2_T(n * n);

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            m2_T[j * n + i] = m2[i * n + j];
        }
    }

    #pragma omp parallel for schedule(static)
    for (int bi = 0; bi < n; bi += blockSize) {
        for (int bj = 0; bj < n; bj += blockSize) {
            for (int bk = 0; bk < n; bk += blockSize) {
                int i_end = min(bi + blockSize, n);
                int j_end = min(bj + blockSize, n);
                int k_end = min(bk + blockSize, n);

                for (int i = bi; i < i_end; i++) {
                    int j = bj;
                    for (; j + 3 < j_end; j += 4) {

                        __m256d c =
                            _mm256_loadu_pd(&result[i * n + j]);

                        for (int k = bk; k < k_end; k++) {

                            __m256d a =
                                _mm256_set1_pd(m1[i * n + k]);

                            __m256d b =
                                _mm256_loadu_pd(
                                    &m2_T[j * n + k]);

                            c = _mm256_fmadd_pd(a, b, c);
                        }

                        _mm256_storeu_pd(
                            &result[i * n + j], c);
                    }

                    for (; j < j_end; j++) {
                        for (int k = bk; k < k_end; k++) {
                            result[i * n + j] +=
                                m1[i * n + k] *
                                m2[k * n + j];
                        }
                    }
                }
            }
        }
    }
    return result;
}

LinearMatrix dgemmOpt6(const LinearMatrix &m1, const LinearMatrix &m2, int n, int blockSize) {
    // Оптимизация за счет распараллеливания вычислений при помощи Intel TBB
    LinearMatrix result(n * n, 0.0);
    LinearMatrix m2_T(n * n);

    tbb::parallel_for(
        tbb::blocked_range2d<int>(
            0, n, 64,
            0, n, 64
        ),

        [&](const tbb::blocked_range2d<int>& r) {

            for (int i = r.rows().begin();
                 i < r.rows().end();
                 ++i) {

                for (int j = r.cols().begin();
                     j < r.cols().end();
                     ++j) {

                    m2_T[j * n + i] =
                        m2[i * n + j];
                }
            }
        }
    );

    tbb::parallel_for(

        tbb::blocked_range<int>(0, n, 256),

        [&](const tbb::blocked_range<int>& r) {

            for (int bi = r.begin();
                 bi < r.end();
                 bi += blockSize) {

                for (int bj = 0;
                     bj < n;
                     bj += blockSize) {

                    for (int bk = 0;
                         bk < n;
                         bk += blockSize) {

                        int i_end =
                            min(bi + blockSize, r.end());

                        int j_end =
                            min(bj + blockSize, n);

                        int k_end =
                            min(bk + blockSize, n);

                        for (int i = bi;
                             i < i_end;
                             ++i) {

                            int j = bj;

                            for (; j + 3 < j_end; j += 4) {

                                __m256d c =
                                    _mm256_loadu_pd(
                                        &result[i * n + j]);

                                for (int k = bk;
                                     k < k_end;
                                     ++k) {

                                    __m256d a =
                                        _mm256_set1_pd(
                                            m1[i * n + k]);

                                    __m256d b =
                                        _mm256_loadu_pd(
                                            &m2_T[j * n + k]);

                                    c =
                                        _mm256_fmadd_pd(a, b, c);
                                }

                                _mm256_storeu_pd(
                                    &result[i * n + j], c);
                            }

                            for (; j < j_end; ++j) {

                                for (int k = bk;
                                     k < k_end;
                                     ++k) {

                                    result[i * n + j] +=
                                        m1[i * n + k] *
                                        m2_T[j * n + k];
                                }
                            }
                        }
                    }
                }
            }
        },

        tbb::auto_partitioner()
    );
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
        
        for (int matrixSize = step < maxValue ? step : 0; matrixSize <= maxValue; matrixSize += step) {
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
            } else {
                LinearMatrix firstLinearMatrix(matrixSize * matrixSize);
                LinearMatrix secondLinearMatrix(matrixSize * matrixSize);
                
                createLinearMatrix(firstLinearMatrix, matrixSize);
                createLinearMatrix(secondLinearMatrix, matrixSize);
                
                LinearMatrix linearResult;
                auto start = chrono::high_resolution_clock::now();
                if (name == "dgemmOpt3") linearResult = dgemmOpt3(firstLinearMatrix, secondLinearMatrix, matrixSize, 32);
                else if (name == "dgemmOpt4") linearResult = dgemmOpt4(firstLinearMatrix, secondLinearMatrix, matrixSize, 32, 16);
                else if (name == "dgemmOpt5") linearResult = dgemmOpt5(firstLinearMatrix, secondLinearMatrix, matrixSize, 32);
                else if (name == "dgemmOpt6") linearResult = dgemmOpt6(firstLinearMatrix, secondLinearMatrix, matrixSize, 32);
                auto end = chrono::high_resolution_clock::now();
                duration = chrono::duration<double>(end - start).count();
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