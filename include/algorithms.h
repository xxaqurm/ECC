#pragma once
#include <vector>
#include <string>

using namespace std;

enum class Algorithms {
    dgemmBlas,
    dgemmOpt1,
    dgemmOpt2,
    dgemmOpt3
};

using Matrix = vector<vector<double>>;
using LinearMatrix = vector<double>;

Matrix dgemmBlas(const Matrix &m1, const Matrix &m2);
Matrix dgemmOpt1(const Matrix &m1, const Matrix &m2);
Matrix dgemmOpt2(const Matrix &m1, const Matrix &m2, int blockSize);
LinearMatrix dgemmOpt3(const LinearMatrix &m1, const LinearMatrix &m2, int n, int blockSize);

void createMatrix(Matrix &m, int mtxSize);
void createLinearMatrix(LinearMatrix &m, int mtxSize);
void displayMatrix(const Matrix &matrix);
void displayLinearMatrix(const LinearMatrix &matrix, int matrixSize);
void checkDgemmxPerformance(int maxValue, int step, const string &name);