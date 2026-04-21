#pragma once
#include <vector>
#include <string>

using namespace std;

enum class Algorithms {
    dgemmBlas,
    dgemmOpt1,
    dgemmOpt2
};

using Matrix = vector<vector<double>>;

Matrix dgemmBlas(const Matrix &m1, const Matrix &m2);
Matrix dgemmOpt1(const Matrix &m1, const Matrix &m2);
Matrix dgemmOpt2(const Matrix &m1, const Matrix &m2, int blockSize);

void createMatrix(Matrix &m, int mtxSize);
void displayMatrix(const Matrix &matrix);
void checkDgemmxPerformance(int maxValue, int step, const string &name);