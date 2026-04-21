#include <iostream>
#include <vector>
#include <algorithm>
#include <stdexcept>

#include "../include/algorithms.h"
#include "../include/menu.h"

using namespace std;

int main() {
    try {
        while (true) {
            MenuMode mode = static_cast<MenuMode>(selectMode());

            if (mode == MenuMode::exit) {
                return 0;
            }

            Algorithms algorithm = static_cast<Algorithms>(selectAlgorithm());
            int matrixSize = getMatrixSize();
            
            switch (mode) {
                case MenuMode::speedTest:
                {
                    int step = getStep();
                    matrixSize++;
                    switch(algorithm) {
                        case Algorithms::dgemmBlas:
                            checkDgemmxPerformance(matrixSize, step, "dgemmBlas");
                            break;
                        case Algorithms::dgemmOpt1:
                            checkDgemmxPerformance(matrixSize, step, "dgemmOpt1");
                            break;
                        case Algorithms::dgemmOpt2:
                            checkDgemmxPerformance(matrixSize, step, "dgemmOpt2");
                            break;
                    }
                    cout << ">>> Plot stored in: ./data/plots\n"
                            ">>> Performance data stored in: ./data/stats\n\n";
                    break;
                }
                case MenuMode::calculation:
                {
                    Matrix firstMatrix(matrixSize, vector<double>(matrixSize));
                    Matrix secondMatrix(matrixSize, vector<double>(matrixSize));

                    createMatrix(firstMatrix, matrixSize);
                    createMatrix(secondMatrix, matrixSize);

                    Matrix result(matrixSize, vector<double>(matrixSize, 0.0));
                    switch(algorithm) {
                        case Algorithms::dgemmBlas:
                            result = dgemmBlas(firstMatrix, secondMatrix);
                            break;
                        case Algorithms::dgemmOpt1:
                            result = dgemmOpt1(firstMatrix, secondMatrix);
                            break;
                        case Algorithms::dgemmOpt2:
                            result = dgemmOpt2(firstMatrix, secondMatrix, 32);
                            break;
                    }
                    cout << "\n>>> First matrix:\n";
                    displayMatrix(firstMatrix);
                    
                    cout << "\n>>> Second matrix:\n";
                    displayMatrix(secondMatrix);

                    cout << "\n>>> Result:\n";
                    displayMatrix(result);
                    cout << "\n";
                    break;
                }
            }
        }
    } catch (const exception &e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
}