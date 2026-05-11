#include <iostream>
#include <vector>
#include <algorithm>
#include <stdexcept>

#include "../include/algorithms.h"
#include "../include/menu.h"

using namespace std;

int main(int argc, char* argv[]) {
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
                        case Algorithms::dgemmOpt3:
                            checkDgemmxPerformance(matrixSize, step, "dgemmOpt3");
                            break;
                        case Algorithms::dgemmOpt4:
                            checkDgemmxPerformance(matrixSize, step, "dgemmOpt4");
                            break;
                        case Algorithms::dgemmOpt5:
                            checkDgemmxPerformance(matrixSize, step, "dgemmOpt5");
                            break;
                        case Algorithms::dgemmOpt6:
                            checkDgemmxPerformance(matrixSize, step, "dgemmOpt6");
                            break;
                        // case Algorithms::dgemmOpt7:
                        //     checkDgemmxPerformance(matrixSize, step, "dgemmOpt7");
                        //     break;
                        // case Algorithms::dgemmOpt8:
                        //     checkDgemmxPerformance(matrixSize, step, "dgemmOpt8");
                        //     break;
                    }
                    cout << ">>> Plot stored in: ./data/plots\n"
                            ">>> Performance data stored in: ./data/stats\n\n";
                    break;
                }
                case MenuMode::calculation:
                {
                    Matrix firstMatrix(matrixSize, vector<double>(matrixSize));
                    Matrix secondMatrix(matrixSize, vector<double>(matrixSize));
                    Matrix result(matrixSize, vector<double>(matrixSize, 0.0));
                    
                    createMatrix(firstMatrix, matrixSize);
                    createMatrix(secondMatrix, matrixSize);

                    LinearMatrix firstLinearMatrix(matrixSize * matrixSize);
                    LinearMatrix secondLinearMatrix(matrixSize * matrixSize);
                    LinearMatrix linResult(matrixSize * matrixSize, 0.0);

                    createLinearMatrix(firstLinearMatrix, matrixSize);
                    createLinearMatrix(secondLinearMatrix, matrixSize);
                    
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
                        case Algorithms::dgemmOpt3:
                            linResult = dgemmOpt3(firstLinearMatrix, secondLinearMatrix, matrixSize, 32);
                            break;
                        case Algorithms::dgemmOpt4:
                            linResult = dgemmOpt4(firstLinearMatrix, secondLinearMatrix, matrixSize, 32, 16);
                            break;
                        case Algorithms::dgemmOpt5:
                            linResult = dgemmOpt5(firstLinearMatrix, secondLinearMatrix, matrixSize, 32);
                            break;
                        case Algorithms::dgemmOpt6:
                            linResult = dgemmOpt6(firstLinearMatrix, secondLinearMatrix, matrixSize, 32);
                            break;
                        // case Algorithms::dgemmOpt7:
                        //     linResult = dgemmOpt7(firstLinearMatrix, secondLinearMatrix, matrixSize, 32, rank, size);
                        //     break;
                        // case Algorithms::dgemmOpt8:
                        //     linResult = dgemmOpt8(firstLinearMatrix, secondLinearMatrix, matrixSize, 32);
                        //     break;
                    }
                    
                    if (algorithm == Algorithms::dgemmOpt3) {
                        cout << "\n>>> First matrix:\n";
                        displayLinearMatrix(firstLinearMatrix, matrixSize);

                        cout << "\n>>> Second matrix:\n";
                        displayLinearMatrix(secondLinearMatrix, matrixSize);

                        cout << "\n>>>Result:\n";
                        displayLinearMatrix(linResult, matrixSize);
                        cout << "\n";
                        break;
                    } else if (algorithm == Algorithms::dgemmOpt4) {
                        cout << "\n>>> First matrix:\n";
                        displayLinearMatrix(firstLinearMatrix, matrixSize);

                        cout << "\n>>> Second matrix:\n";
                        displayLinearMatrix(secondLinearMatrix, matrixSize);

                        cout << "\n>>>Result:\n";
                        displayLinearMatrix(linResult, matrixSize);
                        cout << "\n";
                        break;
                    } else if (algorithm == Algorithms::dgemmOpt5) {
                        cout << "\n>>> First matrix:\n";
                        displayLinearMatrix(firstLinearMatrix, matrixSize);

                        cout << "\n>>> Second matrix:\n";
                        displayLinearMatrix(secondLinearMatrix, matrixSize);

                        cout << "\n>>>Result:\n";
                        displayLinearMatrix(linResult, matrixSize);
                        cout << "\n";
                        break;
                    } else if (algorithm == Algorithms::dgemmOpt6) {
                        cout << "\n>>> First matrix:\n";
                        displayLinearMatrix(firstLinearMatrix, matrixSize);

                        cout << "\n>>> Second matrix:\n";
                        displayLinearMatrix(secondLinearMatrix, matrixSize);

                        cout << "\n>>>Result:\n";
                        displayLinearMatrix(linResult, matrixSize);
                        cout << "\n";
                        break;
                    // } else if (algorithm == Algorithms::dgemmOpt7) {
                    //     cout << "\n>>> First matrix:\n";
                    //     displayLinearMatrix(firstLinearMatrix, matrixSize);

                    //     cout << "\n>>> Second matrix:\n";
                    //     displayLinearMatrix(secondLinearMatrix, matrixSize);

                    //     cout << "\n>>>Result:\n";
                    //     displayLinearMatrix(linResult, matrixSize);
                    //     cout << "\n";
                    //     break;
                    // } else if (algorithm == Algorithms::dgemmOpt8) {
                    //     cout << "\n>>> First matrix:\n";
                    //     displayLinearMatrix(firstLinearMatrix, matrixSize);

                    //     cout << "\n>>> Second matrix:\n";
                    //     displayLinearMatrix(secondLinearMatrix, matrixSize);

                    //     cout << "\n>>>Result:\n";
                    //     displayLinearMatrix(linResult, matrixSize);
                    //     cout << "\n";
                    //     break;
                    } else {
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
        }
    } catch (const exception &e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}