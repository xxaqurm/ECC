#include <iostream>
#include <limits>

using namespace std;

int selectMode(void) {
    int input;
    cout << "Program modes: \n"
            "0. Exit from program\n"
            "1. SpeedTest\n"
            "2. Calculation\n";

    while (true) {
        cout << "Select configuratioin (enter number): ";
        cin >> input;

        if (input < 0 || input > 2) {
            cout << "Invalid input. Please try again." << endl;
            continue;
        }

        return input;
    }
}

int selectAlgorithm(void) {
    int algorithmNumber = 0;
    while (algorithmNumber < 1 || algorithmNumber > 4) {
        cout << "Algorithms:\n";
        cout << "\t1. dgemmBlas\n";
        cout << "\t2. dgemmOpt1\n";
        cout << "\t3. dgemmOpt2\n";
        cout << "\t4. dgemmOpt3\n";
        cout << "Select Algorithm (enter number): ";
        cin >> algorithmNumber;
        
        if (cin.fail()) {
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            algorithmNumber = 0;
            cout << "Invalid input. Please enter a number between 1 and 3.\n";
        }
    }
    return algorithmNumber - 1;
}

int getMatrixSize(void) {
    int matrixSize = 0;
    while (matrixSize < 1) {
        cout << "Enter Matrix Size: ";
        cin >> matrixSize;
        
        if (cin.fail() || matrixSize < 1) {
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            matrixSize = 0;
            cout << "Invalid input. Please enter a positive integer.\n";
        }
    }
    return matrixSize;
}

int getStep(void) {
    int step = 0;
    while (step < 1) {
        cout << "Enter Step: ";
        cin >> step;
        
        if (cin.fail() || step < 1) {
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            step = 0;
            cout << "Invalid input. Please enter a positive integer.\n";
        }
    }
    return step;
}