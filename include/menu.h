#pragma once

enum class MenuMode {
    exit,
    speedTest,
    calculation
};

int selectMode(void);
int selectAlgorithm(void);
int getMatrixSize(void);
int getStep(void);