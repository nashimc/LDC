// #pragma once
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <vector>
#include "matrix.h"
// #include "ldcsolver.h"


Matrix::Matrix(double rows, double columns){
    this->rows = rows;
    this->columns = columns;
    this->data = (rows * columns);
}

double& Matrix::operator(size_t i, size_t j){
    return mData[i * mCols + j];
}

double Matrix::operator(size_t i, size_t j) const{
     return mData[i * mCols + j];
}