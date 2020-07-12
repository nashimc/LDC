// #pragma once
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <vector>


class Matrix{
private:
  size_t rows;
  size_t columns;
  std::vector<double> data;

public:
  Matrix(size_t rows, size_t columns);
  double& operator()(size_t i, size_t j);
  double operator()(size_t i, size_t j) const;
};
