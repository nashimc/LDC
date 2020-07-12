#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <vector>
// #include "matrix.h"
#include "ldcsolver.h"



int main(){
  LidDrivenCavitySolver ldc(41, 0.02, 0);
  ldc.runSimulation();
  // std::cout << ldc.grid  << std::endl;


  system("PAUSE");

  return 0;
}

