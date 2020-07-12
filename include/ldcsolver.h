// #pragma once
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <vector>

class LidDrivenCavitySolver{

// Example indexed 4x4 matrix:
//          [ 0,0 , 0,1 , 0,2 , 0,3 ]
//          [ 1,0 , 1,1 , 1,2 , 1,3 ]
//          [ 2,0 , 2,1 , 2,2 , 2,3 ]
//          [ 3,0 , 3,1 , 3,2 , 3,3 ]

public:
    int grid;
    double lidVelocity;

    std::vector<std::vector<double>> u, uNew, uC;
    std::vector<std::vector<double>> v, vNew, vC;
    std::vector<std::vector<double>> p, pNew, pC;
    std::vector<std::vector<double>> m, mn, mc;

    // double u, uNew, uC;
    // double v, vNew, vC;
    // double p, pNew, pC;
    // double m, mn, mc;
   
    
    double dx;
    double dy;
    double dt;
    double delta = 2;
    double error = 1.0;
    double rho = 1;
    double nu = 0.1;
    double Re = 100;
    int step = 1;
    
public:
    LidDrivenCavitySolver(int grid, double dt, double lidVelocity);
    void initU();
    void initV();
    void initP();  
    void solveUVMomentum();  
    void solveP();       
    void checkError();
    void iterate();
    void runSimulation();
    // void writeData();
};