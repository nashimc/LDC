
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <vector>
// #include "matrix.h"
#include "ldcsolver.h"

LidDrivenCavitySolver::LidDrivenCavitySolver(int grid, double dt, double lidVelocity)
:   
    grid(grid),
    lidVelocity(lidVelocity),

    u((grid), std::vector<double>(grid, 0)),
    uNew((grid), std::vector<double>(grid, 0)),
    uC((grid), std::vector<double>(grid, 0)),

    v((grid), std::vector<double>((grid), 0)),
    vNew((grid), std::vector<double>((grid), 0)),
    vC((grid), std::vector<double>((grid), 0)),

    p((grid), std::vector<double>((grid), 0)),
    pNew((grid), std::vector<double>((grid), 0)),
    pC((grid), std::vector<double>((grid), 0)),

    m((grid), std::vector<double>((grid), 0)),
    mn((grid), std::vector<double>((grid), 0)),
    mc((grid), std::vector<double>((grid), 0)),

    // dx();
    // dy();
    dt(dt)
{
    // this->Re = Re;
    // this->rho = rho;
    // this->nu = nu;
    this->dx = 1.0/(grid - 1);
    this->dy = 1.0/(grid - 1);


    this->initU();
    this->initV();
    this->initP();

}

void LidDrivenCavitySolver::initU(){
    // Init U
    //      initialise 0 everywhere
    //      first and second layer U = lidVelocity to average lidVelocity m/s lid movement at 1st layer nodes.
    for (int i = 0; i <= (this->grid - 1); ++i){            // 0
        for (int j = 0; j <= (this->grid - 1); ++j){        // -1
            this->u[i][j] = 0;
            this->u[0][j] = this->lidVelocity;
            this->u[1][j] = this->lidVelocity;
        }           
    }
}

void LidDrivenCavitySolver::initV(){
    // Init V
    //      initialise 0 everywhere.
    for (int i = 0; i <= (this->grid - 1); ++i){        // -1
        for (int j = 0; j <= (this->grid - 1); ++j){    // 0
            this->v[i][j] = 0.0;
        }           
    }
}

void LidDrivenCavitySolver::initP(){
    // Init P
    //      no flow, initialised as p = 1 everwhere.
    for (int i = 0; i <= (this->grid - 1); ++i){          // 0
        for (int j = 0; j <= (this->grid - 1); ++j){      // 0
            this->p[i][j] = 1;
        }
    } 
}

void LidDrivenCavitySolver::solveUVMomentum(){
    // Solve u-momentum in interior nodes.
    for (int i = 1; i <= (this->grid - 2); ++i){        // -1
        for (int j = 1; j <= (this->grid - 2); ++j){    // -2
            this->uNew[i][j] =  this->u[i][j]   - this->u[i][j] * (this->dt / this->dx) * (this->u[i][j] - this->u[i-1][j]) 
                                        - this->v[i][j] * (this->dt / this->dy) * (this->u[i][j] - this->u[i][j-1])
                                        - (this->dt / this->rho * 2 * this->dx) * (this->p[i+1][j] - this->p[i-1][j])
                                        + this->nu * (this->dt / pow(this->dx, 2) * (this->u[i+1][j] - 2 * this->u[i][j] + this->u[i-1][j])
                                                  + this->dt / pow(this->dy, 2) * (this->u[i][j+1] - 2 * this->u[i][j] + this->u[i][j-1]));
        }   
    }

    for (int i = 1; i <= this->grid - 2; ++i){      // -1
        this->uNew[i][0] = 0;                                                       // left wall
        this->uNew[i][this->grid-1] = 0;                                            // right wall
    }
    for (int j = 0; j <= this->grid - 1; ++j){     
        this->uNew[this->grid-1][j] = -this->uNew[this->grid-1][j];                 // bottom
        this->uNew[0][j] = 2 - this->uNew[1][j];                                    // top
    }
    

    // Solve v-momentum in interior nodes.
    for (int i = 1; i <= (this->grid - 2); ++i){      
        for (int j = 1; j <= (this->grid - 2); ++j){    // -1
            this->vNew[i][j] =  this->v[i][j]   - this->u[i][j] * (this->dt / this->dx) * (this->v[i][j] - this->v[i-1][j]) 
                                        - this->v[i][j] * (this->dt / this->dy) * (this->v[i][j] - this->v[i][j-1])
                                        - (this->dt / this->rho * 2 * this->dx) * (this->p[i+1][j] - this->p[i-1][j])
                                        + this->nu * (this->dt / pow(this->dx, 2) * (this->v[i+1][j] - 2 * this->v[i][j] + this->v[i-1][j])
                                                  + this->dt / pow(this->dy, 2) * (this->v[i][j+1] - 2 * this->v[i][j] + this->v[i][j-1]));
        }                 
    }

    for (int i = 1; i <= this->grid - 2; ++i){
        this->vNew[i][0] = -this->vNew[i][1];                                       // left wall
        this->vNew[i][this->grid] = -this->vNew[i][this->grid-1];                   // right wall
    }
    for (int j = 0; j <= (this->grid - 1); ++j){
        this->vNew[this->grid-1][j] = 0;                                            // bottom
        this->vNew[0][j] = 0;                                                       // top
    }
}

void LidDrivenCavitySolver::solveP(){
    for (int i = 1; i <= (this->grid - 2); ++i){        // -1
        for (int j = 1; j <= (this->grid - 2); ++j){    // -1
            this->p[i][j] = (((this->p[i+1][j] + this->p[i-1][j]) * pow(this->dy, 2) + (this->p[i][j+1] + this->p[i][j-1]) * pow(this->dx, 2)) / (2 * (pow(this->dx, 2) + pow(this->dy, 2))))
                                - ((this->rho * pow(this->dx, 2) * pow(this->dy, 2)) / (2 * (pow(this->dx, 2) + pow(this->dy, 2))))
                                    * ((1 / this->dt) * (((this->u[i+1][j] - this->u[i-1][j]) / (2 * this->dx)) + ((this->v[i][j+1] - this->v[i][j-1]) / (2 * this->dy))) 
                                        - ((this->u[i+1][j] - this->u[i-1][j]) / (2 * this->dx)) * ((this->u[i+1][j] - this->u[i-1][j]) / (2 * this->dx))
                                        - 2 * ((this->u[i][j+1] - this->u[i][j-1]) / (2 * this->dy))
                                        * ((this->v[i+1][j] - this->v[i-1][j]) / (2 * this->dx))
                                        - ((this->v[i][j+1] - this->v[i][j-1]) / (2 * this->dy))
                                        * ((this->v[i][j+1] - this->v[i][j-1]) / (2 * this->dy)));
        }
    }


    for (int i = 0; i <= (this->grid - 1); ++i){      // 0
        this->pNew[i][0] = this->pNew[i][1];                                            // left wall
        this->pNew[i][this->grid-1] = this->pNew[i][this->grid-2];                      // right wall
    }

    for (int j = 1; j <= (this->grid - 2); ++j){
        this->pNew[this->grid-1][j] = this->pNew[this->grid-2][j];                      // bottom
        this->pNew[0][this->grid-1] = this->pNew[1][this->grid-1];                      // top
    }
}  

void LidDrivenCavitySolver::checkError(){
    this->error = 0.0;
    
    for (int i = 1; i <= (this->grid - 2); ++i){
        for (int j = 1; j <= (this->grid - 2); ++j){
            this->m[i][j] = ((this->uNew[i][j] - this->uNew[i-1][j] ) / this->dx + (this->vNew[i][j] - this->vNew[i][j-1] ) / this->dy);
            this->error = this->error + fabs(this->m[i][j]);
        }
    }
    
    if (this->step % 10000 == 1){
        printf("Error is %5.8lf for the step %d\n", this->error, this->step);
    }

    this->iterate();
}

void LidDrivenCavitySolver::iterate(){
    // Iterating u
    for (int i = 0; i <= (this->grid - 1); ++i){
        for (int j = 0; j <= (this->grid - 1); ++j){
            this->u[i][j] = this->uNew[i][j];
        }
    }
    
    // Iterating v
    for (int i = 0; i <= (this->grid - 1); ++i){
        for (int j = 0; j <= (this->grid - 1); ++j){
            this->v[i][j] = this->vNew[i][j];
        }
    }
    
    // Iterating p
    for (int i = 0; i <= (this->grid - 1); ++i){
        for (int j = 0; j <= (this->grid - 1); ++j){
            this->p[i][j] = this->pNew[i][j];
        }
    }

    this->step = this->step + 1;
}

void LidDrivenCavitySolver::runSimulation(){
    while (this->error > 0.0000001){
        this->solveUVMomentum();
        this->solveP();
        this->checkError();
    }
    // this->writeData();
}

// void LidDrivenCavitySolver::writeData(){
//     FILE *fout2, *fout3;
//     fout2 = fopen("UVP.plt","w+t");
//     fout3 = fopen("Central_U.plt","w+t");

//     if (fout2 == NULL){
//         printf("\nERROR when opening file\n");
//         fclose(fout2);
//     }
//     else{
//         fprintf(fout2, "VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"P\"\n");
//         fprintf(fout2, "ZONE  F=POINT\n");
//         fprintf(fout2, "I=%d, J=%d\n", this->grid, this->grid );

//         for (int j = 0 ; j < (this->grid) ; ++j){
//             for (int i = 0 ; i < (this->grid) ; ++i){
//                 double xpos, ypos;
//                 xpos = i * this->dx;
//                 ypos = j * this->dy;

//                 fprintf(fout2, "%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\n", xpos, ypos, this->uC[i][j], this->vC[i][j], this->pC[i][j] );
//             }
//         }
//     }

//     fclose(fout2);
    
//     // CENTRAL --U
//     fprintf(fout3, "VARIABLES=\"U\",\"Y\"\n");
//     fprintf(fout3, "ZONE F=POINT\n");
//     fprintf(fout3, "I=%d\n", grid );

//     for (int j = 0; j < this->grid ; ++j){
//         double ypos;
//         ypos = (double) j * this->dy;

//         fprintf(fout3, "%5.8lf\t%5.8lf\n", (this->uC[grid/2][j] + this->uC[(grid/2)+1][j])/(2.), ypos);
//     }   
// }