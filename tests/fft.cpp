#include <unittest++/UnitTest++.h>
#include <armadillo>
#include <iostream>


#include "stimuli/patchgrating.h"
#include "spatialKernels/dog.h"
#include <lib.h>


using namespace std;
using namespace arma;




//TEST(dogFFT){

//    int nPoints = 128;
//    double dr = 0.1;

//    vec spatialMesh = linspace(0, 1, nPoints);
//    vec freqMesh = linspace(-PI/dr+dr, PI/dr, nPoints);

//    mat spatial = zeros<mat>(nPoints, nPoints);
//    mat freq = zeros<cx_mat>(nPoints, nPoints);

//    DOG dog(1.0, 0.25, 0.85, 0.83);

//    for(int i = 0; i < int(spatialMesh.n_elem); i++){
//        for(int j = 0; j < int(spatialMesh.n_elem); j++){
//            spatial(i,j) = dog.real({m_spatialMesh[i], m_spatialMesh[j]});
//        }
//    }


//    for(int i = 0; i < int(m_freqMesh.n_elem); i++){
//        for(int j = 0; j < int(m_freqMesh.n_elem); j++){

//        freq(i,j) = dog.complex({m_freqMesh[i], m_freqMesh[j]});


//        }
//    }

//    complexResponse = Functions::fftShift(freq);
//    m_response = Functions::fftShift(real(ifft2(complexResponse)));




//}


