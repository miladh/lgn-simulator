#ifndef FFTHELPER_H
#define FFTHELPER_H

#include <armadillo>
#include "special.h"

using namespace arma;

namespace lgnSimulator {

class FFTHelper
{
public:
    FFTHelper();
    ~FFTHelper();

    static vec fftFreq(int windowLength, double sampleSpacing = 1.);

    static cx_vec fftShift(cx_vec x);
    static cx_mat fftShift(cx_mat x);
    static cx_cube fftShift(cx_cube x);

    static cx_cube ifftShift(cx_cube x);

};
}
#endif // FFTHELPER_H
