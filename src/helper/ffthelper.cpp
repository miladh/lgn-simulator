#include "ffthelper.h"

using namespace lgnSimulator;


FFTHelper::FFTHelper()
{

}

FFTHelper::~FFTHelper()
{

}


vec FFTHelper::fftFreq(int windowLength, double sampleSpacing){

    double val = 1. / (windowLength*sampleSpacing);
    int N = int((windowLength-1)/2) + 1;
    vec w1 = linspace(0, N-1, N);
    vec w2 = linspace(-int(windowLength/2), -1, N-Special::isOdd(windowLength));
    vec result = join_cols(w1,w2);

    return result * val;
}


cx_vec FFTHelper::fftShift(cx_vec x)
{
    int n = x.n_rows;
    int n2 = int((n+1)/2);
    cx_vec shifted = join_vert(x.rows(n2, n-1), x.rows(0,n2-1));

    return shifted;

}


cx_mat FFTHelper::fftShift(cx_mat x)
{
    int n = x.n_rows;
    int n2 = int((n+1)/2);

    int m = x.n_cols;
    int m2 = int((m+1)/2);

    cx_mat shifted = join_vert(x.rows(n2, n-1), x.rows(0,n2-1));
    shifted = join_horiz(shifted.cols(m2, m-1), shifted.cols(0,m2-1));

    return shifted;
}

cx_cube FFTHelper::fftShift(cx_cube x)
{
    cx_cube shifted = x;
    int n = x.n_rows;
    int n2 = int((n+1)/2);

    int m = x.n_cols;
    int m2 = int((m+1)/2);

    int o = x.n_slices;
    int o2 = int((o+1)/2);

    for(int i = 0; i < int(x.n_slices); i++){
        shifted.slice(i) = join_vert(x.slice(i).rows(n2, n-1),
                                     x.slice(i).rows(0,n2-1));
        shifted.slice(i) = join_horiz(shifted.slice(i).cols(m2, m-1),
                                      shifted.slice(i).cols(0,m2-1));
    }

    cx_cube tmp = shifted;
    shifted.slices(0, (o-1) - o2) = tmp.slices(o2, o-1);
    shifted.slices((o-1-o2)+1, (o-1-o2)+1+(o2-1)) = tmp.slices(0, o2-1);
    return shifted;
}






cx_cube FFTHelper::ifftShift(cx_cube x)
{
    cx_cube shifted = x;
    int n = x.n_rows;
    int n2 = n-int((n+1)/2);

    int m = x.n_cols;
    int m2 = m-int((m+1)/2);

    int o = x.n_slices;
    int o2 = o-int((o+1)/2);



    for(int i = 0; i < int(x.n_slices); i++){
        shifted.slice(i) = join_vert(x.slice(i).rows(n2, n-1), x.slice(i).rows(0,n2-1));
        shifted.slice(i) = join_horiz(
                    shifted.slice(i).cols(m2, m-1),
                    shifted.slice(i).cols(0,m2-1));
    }

    cx_cube tmp = shifted;
    shifted.slices(0, (o-1) - o2) = tmp.slices(o2, o-1);
    shifted.slices((o-1-o2)+1, (o-1-o2)+1+(o2-1)) = tmp.slices(0, o2-1);


    return shifted;
}


