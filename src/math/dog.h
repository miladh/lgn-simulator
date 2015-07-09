#ifndef DOG_H
#define DOG_H

#define PI 3.14159265359

#include <armadillo>


using namespace std;
using namespace arma;


class DOG
{
public:
    DOG(double A, double a, double B, double b);
    ~DOG();


    double real(vec2 r);
    double complex(vec2 k);

private:
    double m_A, m_a, m_B, m_b = 0.0;

};

#endif // DOG_H
