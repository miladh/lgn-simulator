#ifndef TRAPEZOIDAL_H
#define TRAPEZOIDAL_H

#include "integrator.h"

using namespace std;

class Trapezoidal : public Integrator
{
public:
    Trapezoidal(double upperLim = 1., double lowerLim=-1, double steps = 10);
    ~Trapezoidal();

public:
    void integrate();
    double result() const;

private:
    double m_result = 0;

    double advance(double x) const;

};

#endif // TRAPEZOIDAL_H
