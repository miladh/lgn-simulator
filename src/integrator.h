#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <iostream>
#include <math.h>

class Integrator
{
public:
    Integrator(double upperLim, double lowerLim, double steps);
    ~Integrator();

    virtual void integrate() = 0;
    virtual double result() const = 0;

    double upperLim() const;
    double lowerLim() const;
    double panels() const;
    double stepLength() const;

private:
    double m_upperLim;
    double m_lowerLim;
    double m_panels;
    double m_stepLength;

};

#endif // INTEGRATOR_H
