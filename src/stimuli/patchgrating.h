#ifndef PATCHGRATING_H
#define PATCHGRATING_H

#include "stimuli.h"


class PatchGrating : public Stimuli
{
public:
    PatchGrating(const Config *cfg);
    ~PatchGrating();

    double real(vec2 rVec, double t);
    double complex(vec2 kVec, double w);


private:
    double m_C = 0.0;
    double m_d = 0.0;

};

#endif // PATCHGRATING_H
