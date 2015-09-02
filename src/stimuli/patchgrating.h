#ifndef PATCHGRATING_H
#define PATCHGRATING_H

#include "stimuli.h"


class PatchGrating : public Stimuli
{
public:
    PatchGrating(const Config *cfg);
    ~PatchGrating();

    double spatial(vec2 rVec, double t);
    double frequency(vec2 kVec, double w);


private:
    double m_contrast = 0.0;
    double m_spotDiameter = 0.0;

};

#endif // PATCHGRATING_H
