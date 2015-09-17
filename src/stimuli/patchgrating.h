#ifndef PATCHGRATING_H
#define PATCHGRATING_H

#include "stimuli.h"


class PatchGrating : public Stimuli
{
public:
    PatchGrating(const Config *cfg);
    ~PatchGrating();

private:
    double m_contrast = 0.0;
    double m_spotDiameter = 0.0;


    double valueAtPoint(vec2 rVec, double t);
    double fourierTransformAtFrequency(vec2 kVec, double w);



};

#endif // PATCHGRATING_H
