#ifndef CALCULATIONS_H
#define CALCULATIONS_H

#include "component.h"
#include <cmath>

class Calculations {
public:
    static double wilson_K(const Component& component, double T, double P);
};

#endif // CALCULATIONS_H