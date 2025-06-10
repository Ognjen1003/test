#pragma once
#ifndef WRAPPER_H
#define WRAPPER_H

#include <vector>
#include "Component.h"
#include "SolveResult.h"

std::vector<SolveResult> mainFromPython(
    const std::vector<Component>& components,
    const std::vector<double>& temperatures,
    const std::vector<double>& pressures);

#endif
