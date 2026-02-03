from enum import Enum


class SolveMethod(str, Enum):
    FSOLVE = "FSOLVE"
    ROOT_SCALAR = "ROOT_SCALAR"

class EOSType(str, Enum):
    PR = "PR"
    SRK = "SRK"

class Phase(str, Enum):
    VAPOR = "VAPOR"
    LIQUID = "LIQUID"
    VAPORLIQUID = "VAPORLIQUID"

class CASES(str, Enum):
    CASE1 = "CASE1"
    #PVT = "PVT"
    CO2 = "CO2"
    OXY1 = "OXY1"
    OXY2 = "OXY2"
    OXY3 = "OXY3"

class CONSTANTS(float, Enum):
    R = 8.314462618
    NA = 6.02214076e23

class StepMode(str, Enum):
    FIXED = "FIXED"              # koraci jednake duljine (max_step_m)
    BY_FITTINGS = "BY_FITTINGS"  # rezanje trase po lokacijama koljena


class ViscosityMethod(str, Enum):
    WILKE = "WILKE"              # dilute-gas Wilke
    LBC = "LBC"                  # LBC/JST korekcija (gusti plin)
    AUTO = "AUTO"                # heuristika, mozda nepotrebno ali neka bude