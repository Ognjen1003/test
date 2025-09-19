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
    PVT = "PVT"
    CO2 = "CO2"

class CONSTANTS(float, Enum):
    R = 8.314462618
    NA = 6.02214076e23