from enum import Enum

class SolveMethod(str, Enum):
    FSOLVE = "FSOLVE"
    ROOT_SCALAR = "ROOT_SCALAR"

class EOSType(str, Enum):
    PR = "PR"
    SRK = "SRK"