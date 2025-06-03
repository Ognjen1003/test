from enum import Enum

class SolveMethod(str, Enum):
    FSOLVE = "fsolve"
    ROOT_SCALAR = "root_scalar"

class EOSType(str, Enum):
    PR = "PR"
    SRK = "SRK"