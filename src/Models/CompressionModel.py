from pydantic import BaseModel, ConfigDict
from typing import List

class CompressorInputModel(BaseModel):
    fractions: List[float]
    T1: float
    P1: float
    P2: float
    mass_flow: float   #kg/s
    isentropic_efficiency: float
    polytropic_efficiency: float
    NASA_9: str
    full_report: int
    polytropic_exponent: float

