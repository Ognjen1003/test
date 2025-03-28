from typing import Optional
from pydantic import BaseModel

class InputModel(BaseModel):
    nsteps: Optional[int] = 2
    L: Optional[int] = 40000
    d_in: Optional[float] = 0.315925
    e: Optional[float] = 0.00001
    p: Optional[float] = 4000000
    T: Optional[float] = 293.15
    qm: Optional[float] = 23.75
    case: Optional[str] = "case1"