from typing import Optional, List
from pydantic import BaseModel
from src.EnumsClasses.MethodsAndTypes import CASES, StepMode, ViscosityMethod
from src.Models.Fitting import FittingK

class FlowInputModel(BaseModel):
    nsteps: Optional[int] = 10
    L: Optional[float] = 40000
    d_in: Optional[float] = 0.315925
    e: Optional[float] = 0.00001
    p: Optional[float] = 4000000
    T: Optional[float] = 293.15
    qm: Optional[float] = 23.75
    case: Optional[str] = "case1"
    visual: Optional[int] = 0

    step_mode: Optional[StepMode] = StepMode.BY_FITTINGS
    max_step_m: Optional[float] = 5000.0

    fittings: Optional[List[FittingK]] = []
    viscosity_method: Optional[ViscosityMethod] = ViscosityMethod.AUTO
    virtual_steps: Optional[int] = None

    class Config:
        extra = "forbid"