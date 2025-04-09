from pydantic import BaseModel
from typing import List

class ComponentInput(BaseModel):
    name: str
    Tc: float
    Pc: float
    omega: float
    A: float
    B: float
    C: float
    
class EOSInputModel(BaseModel):
    components: List[ComponentInput]
    T: float
    P: float
    z: List[float]
    eos_type: str  # "PR" or "SRK"

