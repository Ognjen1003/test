from pydantic import BaseModel, ConfigDict
from typing import List
from src.Classes.Component import Component


class ComponentInput(BaseModel):
    Cmp:Component
    class Config:
        extra = "forbid"   # NE dozvoli polja koja nisu definirana u Component


class EOSInputModel(BaseModel):
    components: List[ComponentInput]
    T: float
    P: float
    eos_type: str  # "PR" or "SRK"
    method: str
