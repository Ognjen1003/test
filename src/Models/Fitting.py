from pydantic import BaseModel, Field
from typing import Optional


class FittingK(BaseModel):
    at_m: Optional[float]
    K: Optional[float] = 0.0
    kind: Optional[str] = None
    d_in_new: Optional[float] = None
    K_ref: Optional[str] = None         # opcionalno: "UPSTREAM"/"DOWNSTREAM"/"MIN_DIAMETER"

    class Config:
        extra = "forbid"
