from pydantic import BaseModel, Field
from typing import Optional


class FittingK(BaseModel):
    at_m: Optional[float]
    K: Optional[float]
    kind: Optional[str] = None

    class Config:
        extra = "forbid"
