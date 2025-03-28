from fastapi import FastAPI, HTTPException
from calculations import calculate_steps
from Classes.models import InputModel

#import sys
#import os
#sys.path.append(os.path.abspath(os.path.dirname(__file__)))

app = FastAPI()

@app.post("/calculate", response_model=float)
async def calculate(input_data: InputModel = None):
    try:
        if input_data is None:
            input_data = InputModel()
        result = calculate_steps(input_data.nsteps, input_data.L, input_data.d_in, 
                                input_data.e, input_data.p, input_data.T, input_data.qm, input_data.case)
        return result
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    
