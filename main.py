from fastapi import FastAPI
from datetime import datetime
from pydantic import BaseModel
from calculations import calculate_division

app = FastAPI()

class InputModel(BaseModel):
    nsteps: int
    length: int


@app.get("/calculateBlank")
def calculate():
    result = calculate_division(40, 40000)
    return result

@app.post("/calculate", response_model=float)
async def calculate(input_data: InputModel):
    result = calculate_division(input_data.nsteps, input_data.length)
    return result

@app.get("/time")
def get_current_time():
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return {"Ovo je FastAPI prvi service point - vrijeme": current_time}

@app.get("/hello/{name}")
def say_hello(name: str):
    return {"message": f"Hello, {name}!"}