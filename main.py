from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import HTMLResponse
from calculations import calculate_steps
from Classes.models import InputModel
from Classes.LoggerSingleton import LoggerSingleton
from fastapi.responses import JSONResponse
import pandas as pd

#import sys
#import os
#sys.path.append(os.path.abspath(os.path.dirname(__file__)))
LoggerSingleton()

app = FastAPI()

@app.post("/calculate", response_model=dict)
async def calculate(input_data: InputModel = None, request: Request = None):
    
    LoggerSingleton().log_info(f"Received input data from {request.client.host}: {input_data}")

    try:
        if input_data is None:
            input_data = InputModel()
        result = calculate_steps(input_data.nsteps, input_data.L, input_data.d_in, 
                                input_data.e, input_data.p, input_data.T, input_data.qm, input_data.case)
        if input_data.visual == 0:
            return JSONResponse(content={"result": result})
        else: 
            df = pd.DataFrame(result)
            return HTMLResponse(content=df.to_html(index=False))
    except ValueError as e:
        LoggerSingleton().log_info(f"ValueError: Received input data from {request.client.host}: {input_data}, error: {e}")
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        LoggerSingleton().log_info(f"Exception: Received input data from {request.client.host}: {input_data}, error: {e}")
        raise HTTPException(status_code=400, detail=str("Exception: " + e.__str__()))
    
@app.middleware("http")
async def log_errors_middleware(request: Request, call_next):
    try:
        response = await call_next(request)
        return response
    except Exception as e:
        client_ip = request.client.host  # Dohvati IP adresu
        LoggerSingleton.log_error(f"Uncaught error: {str(e)}", ip_address=client_ip)
        return JSONResponse(
            status_code=500,
            content={"detail": "Internal Server Error"}
        )

@app.get("/error")
async def error_route(request: Request):
    client_ip = request.client.host  # Dohvati IP adresu
    LoggerSingleton.log_error("Route error", ip_address=client_ip)
    raise RuntimeError("Route error!")

@app.get("/", response_class=HTMLResponse)
async def root():
    return """
    <html>
        <body>
            <h2>Welcome to the calculation service!</h2>

            <p>To use this service, send a POST request to /calculate with the following JSON data:</p>

            <pre>
            {
                "nsteps": 2,           # Number of steps for the calculation
                "L": 40000,            # Length (m)
                "d_in": 0.315925,      # Inner diameter (m)
                "e": 0.0001,           # Pipe roughness
                "p": 4000000,          # Pressure (Pa)
                "T": 293.15,           # Temperature (K)
                "qm": 23.75,           # Mass flow rate (kg/m^3)
                "case": "case1"        # Case name (case1, case2, case3, CO2)
            }
            </pre>

            <p><strong>Example:</strong></p>

            <pre>
            POST /calculate
            {
                "nsteps": 2,
                "L": 40000,
                "d_in": 0.315925,
                "e": 0.0001,
                "p": 4000000,
                "T": 293.15,
                "qm": 23.75,
                "case": "case1"
            }
            </pre>

            <p>The response will return the result of the calculation as a float.</p>
        </body>
    </html>
    """