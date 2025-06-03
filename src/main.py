from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import JSONResponse, HTMLResponse
from src.Endpoints.FlowModul import calculate_steps
from src.Endpoints.EOSModul import perform_eos_calculation
from src.Models.FlowModels import FlowInputModel
from src.Models.EOSModels import EOSInputModel
from src.Classes.Logging.LoggerSingleton import LoggerSingleton
from typing import List
from src.Classes.Component import Component
import pandas as pd

#import sys
#import os
#sys.path.append(os.path.abspath(os.path.dirname(__file__)))
LoggerSingleton()

app = FastAPI()



###############################################################
#
#     Endpoints
#
###############################################################

#Flowmodel
@app.post("/calculate_flow", response_model=dict)
async def calculate(input_data: FlowInputModel = None, request: Request = None):
    
    LoggerSingleton().log_info(f"calculate_flow: Received input data from {request.client.host}: {input_data}")

    try:
        if input_data is None:
            input_data = FlowInputModel()
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


#EOS calc
@app.post("/eos_calc")
async def calculate_EOS(input_data: EOSInputModel, request: Request):
    
    LoggerSingleton().log_info(f"eos_calc: Received input data from {request.client.host}: {input_data.model_dump()}")

    try:
        components: List[Component] = [comp_input.Cmp for comp_input in input_data.components]

        result = perform_eos_calculation(
            components=components,
            T=input_data.T,
            P=input_data.P,
            eos_type=input_data.eos_type,
            method=input_data.method
        )

        return result

    except ValueError as e:
        LoggerSingleton().log_info(f"ValueError: Received input data from {request.client.host}: error: {e}")
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        LoggerSingleton().log_info(f"Exception: Received input data from {request.client.host}: error: {e}")
        raise HTTPException(status_code=500, detail=f"Exception: {e}")



###############################################################
#
#     Maintenance and exception handling
#
###############################################################


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
            <h2>Flow calculation</h2>

            <p>To use this service, send a POST request to /calculate_flow with the following JSON data:</p>

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
                "visual": 0            # 0 or 1
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
                "visual": 0
            }
            </pre>

            <p>The response will return the result of the calculation as a JSON or HTML table ("visual": 1).</p>
            
            <pre>
                {"result":
                        [
                            {
                                "step":1,
                                "L":40000.0,
                                "p1":40.0,
                                "t":20.0,
                                "mu":1.66573821741095e-05,
                                "rho_g":92.8530388477516,
                                "u":3.2629440171282442,
                                "Re":5746229.778777943,
                                "ff":0.015240188214917173,
                                "dp":9.537876591534586,
                                "p2":30.462123408465413
                            }
                        ]
                }
            </pre>

            <h3>EOS Calculation</h3>

            <p>To perform an EOS calculation, send a POST request to /eos_calc with the following JSON data:</p>

            <pre>
            {
                "components": [
                    {"name": "Methane", "Tc": 190.6, "Pc": 4599000, "omega": 0.011, "A": 8.07131, "B": 1730.63, "C": 233.426},
                    {"name": "Ethane", "Tc": 305.4, "Pc": 4872000, "omega": 0.099, "A": 8.21201, "B": 1652.57, "C": 229.387}
                ], 
                "T": 300,           # Temperature (K)
                "P": 5e6,           # Pressure (Pa)
                "z": [0.5, 0.5],    # Mole fraction in mixture
                "eos_type": "PG"  # Type of EOS (PG or SRK)
            }
            </pre>

            <p><strong>Example:</strong></p>

            <pre>
            POST /eos_calc
            {
                "components": [
                    {"name": "Methane", "Tc": 190.6, "Pc": 4599000, "omega": 0.011, "A": 8.07131, "B": 1730.63, "C": 233.426},
                    {"name": "Ethane", "Tc": 305.4, "Pc": 4872000, "omega": 0.099, "A": 8.21201, "B": 1652.57, "C": 229.387}
                ],
                "T": 300,
                "P": 5e6,
                "z": [0.5, 0.5],
                "eos_type": "PG"
            }
            </pre>

            <p>The response will return the EOS calculation results as a JSON object with the following structure:</p>

            <pre>
                {
                    'V': 0.5, 
                    'x': [0.5000000197840379, 0.5000001082842513], 
                    'y': [0.4999999802159621, 0.4999998917157487], 
                    'eos_type': 'SRK'
                }
            </pre>

            <p>The result includes the volume, mole fractions in liquid and vapor phases, and a timestamp.</p>
        </body>
    </html>
    """