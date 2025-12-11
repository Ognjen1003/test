from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import JSONResponse, HTMLResponse
from src.Endpoints.FlowModul import calculate_steps
from src.Endpoints.EOSModul import perform_eos_calculation
from src.Endpoints.CompressorModul import compressor_thermodynamics
from src.Models.FlowModels import FlowInputModel
from src.Models.EOSModels import EOSInputModel
from src.Models.CompressionModel import CompressorInputModel
from src.Classes.Logging.LoggerSingleton import LoggerSingleton
from typing import List
from src.Models.Component import Component
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
            T_K=input_data.T,
            P_bar=input_data.P,
            eos_type=input_data.eos_type,
            method=input_data.method,
            calculate_enthalpy=input_data.calculate_enthalpy
        )

        return result

    except ValueError as e:
        LoggerSingleton().log_info(f"ValueError: Received input data from {request.client.host}: error: {e}")
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        LoggerSingleton().log_info(f"Exception: Received input data from {request.client.host}: error: {e}")
        raise HTTPException(status_code=500, detail=f"Exception main.py: {e}")


#Compressor/thermodynamics calc
@app.post("/compressor_calc")
async def calculate_compressor_thermodynamics(input_data: CompressorInputModel, request: Request):
    
    LoggerSingleton().log_info(f"compressor_calc: Received input data from {request.client.host}: {input_data.model_dump()}")

    try:

        result = compressor_thermodynamics(
            fractions=input_data.fractions,
            P1=input_data.P1,
            P2=input_data.P2,
            T1=input_data.T1,
            mass_flow=input_data.mass_flow,
            isentropic_efficiency=input_data.isentropic_efficiency,
            polytropic_efficiency= input_data.polytropic_efficiency,
            full_report= input_data.full_report,
            polytropic_exponent = input_data.polytropic_exponent
        )

        return result

    except ValueError as e:
        LoggerSingleton().log_info(f"ValueError: Received input data from {request.client.host}: error: {e}")
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        LoggerSingleton().log_info(f"Exception: Received input data from {request.client.host}: error: {e}")
        raise HTTPException(status_code=500, detail=f"Exception main.py: {e}")

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
            <h2>Flow calculation (/calculate_flow)</h2>

            <p>Pošalji POST zahtjev na <code>/calculate_flow</code> sa sljedećim JSON-om:</p>

            <pre>
            {
                "nsteps": 2,           # broj koraka diskretizacije
                "L": 40000,            # duljina cjevovoda (m)
                "d_in": 0.315925,      # unutarnji promjer (m)
                "e": 0.0001,           # hrapavost cijevi
                "p": 4000000,          # početni tlak (Pa)
                "T": 293.15,           # početna temperatura (K)
                "qm": 23.75,           # maseni protok (kg/m^3)
                "case": "CASE1",       # CASE1, CO2 ili PVT
                "composition": [],     # samo za case="PVT" (lista kompozicije / podataka)
                "visual": 0            # 0 = JSON, 1 = HTML tablica
            }
            </pre>

            <p><strong>Napomena:</strong></p>
            <ul>
                <li><code>case = "CASE1"</code> – koristi lookup tablicu za case1.</li>
                <li><code>case = "CO2"</code> – pretpostavke za čisti CO₂.</li>
                <li><code>case = "PVT"</code> – koristi dostavljeni <code>composition</code> kao izvor podataka.</li>
            </ul>

            <p><strong>Primjer poziva:</strong></p>

            <pre>
            POST /calculate_flow
            {
                "nsteps": 2,
                "L": 40000,
                "d_in": 0.315925,
                "e": 0.0001,
                "p": 4000000,
                "T": 293.15,
                "qm": 23.75,
                "case": "CASE1",
                "visual": 0
            }
            </pre>

            <p>Povratna vrijednost za <code>"visual": 0</code> je JSON s ključem <code>"result"</code> koji sadrži niz zapisa po koracima (ili JSON string, ovisno o implementaciji):</p>

            <pre>
            {
                "result":
                    [
                        {
                            "step": 1,
                            "L": 40000.0,
                            "p1": 40.0,
                            "t": 20.0,
                            "mu": 1.66573821741095e-05,
                            "rho_g": 92.8530388477516,
                            "u": 3.2629440171282442,
                            "Re": 5746229.778777943,
                            "ff": 0.015240188214917173,
                            "dp": 9.537876591534586,
                            "p2": 30.462123408465413
                        }
                    ]
            }
            </pre>

            <hr/>

            <h3>EOS Calculation (/eos_calc)</h3>

            <p>Za EOS račun, pošalji POST na <code>/eos_calc</code> sa sljedećim JSON-om:</p>

            <pre>
            {
                "components": [
                    {
                        "Cmp": {
                            "name": "Methane",
                            "Tc": 190.6,
                            "Pc": 4599000,
                            "omega": 0.011,
                            "A": 8.07131,
                            "B": 1730.63,
                            "C": 233.426
                        }
                    },
                    {
                        "Cmp": {
                            "name": "Ethane",
                            "Tc": 305.4,
                            "Pc": 4872000,
                            "omega": 0.099,
                            "A": 8.21201,
                            "B": 1652.57,
                            "C": 229.387
                        }
                    }
                ],
                "T": 300,                    # temperatura (K)
                "P": 5e6,                    # tlak (Pa)
                "eos_type": "SRK",           # "PR" ili "SRK"
                "method": "FSOLVE",          # metoda rješavanja (npr. FSOLVE)
                "calculate_enthalpy": false  # ako je implementirano u modelu
            }
            </pre>

            <p><strong>Primjer:</strong></p>

            <pre>
            POST /eos_calc
            {
                "components": [
                    {
                        "Cmp": {
                            "name": "Methane",
                            "Tc": 190.6,
                            "Pc": 4599000,
                            "omega": 0.011,
                            "A": 8.07131,
                            "B": 1730.63,
                            "C": 233.426
                        }
                    },
                    {
                        "Cmp": {
                            "name": "Ethane",
                            "Tc": 305.4,
                            "Pc": 4872000,
                            "omega": 0.099,
                            "A": 8.21201,
                            "B": 1652.57,
                            "C": 229.387
                        }
                    }
                ],
                "T": 300,
                "P": 5e6,
                "eos_type": "PR",
                "method": "FSOLVE",
                "calculate_enthalpy": false
            }
            </pre>

            <p>Tipičan odgovor:</p>

            <pre>
            {
                "V": 0.5,              # udjel pare (vapor fraction)
                "x": [0.50, 0.50],     # sastav tekuće faze
                "y": [0.50, 0.50],     # sastav plinske faze
                "method": "FSOLVE",    # korištena metoda
                "iteration": 7,        # broj iteracija
                "Zl": 0.82,            # kompresibilnost tekućine
                "Zv": 0.95             # kompresibilnost plina
            }
            </pre>

            <hr/>

            <h3>Compressor / Thermodynamics (/compressor_calc)</h3>

            <p>Za termodinamički proračun kompresora pošalji POST na <code>/compressor_calc</code> sa sljedećim JSON-om:</p>

            <pre>
            {
                "fractions": [          # molne frakcije smjese (mora ih biti koliko ima komponenti u gass_data)
                    0.85,               #CO2
                    0.0469,             #O2
                    0.058,              #N2
                    0.0447,             #Ar
                    0.0001,             #H2O
                    0.0001,             #NO
                    0.00005,            #SO2
                    0.00005,            #SO3
                    0.00005,            #CO
                    0.00005             #H2S
                ],
                "P1": 100000,                  # ulazni tlak (Pa)
                "P2": 500000,                  # izlazni tlak (Pa)
                "T1": 293.15,                  # ulazna temperatura (K)
                "mass_flow": 10.0,             # maseni protok (kg/s)
                "isentropic_efficiency": 0.8,  # izentropski stupanj djelovanja
                "polytropic_efficiency": 0.85, # polytropski stupanj djelovanja
                "polytropic_exponent": 1.3,    # polytropski eksponent (n)
                "full_report": 1               # 0, 1 ili 2
            }
            </pre>

            <p><strong>full_report:</strong></p>
            <ul>
                <li><code>0</code> ili bilo koja druga vrijednost osim 1 i 2 – vraća puni rječnik s realnim plinom (ključevi poput <code>"adiabatic_real"</code>, <code>"polytropic_real"</code> itd.).</li>
                <li><code>1</code> – vraća rječnik sa:
                    <ul>
                        <li><code>"ideal"</code> – rezultati idealnog plina (iz <code>calc_ideal_gas_sanity_check</code>)</li>
                        <li><code>"real"</code> – rezultati realnog plina (iz <code>calc_real_gas_thermo</code>)</li>
                    </ul>
                </li>
                <li><code>2</code> – vraća kratki string sažetak snaga i specifičnog rada.</li>
            </ul>

            <p><strong>Primjer odgovora za <code>full_report = 2</code>:</strong></p>

            <pre>
            "Adiabatic P: 5.123 MW, W: 180.345 kJ/kg | Polytropic P: 4.987 MW, W: 175.210 kJ/kg"
            </pre>

            <p><strong>Primjer odgovora za <code>full_report = 1</code> (skraćeno):</p>

            <pre>
            {
                "ideal": {
                    "adiabatic_ideal": { ... },
                    "polytropic_ideal": { ... }
                },
                "real": {
                    "adiabatic_real": { ... },
                    "polytropic_real": { ... }
                }
            }
            </pre>

        </body>
    </html>
    """
