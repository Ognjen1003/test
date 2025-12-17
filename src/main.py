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
            method=input_data.method
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
      <body style="font-family: Arial, sans-serif; line-height: 1.4;">
        <h2>Flow loop web servis - API dokumentacija</h2>
        <p>Servis izlaže tri POST endpointa:</p>
        <ul>
          <li><code>/calculate_flow</code> - hidraulički proračun (pad tlaka po koracima)</li>
          <li><code>/eos_calc</code> - EOS/flash proračun (V, x, y, Z-faktori)</li>
          <li><code>/compressor_calc</code> - termodinamika kompresora (ideal/real, adijabatski/politropski)</li>
        </ul>

        <hr/>

        <h3>1) Flow calculation (<code>/calculate_flow</code>)</h3>
        <p>Računa pad tlaka po diskretiziranim koracima koristeći Darcy-Weisbach uz faktor trenja iz Colebrook-White.</p>

        <h4>Request (JSON)</h4>
        <pre style="background:#f6f6f6;padding:12px;border-radius:8px;">
{
  "nsteps": 1,
  "L": 40000,
  "d_in": 0.315925,
  "e": 0.0001,
  "p": 4000000,
  "T": 293.15,
  "qm": 23.75,
  "case": "OXY3",
  "visual": 0
}
        </pre>

        <p><strong>Opis polja (sažeto):</strong></p>
        <ul>
          <li><code>nsteps</code> - broj diskretnih koraka (npr. 1…50)</li>
          <li><code>L</code> [m] - duljina cjevovoda</li>
          <li><code>d_in</code> [m] - unutarnji promjer</li>
          <li><code>e</code> [m] - hrapavost cijevi</li>
          <li><code>p</code> [Pa] - ulazni tlak</li>
          <li><code>T</code> [K] - ulazna temperatura</li>
          <li><code>qm</code> [kg/s] - maseni protok</li>
          <li><code>case</code> - <code>CO2</code>, <code>CASE1</code>, <code>OXY1</code>, <code>OXY2</code>, <code>OXY3</code></li>
          <li><code>visual</code> - 0: JSON, 1: HTML tablica</li>
        </ul>

        <p><strong>Napomena:</strong> odgovor trenutno vraća <code>result</code> kao JSON-string ili HTML tablicu.</p>

        <h4>Response (primjer - cca, stvarni rezultat)</h4>
        <pre style="background:#f6f6f6;padding:12px;border-radius:8px;">
{
  "result": [
    {
      "step": 1.0,
      "L": 40000.0,
      "p1": 40.0,
      "t": 20.0,
      "mu": 8.106801262939821e-06,
      "rho_g": 91.74428392217587,
      "u": 3.3023775937632505,
      "Re": 11807017.636280673,
      "ff": 0.015173228688780701,
      "dp": 9.610732251173848,
      "p2": 30.38926774882615
    }
  ]
}
        </pre>

        <p><em>Izlaz:</em> <code>dp</code> je pad tlaka [bar], <code>p1</code>/<code>p2</code> su tlakovi [bar], <code>t</code> je temperatura [°C].</p>

        <hr/>

        <h3>2) EOS / flash proračun (<code>/eos_calc</code>)</h3>
        <p>Flash proračun smjese korištenjem Peng-Robinson (PR) ili SRK jednadžbe stanja. Vraća udio pare (V), fazne sastave (x,y) i Z-faktore.</p>

        <h4>Request (JSON - primjer)</h4>
        <p>Svaka komponenta ide kao objekt <code>{"Cmp": {...}}</code>. NASA-9 polja su opcionalna, ali se koriste za idealna svojstva u <code>Component</code> i real-gas termodinamici (npr. residual h u PR).</p>

        <pre style="background:#f6f6f6;padding:12px;border-radius:8px;">
{
  "components": [
    { "Cmp": { "name": "CO2", "formula": "CO2", "Mw": 44.01, "Tc": 304.7, "Pc": 73.866, "omega": 0.225, "fraction": 0.85,
               "nasa_mid": 1000.0, "nasa_low": [ ...9... ], "nasa_high": [ ...9... ] } },
    { "Cmp": { "name": "O2",  "formula": "O2",  "Mw": 32.0,  "Tc": 154.6, "Pc": 50.43,  "omega": 0.0222,"fraction": 0.0469,
               "nasa_mid": 1000.0, "nasa_low": [ ...9... ], "nasa_high": [ ...9... ] } }
    // ...
  ],
  "T": 245.15,
  "P": 40.0,
  "eos_type": "PR",
  "method": "FSOLVE"
}
        </pre>

        <h4>Response (primjer - cca, stvarni rezultat)</h4>
        <pre style="background:#f6f6f6;padding:12px;border-radius:8px;">
{
  "V": 0.17858923988745115,
  "x": [0.9336537899, 0.0241024354, 0.0192631432, 0.0226127973, 0.0001215454,
        0.0000488753, 0.0000603434, 0.0000607323, 0.0000180098, 0.0000583279],
  "y": [0.4652392876, 0.1517560644, 0.2361679063, 0.1462887964, 0.0000009031,
        0.0003351451, 0.0000024259, 0.0000006375, 0.0001971371, 0.0000116965],
  "method": "FSOLVE",
  "iteration": 8,
  "Zl": 0.07987173641093182,
  "Zv": 0.8000073759944746
}
        </pre>

        <p><em>Napomena:</em> <code>V</code> ~ 0.18 znači oko 18% pare i 82% tekuće faze; <code>x</code>/<code>y</code> su sastavi tekuće/parne faze; <code>Zl</code>/<code>Zv</code> su kompresibilnosti faza.</p>
        <p><em>Napomena:</em> Frakcije moraju biti normalizirane.</p>
        <hr/>

        <h3>3) Compressor / thermodynamics (<code>/compressor_calc</code>)</h3>
        <p>Računa adijabatsku i politropsku kompresiju oxyfuel smjese uz zadane učinkovitosti. Podržava “full report” i sažeti izlaz.</p>

        <h4>Request (JSON - primjer)</h4>
        <pre style="background:#f6f6f6;padding:12px;border-radius:8px;">
{
  "fractions": [0.85, 0.0469, 0.058, 0.0447, 0.0001, 0.0001, 0.00005, 0.00005, 0.00005, 0.00005],
  "T1": 300.0,
  "P1": 1.0,
  "P2": 10.0,
  "mass_flow": 10.0,
  "isentropic_efficiency": 0.8,
  "polytropic_efficiency": 0.9,
  "NASA_9": "true",
  "full_report": 1,
  "polytropic_exponent": 1.25
}
        </pre>

        <p><strong>Važno:</strong> <code>fractions</code> mora imati isti broj elemenata i isti redoslijed kao interna lista komponenti: <code>CO2, O2, N2, Ar, H2O, NO, SO2, SO3, CO, H2S </code>.</p>

        <h4>Response (primjeri - cca)</h4>

        <p><strong>A) <code>full_report = 2</code> (sažeti string):</strong></p>
        <pre style="background:#f6f6f6;padding:12px;border-radius:8px;">
"Adiabatic P: ~X.XXX MW, W: ~YYY.YYY kJ/kg | Polytropic P: ~X.XXX MW, W: ~YYY.YYY kJ/kg"
        </pre>

        <p><strong>B) <code>full_report = 0</code> (real-gas dict - skraćeno):</strong></p>
        <pre style="background:#f6f6f6;padding:12px;border-radius:8px;">
{
  "adiabatic_real": {
    "p_ratio": 10.0,
    "eta_s": 0.8,
    "w_s": "... kJ/kg",
    "w_actual": "... kJ/kg",
    "P_s_MW": "... MW",
    "P_MW": "... MW",
    "state1": { "p": 1.0,  "T": 300.0, "h": "...", "s": "...", "v": "...", "z": [ ... ] },
    "state2s":{ "p": 10.0, "T": "...", "h": "...", "s": "...", "v": "...", "z": [ ... ] },
    "state2": { "p": 10.0, "T": "...", "h": "...", "s": "...", "v": "...", "z": [ ... ] }
  },
  "polytropic_real": {
    "p_ratio": 10.0,
    "eta_p": 0.9,
    "n_steps": "...",
    "w": "... kJ/kg",
    "P_MW": "... MW",
    "state1": { ... },
    "state2": { ... }
  }
}
        </pre>

        <p><strong>C) <code>full_report = 1</code> (ideal + real):</strong></p>
        <pre style="background:#f6f6f6;padding:12px;border-radius:8px;">
{
  "ideal": { "adiabatic_ideal": { ... }, "polytropic_ideal": { ... } },
  "real":  { "adiabatic_real":  { ... }, "polytropic_real":  { ... } }
}
        </pre>

        <hr/>

        <h3>Općenito - greške</h3>
        <ul>
          <li>HTTP 400 - neispravan ulaz (npr. krivi <code>case</code>, kriva duljina <code>fractions</code>, zbroj udjela ≠ 1.0, itd.)</li>
          <li>HTTP 500 - neočekivana greška u izvođenju (numerički problem ili bug)</li>
        </ul>

      </body>
    </html>
    """
