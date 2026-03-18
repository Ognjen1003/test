import requests
import datetime
import pandas as pd


url1 = "http://127.0.0.1:8000/calculate_flow2"
# url2 = "https://electroacoustic-junctional-perry.ngrok-free.dev/calculate_flow"

def print_res(result):
    cols_round = ["p1", "p2", "dp_fric", "dp_minor", "dp_total", "u", "Re", "ff", "rho_g", "mu"]
    df = pd.DataFrame(result)  # result = ona lista dictova

    for c in cols_round:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.round({
    "L": 0,
    "p1": 3,
    "p2": 3,
    "dp_fric": 3,
    "dp_minor": 3,
    "dp_total": 3,
    "u": 3,
    "Re": 0,
    "ff": 6,
    "rho_g": 3,
    "mu": 10,
    "K_step": 3 })

    print(df.to_string(index=False))


begin = datetime.datetime.now()

input_data1 = {
    "L": 40000,
    "d_in": 0.315925,
    "e": 0.0001,
    "p": 4100000,
    "T": 275.0,
    "qm": 23.75,
    "case": "case1",
    "visual": 1,

    "step_mode": "BY_FITTINGS",
    "virtual_steps": 8,

    "fittings": [],
    "viscosity_method": "AUTO"
}

input_data2 = {
    "L": 40000,
    "d_in": 0.315925,
    "e": 0.0001,
    "p": 4100000,
    "T": 275.0,
    "qm": 23.75,
    "case": "case1",
    "visual": 1,

    "step_mode": "BY_FITTINGS",
    "max_step_m": 5000.0,

    "fittings": [
        {"at_m": 10000.0, "K": 0.9, "kind": "elbow"},
        {"at_m": 17000.0, "K": 1.2, "kind": "elbow"}
    ],

    "viscosity_method": "AUTO"
}

response = requests.post(url1, json=input_data2)

if response.status_code == 200:
    if input_data2["visual"] == 0:
        result = response.json()
        print(result)
    else:
        html = response.text
        print(html)
else:
    print(f"Greska: {response.status_code}\n{response.text}")
