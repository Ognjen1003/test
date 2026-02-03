import requests
import datetime

url1 = "http://127.0.0.1:8000/calculate_flow"
# url2 = "https://electroacoustic-junctional-perry.ngrok-free.dev/calculate_flow"

begin = datetime.datetime.now()

input_data1 = {
    "L": 40000,
    "d_in": 0.315925,
    "e": 0.0001,
    "p": 4100000,
    "T": 275.0,
    "qm": 23.75,
    "case": "case1",
    "visual": 0,

    "step_mode": "VIRTUAL_STEPS",
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
    "visual": 0,

    "step_mode": "BY_FITTINGS",
    "max_step_m": 5000.0,

    "fittings": [
        {"at_m": 10000.0, "K": 0.9, "kind": "elbow"},
        {"at_m": 17000.0, "K": 1.2, "kind": "elbow"}
        {"at_m": 21000, "K": 5.0, "kind": "valve"},
        {"at_m": 24000, "K": 1.2, "d_in_new": 0.250, "kind": "contraction"},
        {"at_m": 33000, "K": 0.3, "d_in_new": 0.315925, "kind": "expansion"}
    ],

    "viscosity_method": "AUTO"
}

response = requests.post(url1, json=input_data1)

if response.status_code == 200:
    result = response.json()  # sad je bolje json
    end = datetime.datetime.now()
    print(f"OK in {end - begin}")
    print(result)
else:
    print(f"Greska: {response.status_code}\n{response.text}")
