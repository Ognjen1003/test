import requests
import json

url1 = (
    "http://127.0.0.1:8000/compressor_calc"
)

url2 = (
    "https://test-production-f873.up.railway.app/compressor_calc"
)


# molar fractions CO2, O2, N2, Ar, H2O, NO, SO2, SO3, CO, H2S 

data = {
    "fractions": [0.85, 0.0469, 0.058, 0.0447, 0.0001, 0.0001, 0.00005, 0.00005, 0.00005, 0.00005],   # molni udjeli smjese
    "T1": 300.0,                      # K
    "P1": 1,                          # bar
    "P2": 10,                         # bar
    "mass_flow": 10,                  # kg/s
    "isentropic_efficiency": 0.9,    # ηₛ
    "polytropic_efficiency": 1.3,     # ηₚ
    "NASA_9": "true",                 
    "full_report": 1,
    "polytropic_exponent": 1.2                  # 1 = vrati puni report, 0 = sažetak
}

response = requests.post(url1, json=data)

if response.status_code == 200:
    print("Proslo je")
    print("Response:", response.json())
else:
    print(f"Greska: {response.status_code}")
    print("Response:", response.text)

