import requests
import json

url1 = (
    "http://127.0.0.1:8000/eos_calc"
)

url2 = (
    "https://test-production-f873.up.railway.app/eos_calc"
)




data = {
    "components": [
        {
            "Cmp": {
                "name": "Methane",
                "formula": "CH4",
                "Mw": 16.04,
                "Tc": 190.6,
                "Pc": 4.5992,
                "omega": 0.008,
                "fraction": 0.3
            }
        },
        {
            "Cmp": {
                "name": "Ethane",
                "formula": "C2H6",
                "Mw": 30.07,
                "Tc": 305.3,
                "Pc": 4.872,
                "omega": 0.1,
                "fraction": 0.3
            }
        },
        {
            "Cmp": {
                "name": "n-Butane",
                "formula": "C4H10",
                "Mw": 58.12,
                "Tc": 425.2,
                "Pc": 3.796,
                "omega": 0.2,
                "fraction": 0.4
            }
        },
    ],
    "T": 320.0,
    "P": 50.0,
    "eos_type": "PR",  # PR ili SRK
    "method": "FSOLVE", #F SOLVE ili ROOT_SCALAR
    "calculate_enthalpy": False

}

response = requests.post(url1, json=data)

if response.status_code == 200:
    print("Proslo je")
    print("Response:", response.json())
else:
    print(f"Greska: {response.status_code}")
    print("Response:", response.text)

