import requests
import json


url = "http://127.0.0.1:8000/eos_calc"
#url = "https://test-production-f873.up.railway.app/eos_calc"

data = {
    "components": [
        {
            "name": "Methane",
            "Tc": 190.6,
            "Pc": 4.5992,
            "omega": 0.008,
            "A": 8.07131,
            "B": 1730.63,
            "C": 233.426
        },
        {
            "name": "Ethane",
            "Tc": 305.3,
            "Pc": 4.872,
            "omega": 0.1,
            "A": 8.21201,
            "B": 1652.57,
            "C": 229.387
        }
    ],
    "z": [0.5, 0.5],
    "T": 320.0,
    "P": 50.0,
    "eos_type": "SRK"
}


response = requests.post(url, json=data)

if response.status_code == 200:
    print("Proslo je")
    print("Response:", response.json())
else:
    print(f"Greska: {response.status_code}")
    print("Response:", response.text)