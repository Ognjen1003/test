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
            "AntoineA": 8.07131,
            "AntoineB": 1730.63,
            "AntoineC": 233.426,
            "CpA": 19.89,
            "CpB": 5.024e-2,
            "CpC": 1.269e-5,
            "CpD": -11.01e-9
        },
        {
            "name": "Ethane",
            "Tc": 305.3,
            "Pc": 4.872,
            "omega": 0.1,
            "AntoineA": 8.21201,
            "AntoineB": 1652.57,
            "AntoineC": 229.387,
            "CpA": 21.13,
            "CpB": 78.60e-3,
            "CpC": -11.85e-6,
            "CpD": 18.99e-9,
        },
        {
            "name": "n-Butane",
            "Tc": 425.2,
            "Pc": 3.796,
            "omega": 0.2,
            "AntoineA": 6.80896,
            "AntoineB": 935.86,
            "AntoineC": 238.73,
            "CpA": 24.86,
            "CpB": 133.4e-3,
            "CpC": -30.1e-6,
            "CpD": 34.1e-9,
        }
    ],
    "z": [0.3, 0.3, 0.4],
    "T": 320.0,
    "P": 50.0,
    "eos_type": "PG"
}


response = requests.post(url, json=data)

if response.status_code == 200:
    print("Proslo je")
    print("Response:", response.json())
else:
    print(f"Greska: {response.status_code}")
    print("Response:", response.text)