import requests
import datetime

url = "http://127.0.0.1:8000/calculate_flow"
#url = "https://test-production-f873.up.railway.app/calculate_flow"


begin = datetime.datetime.now()

input_data = {
    "nsteps": 1,
    "L": 40000, 
    "d_in": 0.315925,
    "e": 0.0001,
    "p": 4000000,
    "T": 293.15,
    "qm": 23.75,
    "case": "case1",   
    "visual": 0
}

response = requests.post(url, json=input_data)

if response.status_code == 200:
    #result = response.json()
    result = response.text
    end = datetime.datetime.now()
    runtime = end - begin
    print(f"{result}  in {runtime} seconds")
else:
    print(f"Greska: {response.status_code}, {response.text}")