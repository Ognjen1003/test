import requests
import datetime


url1 = "http://127.0.0.1:8000/calculate_flow"
url2 = "https://electroacoustic-junctional-perry.ngrok-free.dev/calculate_flow"


begin = datetime.datetime.now()

input_data = {
    "nsteps": 1,
    "L": 40000, 
    "d_in": 0.315925,
    "e": 0.0001,
    "p": 4000000,
    "T": 293.15,
    "qm": 23.75,
    "case": "oxy3",   
    "visual": 0
}

response = requests.post(url1, json=input_data)

if response.status_code == 200:
    #result = response.json()
    result = response.text
    end = datetime.datetime.now()
    runtime = end - begin
    print(f"{result}  in {runtime} seconds")
else:
    print(f"Greska: {response.status_code}, {response.text}")