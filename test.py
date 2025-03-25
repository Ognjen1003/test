import requests

# Define the URL of the FastAPI endpoint
url = "http://127.0.0.1:8000/calculate"

# Create a dictionary for the input data
input_data = {
    "nsteps": 40,
    "length": 40000 # example values for nsteps and length
}

# Send a POST request to the FastAPI server
response = requests.post(url, json=input_data)

# Check if the request was successful
if response.status_code == 200:
    result = response.json()
    print(f"Calculation result: {result}")
else:
    print(f"Error: {response.status_code}, {response.text}")