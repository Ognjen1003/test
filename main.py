from fastapi import FastAPI
from datetime import datetime

app = FastAPI()

# Endpoint koji vraća trenutno vreme i datum
@app.get("/time")
def get_current_time():
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return {"Ovo je FastAPI prvi service point - vrijeme": current_time}

# Endpoint koji prima parametar 'ime' i vraća pozdravnu poruku
@app.get("/hello/{name}")
def say_hello(name: str):
    return {"message": f"Hello, {name}!"}