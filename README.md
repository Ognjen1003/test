commands.....


git clone https://github.com/your_username/railway-test.git
cd railway-test


virtualno okruzenje 
    python -m venv venv
    venv\Scripts\activate
    deactivate     

pip
    pip show _library_
    pip list
    pip install -r requirements.txt

pokretanje servera
    uvicorn src.main:app --reload

pokretanje iz testova lokalno
    set PYTHONPATH=%cd%       # da ide na src folder, absolute path
    python -m tests.local.test_main_local
    ----> ali koristi .vscode launch.json!!!

CMake 
    otvoriti x64 Native Tools Command Prompt
    rm -rf build
    mkdir build
    cd build
    cmake .. -G "NMake Makefiles" -DCMAKE_BUILD_TYPE=Release
    nmake
    ista verzija pythona kao i pyd....
    ------> ali sve to je u automat.bat-u u C++!!!


Svi importi su apsolutni, npr.:
    from src.Endpoints.EOSModul import perform_eos_calculation
    from src.Classes.Component import Component



---

Setup za razvoj (mozda...)

    `setup.py`:
    from setuptools import setup, find_packages
        setup(
            name="railway_test",
            version="0.1",
            packages=find_packages(where="src"),
            package_dir={"": "src"},
        )

ngrok 
url: https://electroacoustic-junctional-perry.ngrok-free.dev
config: ngrok config check
start: ngrok http 8000 

config yaml: usual places, refer code sa ngrok.com
    Linux: ~/.config/ngrok/ngrok.yml
    MacOS (Darwin): ~/Library/Application Support/ngrok/ngrok.yml
    Windows: "%HOMEPATH%\AppData\Local\ngrok\ngrok.yml"



