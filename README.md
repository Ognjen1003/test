commands.....


git clone https://github.com/your_username/railway-test.git
cd railway-test


virtualno okruzenje 
    python -m venv venv
    venv\Scripts\activate
    deactivate     

pip
    pip show library
    pip list
    pip install -r requirements.txt

pokretanje servera
    uvicorn src.main:app --reload

pokretanje iz testova lokalno
    set PYTHONPATH=%cd%       # da ide na src folder, absolute path
    python -m tests.local.test_main_local


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


