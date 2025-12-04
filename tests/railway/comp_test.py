import requests
import json

url1 = (
    "http://127.0.0.1:8000/compressor_calc"
)

url2 = (
    "https://electroacoustic-junctional-perry.ngrok-free.dev/compressor_calc"
)


# molar fractions CO2, O2, N2, Ar, H2O, NO, SO2, SO3, CO, H2S 
mass_flow = 10
data = {
    "fractions": [0.85, 0.0469, 0.058, 0.0447, 0.0001, 0.0001, 0.00005, 0.00005, 0.00005, 0.00005],   # molni udjeli smjese
    "T1": 300.0,                      # K
    "P1": 1,                          # bar
    "P2": 10,                         # bar
    "mass_flow": mass_flow,            # kg/s
    "isentropic_efficiency": 0.9,     # ηₛ
    "polytropic_efficiency": 0.9,     # ηₚ
    "NASA_9": "true",                 
    "full_report": 1,
    "polytropic_exponent": 1.1        # 1 = vrati puni report, 0 = sažetak
}

response = requests.post(url1, json=data)


def print_full_compression_results(res_all: dict, m_dot: float) -> None:
    """
    Lijep ispis rezultata iz odgovora API-ja, s jasnim razdvajanjem:
      - TERMODINAMIČKI MODEL (idealni plin / realni plin PR EOS)
      - STROJ (korekcija preko izentropskog / politropskog stupnja)

    Struktura res_all:
      res_all["ideal"]["adiabatic_ideal"]      -> T2s, T2, h1, h2s, h2, w_s, w_actual, P_s_MW, P_MW, p_ratio
      res_all["ideal"]["polytropic_ideal"]     -> T2, h1, h2, w, P_MW, p_ratio

      res_all["real"]["adiabatic_ideal"]       -> state1, state2s, state2, w_s, w_actual, P_s_MW, P_MW, p_ratio, eta_s
      res_all["real"]["polytropic_ideal"]      -> state1, state2, w, P_MW, p_ratio, eta_p, n_steps
    """

    def perc_diff(a, b):
        """Relativna razlika b u odnosu na a, u %."""
        return (b - a) / a * 100.0 if a != 0 else float("nan")

    # ==========================
    # 1) Parsiranje podataka
    # ==========================

    ideal_ad   = res_all["ideal"]["adiabatic_ideal"]
    ideal_poly = res_all["ideal"]["polytropic_ideal"]

    real_ad    = res_all["real"]["adiabatic_ideal"]
    real_poly  = res_all["real"]["polytropic_ideal"]

    # Realni adijabatski slučaj ima kompletna stanja -> koristimo ih za p1, p2, T1
    st1_r   = real_ad["state1"]
    st2s_r  = real_ad["state2s"]
    st2_r   = real_ad["state2"]

    st1_pr  = real_poly["state1"]
    st2_pr  = real_poly["state2"]

    p1 = st1_r["p"]
    p2 = st2_r["p"]
    T1 = st1_r["T"]

    # Efektivna izentropska učinkovitost iz idealnog i realnog (samo info)
    eta_s_ideal = ideal_ad["w_s"] / ideal_ad["w_actual"] if ideal_ad["w_actual"] != 0 else float("nan")
    eta_s_real  = real_ad["w_s"]  / real_ad["w_actual"]  if real_ad["w_actual"]  != 0 else float("nan")
    polytropic_exponent = ideal_ad["n"]

    # ==========================
    # 2) ADIJABATSKA – IDEALNI PLIN
    # ==========================

    print("\n====================================================")
    print(" ADIJABATSKA KOMPRESIJA – IDEALNI PLIN")
    print("====================================================")
    print(f"Ulazni tlak p1:                 {p1:.3f} bar")
    print(f"Izlazni tlak p2:                {p2:.3f} bar")
    print(f"Omjer tlakova p2/p1:            {ideal_ad['p_ratio']:.3f} [-]\n")

    print(f"Ulazna temperatura T1:          {T1:.2f} K ({T1 - 273.15:.2f} °C)\n")

    # --- A) TERMODINAMIČKI MODEL – idealni plin (η_s = 1.0) ---
    print("  [A] TERMODINAMIČKI MODEL – idealni plin (izentropski, η_s = 1.0)")
    print(f"      Izentropska izlazna T2s:  {ideal_ad['T2s']:.2f} K ({ideal_ad['T2s'] - 273.15:.2f} °C)")
    print("      Entalpije (h = cp·T, referentna nula proizvoljna):")
    print(f"        h1  (ulaz, idealno):    {ideal_ad['h1']:.3f} kJ/kg")
    print(f"        h2s (izentropski):      {ideal_ad['h2s']:.3f} kJ/kg")
    print(f"      Specifični rad idealnog izentropskog procesa:")
    print(f"        w_s,ideal (iz h2s - h1): {ideal_ad['w_s']:.3f} kJ/kg")
    print(f"      Snaga idealnog izentropskog procesa:")
    print(f"        P_s,ideal:              {ideal_ad['P_s_MW']:.4f} MW\n")

    # --- B) STROJ – korekcija na zadani izentropski stupanj (realni stroj, idealni model plina) ---
    print("  [B] STROJ – učinak izentropske učinkovitosti na idealni model plina")
    print("      (ovdje je TERMODINAMIKA još uvijek idealni plin,")
    print("       ali uvažava se realnost KOMESORA preko zadane η_s)\n")

    print(f"      Zadana/efektivna izentropska učinkovitost η_s: {eta_s_ideal:.3f}")
    print(f"      'Stvarna' izlazna temperatura T2 (iz η_s):      {ideal_ad['T2']:.2f} K ({ideal_ad['T2'] - 273.15:.2f} °C)")
    print("      Entalpija stvarnog izlaza (idealni plin + η_s):")
    print(f"        h2,stroj:               {ideal_ad['h2']:.3f} kJ/kg")
    print("      Specifični rad realnog stroja (na temelju idealnog plina):")
    print(f"        w_actual,ideal:         {ideal_ad['w_actual']:.3f} kJ/kg")
    print("      Snaga realnog stroja (na temelju idealnog plina):")
    print(f"        P_actual,ideal:         {ideal_ad['P_MW']:.4f} MW")
    print(f"      Maseni protok ṁ:         {m_dot:.3f} kg/s")
    print("====================================================\n")

    # ==========================
    # 3) ADIJABATSKA – REALNI PLIN (PR EOS)
    # ==========================

    print("====================================================")
    print(" ADIJABATSKA KOMPRESIJA – REALNI PLIN (PR EOS)")
    print("====================================================")
    print(f"Ulazni tlak p1:                 {st1_r['p']:.3f} bar")
    print(f"Izlazni tlak p2:                {st2_r['p']:.3f} bar")
    print(f"Omjer tlakova p2/p1:            {real_ad['p_ratio']:.3f} [-]\n")

    print(f"Ulazna temperatura T1:          {st1_r['T']:.2f} K ({st1_r['T'] - 273.15:.2f} °C)\n")

    # --- C) TERMODINAMIČKI MODEL – realni plin (EOS, η_s = 1.0) ---
    print("  [C] TERMODINAMIČKI MODEL – realni plin (PR EOS, izentropski, η_s = 1.0)")
    print(f"      Izentropska izlazna T2s:  {st2s_r['T']:.2f} K ({st2s_r['T'] - 273.15:.2f} °C)")
    print("      Entalpije realnog plina (apsolutna razina proizvoljna):")
    print(f"        h1,real  (ulaz):        {st1_r['h']:.3f} kJ/kg")
    print(f"        h2s,real (izentropski): {st2s_r['h']:.3f} kJ/kg")
    print("      Specifični rad idealnog izentropskog procesa realnog plina:")
    print(f"        w_s,real (iz h2s - h1): {real_ad['w_s']:.3f} kJ/kg")
    print("      Snaga idealnog izentropskog procesa (realni plin):")
    print(f"        P_s,real:               {real_ad['P_s_MW']:.4f} MW\n")

    # --- D) STROJ – učinak η_s na realni plin ---
    print("  [D] STROJ – učinak izentropske učinkovitosti na realni plin (PR EOS)")
    print(f"      Zadani izentropski stupanj kompresora η_s: {real_ad.get('eta_s', eta_s_real):.3f}")
    print(f"      Stvarna izlazna T2,real (EOS + η_s):       {st2_r['T']:.2f} K ({st2_r['T'] - 273.15:.2f} °C)")
    print("      Entalpija stvarnog izlaza (realni plin + η_s):")
    print(f"        h2,stroj,real:          {st2_r['h']:.3f} kJ/kg")
    print("      Specifični rad realnog stroja (na temelju realnog plina):")
    print(f"        w_actual,real:          {real_ad['w_actual']:.3f} kJ/kg")
    print("      Snaga realnog stroja (na temelju realnog plina):")
    print(f"        P_actual,real:          {real_ad['P_MW']:.4f} MW")
    print(f"      Maseni protok ṁ:         {m_dot:.3f} kg/s")
    print("====================================================\n")


    # ==========================
    # 4) POLITROPSKA – IDEALNI PLIN
    # ==========================

    print("====================================================")
    print(" POLITROPSKA KOMPRESIJA – IDEALNI PLIN")
    print("====================================================")
    print(f"Ulazni tlak p1:                 {st1_pr['p']:.3f} bar")
    print(f"Izlazni tlak p2:                {st2_pr['p']:.3f} bar")
    print(f"Omjer tlakova p2/p1:            {ideal_poly['p_ratio']:.3f} [-]\n")

    # Napomena: ovdje je već uračunata zadana politropska učinkovitost η_p (iz ulaza)
    eta_p = real_poly.get("eta_p", float("nan"))
    print("  [E] TERMODINAMIČKO–STROJNI MODEL – idealni plin + zadani η_p")
    print("      (kod politropske kompresije već imamo miješano: ")
    print("       termodinamika idealnog plina + učinkovitost procesa η_p)\n")

    print(f"      Polytropic exponent (zadan): {polytropic_exponent:.3f}")
    print(f"      Ulazna temperatura T1:                 {st1_pr['T']:.2f} K ({st1_pr['T'] - 273.15:.2f} °C)")
    print(f"      Izlazna temperatura T2 (idealni plin): {ideal_poly['T2']:.2f} K ({ideal_poly['T2'] - 273.15:.2f} °C)\n")

    print("      Entalpije (idealni plin):")
    print(f"        h1,ideal (ulaz):       {ideal_poly['h1']:.3f} kJ/kg")
    print(f"        h2,ideal (izlaz):      {ideal_poly['h2']:.3f} kJ/kg\n")

    print("      Specifični rad politropske kompresije – idealni plin:")
    print(f"        w,ideal:               {ideal_poly['w']:.3f} kJ/kg")
    print(f"      Snaga P,ideal:           {ideal_poly['P_MW']:.4f} MW")
    print(f"      Maseni protok ṁ:        {m_dot:.3f} kg/s")
    print("====================================================\n")

    # ==========================
    # 5) POLITROPSKA – REALNI PLIN (PR EOS)
    # ==========================

    print("====================================================")
    print(" POLITROPSKA KOMPRESIJA – REALNI PLIN (PR EOS)")
    print("====================================================")
    print(f"Ulazni tlak p1:                 {st1_pr['p']:.3f} bar")
    print(f"Izlazni tlak p2:                {st2_pr['p']:.3f} bar")
    print(f"Omjer tlakova p2/p1:            {real_poly['p_ratio']:.3f} [-]\n")

    print("  [F] TERMODINAMIČKO–STROJNI MODEL – realni plin (PR EOS) + zadani η_p")
    print("      (opet miješano: EOS za plin, ali proces je politropski sa zadanom η_p)\n")

    print(f"      Politropska učinkovitost η_p (zadana): {real_poly['eta_p']:.3f}")
    print(f"      Broj diskretizacijskih koraka:         {real_poly['n_steps']}\n")

    print(f"      Ulazna temperatura T1:                 {st1_pr['T']:.2f} K ({st1_pr['T'] - 273.15:.2f} °C)")
    print(f"      Izlazna temperatura T2 (realni plin):  {st2_pr['T']:.2f} K ({st2_pr['T'] - 273.15:.2f} °C)\n")

    print("      Entalpije (realni politropski slučaj):")
    print(f"        h1,real (ulaz):        {st1_pr['h']:.3f} kJ/kg")
    print(f"        h2,real (izlaz):       {st2_pr['h']:.3f} kJ/kg\n")

    print("      Specifični rad politropske kompresije – realni plin:")
    print(f"        w,real:                {real_poly['w']:.3f} kJ/kg")
    print(f"      Snaga P,real:            {real_poly['P_MW']:.4f} MW")
    print(f"      Maseni protok ṁ:        {m_dot:.3f} kg/s")
    print("====================================================\n")


    # Kratka usporedba adijabatske kompresije – strojevi (idealni model vs realni model plina)
    print("******** USPOREDBA STROJA – ADIJABATSKA, IDEALNI MODEL vs REALNI MODEL ********")
    print(" (ovdje uspoređujemo samo učinak stroja – w_actual i P – za dva različita")
    print("  termodinamička modela plina: idealni plin vs. realni plin PR EOS)\n")
    print(f"  w_actual (idealni plin): {ideal_ad['w_actual']:.3f} kJ/kg")
    print(f"  w_actual (realni plin):  {real_ad['w_actual']:.3f} kJ/kg  "
          f"(razlika {perc_diff(ideal_ad['w_actual'], real_ad['w_actual']):.2f} %)")
    print(f"  P_actual (idealni plin): {ideal_ad['P_MW']:.4f} MW")
    print(f"  P_actual (realni plin):  {real_ad['P_MW']:.4f} MW  "
          f"(razlika {perc_diff(ideal_ad['P_MW'], real_ad['P_MW']):.2f} %)")
    print("*********************************************************************\n")

    # Kratka usporedba politropske kompresije – idealni vs realni model plina
    print("******** USPOREDBA – POLITROPSKA, IDEALNI vs REALNI MODEL PLINA ********")
    print(" (ovdje su u oba slučaja uključeni i proces (η_p) i stroj, ali se razlikuje")
    print("  TERMODINAMIČKI model plina: idealni plin vs realni plin PR EOS)\n")
    print(f"  w,ideal: {ideal_poly['w']:.3f} kJ/kg")
    print(f"  w,real : {real_poly['w']:.3f} kJ/kg  "
          f"(razlika {perc_diff(ideal_poly['w'], real_poly['w']):.2f} %)")
    print(f"  P,ideal: {ideal_poly['P_MW']:.4f} MW")
    print(f"  P,real : {real_poly['P_MW']:.4f} MW  "
          f"(razlika {perc_diff(ideal_poly['P_MW'], real_poly['P_MW']):.2f} %)")
    print("*************************************************************************\n")







if response.status_code == 200:
    #print("Proslo je")
    #print("Response:", response.json())
    print_full_compression_results(response.json(), mass_flow)
else:
    print(f"Greska: {response.status_code}")
    print("Response:", response.text)





