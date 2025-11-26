from typing import Callable
from data.testData import ComponentData
from src.Classes.EOS import PengRobinsonEOS 
from src.Classes.UtilClass import State






def build_state(p: float, T: float, eos: PengRobinsonEOS) -> State:
    """Iz tlaka, temperature i kompozicije izgradi kompletno stanje."""
    h = eos.h(p, T)
    s = eos.s(p, T)
    v = eos.v(p, T)
    return State(p=p, T=T, h=h, s=s, v=v, z=eos.components)


# ================================================================
# 3) Jednostavan 1D root-finder po temperaturi
#    (brutalno jednostavno, ali dovoljno za strukturu)
# ================================================================

def find_T_for_target(
    func_T: Callable[[float], float],
    T_min: float,
    T_max: float,
    tol: float = 1e-4,
    max_iter: int = 50,
) -> float:
    """
    Bisection za rješavanje func_T(T) = 0 u [T_min, T_max].
    Koristi se za:
      - s(p, T) - s_target = 0
      - h(p, T) - h_target = 0
    Pretpostavka: funkcija mijenja predznak na krajevima intervala.
    """
    f_min = func_T(T_min)
    f_max = func_T(T_max)

    if f_min * f_max > 0:
        raise ValueError("Root-finder: funkcija nema promjenu predznaka u zadanom intervalu.")

    for _ in range(max_iter):
        T_mid = 0.5 * (T_min + T_max)
        f_mid = func_T(T_mid)

        if abs(f_mid) < tol:
            return T_mid

        if f_min * f_mid < 0:
            T_max = T_mid
            f_max = f_mid
        else:
            T_min = T_mid
            f_min = f_mid

    # Ako nismo konvergirali, svejedno vraćamo sredinu
    return 0.5 * (T_min + T_max)


# Helperi specifični za EOS:

def solve_T_at_const_s(p: float, s_target: float, eos: PengRobinsonEOS,
                       T_min: float, T_max: float) -> float:
    """Nađi T tako da s(p, T, z) = s_target."""
    def f(T):
        return eos.s(p, T) - s_target
    return find_T_for_target(f, T_min, T_max)


def solve_T_at_const_h(p: float, h_target: float, eos: PengRobinsonEOS,
                       T_min: float, T_max: float) -> float:
    """Nađi T tako da h(p, T, z) = h_target."""
    def f(T):
        return eos.h(p, T) - h_target
    return find_T_for_target(f, T_min, T_max)


# ================================================================
# 4) Adijabatska kompresija (idealno izentropska + stvarna s η_s)
# ================================================================

def adiabatic_compression_real(
    p1: float,
    T1: float,
    p2: float,
    m_dot: float,
    eos: PengRobinsonEOS,
    eta_s: float = 1.0,
    T_bracket: tuple = (250.0, 1500.0),
) -> dict:
    """
    Adijabatska kompresija realne smjese (PR/GERG EOS).

    p1, p2   – tlak (npr. bar ili Pa, ali konzistentno u cijelom kodu)
    T1       – ulazna temperatura [K]
    z        – kompozicija smjese (molni ili maseni udjeli, konzistentno s EOS-om)
    m_dot    – maseni protok [kg/s]
    eos      – objekt koji implementira EoSModel (tvoj PR kod)
    eta_s    – izentropski stupanj djelovanja kompresora (1.0 = idealno)
    T_bracket – interval [T_min, T_max] unutar kojeg tražimo rješenja za T [K].

    Vraća dict s ključnim veličinama.
    """

    if not (0 < eta_s <= 1.0):
        raise ValueError("eta_s treba biti u intervalu (0, 1].")

    # 1) Početno stanje
    st1 = build_state(p1, T1, eos)

    # 2) Idealno izentropsko izlazno stanje: s2s = s1, p2 zadano
    s1 = st1.s
    T2s = solve_T_at_const_s(
        p=p2,
        s_target=s1,
        eos=eos,
        T_min=T_bracket[0],
        T_max=T_bracket[1],
    )
    st2s = build_state(p2, T2s, eos)

    # 3) Idealni izentropski specifični rad
    w_s = st2s.h - st1.h   # [kJ/kg]

    # 4) Stvarni izlaz s eta_s
    #    w_actual = w_s / eta_s → h2 = h1 + w_actual
    w_actual = w_s / eta_s
    h2 = st1.h + w_actual

    T2 = solve_T_at_const_h(
        p=p2,
        h_target=h2,
        eos=eos,
        T_min=T_bracket[0],
        T_max=T_bracket[1],
    )
    st2 = build_state(p2, T2, z, eos)

    # 5) Snaga
    P_s_kW = m_dot * w_s          # [kW]
    P_kW = m_dot * w_actual       # [kW]
    P_s_MW = P_s_kW / 1e3
    P_MW = P_kW / 1e3

    return {
        "state1": st1,
        "state2s": st2s,
        "state2": st2,
        "w_s": w_s,
        "w_actual": w_actual,
        "P_s_MW": P_s_MW,
        "P_MW": P_MW,
        "p_ratio": p2 / p1,
        "eta_s": eta_s,
    }


# ================================================================
# 5) Politropska kompresija realne smjese (diskretizacija)
# ================================================================

def polytropic_compression_real(
    p1: float,
    T1: float,
    p2: float,
    m_dot: float,
    eos: PengRobinsonEOS,
    eta_p: float,
    n_steps: int = 50,
    T_bracket: tuple = (250.0, 1500.0),
) -> dict:
    """
    Politropska kompresija realne smjese korištenjem diskretizacije po tlaku.

    Ideja:
      - Podijelimo raspon tlaka [p1, p2] na n_steps
      - U svakom koraku:
          * idealni izentropski skok pri Δp
          * korekcija na stvarni politropski skok preko η_p
      - w_total = h2 - h1 (zbroj Δh_actual)

    Parametri:
    p1, p2   – početni i krajnji tlak
    T1       – početna temperatura [K]
    z        – kompozicija smjese
    m_dot    – maseni protok [kg/s]
    eos      – PR/GERG model
    eta_p    – politropska učinkovitost kompresora (0–1)
    n_steps  – broj diskretnih koraka po tlaku
    T_bracket – interval za traženje T po koracima

    Vraća dict s ukupnim radom, snagom i početnim/završnim stanjem.
    """

    if not (0 < eta_p <= 1.0):
        raise ValueError("eta_p treba biti u intervalu (0, 1].")

    # 1) Početno stanje
    st = build_state(p1, T1, eos)
    st1 = st
    p_start = p1
    p_end = p2

    # Pretpostavimo p2 > p1 (kompresija). Ako nije, može se dodati provjera.
    dp = (p_end - p_start) / n_steps

    for i in range(n_steps):
        p_old = p_start + i * dp
        p_new = p_old + dp

        # 2) Idealni izentropski skok u ovom malom koraku:
        s_old = st.s
        T_iso = solve_T_at_const_s(
            p=p_new,
            s_target=s_old,
            eos=eos,
            T_min=T_bracket[0],
            T_max=T_bracket[1],
        )
        st_iso = build_state(p_new, T_iso, eos)

        Δh_s = st_iso.h - st.h    # idealni Δh u tom koraku

        # 3) Politropska korekcija sa η_p:
        #    Δh_actual = Δh_s / η_p
        Δh_actual = Δh_s / eta_p
        h_new = st.h + Δh_actual

        # 4) Nađi T_new iz uvjeta h(p_new, T_new, z) = h_new
        T_new = solve_T_at_const_h(
            p=p_new,
            h_target=h_new,
            eos=eos,
            T_min=T_bracket[0],
            T_max=T_bracket[1],
        )
        st = build_state(p_new, T_new, eos)

    st2 = st

    # 5) Ukupni specifični rad i snaga
    w_total = st2.h - st1.h   # [kJ/kg]
    P_kW = m_dot * w_total
    P_MW = P_kW / 1e3

    return {
        "state1": st1,
        "state2": st2,
        "w": w_total,
        "P_MW": P_MW,
        "p_ratio": p2 / p1,
        "eta_p": eta_p,
        "n_steps": n_steps,
    }


# ================================================================
# 6) Primjer kako bi se to koristilo (pseudo, jer EOS ti već imaš)
# ================================================================

if __name__ == "__main__":
    # Pretpostavimo da imaš klasu MyPengRobinson koja implementira EoSModel:
    #
    #   class MyPengRobinson:
    #       def h(self, p, T, z): ...
    #       def s(self, p, T, z): ...
    #       def v(self, p, T, z): ...
    #       def cp(self, p, T, z): ...
    #
    # eos = MyPengRobinson(...)
    #
    # Ovdje samo stavljam placeholder:

    components = ComponentData.oxyfuel_comp1
    P1 = 1.0      # npr. bar
    T1 = 300.0    # K
    p2 = 100.0    # bar
    m_dot = 10.0  # kg/s

    eos = PengRobinsonEOS(components, T1, P1)

    # Adijabatska kompresija (izentropski stupanj 0.8)
    res_ad = adiabatic_compression_real(P1, T1, p2, m_dot, eos, eta_s=0.8)

    # Politropska kompresija (politropska učinkovitost 0.75)
    res_poly = polytropic_compression_real(P1, T1, p2, m_dot, eos, eta_p=0.75)
