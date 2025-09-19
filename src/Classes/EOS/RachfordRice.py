from typing import List
from scipy.optimize import root_scalar, fsolve
from src.Models.Component import Component
from src.Classes.EOS.EOSUtil import Calculations
import EnumsClasses.MethodsAndTypes as MT
import numpy as np

class RachfordRice:

    z = None

    @staticmethod
    def solve(eos_class, components: List[Component], T: float, P: float,
              method: MT.SolveMethod = MT.SolveMethod.FSOLVE, phase_detect: bool = False, max_iter: int = 100, tol: float = 1e-6):
        
        # if T == 250 and P == 41:
        #     print("unutra")

        if RachfordRice.z is None:
            RachfordRice.z = np.array([comp.fraction for comp in components])

        K = Calculations.wilson_K_vectorized(components, T, P)
        #print(f"[DEBUG] Inicijalni K: {K}")

        last_valid_solution = None
        eos = eos_class(components, T, P)
        for iteration in range(max_iter):

            ####################################   start osnovne provjere   #####################################    
            if np.all(K <= 1.0):
                 if phase_detect:
                    Zl = eos.ensure_Z_for_phase(RachfordRice.z, MT.Phase.LIQUID)
                    return  (-3.0, None, None, method, iteration + 1, Zl, None ) #tekuce
                 return  (-1.0, None, None, method, iteration + 1, None, None )
            if np.all(K >= 1.0):
                 if phase_detect:
                    Zv = eos.ensure_Z_for_phase(RachfordRice.z, MT.Phase.VAPOR)
                    return  (-9.0, None, None, method, iteration + 1, None, Zv ) #plinovito
                 return  (-1.0, None, None, method, iteration + 1, None, None ) 
            

            K_min, K_max = min(K), max(K)
            V_min = max(0.0, 1.0 / (1.0 - K_max))
            V_max = min(1.0, 1.0 / (1.0 - K_min))

            f_min = RachfordRice.rachford_rice_function(V_min, RachfordRice.z, K)
            f_max = RachfordRice.rachford_rice_function(V_max, RachfordRice.z, K)

            if V_min > V_max or f_min * f_max > 0:

                if not phase_detect:
                    return  (-1.0, None, None, method, iteration + 1, None, None )



                # Negative flash - Whitson ponovo
                #https://www.researchgate.net/publication/223072443_The_Negative_Flash
                phase_flag, V_nf = RachfordRice.classify_negative_flash(RachfordRice.z, K, V_guess=0.5)

                # V_nf<0 ili >1
                if phase_flag is not None:
                    if phase_flag == MT.Phase.LIQUID:
                        Zl = eos.ensure_Z_for_phase(RachfordRice.z, MT.Phase.LIQUID)
                        return (-2, None, None, method, iteration + 1, Zl, None)
                    else:
                        Zv = eos.ensure_Z_for_phase(RachfordRice.z, MT.Phase.VAPOR)
                        return (-8, None, None, method, iteration + 1, None, Zv)


                # Ako V_nf ∈ [0,1], to bi impliciralo dvofazno, ali do ovdje nema presjeka.
                # Fallback (heuristika preko K ili Z), K je gore na pocetku jer imam koristi od nje i inace
                # Posljednji fallback: heuristika preko jednofaznog Z (manji Z ≈ liquid-like)
                a_mix, b_mix, A, B = eos.get_mixture_parameters(RachfordRice.z)
                Z_single = eos.calc_Z_factor(A, B, phase='vapor')
                code = 2.0 if Z_single < 0.8 else 3.0
                return (code, None, None, method, iteration + 1, eos.Zl, eos.Zv)
                #ATTENTION!!!! Ovdje nikad vise ne ulazi!

            ####################################   end osnovne provjere   #####################################

        
            V = None

            if method == MT.SolveMethod.ROOT_SCALAR:
                result = root_scalar(RachfordRice.rachford_rice_function, args=(RachfordRice.z, K),
                                         method='bisect', bracket=[V_min, V_max], xtol=tol)
                if result.converged:
                    V = result.root
                    print(f"  V (root_scalar) = {V:.6f}")
                else:
                    print("[DEBUG] root_scalar nije konvergirao - fallback")
                    return last_valid_solution if last_valid_solution else (-1, RachfordRice.z, [K[i] * RachfordRice.z[i] for i in range(len(RachfordRice.z))], method, iteration + 1, eos.Zl, eos.Zv)
               
            elif method == MT.SolveMethod.FSOLVE:
                V = fsolve(lambda V: RachfordRice.rachford_rice_function(V, RachfordRice.z, K), 0.5)[0]
                if not (0.0 <= V <= 1.0):
                    if phase_detect:
                        if V < 0:
                            code = 2
                            if eos.Zl is None:
                                eos.ensure_Z_for_phase(RachfordRice.z, MT.Phase.LIQUID)
                        else:
                            code = 3
                            if eos.Zv is None:
                                eos.ensure_Z_for_phase(RachfordRice.z, MT.Phase.VAPOR)

                        return (code, None, None, method, iteration + 1, eos.Zl, eos.Zv)

                    return (-1.0, None, None, method, iteration + 1, None, None)
            else:
                raise ValueError(f"Nepoznata metoda: {method}")

            denom = 1 + V * (K - 1)
            x = RachfordRice.z / denom
            y = K * x


            phi_v = eos.fugacity_coeff(y, MT.Phase.VAPOR)
            phi_l = eos.fugacity_coeff(x, MT.Phase.LIQUID)

            #=================== prva verzija , ima jos tri za sada ovo ide ====================
            new_K = phi_l / phi_v

            last_valid_solution = (V, x, y, method, iteration + 1, eos.Zl, eos.Zv)

            if np.all(np.abs(new_K - 1.0) < 1e-3):
                return 0.5, x.tolist(), y.tolist(), method, iteration + 1, eos.Zl, eos.Zv

            if np.all(np.abs((new_K - K) / K) < tol):
                return V, x.tolist(), y.tolist(), method, iteration + 1, eos.Zl, eos.Zv # većina izlazi ovdje

            K = new_K
            #=================== pva verzija ====================


        #print("[DEBUG] Dosegnut max_iter - vraćam zadnje poznato rješenje") - tesko tu dolazi gotovo nikad sta je u redu
        return last_valid_solution 


    @staticmethod
    def solve_single_optimized(eos_class, components: List[Component], T: float, P: float, method: MT.SolveMethod = MT.SolveMethod.FSOLVE, 
                               max_iter: int = 100, tol: float = 1e-6,):

        if RachfordRice.z is None:
            RachfordRice.z = np.array([comp.fraction for comp in components], dtype=float)
        z = np.asarray(RachfordRice.z, dtype=float)

        K = np.asarray(Calculations.wilson_K_vectorized(components, T, P), dtype=float)

        eos = eos_class(components, T, P)

        if np.all(K <= 1.0):  # tekućina stabilna
            x = z.copy()
            eos.fugacity_coeff(x, 'liquid')   # postavlja eos.Zl
            
            return 0.0, x.tolist(), (K * x).tolist(), method, 0, eos.Zl, None

        if np.all(K >= 1.0):  # para stabilna
            y = z.copy()
            eos.fugacity_coeff(y, 'vapor')    # postavlja eos.Zv
            return 1.0, (y / np.maximum(K, 1e-30)).tolist(), y.tolist(), method, 0, None, eos.Zv

        last_valid_solution = None

        # --- glavni iterativni loop ---
        for it in range(1, max_iter + 1):

            # 1) odredi bracket za RR
            eps = 1e-12
            V_lo = 0.0 + eps
            V_hi = 1.0 - eps
            if np.any(K < 1.0):
                # singularnost: 1 + V*(K_i-1) = 0  => V = 1/(1-K_i) za K_i < 1 (pozitivno)
                V_hi = min(V_hi, float(np.min(1.0 / (1.0 - K[K < 1.0])) - eps))

            f_lo = rr(V_lo, z, K)
            f_hi = rr(V_hi, z, K)

            # 2) ako nema promjene znaka, tretiraj kao jednofazno (fallback heuristika)
            if not (np.isfinite(f_lo) and np.isfinite(f_hi)) or f_lo * f_hi > 0:
                if np.all(K <= 1.0):
                    x = z.copy()
                    eos.fugacity_coeff(x, 'liquid')
                    return 0.0, x.tolist(), (K * x).tolist(), method, it, eos.Zl, None
                if np.all(K >= 1.0):
                    y = z.copy()
                    eos.fugacity_coeff(y, 'vapor')
                    return 1.0, (y / np.maximum(K, 1e-30)).tolist(), y.tolist(), method, it, None, eos.Zv
                # ako je i dalje neodlučivo, tretiraj kao jednofaznu s manjim Z korijenom (često tekućina)
                x = z.copy()
                eos.fugacity_coeff(x, 'liquid')
                return 0.0, x.tolist(), (K * x).tolist(), method, it, eos.Zl, None

            # 3) riješi RR za V
            V = None
            try:
                if method == SolveMethod.ROOT_SCALAR:
                    from scipy.optimize import root_scalar
                    result = root_scalar(rr, args=(z, K), method='bisect', bracket=[V_lo, V_hi], xtol=tol)
                    if result.converged:
                        V = float(result.root)
                    else:
                        raise RuntimeError("root_scalar nije konvergirao")
                elif method == SolveMethod.FSOLVE:
                    from scipy.optimize import fsolve
                    V_guess = 0.5
                    V = float(fsolve(lambda V_: rr(V_, z, K), V_guess)[0])
                else:
                    raise ValueError(f"Nepoznata metoda: {method}")
            except Exception:
                # fallback na jednofazno prema K heuristici
                if np.all(K <= 1.0):
                    x = z.copy()
                    eos.fugacity_coeff(x, 'liquid')
                    return 0.0, x.tolist(), (K * x).tolist(), method, it, eos.Zl, None
                if np.all(K >= 1.0):
                    y = z.copy()
                    eos.fugacity_coeff(y, 'vapor')
                    return 1.0, (y / np.maximum(K, 1e-30)).tolist(), y.tolist(), method, it, None, eos.Zv
                # konzervativno:
                x = z.copy()
                eos.fugacity_coeff(x, 'liquid')
                return 0.0, x.tolist(), (K * x).tolist(), method, it, eos.Zl, None

            # ako V izvan [0,1], nije dvjefazno → jednofazno
            if not (0.0 <= V <= 1.0) or not np.isfinite(V):
                if V < 0.5:
                    x = z.copy()
                    eos.fugacity_coeff(x, 'liquid')
                    return 0.0, x.tolist(), (K * x).tolist(), method, it, eos.Zl, None
                else:
                    y = z.copy()
                    eos.fugacity_coeff(y, 'vapor')
                    return 1.0, (y / np.maximum(K, 1e-30)).tolist(), y.tolist(), method, it, None, eos.Zv

            # 4) izračun x,y
            den = 1.0 + V * (K - 1.0)
            if np.any(np.isclose(den, 0.0, atol=1e-12)):
                # lagano pomakni V prema unutra
                V = min(max(V, V_lo + 1e-9), V_hi - 1e-9)
                den = 1.0 + V * (K - 1.0)

            x = z / den
            y = K * x

            # negativne/naN vrijednosti -> fallback
            if np.any(~np.isfinite(x)) or np.any(~np.isfinite(y)) or np.any(x < 0) or np.any(y < 0):
                # probaj opet s V na sredini bracketa
                V = 0.5 * (V_lo + V_hi)
                den = 1.0 + V * (K - 1.0)
                x = z / den
                y = K * x

            x = normalize_safe(x)
            y = normalize_safe(y)

            # 5) fugacitetski koeficijenti i update K (damping u log-prostoru)
            phi_v = np.asarray(eos.fugacity_coeff(y, 'vapor'), dtype=float)
            phi_l = np.asarray(eos.fugacity_coeff(x, 'liquid'), dtype=float)

            # mjera ravnoteže fugaciteta
            lnF_gap = np.log(phi_l) - np.log(phi_v)   # treba → 0
            gap_inf = float(np.max(np.abs(lnF_gap)))

            # damping
            lam = 0.5
            new_K = K * np.exp(lam * lnF_gap)

            last_valid_solution = (float(V), x.tolist(), y.tolist(), method, it, eos.Zl, eos.Zv)

            # 6) kriteriji konvergencije
            # (a) fugacitetska ravnoteža
            tol_f = 1e-8
            if gap_inf < tol_f and 0.0 < V < 1.0:
                return float(V), x.tolist(), y.tolist(), method, it, eos.Zl, eos.Zv

            # (b) promjena K
            dlnK = np.max(np.abs(np.log(new_K / np.maximum(K, 1e-300))))
            if dlnK < tol and 0.0 <= V <= 1.0:
                return float(V), x.tolist(), y.tolist(), method, it, eos.Zl, eos.Zv

            K = new_K

        # --- max_iter dosegnut: vrati zadnje poznato ---
        return last_valid_solution


    @staticmethod
    def rachford_rice_function(V, z, K) -> float:
        denom = 1 + V * (K - 1)
        mask = np.abs(denom) > 1e-12
        return np.sum(z[mask] * (K[mask] - 1) / denom[mask])

    @staticmethod
    def rr_unrestricted(V, z, K):
        """RR funkcija bez ograničenja na V (može biti <0 ili >1)."""
        den = 1.0 + V * (K - 1.0)
        # zaštita od singularnosti: ako je koji nazivnik ~0, vrati veliku vrijednost s ispravnim predznakom
        if np.any(np.isclose(den, 0.0, atol=1e-12)):
            sgn = np.sign(np.sum(z * (K - 1.0)))
            return sgn * 1e30
        return float(np.sum(z * (K - 1.0) / den))

    @staticmethod
    def classify_negative_flash(z, K, V_guess=0.5):
        """
        Vrati (phase_flag, V_nf) koristeći negative-flash:
        - 'liquid-like' ako V_nf < 0
        - 'vapor-like'  ako V_nf > 1
        - None          ako V_nf ∈ [0,1] (tj. dvjefazno po ovom testu)
        """
        try:
            V_nf = float(fsolve(lambda V_: RachfordRice.rr_unrestricted(V_, z, K), V_guess, xtol=1e-12)[0])
            if not np.isfinite(V_nf):
                return None, V_nf
        except Exception:
            return None, np.nan

        if V_nf < 0.0:
            return MT.Phase.LIQUID, V_nf
        if V_nf > 1.0:
            return MT.Phase.VAPOR, V_nf
        return None, V_nf

    @staticmethod
    def rr(V: float, z: np.ndarray, K: np.ndarray) -> float:
        den = 1.0 + V * (K - 1.0)
        if np.any(np.isclose(den, 0.0, atol=1e-12)):
                # vrati "veliku" vrijednost s potpisom – da bi root-finder znao da je s druge strane
            return np.sign(np.sum(z * (K - 1.0))) * np.inf
        return float(np.sum(z * (K - 1.0) / den))
    
    @staticmethod
    def normalize_safe(x: np.ndarray) -> np.ndarray:
        s = x.sum()
        if s <= 0:
            return x
        return x / s
    

def saturation_pressure(eos_class, components, T,
                        kind='bubble', P_init=None, tol=1e-8, max_iter=50,
                        inner_K_iters=20, inner_tol=1e-8):
    """
    Vrati (psat, kind) gdje je kind ∈ {'bubble','dew'} za zadane T i z.
    Ako ne uspije (npr. superkritično), vrati (None, None).

    - Koristi tvoje: Calculations.wilson_K_vectorized, eos_class(...).fugacity_coeff(...)
    - K-update: K <- K * (phi_L / phi_V) (bubble) ili recipročni oblik (dew)
    - Vanjska petlja: Newton po tlaku na F(P) = S(P) - 1,
        s bubble S = sum(z*K), dew S = sum(z / K)
    """
    z = np.array([comp.fraction for comp in components])
    z = np.asarray(z, dtype=float)
    assert abs(z.sum() - 1.0) < 1e-10, "z mora biti normiran"

    # --- početna pretpostavka za tlak
    P = float(P_init) if P_init is not None else 1.0e6  # 10 bar kao default; slobodno promijeni

    # --- inicijalni K s Wilsonom (na početnom P)
    K = Calculations.wilson_K_vectorized(components, T, P).astype(float)
    # lagana zaštita
    K = np.clip(K, 1e-12, 1e12)

    def refine_K_and_S(P_local, K_guess):
        """
        Na fiksnom P_local (i T,z) odradi nekoliko SS K-update-a s fugacitetima,
        vrati (S, K_new) gdje je S=∑z*K (bubble) ili S=∑z/K (dew) koristeći konvergirane K.
        """
        K_loc = K_guess.copy()
        for _ in range(inner_K_iters):
            if kind == 'bubble':
                # incipient vapour: y ∝ z*K
                y = z * K_loc
                y_sum = y.sum()
                if y_sum == 0.0 or not np.isfinite(y_sum):
                    break
                y /= y_sum

                eos_mix_L = eos_class(components, T, P_local)  # za z kao (liquid-like) stranu
                eos_inc_V = eos_class(components, T, P_local)  # za y kao (vapor-like) stranu
                phi_L = eos_mix_L.fugacity_coeff(z, 'liquid')
                phi_V = eos_inc_V.fugacity_coeff(y, 'vapor')

                K_new = K_loc * (phi_L / phi_V)
                # kriterij promjene K
                if np.all(np.abs((K_new - K_loc) / np.maximum(K_loc, 1e-16)) < inner_tol):
                    K_loc = K_new
                    break
                K_loc = K_new

                S = float(np.sum(z * K_loc))  # bubble: ∑ z_i K_i  -> treba biti 1
            else:
                # kind == 'dew'
                x = z / np.maximum(K_loc, 1e-16)
                x_sum = x.sum()
                if x_sum == 0.0 or not np.isfinite(x_sum):
                    break
                x /= x_sum

                eos_inc_L = eos_class(components, T, P_local)  # za x kao (liquid-like)
                eos_mix_V = eos_class(components, T, P_local)  # za z kao (vapor-like)
                phi_L = eos_inc_L.fugacity_coeff(x, 'liquid')
                phi_V = eos_mix_V.fugacity_coeff(z, 'vapor')

                # za dew je zgodno ažurirati 1/K, ali ostanimo s K:
                K_new = K_loc * (phi_L / phi_V)  # ista forma radi; S i normalizacija to "poslože"
                if np.all(np.abs((K_new - K_loc) / np.maximum(K_loc, 1e-16)) < inner_tol):
                    K_loc = K_new
                    break
                K_loc = K_new

                S = float(np.sum(z / np.maximum(K_loc, 1e-16)))  # dew: ∑ z_i / K_i -> 1

        return S, np.clip(K_loc, 1e-12, 1e12)

    # --- Newton po tlaku: F(P)=S(P)-1
    for _ in range(max_iter):
        S, K = refine_K_and_S(P, K)
        F = S - 1.0
        if not np.isfinite(F):
            break
        if abs(F) < tol:
            return P, kind

        # numerička derivacija po P
        dP = max(1e-3 * P, 10.0)  # relativni ili min 10 Pa
        S2, _ = refine_K_and_S(P + dP, K)  # koristi trenutni K kao dobar guess
        F2 = S2 - 1.0
        dFdP = (F2 - F) / dP if np.isfinite(F2) else None
        if (dFdP is None) or (dFdP == 0.0) or (not np.isfinite(dFdP)):
            break

        P_new = P - F / dFdP
        if P_new <= 0 or not np.isfinite(P_new):
            break
        # damping ako želiš: P_new = 0.5*P_new + 0.5*P
        P = P_new

    # probaj alternativni tip (bubble<->dew)
    if kind == 'bubble':
        return saturation_pressure(eos_class, components, T, kind='dew',
                                   P_init=P_init, tol=tol, max_iter=max_iter,
                                   inner_K_iters=inner_K_iters, inner_tol=inner_tol)
    else:
        return (None, None)


def classify_single_phase(eos_class, components, T, P, z, P_hint=None):
    """
    Vraća 'liquid-like' ili 'vapor-like' za jednofazno stanje.
    Logika:
      - Izračunaj bubble psat_b i/ili dew psat_d na istoj T i z.
      - Ako psat_b postoji i P > psat_b -> liquid-like (iznad bubble linije).
      - Ako psat_d postoji i P < psat_d -> vapor-like (ispod dew linije).
      - Ako nijedan ne postoji (superkritično) -> default 'vapor-like'.
      - Ako oba postoje (u pravilu bubble < dew), izvan omotača biramo po odnosu P naspram oba.
    """
    z = np.asarray(z, dtype=float)

    psat_b, kind_b = saturation_pressure(eos_class, components, T, z, kind='bubble', P_init=P_hint or P)
    psat_d, kind_d = saturation_pressure(eos_class, components, T, z, kind='dew',    P_init=P_hint or P)

    # Ako oba None -> superkritično na toj T: u praksi vapor-like
    if psat_b is None and psat_d is None:
        return 'vapor-like', {'psat_bubble': None, 'psat_dew': None}

    # Ako bubble postoji i smo iznad njega, to je liquid-like
    if psat_b is not None and P > psat_b:
        return 'liquid-like', {'psat_bubble': psat_b, 'psat_dew': psat_d}

    # Ako dew postoji i smo ispod njega, to je vapor-like
    if psat_d is not None and P < psat_d:
        return 'vapor-like', {'psat_bubble': psat_b, 'psat_dew': psat_d}

    # Ako smo između (psat_b < P < psat_d) to je dvofazno – ali ti kažeš da već znaš da je jednofazno,
    # pa ovo vjerojatno nećeš pogoditi. Svejedno, dodaj fallback:
    # heuristika (Z-smjer) – jednofazni z:
    eos_single = eos_class(components, T, P)
    a_mix, b_mix, A, B = eos_single.get_mixture_parameters(z)
    Z_try = eos_single.calc_Z_factor(A, B, phase='vapor')  # jedini realni korijen
    return ('liquid-like' if Z_try < 0.8 else 'vapor-like',
            {'psat_bubble': psat_b, 'psat_dew': psat_d, 'Z': Z_try})