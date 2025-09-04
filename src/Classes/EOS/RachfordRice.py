from typing import List
from scipy.optimize import root_scalar, fsolve
from src.Models.Component import Component
from src.Classes.EOS.EOSUtil import Calculations
from src.EnumsClasses import SolveMethod
import numpy as np

class RachfordRice:

    z = None

    @staticmethod
    def solve(eos_class, components: List[Component], T: float, P: float,
              method: SolveMethod = SolveMethod.FSOLVE, max_iter: int = 100, tol: float = 1e-6):
        
        if RachfordRice.z is None:
            RachfordRice.z = np.array([comp.fraction for comp in components])

        K = Calculations.wilson_K_vectorized(components, T, P)
        #print(f"[DEBUG] Inicijalni K: {K}")

        last_valid_solution = None
        eos = eos_class(components, T, P)

        for iteration in range(max_iter):

            ####################################   start osnovne provjere   #####################################    
            """if all(Ki < 1.0 for Ki in K):
                 print("[DEBUG] Svi K < 1 → tekuća faza")
                 return 0.0, z, [K[i] * z[i] for i in range(len(z))], method, 0
            if all(Ki > 1.0 for Ki in K):
                 print("[DEBUG] Svi K > 1 → parna faza")
                 return 1.0, [z[i] / K[i] for i in range(len(z))], z, method, 0 
            """

            K_min, K_max = min(K), max(K)
            V_min = max(0.0, 1.0 / (1.0 - K_max))
            V_max = min(1.0, 1.0 / (1.0 - K_min))
            #print(f"  V_min = {V_min:.6f}, V_max = {V_max:.6f}")

            f_min = RachfordRice.rachford_rice_function(V_min, RachfordRice.z, K)
            f_max = RachfordRice.rachford_rice_function(V_max, RachfordRice.z, K)
            #print(f"  f(V_min) = {f_min:.6f}, f(V_max) = {f_max:.6f}")

            if V_min > V_max or f_min * f_max > 0:
                #print("[DEBUG] Nema valjanog presjeka - fallback")
                return  (-1.0, RachfordRice.z, [K[i] * RachfordRice.z[i] for i in range(len(RachfordRice.z))], method, iteration + 1, eos.Zl, eos.Zv )
                # -1 uglavnom odavde izlazi

            ####################################   end osnovne provjere   #####################################

            
            #print(f"[DEBUG] Iteracija {iteration + 1}")
            V = None

            if method == SolveMethod.ROOT_SCALAR:
                result = root_scalar(RachfordRice.rachford_rice_function, args=(RachfordRice.z, K),
                                         method='bisect', bracket=[V_min, V_max], xtol=tol)
                if result.converged:
                    V = result.root
                    print(f"  V (root_scalar) = {V:.6f}")
                else:
                    print("[DEBUG] root_scalar nije konvergirao - fallback")
                    return last_valid_solution if last_valid_solution else (0.0, RachfordRice.z, [K[i] * RachfordRice.z[i] for i in range(len(RachfordRice.z))], method, iteration + 1, eos.Zl, eos.Zv)
               
            elif method == SolveMethod.FSOLVE:
                V = fsolve(lambda V: RachfordRice.rachford_rice_function(V, RachfordRice.z, K), 0.5)[0]
                #print(f"  V (fsolve) = {V:.6f}")
                if not (0.0 <= V <= 1.0):
                    #print("[DEBUG] fsolve dao nerealno V - fallback")
                    return (-1.0, RachfordRice.z, [K[i] * RachfordRice.z[i] for i in range(len(RachfordRice.z))], method, iteration + 1, eos.Zl, eos.Zv)
            else:
                raise ValueError(f"Nepoznata metoda: {method}")

            denom = 1 + V * (K - 1)
            x = RachfordRice.z / denom
            y = K * x


            phi_v = eos.fugacity_coeff(y, 'vapor')
            phi_l = eos.fugacity_coeff(x, 'liquid')

            
            new_K = phi_l / phi_v
            #print(f"  Novi K: {new_K}")

            last_valid_solution = (V, x, y, method, iteration + 1, eos.Zl, eos.Zv)

            #if V < 1e-3:
            #    print("[DEBUG] V vrlo mali - tekuća faza")
            #   return 0.0, z, [K[i] * z[i] for i in range(len(z))], method, iteration + 1
            if np.all(np.abs(new_K - 1.0) < 1e-3):
                return 0.5, x.tolist(), y.tolist(), method, iteration + 1, eos.Zl, eos.Zv

            if np.all(np.abs((new_K - K) / K) < tol):
                return V, x.tolist(), y.tolist(), method, iteration + 1, eos.Zl, eos.Zv # većina izlazi ovdje

            K = new_K

        #print("[DEBUG] Dosegnut max_iter - vraćam zadnje poznato rješenje") - tesko tu dolazi
        return last_valid_solution 


    @staticmethod
    def solve_single_optimized(eos_class, components: List[Component], T: float, P: float, method: SolveMethod = SolveMethod.FSOLVE, 
                               max_iter: int = 100, tol: float = 1e-6,):

        if RachfordRice.z is None:
            RachfordRice.z = np.array([comp.fraction for comp in components], dtype=float)
        z = np.asarray(RachfordRice.z, dtype=float)

        K = np.asarray(Calculations.wilson_K_vectorized(components, T, P), dtype=float)

        eos = eos_class(components, T, P)

        # --- helperi ---
        def rr(V: float, z: np.ndarray, K: np.ndarray) -> float:
            den = 1.0 + V * (K - 1.0)
            if np.any(np.isclose(den, 0.0, atol=1e-12)):
                # vrati "veliku" vrijednost s potpisom – da bi root-finder znao da je s druge strane
                return np.sign(np.sum(z * (K - 1.0))) * np.inf
            return float(np.sum(z * (K - 1.0) / den))

        def normalize_safe(x: np.ndarray) -> np.ndarray:
            s = x.sum()
            if s <= 0:
                return x
            return x / s

        # --- brza jednofazna heuristika prije RR ---
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
    def rachford_rice_function(V, z, K):
        denom = 1 + V * (K - 1)
        mask = np.abs(denom) > 1e-12
        return np.sum(z[mask] * (K[mask] - 1) / denom[mask])
