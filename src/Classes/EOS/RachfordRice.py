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
              method: MT.SolveMethod = MT.SolveMethod.FSOLVE, phase_detect: bool = False, 
              max_iter: int = 100, tol: float = 1e-6, k_ij: np.ndarray | None = None):
        

        if RachfordRice.z is None:
            RachfordRice.z = np.array([comp.fraction for comp in components])

        K = Calculations.wilson_K_vectorized(components, T, P)
        #print(f"[DEBUG] Inicijalni K: {K}")

        last_valid_solution = None
        eos = eos_class(components, T, P, k_ij=k_ij)
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
            

            mask = RachfordRice.z > 0.0 #potrebno ako zbog BIC-a stavim neke frakcije na nulu da ne mijenjam cijelu matricu
            if np.any(mask):
                K_min = np.min(K[mask])
                K_max = np.max(K[mask])
            else:
                K_min, K_max = np.min(K), np.max(K)  # fallback


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
                if not (0.0 < V < 1.0):
                    if phase_detect:
                        if V < 0:
                            code = -10
                            if eos.Zl is None:
                                eos.ensure_Z_for_phase(RachfordRice.z, MT.Phase.LIQUID)
                        else:
                            code = 5
                            if eos.Zv is None:
                                eos.ensure_Z_for_phase(RachfordRice.z, MT.Phase.VAPOR)

                        return (code, None, None, method, iteration + 1, eos.Zl, eos.Zv)

                    return (-1, None, None, method, iteration + 1, None, None)
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