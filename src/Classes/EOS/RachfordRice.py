from typing import List, Tuple
from scipy.optimize import root_scalar, fsolve
from src.Classes.Component import Component
from src.Classes.EOS.EOSUtil import Calculations
from src.EnumsClasses import SolveMethod

class RachfordRice:

    @staticmethod
    def solve(eos_class, components: List[Component], T: float, P: float,
              method: SolveMethod = SolveMethod.FSOLVE, max_iter: int = 100, tol: float = 1e-6):
        
        z = [comp.fraction for comp in components]
        K = [Calculations.wilson_K(comp, T, P) for comp in components]
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

            f_min = RachfordRice.rachford_rice_function(V_min, z, K)
            f_max = RachfordRice.rachford_rice_function(V_max, z, K)
            #print(f"  f(V_min) = {f_min:.6f}, f(V_max) = {f_max:.6f}")

            if V_min > V_max or f_min * f_max > 0:
                #print("[DEBUG] Nema valjanog presjeka - fallback")
                return  (-1.0, z, [K[i] * z[i] for i in range(len(z))], method, iteration + 1, eos.Zl, eos.Zv )

            ####################################   end osnovne provjere   #####################################


            #print(f"[DEBUG] Iteracija {iteration + 1}")
            V = None

            if method == SolveMethod.ROOT_SCALAR:
                result = root_scalar(RachfordRice.rachford_rice_function, args=(z, K),
                                         method='bisect', bracket=[V_min, V_max], xtol=tol)
                if result.converged:
                    V = result.root
                    print(f"  V (root_scalar) = {V:.6f}")
                else:
                    print("[DEBUG] root_scalar nije konvergirao - fallback")
                    return last_valid_solution if last_valid_solution else (0.0, z, [K[i] * z[i] for i in range(len(z))], method, iteration + 1, eos.Zl, eos.Zv)
               
            elif method == SolveMethod.FSOLVE:
                V = fsolve(lambda V: RachfordRice.rachford_rice_function(V, z, K), 0.5)[0]
                #print(f"  V (fsolve) = {V:.6f}")
                if not (0.0 <= V <= 1.0):
                    #print("[DEBUG] fsolve dao nerealno V - fallback")
                    return (-1.0, z, [K[i] * z[i] for i in range(len(z))], method, iteration + 1, eos.Zl, eos.Zv)
            else:
                raise ValueError(f"Nepoznata metoda: {method}")

            x = [z[i] / (1 + V * (K[i] - 1)) for i in range(len(z))]
            y = [K[i] * x[i] for i in range(len(z))]


            phi_v = eos.fugacity_coeff(y, 'vapor')
            phi_l = eos.fugacity_coeff(x, 'liquid')

            new_K = [phi_l[i] / phi_v[i] for i in range(len(z))]
            #print(f"  Novi K: {new_K}")

            last_valid_solution = (V, x, y, method, iteration + 1, eos.Zl, eos.Zv)

            #if V < 1e-3:
            #    print("[DEBUG] V vrlo mali - tekuća faza")
            #   return 0.0, z, [K[i] * z[i] for i in range(len(z))], method, iteration + 1
            if all(abs(Ki - 1.0) < 1e-3 for Ki in new_K):
                #print("[DEBUG] Svi K ≈ 1 - ravnoteža")
                return 0.5, z, z, method, iteration + 1, eos.Zl, eos.Zv
            if all(abs((new_K[i] - K[i]) / K[i]) < tol for i in range(len(K))):
                #print("[DEBUG] K konvergirao")
                return V, x, y, method, iteration + 1, eos.Zl, eos.Zv

            K = new_K

        #print("[DEBUG] Dosegnut max_iter - vraćam zadnje poznato rješenje")
        return last_valid_solution 

    @staticmethod
    def rachford_rice_function(V, z, K):
        eps = 1e-12  # tolerancija da izbjegnemo dijeljenje s nulom
        return sum(
            z[i] * (K[i] - 1) / (1 + V * (K[i] - 1)) 
            for i in range(len(z)) 
            if abs(1 + V * (K[i] - 1)) > eps  # samo ako nije blizu nule
        )
