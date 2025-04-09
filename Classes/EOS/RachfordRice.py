from Classes.Component import Component
from Classes.EOS.EOSUtil import Calculations
from typing import List, Tuple

class RachfordRice:
    """
    Iterativna RR metoda s ažuriranjem K iz fugaciteta.
    """
    @staticmethod
    def solve(z: List[float], eos_class, components: List[Component], T: float, P: float, max_iter: int = 100, tol: float = 1e-6) -> Tuple[float, List[float], List[float]]:
        K = [Calculations.wilson_K(comp, T, P) for comp in components]

        for iteration in range(max_iter):
            V_low, V_high = 0.0, 1.0
            for _ in range(100):  # Bisection
                V = (V_low + V_high) / 2
                f_V = sum(z[i] * (K[i] - 1) / (1 + V * (K[i] - 1)) for i in range(len(z)))
                if abs(f_V) < tol:
                    break
                if f_V > 0:
                    V_low = V
                else:
                    V_high = V

            x = [z[i] / (1 + V * (K[i] - 1)) for i in range(len(z))]
            y = [K[i] * x[i] for i in range(len(z))]
            eos_vap = eos_class(components, T, P)
            eos_liq = eos_class(components, T, P)
            phi_v = eos_vap.fugacity_coeff(y, 'vapor')
            phi_l = eos_liq.fugacity_coeff(x, 'liquid')

            new_K = [phi_l[i] / phi_v[i] for i in range(len(z))]

            if all(abs((new_K[i] - K[i]) / K[i]) < tol for i in range(len(K))):
                break
            K = new_K

        return V, x, y