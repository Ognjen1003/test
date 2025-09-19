import math
from typing import List, Tuple
import EnumsClasses.MethodsAndTypes
from Models.Component import Component
import numpy as np
import EnumsClasses.MethodsAndTypes as MT


class EOSBase:

    R = MT.CONSTANTS.R
    
    def __init__(self, components: List[Component], T: float, P: float):
        self.components = components  
        self.T = T                    
        self.P = P                   
        self.a_i = []                 
        self.b_i = []
        self.Zv = None
        self.Zl = None
        self.sqrt_ai = None
        self.calc_parameters()
        self.k_ij = None

    def calc_parameters(self):
        alpha_vals = np.array([self.alpha(comp) for comp in self.components])
        self.a_i = np.array([self.a_formula(comp) for comp in self.components]) * alpha_vals
        self.b_i = np.array([self.b_formula(comp) for comp in self.components])
        self.sqrt_ai = np.sqrt(self.a_i)

    def get_Z_factors(self, x: List[float], y: List[float]) -> Tuple[float, float]:
        
        x = np.array(x)
        a_mix, b_mix, A, B = self.get_mixture_parameters(x)
        Z_liq= self.calc_Z_factor(A, B, 'liquid')

        y = np.array(y)
        a_mix, b_mix, A, B = self.get_mixture_parameters(y)
        Z_vap = self.calc_Z_factor(A, B, 'vapor')

        self.Zl = Z_liq
        self.Zv = Z_vap
        return Z_liq, Z_vap
    
    def alpha(self, comp: Component) -> float:
        """ 
        Placeholder method alpha correction
        """
        raise NotImplementedError

    def a_formula(self, comp: Component) -> float:
        raise NotImplementedError

    def b_formula(self, comp: Component) -> float:
        raise NotImplementedError

    def get_coefficients(self, A: float, B: float) -> List[float]:
        raise NotImplementedError

    def solve_cubic(self, coeffs: List[float]) -> List[float]:
        """
        Cubic eq. - Cardano.
        """

        a, b, c, d = coeffs
        f = ((3 * c / a) - (b ** 2 / a ** 2)) / 3
        g = ((2 * b ** 3 / a ** 3) - (9 * b * c / a ** 2) + (27 * d / a)) / 27
        h = g ** 2 / 4 + f ** 3 / 27

        roots = []
        if h > 0:
            R_ = -(g / 2) + math.sqrt(h)
            S = math.copysign(abs(R_) ** (1 / 3), R_)
            T_ = -(g / 2) - math.sqrt(h)
            T = math.copysign(abs(T_) ** (1 / 3), T_)
            root1 = (S + T) - (b / (3 * a))
            roots.append(root1)
        else:
            i = math.sqrt((g ** 2 / 4) - h)
            j = i ** (1 / 3)
            k = math.acos(-(g / (2 * i)))
            M = math.cos(k / 3)
            N = math.sqrt(3) * math.sin(k / 3)
            root1 = 2 * j * M - (b / (3 * a))
            root2 = -j * (M + N) - (b / (3 * a))
            root3 = -j * (M - N) - (b / (3 * a))
            roots.extend([root1, root2, root3])

        return [r for r in roots if isinstance(r, float) and not math.isnan(r)]

    def calc_Z_factor(self, A: float, B: float, phase: MT.Phase) -> float:

        coeffs = self.get_coefficients(A, B)
        roots = self.solve_cubic(coeffs)

        if not roots:
            raise ValueError(
                f"No real roots found in calc_Z_factor. "
                f"A={A:.5e}, B={B:.5e}, coeffs={coeffs}, phase='{phase}', "
                f"T={self.T}, P={self.P}. Check EOS parameters or system conditions."
            )
        
        if phase == MT.Phase.VAPOR:
            self.Zv = max(roots)
        else:
            self.Zl = min(roots)

        return self.Zv if phase == MT.Phase.VAPOR else self.Zl

    def ensure_Z_for_phase(self, composition, phase: MT.Phase) -> float:
        x = np.asarray(composition, dtype=float)
        _, _, A, B = self.get_mixture_parameters(x)
        return self.calc_Z_factor(A, B, phase)

    def fugacity_coeff(self, x: List[float], phase: MT.Phase) -> np.ndarray:
        x = np.array(x, dtype=float)
        a_mix, b_mix, A, B = self.get_mixture_parameters(x)
        Z = self.calc_Z_factor(A, B, phase)

        # --- preduvjeti i konstante ---
        sqrt2 = np.sqrt(2.0)
        eps = 1e-14

        # vektorizacija: sum_a_i = sqrt(a_i) * sum_j x_j sqrt(a_j)
        sqrt_ai = self.sqrt_ai                    # (nc,)
        S = float(x @ sqrt_ai)                        # skalar
        sum_a_vec = sqrt_ai * S                       # (nc,)

        bi = self.b_i                                 # (nc,)

        # skalari
        term_log = np.log(np.clip(Z - B, eps, None))
        C = A / (2.0 * sqrt2 * B)
        L = np.log(
            np.clip(Z + (1.0 + sqrt2) * B, eps, None) /
            np.clip(Z + (1.0 - sqrt2) * B, eps, None)
        )

        # vektori
        term1 = (bi / b_mix) * (Z - 1.0) - term_log
        term3 = 2.0 * (sum_a_vec / a_mix) - (bi / b_mix)

        ln_phi = term1 - C * term3 * L                # (nc,)
        return np.exp(ln_phi)

    def get_pressure(self, v_molar: float) -> float: 
        a_mix, b_mix, _, _ = self.get_mixture_parameters()

        T = self.T
        R = self.R

        num1 = R * T
        denom1 = v_molar - b_mix

        num2 = a_mix
        denom2 = v_molar ** 2 + 2 * v_molar * b_mix - b_mix ** 2

        P = num1 / denom1 - num2 / denom2

        return P
  
    def get_mixture_parameters(self, composition: np.ndarray = None) -> tuple[float, float, float, float]:

        if composition is None:
            x = np.array([c.fraction for c in self.components], dtype=float)
        else:
            x = np.asarray(composition, dtype=float)
        s = x.sum()
        if s <= 0.0:
            raise ValueError("Suma sastava mora biti > 0.")
        x = x / s

        # pomocni vektor u_i = x_i * sqrt(a_i)
        u = x * self.sqrt_ai

        if self.k_ij is None or not np.any(self.k_ij):
            S = float(u.sum())
            a_mix = S * S
        else:
            K = np.asarray(self.k_ij, dtype=float)
            a_mix = float(u @ ((1.0 - K) @ u))

        # 3) b_mix
        b_mix = float(x @ self.b_i)

        R, T, P = self.R, self.T, self.P
        A = a_mix * P / (R * R * T * T)
        B = b_mix * P / (R * T)

        return a_mix, b_mix, A, B
