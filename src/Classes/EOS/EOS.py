import math
from typing import List, Tuple
from Models.Component import Component
import numpy as np

class EOSBase:
    """
    Base class
    """

    R = 8.314
    
    def __init__(self, components: List[Component], T: float, P: float):
        self.components = components  
        self.T = T                    
        self.P = P                   
        self.a_i = []                 
        self.b_i = []
        self.Zv = None
        self.Zl = None                 
        self.calc_parameters()   


    def calc_parameters(self):
        alpha_vals = np.array([self.alpha(comp) for comp in self.components])
        self.a_i = np.array([self.a_formula(comp) for comp in self.components]) * alpha_vals
        self.b_i = np.array([self.b_formula(comp) for comp in self.components])


    def get_Z_factors(self, x: List[float], y: List[float]) -> Tuple[float, float]:
        x = np.array(x)
        y = np.array(y)

        def compute_A_B(comp_frac: np.ndarray) -> Tuple[float, float]:
            b_mix = np.dot(comp_frac, self.b_i)
            a_mix = np.sum(np.outer(comp_frac, comp_frac) * np.sqrt(np.outer(self.a_i, self.a_i)))
            A = a_mix * self.P / (self.R ** 2 * self.T ** 2)
            B = b_mix * self.P / (self.R * self.T)
            return A, B

        A_liq, B_liq = compute_A_B(x)
        A_vap, B_vap = compute_A_B(y)

        Z_liq = self.calc_Z_factor(A_liq, B_liq, phase='liquid')
        Z_vap = self.calc_Z_factor(A_vap, B_vap, phase='vapor')

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

    def calc_Z_factor(self, A: float, B: float, phase: str = 'vapor') -> float:
        """
        Z factor
        """
        coeffs = self.get_coefficients(A, B)
        roots = self.solve_cubic(coeffs)

        if not roots:
            raise ValueError(
                f"No real roots found in calc_Z_factor. "
                f"A={A:.5e}, B={B:.5e}, coeffs={coeffs}, phase='{phase}', "
                f"T={self.T}, P={self.P}. Check EOS parameters or system conditions."
            )
        
        if phase == 'vapor':
            self.Zv = max(roots)
        else:
            self.Zl = min(roots)

        return self.Zv if phase == 'vapor' else self.Zl


    def fugacity_coeff(self, x: List[float], phase: str = 'vapor') -> list:
        
        x = np.array(x)
        a_mix, b_mix, A, B = self.get_mixture_parameters(x)
        Z = self.calc_Z_factor(A, B, phase)

        #sqrt_ai = np.sqrt(self.a_i)
        phi = []

        for i in range(len(x)):
            sum_a = np.sum(x * np.sqrt(self.a_i[i] * self.a_i))
            bi = self.b_i[i]
            term1 = bi / b_mix * (Z - 1) - np.log(Z - B)
            term2 = A / (2 * np.sqrt(2) * B)
            term3 = 2 * sum_a / a_mix - bi / b_mix
            term4 = np.log((Z + (1 + np.sqrt(2)) * B) / (Z + (1 - np.sqrt(2)) * B))
            ln_phi = term1 - term2 * term3 * term4
            phi.append(np.exp(ln_phi))

        return np.array(phi)

    
    def get_pressure(self, v_molar: float, phase: str = 'vapor') -> float: 
        a_mix, b_mix, _, _ = self.get_mixture_parameters()

        T = self.T
        R = self.R

        num1 = R * T
        denom1 = v_molar - b_mix

        num2 = a_mix
        denom2 = v_molar ** 2 + 2 * v_molar * b_mix - b_mix ** 2

        P = num1 / denom1 - num2 / denom2
        return P

    
    def get_mixture_parameters(self, composition: np.ndarray = None) -> Tuple[float, float, float, float]:
        if composition is None:
            composition = np.array([comp.fraction for comp in self.components])
        
        a_mix = np.sum(np.outer(composition, composition) * np.sqrt(np.outer(self.a_i, self.a_i)))
        b_mix = np.dot(composition, self.b_i)

        A = a_mix * self.P / (self.R ** 2 * self.T ** 2)
        B = b_mix * self.P / (self.R * self.T)

        return a_mix, b_mix, A, B