import math
from typing import List, Tuple
from src.Classes.Component import Component

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


    def get_Z_factors(self, x: List[float], y: List[float]) -> Tuple[float, float]: 

        def compute_A_B(comp_frac: List[float]) -> Tuple[float, float]:
            a_mix = 0.0
            b_mix = 0.0
            N = len(comp_frac)

            for i in range(N):
                b_mix += comp_frac[i] * self.b_i[i]
                for j in range(N):
                    a_mix += comp_frac[i] * comp_frac[j] * math.sqrt(self.a_i[i] * self.a_i[j])

            A = a_mix * self.P / (self.R**2 * self.T**2)
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

    def calc_parameters(self):
        """
        ai, bi for all components
        """
        for comp in self.components:
            alpha_i = self.alpha(comp)
            a = self.a_formula(comp) * alpha_i
            b = self.b_formula(comp)
            self.a_i.append(a)
            self.b_i.append(b)

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
        
        if phase == 'vapor':
            self.Zv = max(roots)
        else:
            self.Zl = min(roots)

        return max(roots) if phase == 'vapor' else min(roots)

    def fugacity_coeff(self, x: List[float], phase: str = 'vapor') -> List[float]:
        """
        fugacity
        """
        a_mix = 0.0
        b_mix = 0.0
        N = len(x)

        for i in range(N):
            b_mix += x[i] * self.b_i[i]
            for j in range(N):
                a_mix += x[i] * x[j] * math.sqrt(self.a_i[i] * self.a_i[j])

        A = a_mix * self.P / (self.R**2 * self.T**2)
        B = b_mix * self.P / (self.R * self.T)
        Z = self.calc_Z_factor(A, B, phase)

        phi = []
        for i in range(N):
            bi = self.b_i[i]
            ai = self.a_i[i]
            sum_a = sum(x[j] * math.sqrt(ai * self.a_i[j]) for j in range(N))
            term1 = bi / b_mix * (Z - 1) - math.log(Z - B)
            term2 = A / (2 * math.sqrt(2) * B)
            term3 = 2 * sum_a / a_mix - bi / b_mix
            term4 = math.log((Z + (1 + math.sqrt(2)) * B) / (Z + (1 - math.sqrt(2)) * B))
            ln_phi = term1 - term2 * term3 * term4
            phi.append(math.exp(ln_phi))

        return phi
    
    def get_pressure_with_Z(self, v_molar: float, phase: str = 'vapor') -> float:
        z = [comp.fraction for comp in self.components]

        a_mix = 0.0
        b_mix = 0.0
        for i in range(len(z)):
            b_mix += z[i] * self.b_i[i]
            for j in range(len(z)):
                a_mix += z[i] * z[j] * math.sqrt(self.a_i[i] * self.a_i[j])

        A = a_mix * self.P / (self.R ** 2 * self.T ** 2)
        B = b_mix * self.P / (self.R * self.T)

        Z = self.calc_Z_factor(A, B, phase)
        return Z * self.R * self.T / v_molar
