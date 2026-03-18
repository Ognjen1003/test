import src.EnumsClasses.MethodsAndTypes as MT
from src.Models.Component import Component
import math
from typing import List, Tuple
import numpy as np


class EOSBase:
    R = MT.CONSTANTS.R

    def __init__(self, components: List[Component], T: float, P: float, k_ij: np.ndarray | None = None):
        self.components = components
        self.T = T
        self.P = P

        self.nc = len(components)

        # cache pure-component arrays
        self.z_default = np.array([c.fraction for c in components], dtype=float)
        self.Tc = np.array([c.Tc for c in components], dtype=float)
        self.Pc = np.array([c.Pc for c in components], dtype=float)
        self.omega = np.array([c.omega for c in components], dtype=float)
        self.Mw = np.array([c.Mw for c in components], dtype=float)

        self.a_i = None
        self.b_i = None
        self.a0_i = None
        self.sqrt_ai = None

        self.Zv = None
        self.Zl = None

        self.k_ij = None
        self.one_minus_kij = None
        self.has_kij = False

        # pure component constants that do not depend on T
        self.a0_i = self._build_a0_array()
        self.b_i = self._build_b_array()

        if k_ij is not None:
            self.set_kij(k_ij)

        self.calc_parameters()

    def _build_a0_array(self) -> np.ndarray:
        return np.array([self.a_formula(comp) for comp in self.components], dtype=float)

    def _build_b_array(self) -> np.ndarray:
        return np.array([self.b_formula(comp) for comp in self.components], dtype=float)

    def calc_parameters(self):
        alpha_vals = self.alpha_array()
        self.a_i = self.a0_i * alpha_vals
        self.sqrt_ai = np.sqrt(self.a_i)

    def set_kij(self, k_ij: np.ndarray):
        k_ij = np.asarray(k_ij, dtype=float)
        if k_ij.shape != (self.nc, self.nc):
            raise ValueError(f"k_ij mora biti oblika {(self.nc, self.nc)}, dobiveno {k_ij.shape}")
        if not np.allclose(k_ij, k_ij.T, atol=1e-12):
            raise ValueError("k_ij mora biti simetrična matrica.")
        if not np.allclose(np.diag(k_ij), 0.0, atol=1e-12):
            raise ValueError("Dijagonala k_ij mora biti nula.")

        self.k_ij = k_ij
        self.one_minus_kij = 1.0 - k_ij
        self.has_kij = np.any(np.abs(k_ij) > 0.0)

    def update_state(self, T: float | None = None, P: float | None = None):
        recalc = False
        if T is not None and T != self.T:
            self.T = T
            recalc = True
        if P is not None:
            self.P = P
        if recalc:
            self.calc_parameters()

    def get_Z_factors(self, x: List[float], y: List[float]) -> Tuple[float, float]:
        x = np.asarray(x, dtype=float)
        _, _, A, B = self.get_mixture_parameters(x)
        Z_liq = self.calc_Z_factor(A, B, MT.Phase.LIQUID)

        y = np.asarray(y, dtype=float)
        _, _, A, B = self.get_mixture_parameters(y)
        Z_vap = self.calc_Z_factor(A, B, MT.Phase.VAPOR)

        self.Zl = Z_liq
        self.Zv = Z_vap
        return Z_liq, Z_vap

    def solve_cubic(self, coeffs: List[float]) -> List[float]:
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
            roots.append((S + T) - (b / (3 * a)))
        else:
            i = math.sqrt((g ** 2 / 4) - h)
            j = i ** (1 / 3)
            k = math.acos(-(g / (2 * i)))
            M = math.cos(k / 3)
            N = math.sqrt(3) * math.sin(k / 3)
            roots.extend([
                2 * j * M - (b / (3 * a)),
                -j * (M + N) - (b / (3 * a)),
                -j * (M - N) - (b / (3 * a)),
            ])

        return [r for r in roots if isinstance(r, float) and not math.isnan(r)]

    def calc_Z_factor(self, A: float, B: float, phase: MT.Phase) -> float:
        coeffs = self.get_coefficients(A, B)
        roots = self.solve_cubic(coeffs)

        if not roots:
            raise ValueError(
                f"No real roots found in calc_Z_factor. "
                f"A={A:.5e}, B={B:.5e}, coeffs={coeffs}, phase='{phase}', "
                f"T={self.T}, P={self.P}."
            )

        if phase == MT.Phase.VAPOR:
            self.Zv = max(roots)
            return self.Zv
        else:
            self.Zl = min(roots)
            return self.Zl

    def ensure_Z_for_phase(self, composition, phase: MT.Phase) -> float:
        x = np.asarray(composition, dtype=float)
        _, _, A, B = self.get_mixture_parameters(x)
        return self.calc_Z_factor(A, B, phase)

    def get_mixture_parameters(self, composition: np.ndarray = None) -> tuple[float, float, float, float]:
        x = self.z_default if composition is None else np.asarray(composition, dtype=float)

        s = x.sum()
        if s <= 0.0:
            raise ValueError("Suma sastava mora biti > 0.")
        if abs(s - 1.0) > 1e-14:
            x = x / s

        u_vec = x * self.sqrt_ai

        if not self.has_kij:
            S = float(u_vec.sum())
            a_mix = S * S
        else:
            a_mix = float(u_vec @ (self.one_minus_kij @ u_vec))

        b_mix = float(x @ self.b_i)

        R, T, P = self.R, self.T, self.P
        A = a_mix * P / (R * R * T * T)
        B = b_mix * P / (R * T)

        return a_mix, b_mix, A, B

    def fugacity_coeff(self, x: List[float], phase: MT.Phase) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        a_mix, b_mix, A, B = self.get_mixture_parameters(x)
        Z = self.calc_Z_factor(A, B, phase)

        eps = 1e-14
        bi = self.b_i

        u = self.eos_u
        w = self.eos_w
        delta = math.sqrt(u * u - 4.0 * w)

        term_log = np.log(np.clip(Z - B, eps, None))

        vec = x * self.sqrt_ai
        if not self.has_kij:
            S = float(vec.sum())
            sum_a_vec = self.sqrt_ai * S
        else:
            sum_a_vec = self.sqrt_ai * (self.one_minus_kij @ vec)

        term1 = (bi / b_mix) * (Z - 1.0) - term_log
        term3 = 2.0 * (sum_a_vec / a_mix) - (bi / b_mix)

        log_arg_num = np.clip(2.0 * Z + B * (u + delta), eps, None)
        log_arg_den = np.clip(2.0 * Z + B * (u - delta), eps, None)
        L = np.log(log_arg_num / log_arg_den)

        coeff = A / (B * delta)
        ln_phi = term1 - coeff * term3 * L
        return np.exp(ln_phi)

    def get_pressure(self, v_molar: float, composition: np.ndarray = None) -> float:
        a_mix, b_mix, _, _ = self.get_mixture_parameters(composition)

        u = self.eos_u
        w = self.eos_w

        return (
            self.R * self.T / (v_molar - b_mix)
            - a_mix / (v_molar ** 2 + u * v_molar * b_mix + w * b_mix ** 2)
        )

    def da_mix_dT(self, composition: np.ndarray = None) -> float:
        x = self.z_default if composition is None else np.asarray(composition, dtype=float)

        s = x.sum()
        if s <= 0.0:
            raise ValueError("Suma sastava mora biti > 0.")
        if abs(s - 1.0) > 1e-14:
            x = x / s

        da_i_dT = self.a0_i * self.dalpha_dT_array()

        u_vec = x * self.sqrt_ai
        du_dT = x * 0.5 * da_i_dT / self.sqrt_ai

        if not self.has_kij:
            S = float(u_vec.sum())
            return 2.0 * S * float(du_dT.sum())
        else:
            return 2.0 * float(du_dT @ (self.one_minus_kij @ u_vec))

    def _departure_log_term(self, Z: float, B: float) -> float:
        eps = 1e-14
        u = self.eos_u
        w = self.eos_w
        delta = math.sqrt(u * u - 4.0 * w)

        num = max(2.0 * Z + B * (u + delta), eps)
        den = max(2.0 * Z + B * (u - delta), eps)
        return math.log(num / den)

    def h_residual(self, p: float, T: float) -> float:
        T_old, P_old = self.T, self.P
        self.update_state(T=T, P=p)

        try:
            x = self.z_default
            a_mix, b_mix, A, B = self.get_mixture_parameters(x)
            Z = self.calc_Z_factor(A, B, MT.Phase.VAPOR)

            delta = math.sqrt(self.eos_u * self.eos_u - 4.0 * self.eos_w)
            L = self._departure_log_term(Z, B)
            da_mix_dT = self.da_mix_dT(x)

            hR_molar = self.R * T * (Z - 1.0) + ((T * da_mix_dT - a_mix) / (b_mix * delta)) * L
            return hR_molar / self.mixture_molar_mass()

        finally:
            self.update_state(T=T_old, P=P_old)

    def s_residual(self, p: float, T: float) -> float:
        T_old, P_old = self.T, self.P
        self.update_state(T=T, P=p)

        try:
            x = self.z_default
            _, b_mix, A, B = self.get_mixture_parameters(x)
            Z = self.calc_Z_factor(A, B, MT.Phase.VAPOR)

            delta = math.sqrt(self.eos_u * self.eos_u - 4.0 * self.eos_w)
            L = self._departure_log_term(Z, B)
            da_mix_dT = self.da_mix_dT(x)

            eps = 1e-14
            sR_molar = self.R * math.log(max(Z - B, eps)) + (da_mix_dT / (b_mix * delta)) * L
            return sR_molar / self.mixture_molar_mass()

        finally:
            self.update_state(T=T_old, P=P_old)

    def cp_ideal(self, T: float) -> float:
        cp_molar = 0.0
        for comp in self.components:
            cp_molar += comp.fraction * comp.cp_ideal_molar(T)
        return cp_molar / self.mixture_molar_mass()

    def h_ideal(self, T: float) -> float:
        h_molar = 0.0
        for comp in self.components:
            h_molar += comp.fraction * comp.h_ideal_molar(T)
        return h_molar / self.mixture_molar_mass()

    def s_ideal(self, T: float, p: float) -> float:
        s_molar = 0.0
        for comp in self.components:
            s_molar += comp.fraction * comp.s_ideal_molar(T, p)
        return s_molar / self.mixture_molar_mass()

    def h(self, p: float, T: float) -> float:
        return self.h_ideal(T) + self.h_residual(p, T)

    def s(self, p: float, T: float) -> float:
        return self.s_ideal(T, p) + self.s_residual(p, T)

    def mixture_molar_mass(self) -> float:
        return float(self.z_default @ self.Mw)

    @property
    def eos_u(self) -> float:
        raise NotImplementedError

    @property
    def eos_w(self) -> float:
        raise NotImplementedError

    def alpha(self, comp: Component) -> float:
        raise NotImplementedError

    def dalpha_dT(self, comp: Component) -> float:
        raise NotImplementedError

    def alpha_array(self) -> np.ndarray:
        return np.array([self.alpha(comp) for comp in self.components], dtype=float)

    def dalpha_dT_array(self) -> np.ndarray:
        return np.array([self.dalpha_dT(comp) for comp in self.components], dtype=float)

    def a_formula(self, comp: Component) -> float:
        raise NotImplementedError

    def b_formula(self, comp: Component) -> float:
        raise NotImplementedError

    def get_coefficients(self, A: float, B: float) -> List[float]:
        raise NotImplementedError