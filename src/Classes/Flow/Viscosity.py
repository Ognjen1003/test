from __future__ import annotations
import math
import numpy as np
from typing import List, Literal, Optional, Dict
import src.EnumsClasses.MethodsAndTypes as MT
from math import sqrt, exp
from typing import List


R = 8.314462618  # J/(mol·K)
BAR_TO_PA = 1e5
CM3_PER_M3 = 1e6  # 1 m^3 = 1e6 cm^3


class ViscosityClass:

    @staticmethod
    def omega_mu_neufeld(T_star: float) -> float:
        """Collision integral for viscosity, Neufeld et al. (dimensionless)."""
        # Valid for roughly 0.3 <= T* <= 100
        return (1.16145 * T_star**-0.14874
                + 0.52487 * exp(-0.77320 * T_star)
                + 2.16178 * exp(-2.43787 * T_star))

    @staticmethod
    def estimate_sigma_epsilon(component) -> tuple[float, float]:
        """
        Procijeni (sigma [Å], epsilon_over_k [K]) iz Tc, Pc, omega.
        - Zc ~ 0.291 - 0.08*omega (Edmister)
        - Vc = Zc*R*Tc/Pc (m^3/mol)  -> u cm^3/mol
        - sigma = 0.809 * Vc^(1/3)  (Å), kad je Vc u cm^3/mol
        - epsilon/k ~ 1.2593 * Tc / Zc^(2/3)  (K)
        """
        Tc = component.Tc          # K
        Pc_bar = component.Pc      # bar 
        Pc = Pc_bar * BAR_TO_PA    # Pa
        
        omega = component.omega
        # ograniči Zc na razuman raspon
        Zc = max(0.2, min(0.35, 0.291 - 0.08 * omega))
        Vc_m3_per_mol = Zc * R * Tc / Pc
        Vc_cm3_per_mol = Vc_m3_per_mol * CM3_PER_M3

        sigma = 0.809 * (Vc_cm3_per_mol ** (1.0 / 3.0))  # Å
        epsilon_over_k = 1.2593 * Tc / (Zc ** (2.0 / 3.0))  # K
        return sigma, epsilon_over_k


    @staticmethod
    def mu_component_dilute(component, T: float) -> float:
        """
        Dilute-gas viscosity of a pure component at temperature T [K], in Pa·s.
        Uses Chapman–Enskog with Neufeld collision integral and
        sigma/epsilon derived from critical data.
        """
        M_g_per_mol = component.Mw  # g/mol (tako je u tvojim podacima)
        sigma, eps_over_k = ViscosityClass.estimate_sigma_epsilon(component)
        T_star = T / eps_over_k
        omega_mu = ViscosityClass.omega_mu_neufeld(T_star)

        # Chapman–Enskog: mu(μP) = 26.693 * sqrt(M*T) / (sigma^2 * omega_mu)
        # gdje je M u g/mol, T u K, sigma u Å. Pretvori u Pa·s: 1 μP = 1e-7 Pa·s.
        mu_micropoise = 26.693 * sqrt(M_g_per_mol * T) / (sigma * sigma * omega_mu)
        mu_Pa_s = mu_micropoise * 1e-7
        return mu_Pa_s


    @staticmethod
    def normalize_mole_fractions(ys: List[float]) -> List[float]:
        s = sum(ys)
        if s <= 0.0:
            raise ValueError("Zbroj molnih udjela je nula.")
        return [y / s for y in ys]

    @staticmethod
    def gas_viscosity_wilke(components: List, T: float) -> float:
        """
        Viskoznost plinovite smjese (Pa·s) pri T [K] u dilute-gas limitu.
        - components: lista Component objekata (ima .Mw, .omega, .Tc, .Pc, .fraction)
        - T: temperatura u K
        Koraci:
        1) μ_i(T) za svaku komponentu (Chapman–Enskog + Neufeld; σ, ε/k iz Tc,Pc,ω).
        2) Wilke miješanje s molnim udjelima (fraction).
        """
        if not components:
            raise ValueError("Prazna lista komponenti.")

        y =  ViscosityClass.normalize_mole_fractions([c.fraction for c in components])
        M = [c.Mw for c in components]
        mu_i = [ViscosityClass.mu_component_dilute(c, T) for c in components]

        # Wilke Φ_ij
        # Φ_ij = [1 + (μ_i/μ_j)^{1/2} (M_j/M_i)^{1/4}]^2 / sqrt(8(1 + M_i/M_j))
        n = len(components)
        phi = [[0.0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i == j:
                    phi[i][j] = 1.0
                else:
                    num = (1.0 + (mu_i[i]/mu_i[j])**0.5 * (M[j]/M[i])**0.25)**2
                    den = (8.0 * (1.0 + M[i]/M[j]))**0.5
                    phi[i][j] = num / den

        # μ_mix = sum_i [ y_i μ_i / sum_j ( y_j Φ_ij ) ]
        mu_mix = 0.0
        for i in range(n):
            denom = 0.0
            for j in range(n):
                denom += y[j] * phi[i][j]
            mu_mix += y[i] * mu_i[i] / denom

        return mu_mix
    
    @staticmethod
    def viscosity_lbc(components: List, T: float, rho: float) -> float:
        """
        LBC/JST viskoznost smjese (Pa·s) pri T [K] i rho [kg/m3].
        - μ = μ0(T) + Δμ_LBC(rho_r), gdje je μ0 tvoje Chapman-Enskog+Wilke razrjeđeno-PLIN.
        - rho_r = rho / rho_c, a rho_c iz Kay pseudo-kritičnih (Tc_mix, Pc_mix, Vc_mix).
        - Δμ_LBC (u cP) = ξ * poly(rho_r), pa se zbroji s μ0 u cP.
        
        Napomene:
        * μ0 je iz gas_viscosity_wilke (Pa·s).
        * ξ = Tc_mix^(1/6) * MW_mix^(-1/2) * (Pc_mix[atm])^(-2/3)
        * poly(rho_r) = 0.1023 + 0.023364 rho_r + 0.058533 rho_r^2 - 0.040758 rho_r^3 + 0.0093724 rho_r^4
          (standardna JST/LBC varijanta u praksi)
        """

        y = [c.fraction for c in components]

        # 1) μ0(T) iz tvoje implementacije (Pa·s) -> u cP
        mu0_Pa_s = ViscosityClass.gas_viscosity_wilke(components, T)
        mu0_cP = mu0_Pa_s * 1000.0  # 1 cP = 1e-3 Pa·s

        # 2) Kay pseudo-kritične: Tc_mix, Pc_mix (Pa), Vc_mix (m3/mol)
        Tc_mix = sum(yi * c.Tc for yi, c in zip(y, components))
        # Pc_mix: 1/Pc_mix = sum(y_i / Pc_i)
        Pc_inv_bar = sum(yi / c.Pc for yi, c in zip(y, components))  # 1/bar
        if Pc_inv_bar <= 0.0:
            raise ValueError("Neispravne kritične tlakove.")
        Pc_mix_Pa = (1.0 / Pc_inv_bar) * BAR_TO_PA

        # Vc_i iz Tc,Pc,omega: Vc = Zc R Tc / Pc, Zc ~ 0.291 - 0.08*omega (ograničeno 0.2–0.35)
        Vc_i = []
        for yi, c in zip(y, components):
            Zc_i = max(0.2, min(0.35, 0.291 - 0.08 * c.omega))
            Pc_i_Pa = c.Pc * BAR_TO_PA
            Vc_i.append(Zc_i * R * c.Tc / Pc_i_Pa)  # m3/mol
        Vc_mix = sum(yi * v for yi, v in zip(y, Vc_i))  # m3/mol

        # 3) ρ_c,mix = M_mix / Vc_mix (kg/m3)
        Mw_mix_g_per_mol = sum(yi * c.Mw for yi, c in zip(y, components))  # g/mol
        Mw_mix_kg_per_mol = Mw_mix_g_per_mol / 1000.0
        if Vc_mix <= 0.0:
            raise ValueError("Pseudo-kritični Vc_mix <= 0.")
        rho_c_mix = Mw_mix_kg_per_mol / Vc_mix  # kg/m3

        # 4) Reducirana gustoća
        rho_r = rho / rho_c_mix

        # 5) LBC/JST korekcija: Δμ (cP) = ξ * poly(ρ_r)
        Pc_mix_atm = Pc_mix_Pa / 101325.0
        if Pc_mix_atm <= 0.0:
            raise ValueError("Pc_mix_atm <= 0.")
        xi = (Tc_mix ** (1.0/6.0)) * (Mw_mix_g_per_mol ** -0.5) * (Pc_mix_atm ** (-2.0/3.0))

        # polinom (JST/LBC)
        poly = (0.1023
                + 0.023364 * rho_r
                + 0.058533 * rho_r**2
                - 0.040758 * rho_r**3
                + 0.0093724 * rho_r**4)

        delta_mu_cP = xi * poly

        # 6) μ = μ0 + Δμ  (u cP), pa natrag u Pa·s
        mu_cP = mu0_cP + delta_mu_cP
        mu_Pa_s = max(mu_cP, 0.0) / 1000.0
        return mu_Pa_s

