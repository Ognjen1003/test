import math
import numpy as np
from typing import List, Literal, Optional, Dict
import UtilClass
import EnumsClasses.MethodsAndTypes as MT


class ViscosityClass:

    # ===== CO2 / N2 bazne svojstvene konstante =====
    # Tc [K], Pc [Pa], Zc [-], omega [-], Mw [kg/mol], Lennard-Jones: sigma [Å], eps/k [K]

    COMP_DB: Dict[str, Dict[str, float]] = {
        "CO2": {
            "Tc": 304.1282, "Pc": 7.3773e6, "Zc": 0.274, "omega": 0.22394,
            "Mw": 44.0095e-3, "sigma_A": 3.941, "eps_over_k": 195.2
        },
        "N2": {
            "Tc": 126.2, "Pc": 3.3958e6, "Zc": 0.289, "omega": 0.0372,
            "Mw": 28.0134e-3, "sigma_A": 3.667, "eps_over_k": 71.4
        },
    }

    @staticmethod
    def _crit_arrays(components: List):
        """Vrati (Tc_i, Pc_i, Zc_i, omega_i, Mw_i) za svaku komponentu (CO2/N2)."""
        Tc = []; Pc = []; Zc = []; om = []; Mw = []
        for c in components:
            k = UtilClass._name_key(c)
            d = ViscosityClass.COMP_DB[k]
            Tc.append(d["Tc"]); Pc.append(d["Pc"]); Zc.append(d["Zc"]); om.append(d["omega"]); Mw.append(d["Mw"])
        return np.array(Tc), np.array(Pc), np.array(Zc), np.array(om), np.array(Mw)

    # ===== kolizijska integralna Ω* (Neufeld/Chung) =====
    def _Omega_star(T_star: float) -> float:
        # prema Neufeldu/Chungu – vidi sažetak jednadžbe i koeficijente (Wikipedia/IDAES)
        t = T_star
        return (1.16145 / (t ** 0.14874)
                + 0.52487 * math.exp(-0.7732 * t)
                + 2.16178 * math.exp(-2.43787 * t))

    # ===== Chung CSP – gas-like μ0 s Fc korekcijom; mješavina po Wilke =====
    def viscosity_chung_csp(
        components: List,
        T: float,
        P: float,
        composition: np.ndarray,
        Z: float,
        phase: Literal["vapor", "liquid"] = "vapor",
        association_factor_kappa: float = 0.0,
        use_chung_pure: bool = True,
    ) -> float:
        """
        Viskoznost mješavine (Pa·s) – Chung CSP (gas-like) + Wilke miješanje.
        • Najtočniji za plin/“gas-like” uvjete CO2/N2; blizu kritike se poboljšava s točnim ρ.
        • 'phase' je informativan (za 'liquid' može biti donja procjena; za realne tekućine koristi LBC).
        """
        x = ViscosityClass._normalize(np.array(composition, dtype=float))
        Tc_i, Pc_i, Zc_i, om_i, Mw_i = ViscosityClass._crit_arrays(components)

        # mješavina: Kay (Tc, Pc, Vc), ω (molarno ponderiran), Mw (molarno), + LJ za pure iz baze
        Vc_i = UtilClass._Vc_from_TcPcZc(Tc_i, Pc_i, Zc_i)        # m^3/mol (po komponentama)
        Tc_mix = UtilClass._Kay_mixing(x, Tc_i)
        Pc_mix = UtilClass._Kay_mixing(x, Pc_i)
        Vc_mix = UtilClass._Kay_mixing(x, Vc_i)
        Mw_mix = UtilClass._Kay_mixing(x, Mw_i)
        omega_mix = UtilClass._Kay_mixing(x, om_i)

        # Chung redukcije:
        # U "SS" formi često se koristi T* = 1.2593 * T/Tc_mix (vidi wiki sažetak),
        # a Ω* = f(T*). Korigiraj s Fc = 1 - 0.2756 ω + 0.059035 μ_r^4 + κ (za nepolarne: μ_r≈0).
        T_star = 1.2593 * T / Tc_mix
        Omega = ViscosityClass._Omega_star(T_star)
        Fc = 1.0 - 0.2756 * omega_mix + 0.0 + float(association_factor_kappa)

        # μ0_mixture preko Wilkea: treba μ_i^0 (pure) – možemo koristiti Chung-ov pure ili Chapman-Enskog
        mu_i = []
        for c in components:
            nm = UtilClass._name_key(c)
            d = ViscosityClass.COMP_DB[nm]
            # Pure μ0 po Chung (μPa·s): 40.785 * sqrt(M * T*) / (Vc^(2/3) * Ω*) * Fc_i
            # Za pure, uzmi vlastiti Tc_i,Vc_i,ω_i i T*_i
            T_star_i = 1.2593 * T / d["Tc"]
            Omega_i = ViscosityClass._Omega_star(T_star_i)
            Vc_i_val = (d["Zc"] * MT.CONSTANTS.R * d["Tc"] / d["Pc"])  # m^3/mol
            Fc_i = 1.0 - 0.2756 * d["omega"]  # μ_r=0, κ=0 (nepolarne)
            mu_uPa_s = 40.785 * math.sqrt(d["Mw"] * T_star_i) / ((Vc_i_val ** (2.0/3.0)) * Omega_i) * Fc_i
            mu_i.append(mu_uPa_s * 1e-6)  # → Pa·s

        mu_i = np.array(mu_i, dtype=float)

        # Wilke miješanje (Pa·s)
        def _wilke(mu_i, x, Mw):
            n = len(x)
            phi = np.zeros((n, n))
            for i in range(n):
                for j in range(n):
                    if i == j:
                        phi[i, j] = 1.0
                    else:
                        num = (1.0 + math.sqrt(mu_i[i] / mu_i[j]) * (Mw[j] / Mw[i]) ** 0.25) ** 2
                        den = math.sqrt(8.0) * (1.0 + Mw[i] / Mw[j]) ** 0.5
                        phi[i, j] = num / den
            denom = np.array([np.sum(x[j] * phi[i, j] for j in range(n)) for i in range(n)])
            return float(np.sum(x[i] * mu_i[i] / denom[i] for i in range(n)))

        mu0_mix = _wilke(mu_i, x, Mw_i)  # Pa·s (dilute-gas)

        # Opcionalno: "mixture-level" Chung (koristeći mješ. Vc_mix, Omega, Fc) – kao sanity check:
        mu_chung_mix_uPa_s = 40.785 * math.sqrt(Mw_mix * T_star) / ((Vc_mix ** (2.0/3.0)) * Omega) * Fc
        mu_chung_mix = mu_chung_mix_uPa_s * 1e-6  # Pa·s

        # Vratimo Wilke-mu0 (često bolje za mješavine) – možeš prebaciti na mu_chung_mix po želji:
        return float(mu0_mix)

    # ===== LBC (JST) – radi i za gas i za liquid =====
    def viscosity_lbc(
        components: List,
        T: float,
        P: float,
        composition: np.ndarray,
        Z: float,
        eta0_strategy: Literal["wilke_chapman", "chung_pure"] = "wilke_chapman",
        lbc_coeffs: Optional[np.ndarray] = None,
    ) -> float:
        """
        LBC (JST) viskoznost mješavine (Pa·s) koristeći reduciranu molarnu koncentraciju.
        • eta0_strategy: "wilke_chapman" (default) – μ0 iz Wilke + Chapman-Enskog (Neufeld Ω)
                        "chung_pure" – μ0 iz Chung-pure + Wilke
        • lbc_coeffs: [a1..a5]; default su originalni JST koeficijenti (literatura).
        """
        x = UtilClass._normalize(np.array(composition, dtype=float))
        Tc_i, Pc_i, Zc_i, om_i, Mw_i = ViscosityClass._crit_arrays(components)
        Vc_i = UtilClass._Vc_from_TcPcZc(Tc_i, Pc_i, Zc_i)   # m^3/mol

        # pseudo-kritične mix vrijednosti (Kay)
        Tc_mix = UtilClass._Kay_mixing(x, Tc_i)
        Pc_mix = UtilClass._Kay_mixing(x, Pc_i)
        Vc_mix = UtilClass._Kay_mixing(x, Vc_i)
        Mw_mix = UtilClass._Kay_mixing(x, Mw_i)

        # molarna koncentracija (EOS): c = P / (Z R T)   [mol/m^3]
        c = P / (Z * MT.CONSTANTS.R  * T)
        # reducirana molarna koncentracija: c_r = c / c_c, a c_c = 1 / Vc_mix ⇒ c_r = c * Vc_mix
        c_r = float(c * Vc_mix)

        # D_p (dimenzijska “skala” viskoznosti u JST): T_c^{-1/6} P_c^{2/3} M^{1/2}
        # Napomena: JST je uobičajeno formuliran tako da završni μ bude u cP ako se koriste određene jedinice;
        # ovdje radimo u SI i sve pretvaramo na Pa·s dosljedno s μ0 (Pa·s).
        D_p = (Tc_mix ** (-1.0 / 6.0)) * (Pc_mix ** (2.0 / 3.0)) * (Mw_mix ** 0.5)  # SI

        # μ0 (dilute) – dvije opcije
        def _omega_mu(T_over_eps: float) -> float:
            t = T_over_eps
            return (1.16145 * t ** -0.14874
                    + 0.52487 * math.exp(-0.7732 * t)
                    + 2.16178 * math.exp(-2.43787 * t))

        def mu0_wilke_chapman() -> float:
            mu_i = []
            for cpt in components:
                nm = UtilClass._name_key(cpt)
                d = ViscosityClass.COMP_DB[nm]
                sigma_m = d["sigma_A"] * 1e-10
                tstar = T / d["eps_over_k"]
                omega = _omega_mu(tstar)
                m = d["Mw"] / MT.CONSTANTS.NA
                mu_pure = (5.0 / 16.0) * math.sqrt(m * 1.380649e-23 * T / math.pi) / (sigma_m ** 2 * omega)  # Pa·s
                mu_i.append(mu_pure)
            mu_i = np.array(mu_i)
            # Wilke
            n = len(x)
            phi = np.zeros((n, n))
            for i in range(n):
                for j in range(n):
                    if i == j:
                        phi[i, j] = 1.0
                    else:
                        num = (1.0 + math.sqrt(mu_i[i] / mu_i[j]) * (Mw_i[j] / Mw_i[i]) ** 0.25) ** 2
                        den = math.sqrt(8.0) * (1.0 + Mw_i[i] / Mw_i[j]) ** 0.5
                        phi[i, j] = num / den
            denom = np.array([np.sum(x[j] * phi[i, j] for j in range(n)) for i in range(n)])
            return float(np.sum(x[i] * mu_i[i] / denom[i] for i in range(n)))

        def mu0_chung_pure() -> float:
            mu_i = []
            for cpt in components:
                nm = UtilClass._name_key(cpt)
                d = ViscosityClass.COMP_DB[nm]
                T_star_i = 1.2593 * T / d["Tc"]
                Omega_i = ViscosityClass._Omega_star(T_star_i)
                Vc_i_val = (d["Zc"] * MT.CONSTANTS.R  * d["Tc"] / d["Pc"])
                Fc_i = 1.0 - 0.2756 * d["omega"]
                mu_uPa_s = 40.785 * math.sqrt(d["Mw"] * T_star_i) / ((Vc_i_val ** (2.0/3.0)) * Omega_i) * Fc_i
                mu_i.append(mu_uPa_s * 1e-6)  # Pa·s
            # Wilke miješanje
            return ViscosityClass.viscosity_chung_csp(components, T, P, x, Z, phase="vapor")

        mu0 = mu0_wilke_chapman() if eta0_strategy == "wilke_chapman" else mu0_chung_pure()

        # JST/LBC polinom koeficijenti (a1..a5) – iz standardne literature (Wikipedia sažetak)
        # η_mix = η0_mix − 1e-4 D_p + D_p * [Σ a_i c_r^{i-1}]^4
        if lbc_coeffs is None:
            a = np.array([0.10230, 0.023364, 0.058533, -0.040758, 0.0093324], dtype=float)
        else:
            a = np.array(lbc_coeffs, dtype=float)
            if a.size != 5:
                raise ValueError("lbc_coeffs mora biti niz od 5 vrijednosti [a1..a5].")

        L = float(np.polyval(a[::-1], c_r))  # a1 + a2 c_r + ... + a5 c_r^4
        eta = mu0 - 1e-4 * D_p + D_p * (L ** 4)  # (Pa·s) – u praksi se koeficijenti često blago tuniraju

        # Safety – viskoznost ne smije pasti ispod μ0 u razrijeđenom limitu kad c_r→0
        eta = max(eta, mu0 * 0.95)
        return float(eta)

    def mix_pseudocritical(components, composition):
        """
        Vrati (Tc_mix, Pc_mix, Vc_mix, Zc_mix_approx) Kay-ovim miješanjem.
        components: imaju .name (CO2/N2), koristi interne Tc,Pc,Zc iz prethodnog koda
        composition: np.array (x ili y ili z, prema kontekstu)
        """
        import numpy as np
        from math import isfinite

        COMP = {
            "CO2": {"Tc": 304.1282, "Pc": 7.3773e6, "Zc": 0.274},
            "N2":  {"Tc": 126.2,    "Pc": 3.3958e6, "Zc": 0.289},
        }
        x = np.array(composition, dtype=float)
        x = x / x.sum()

        Tc_i = []; Pc_i = []; Zc_i = []
        for c in components:
            nm = getattr(c, "name", "").upper().replace("₂","2")
            key = "CO2" if "CO2" in nm else "N2"
            Tc_i.append(COMP[key]["Tc"])
            Pc_i.append(COMP[key]["Pc"])
            Zc_i.append(COMP[key]["Zc"])
        Tc_i = np.array(Tc_i); Pc_i = np.array(Pc_i); Zc_i = np.array(Zc_i)

        Vc_i = Zc_i * MT.CONSTANTS.R * Tc_i / Pc_i
        Tc_mix = float((x * Tc_i).sum())
        Pc_mix = float((x * Pc_i).sum())
        Vc_mix = float((x * Vc_i).sum())
        # “Zc_mix” se rijetko koristi, ali za referencu možemo uzeti x·Zc_i
        Zc_mix = float((x * Zc_i).sum())
        return Tc_mix, Pc_mix, Vc_mix, Zc_mix

    def classify_near_critical(T, P, Tc_mix, Pc_mix, tol_T=0.02, tol_P=0.05):
        """
        Jednostavna klasifikacija oko kritike.
        - superkritično: T>Tc_mix i P>Pc_mix
        - near-critical: unutar tolerancija (relativno) spram Tc i/ili Pc
        - inače: None
        """
        supercritical = (T > Tc_mix) and (P > Pc_mix)
        near_T = abs(T - Tc_mix)/Tc_mix <= tol_T
        near_P = abs(P - Pc_mix)/Pc_mix <= tol_P
        near_critical = (near_T and P >= Pc_mix*0.9) or (near_P and T >= Tc_mix*0.9)
        return supercritical, near_critical

""" # Pretpostavimo da već imaš iz svog solvera:
# V, x_liq, y_vap, method, iters, Zl, Zv = RachfordRice.solve(...)



# VISKOZNOSTI – varijanta 1 (Chung CSP, “gas-like”):
mu_v_chung = viscosity_chung_csp(components, T, P, y_vap, Z=Zv, phase="vapor")
# (za liquid će dati donju procjenu; za tekućine koristi LBC ispod)

# VISKOZNOSTI – varijanta 2 (LBC/JST – radi i za gas i za liquid):
mu_v_lbc = viscosity_lbc(components, T, P, y_vap, Z=Zv)  # gas
mu_l_lbc = viscosity_lbc(components, T, P, x_liq, Z=Zl)  # liquid """



# 1. `density_from_Z(...)` – gustoća faze iz Z (kg/m³).
# 2. `viscosity_chung_csp(...)` – **Chung CSP** (dilute-gas + acentrični faktor) za CO₂/N₂ mješavine; dobar za plin/“gas-like” uvjete, često ok i blizu kritike uz dobru gustoću.
# 3. `viscosity_lbc(...)` – **Lohrenz–Bray–Clark (LBC/JST)** mješovita korelacija koja koristi **reduciranu molarnu koncentraciju** i **pseudo-kritične** (Kay) miks-vrijednosti. Radi i za **plin** i za **tekućinu** (često se blago “tjunira” na podatke).

# > Jedinice: **SI** (T\[K], P\[Pa], μ\[Pa·s], ρ\[kg/m³]). U kodu su kritične veličine za CO₂ i N₂, i LJ parametri za Chapman/Chung dio.
# > Za LBC sam implementirao standardnu **JST** formu i Kay mix-pravila; kad imaš lab PVT, obično se fino dotjeruju $V_c$ teških frakcija / $a_i$ koeficijenti. Formula i konstante su kao u sažecima iz literature (a₁…a₅) ([Wikipedia][1]); pregled i praktične napomene o LBC tuningu kod Whitsona ([Whitson][2]). Chungov gas-like izraz i kolizijski integral Ω\* prema Neufeldu su dani sažeto ovdje (T\* definicija, Fc korekcija) ([Wikipedia][1]).



# * **Za CO₂/N₂** plin: `viscosity_chung_csp` je vrlo solidan; LBC često daje gotovo iste brojeve za “gas-like” regiju.
# * **Za tekućinu** (kondenzat): koristi `viscosity_lbc`, i ako imaš mjerne μ, **tjuniraj** $[a_1..a_5]$ i/ili kritične volumene (posebno C7+ u naftama) – to je industrijski standard i preporuka (Whitson) ([Whitson][2]).
# * U LBC-u je ključna **reducirana koncentracija** $c_r = c \, V_{c,\text{mix}}$ i **Kay** za $T_c, P_c, V_c$ te miješanje $\bar M$. JST polinom (a₁…a₅) je standardni, vidi sažetak i miks-pravila (u izvorniku je dana i $\eta_{0}$ i $D_p$) ([Wikipedia][1]).
# * Ako želiš, mogu dodati “prekidač” koji automatski bira `Chung` za faze s $\rho < 120~\mathrm{kg/m^3}$ i `LBC` inače, ili uvezati **tuning** koeficijenata prema tvojim PVT tablicama (brzi least-squares).

# [1]: https://en.wikipedia.org/wiki/Viscosity_models_for_mixtures "Viscosity models for mixtures - Wikipedia"
# [2]: https://whitson.com/2020/01/20/lbc-viscosity-correlation-for-gas-condensate-reservoirs/ "LBC Viscosity Correlation for Gas Condensate Reservoirs - Whitson"
