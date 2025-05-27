class Component:
    """
    Pure component
    """
    def __init__(self, name: str, Tc: float, Pc: float, omega: float, AntoineA: float, AntoineB: float, AntoineC: float,
                 CpA: float, CpB: float, CpC: float, CpD: float):
        self.name = name                    # Name (npr. Methane)
        self.Tc = Tc                        # (K)
        self.Pc = Pc                        # (Pa)
        self.omega = omega                  # Accentric
        self.AntoineA = AntoineA            # Antoine A
        self.AntoineB = AntoineB            # Antoine B
        self.AntoineC = AntoineC            # Antoine C
        self.CpA = CpA                      # Cp A  J/mol·K
        self.CpB = CpB                      # Cp B  J/mol·K
        self.CpC = CpC                      # Cp C  J/mol·K
        self.CpD = CpD                      # Cp D  J/mol·K

    def saturation_pressure(self, T: float) -> float:
        """
        Antoine eq.
        """
        return 10**(self.A - self.B / (T + self.C))
    
class ComponentNew:
    """
    Pure component with extended properties for petroleum applications.
    """
    def __init__(self, name: str, formula: str, Mw: float, Tc: float, Pc: float, omega: float, fraction: float, CpA: float, CpB: float, CpC: float, CpD: float):
        self.name = name              # e.g., "Propane"
        self.formula = formula        # e.g., "C3"
        self.Mw = Mw                  # Molecular weight [g/mol]
        self.Tc = Tc                  # Critical temperature [K]
        self.Pc = Pc                  # Critical pressure [bar]
        self.omega = omega            # Acentric factor
        self.fraction = fraction      # Mole fraction
        self.CpA = CpA                # Cp constant A [J/mol·K]
        self.CpB = CpB                # Cp constant B [J/mol·K²]
        self.CpC = CpC                # Cp constant C [J/mol·K³]
        self.CpD = CpD                # Cp constant D [J·K/mol]    