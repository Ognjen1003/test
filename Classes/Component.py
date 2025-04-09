class Component:
    """
    Pure component
    """
    def __init__(self, name: str, Tc: float, Pc: float, omega: float, A: float, B: float, C: float):
        self.name = name      # Name (npr. Methane)
        self.Tc = Tc          # (K)
        self.Pc = Pc          # (Pa)
        self.omega = omega    # Accentric
        self.A = A            # Antoine A
        self.B = B            # Antoine B
        self.C = C            # Antoine C

    def saturation_pressure(self, T: float) -> float:
        """
        Antoine eq.
        """
        return 10**(self.A - self.B / (T + self.C))