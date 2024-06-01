import numpy as np

freq = np.arange(50, 10e3, 10)
omega = 2 * np.pi * freq
c0 = 343.0
rho0 = 1.213
z0 = rho0 * c0
k0 = omega / c0


def alpha_from_dkz(d, k, z):
    zs = -1j * z / np.tan(k * d)
    R = (zs - z0) / (zs + z0)
    alpha = 1 - np.abs(R) ** 2
    return alpha


class Fluid:
    def __init__(
        self,
        T0=20.0,
        P0=101325.0,
        R=287.05,
        adia=1.4,
        Tref=291.15,
        C=120.0,
        Eta0=18.07e-6,
    ):
        """Instantiate fluid material with given properties

        Args:
            T0 (float, optional): temperature. Defaults to 20.0.
            P0 (float, optional): static pressure. Defaults to 101325.0.
            R (float, optional): gas constant. Defaults to 287.05.
            adia (float, optional): adiabatic index. Defaults to 1.4.
            Tref (float, optional): reference temperature (for dyn viscosity). Defaults to 291.15.
            C (float, optional):   Sunderland constant. Defaults to 120.0.
            Eta0 (float, optional): Reference viscosity. Defaults to 18.07e-6.
        """        
        self.TK = T0 + 273.15
        self.rho = P0 / (R * (self.TK))
        self.Eta = (
            Eta0 * ((Tref + C) / (self.TK + C)) * ((self.TK) / Tref) ** (3.0 / 2.0)
        )
        self.Cp = 1e3 * (
            1.0062 + 3.6028e-5 * T0 - 1.0855e-6 * T0**2 + 1.3791e-8 * T0**3
        )
        self.kappa = (
            1.5207e-11 * (self.TK) ** 3
            - 4.8574e-8 * (self.TK) ** 2
            + 1.0184e-4 * (self.TK)
            - 3.9333e-4
        )
        self.co = np.sqrt(adia * (R * (self.TK)))
        self.c = self.co
        self.speedOfSound = self.co
        self.density = self.rho
        self.adia = adia
        self.gamma = adia
        self.P0 = P0
        self.Pr = self.Cp * self.Eta / self.kappa


def rho_c_JCA(freq, phi, sigma, alpha, v_length, t_length, *args, **kw):
    """Effective density and speed of sound for a porous material using the Johnson-Allard-Champoux model

    Args:
        freq (float): frequency (Hz)
        phi (float): porosity 
        sigma (float): flow resistivity (Rayl/m)
        alpha (float): tortuosity
        v_length (float): viscous characteristic length (m)
        t_length (float): thermal characteristic length (m)

    Returns:
        tuple: rho, c
    """    
    f = Fluid(*args, **kw)
    omega = 2 * np.pi * freq

    Pr = f.Cp * f.Eta / f.kappa
    G_J = (
        1
        + (4j * alpha**2.0 * f.Eta * f.rho * omega)
        / (sigma**2.0 * v_length**2.0 * phi**2.0)
    ) ** 0.5
    rho = alpha * f.rho / phi * (1 + sigma * phi / (1j * f.rho * omega * alpha) * G_J)
    K = (
        (f.gamma * f.P0)
        / phi
        / (
            f.gamma
            - (f.gamma - 1)
            * (
                1
                + (8 * f.Eta)
                / (1j * t_length**2.0 * Pr**2.0 * omega * f.rho)
                * (1 + 1j * f.rho * (omega * Pr**2.0 * t_length**2.0) / (16 * f.Eta))
                ** 0.5
            )
            ** -1
        )
    )
    c = np.sqrt(K / rho)

    return rho, c


def DBM(f, d, sigma):
    """Sound ansorption from the Johnson-Allard-Champoux model

    Args:
        freq (float): frequency (Hz)
        phi (float): porosity 
        sigma (float): flow resistivity (Rayl/m)
        alpha (float): tortuosity
        v_length (float): viscous characteristic length (m)
        t_length (float): thermal characteristic length (m)

    Returns:
        float[]: sound absorption coefficient by frequency
    """    
    z = z0 * (
        1
        + 5.5 * (1e3 * f / sigma) ** (-0.632)
        - 1j * 8.43 * (1e3 * f / sigma) ** (-0.632)
    )
    k = k0 * (
        1
        + 7.81 * (1e3 * f / sigma) ** (-0.618)
        - 1j * 11.41 * (1e3 * f / sigma) ** (-0.618)
    )

    alpha = alpha_from_dkz(d, k, z)

    return alpha


def JCA(f, d, phi, sigma, tort, v_length, t_length):
    """_summary_

    Args:
        f (_type_): _description_
        d (_type_): _description_
        phi (_type_): _description_
        sigma (_type_): _description_
        tort (_type_): _description_
        v_length (_type_): _description_
        t_length (_type_): _description_

    Returns:
        _type_: _description_
    """    
    rho, c = rho_c_JCA(freq, phi, sigma, tort, v_length, t_length)
    k = f * np.pi * 2 / c
    z = rho * c
    alpha = alpha_from_dkz(d, k, z)
    return alpha
