"""This code is made to be used in undergraduation heat transfer classes, so
     in this context there is a lot of simplification assumptions often made 
     (e.g. negleceting fouling factors), and due to this, some arguments may
     take unrealistic default values (e.g. Rf = 0)"""
     
import math
def heatRate(m:float,cp:float,dT:float):
    
    """Calculates heat transfer rate of heat transfer [w]

    Args:
        m (float): Mass flow rate of the fluid [Kg/s].
        cp (float): Specific heat of the fluid [J/Kg.K].
        dT (float): Temerature differance (dT = T1 - T2) [C or K].

    Returns:
        float: Heat Transfer Rate (Q).
    """
    Q = m*cp*dT
    return Q

def LMTD(Th1:float, Th2:float, Tc1:float, Tc2:float, flowType:str):

    """Calculates Logarithmic Mean Temparature Differanse [C].

    Args:
        Th1 (float): Hot fluid inlet temperature [C].
        Th2 (float): Hot fluid outlet temperature [C].
        Tc1 (float): Cold fluid inlet temperature [C].
        Tc2 (float): Cold fluid outlet temperature [C].
        flowType (str): Parallel or counter.

    Returns:
        float : Logarithmic Mean Temparature Differanse
    """
    flowType = flowType.lower()
    if flowType == 'parallel':
        T1 = Th1 - Tc2
        T2 = Th1 - Tc1
    elif flowType == 'counter':
        T1 = Th1 - Tc2
        T2 = Th2 - Tc1
    else:
        print("error: Invalid option!")
        return None

    if T1 == T2:
        return T1
    try:
        return (T1 - T2) / math.log(T1 / T2)
    except (ZeroDivisionError, ValueError):
        print("error: Invalid temperature values!")
        return None
    
def Reynold(rho:float, V:float, D:float, mio:float):
    """Calculates Reynold's number.

    Args:
        rho (float): Density of the fliud [Kg/m^3].
        V (float): Velosity of the fluid [m/s].
        D (float): Tube(s) diameter [m].
        mio (float): Dynaminc viscosity of the fluid [N/m.s]

    Returns:
        float: Reynold Number
    """
    try:
      return (rho * V * D) / mio
    except (ZeroDivisionError):
      print('error: Vescosity cant be zero!')
      return None

def Prandtl(mio, Cp, k):
    """Calculates Prandtl Number.

    Args:
        mio (float): Dynaminc viscosity of the fluid [N/m.s].
        Cp (float): Specific heat of the fluid [J/Kg.K].
        k (float): Thermal conductivity of tube wall [w/m.K]

    Returns:
        float: Prandtl Number
    """
    try:
       Pr =  (mio * Cp) / k
       return Pr
    except (ZeroDivisionError):
       print('error: Conductivity cant be zero!')
       return None

def dittus(Re:float, Pr:float):
    """Calculates nusselt number using Dittus-Boelter equation.

    Args:
        Re (float): Reynold Number
        Pr (float): Prantdl Number

    Returns:
        float : Nusselt Number
    """
    Nu = 0.023 * Re**0.8 * Pr**0.4
    return Nu

def ConvCoeff(Nu:float, k:float, D:float):
    """Calculates convection heat transfer coefficient.
    
    Args:
        Nu (float): Nusselt number.
        k (float): Thermal Conductivity [w/m.K]
        D (float): Tube Diameter [m].

    Returns:
        float: Convection heat transfer coefficient
    """
    if D == 0:
        raise ZeroDivisionError('Diameter can not be zero.')
    else:
        h = (Nu * k) / D
        return h


def TransferArea(D:float, L:float, N=1):
    """Calculates heat transfer area.

    Args:
        D (float): Tube(s) diameter [m].
        L (float): Length of the tube(s) [m].
        N (int, optional): Number of tubes. Defaults to 1.

    Returns:
        float: Heat Transfer Area
    """
    A = N * math.pi * D * L
    return A

    
def OverallCoeff(hi:float, ho:float, Di:float=1, Do:float=1, k:float=1, Rfi:float=0, Rfo:float=0):
    """Calculates overall heat transfer coeffecient based on 
        inner and outer surface area [w/m^2*k].

    Args:
        hi (float): Convection coefficient in tube side  [w/m^2.K].
        ho (float): Convection coefficient in the annulus [w/m^2.K].
        Di (float, optional): Inner dimeter [m]. Defaults to 1.
        Do (float, optional): Outer diameter [m]. Defaults to 1.
        k (float, optional): Conduction coefficient of tube walls [w/m.K]. Defaults to 1.
        Rfi (float, optional): fouling factor on the inner surface [m^2.K/w]. Defaults to 0.
        Rfo (float, optional): fouling factor on the outer surface [m^2.K/w]. Defaults to 0.

    Returns:
        float: Ui, Uo
    """
    Ui = 1/((1/hi)+Rfi+((Di/k)*(math.log(Do/Di)))+((Di/Do)*Rfo)+((Di/Do)*(1/ho)))
    Uo = 1/(((Do/Di)*(1/hi))+((Do/Di)*Rfo)+((Do/k)*math.log(Do/Di))+Rfo+(1/ho))
    return Ui, Uo


def CorrectionFactor(Th1:float, Th2:float, Tc1:float, Tc2:float, N:int=1):
    
    """Calculates TEMA correction factor for shell and tube heat exchanger.

    Args:
        Th1 (float): Hot fluid inlet temperature [C].
        Th2 (float): Hot fluid outlet temperature [C].
        Tc1 (float): Cold fluid inlet temperature [C].
        Tc2 (float): old fluid outlet temperature [C].
        N (int, optional): Number of shell side passes. Defaults to 1.

    Returns:
        float : Correction factor (F)
    """
    P = (Tc2-Tc1)/(Th1-Tc1) #The temperature ratio
    R = (Th1-Th2)/(Tc2-Tc1) #The capacity ratio
    if R != 1:
       alpha = (( 1- P*R)/(1 - P))**(1/N)
       S = (alpha - 1)/(alpha - R)
       numerator = math.sqrt(R**2 + 1)*math.log((1 - S)/(1 - R*S))
       denomerator = (R-1)*math.log((2-S*(R + 1 -math.sqrt(R**2+1)))/(2-S*(R+1+math.sqrt(R**2 + 1))))
       F = numerator / denomerator
    else:
       S = P/(N*(N - 1)*P)
       F = (S*math.sqrt(2))/(1 - S)*math.log((2-S*(2 - math.sqrt(2))/(2-S+math.sqrt(2))))
    return F

  
def effectiveness(mh:float, mc:float, cph:float, cpc:float, U:float, A:float, ExchangerType:str, n=1):
    """Calculates effectiveness of heat exchanger.

    Args:
        mh (float): Mass flow rate of the hot fluid [Kg/s].
        mc (float): Mass flow rate of the cold fluid [Kg/s].
        cph (float): Specific heat of the hot fluid [J/kg.K].
        cpc (float): specific heat of the cold fluid [J/Kg.K].
        U (float): Overall heat transfer coefficient [w/m^2.K].
        A (float): Heat transfer area [m^2].
        ExchangerType (str): The type of the heat exchanger
           (parallel, counter, 1-shell, n-shell, cross_unmix, cross_mix,
           evaporator, condensor).
        n (int, optional, Defaults to 1): n have two values:
           1. n = number of shell passes in shell and tube heat exchanger.
           2. n = 0 or 1 in cross flow heat exchangers (one fluid mixed, one unmixed).
                = 1 if Cmax mixed, Cmin unmixed.
                = 0 if Cmax unmixed, Cmin mixed.
             Other than this two conditions it has no effect.

    Returns:
        float: Effectiveness.
    """
    Ch = mh * cph
    Cc = mc * cpc
    Cmin = min(Ch, Cc)
    Cmax = max(Ch, Cc)
    R = Cmin / Cmax
    NTU = (U * A) / Cmin;

    ExchangerType = ExchangerType.lower()
    if ExchangerType == 'parallel':
        if R == 1:
           e = (1 - math.exp(-2 * NTU)) / 2
        else:
            e = (1 - math.exp(-NTU * (1 + R))) / (1 + R)
        return e
    elif ExchangerType == 'counter':
        if R == 1:
            e = NTU / (1 + NTU)
        else:
            e = (1 - math.exp(-NTU * (1 - R))) / (1 - R * math.exp(-NTU * (1 - R)))
        return e
    elif ExchangerType == '1-shell': # Shell and tube with one shell pass
        sqrt_term = math.sqrt(1 + R**2)
        exp_term = math.exp(-NTU * sqrt_term);
        e = 2 / (1 + R + sqrt_term * (1 + exp_term) / (1 - exp_term))
        return e
    elif ExchangerType == 'n-shell': # Shell and tube with n-shell passes
        sqrt_term = math.sqrt(1 + R**2)
        exp_term = math.exp(-NTU * sqrt_term / n)
        e1 = 2 / (1 + R + sqrt_term * (1 + exp_term) / (1 - exp_term))
        term = (1 - e1 * R) / (1 - e1)
        e = (term**n - 1) / (term**n - R)
        return e
    elif ExchangerType == 'cross_unmix': # Cross-flow (both fluids unmixed)
         e = 1 - math.exp((NTU**0.22 / R) * (math.exp(-R * NTU**0.78) - 1))
         return e
    elif ExchangerType == 'cross_mix': # Cross-flow (one mixed, one unmixed)
         if n == 1 : # Cmax mixed, Cmin unmixed
            e = (1 - math.exp(-R * (1 - math.exp(-NTU)))) / R
         else:     # Cmin mixed, Cmax unmixed'
            e = 1 - math.exp(-(1 - math.exp(-R * NTU)) / R)
         return e
    elif ExchangerType in ('condensor', 'evaporator'): 
         e = 1 - math.exp(-NTU)
         return e
    else:
        print('error: Exchanger type not found.')
        return None
        