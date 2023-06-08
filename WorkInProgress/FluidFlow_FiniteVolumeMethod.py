import numpy as np
from matplotlib import pyplot

def getConserved(rho, vx, vy, P, gamma, vol):
    Mass = rho * vol
    Momx = rho * vx * vol
    Momy = rho * vy * vol
    Energy = (P/(gamma-1) + 0.5*rho*(vx**2+vy**2))*vol

    return Mass, Momx, Momy, Energy

def getPrimitive(Mass, Momx, Momy, Energy, gamma, vol):
    rho = Mass / vol
    vx = Momx / vol / rho
    vy = Momy / vol / rho
    P = (Energy/vol - 0.5*rho * (vx**2+vy**2)) * (gamma-1)

    return rho, vx, vy, P

def getGradient(f, dx):
    f_dx = (f[2:, 1:-1] - f[:-2, 1:-1])/(2*dx)
    f_dy = (f[1:-1, 2:] - f[1:-1, :-2])/(2*dx)

    return f_dx, f_dy

def extrapolateInSpaceToFace(f, f_dx, f_dy, dx):

  R = -1 
  
  f_XL = f - f_dx * dx/2
  f_XL = np.roll(f_XL,R,axis=0)
  f_XR = f + f_dx * dx/2
  
  f_YL = f - f_dy * dx/2
  f_YL = np.roll(f_YL,R,axis=1)
  f_YR = f + f_dy * dx/2
  
  return f_XL, f_XR, f_YL, f_YR

def getFlux(rho_L, rho_R, vx_L, vx_R, vy_L, vy_R, P_L, P_R, gamma):
    en_L = P_L/(gamma-1)+0.5*rho_L * (vx_L**2+vy_L**2)
    en_R = P_R/(gamma-1)+0.5*rho_R * (vx_R**2+vy_R**2)

    rho_star  = 0.5*(rho_L + rho_R)
    momx_star = 0.5*(rho_L * vx_L + rho_R * vx_R)
    momy_star = 0.5*(rho_L * vy_L + rho_R * vy_R)
    en_star   = 0.5*(en_L + en_R)
    P_star = (gamma-1)*(en_star-0.5*(momx_star**2+momy_star**2)/rho_star)

    flux_Mass   = momx_star
    flux_Momx   = momx_star**2/rho_star + P_star
    flux_Momy   = momx_star * momy_star/rho_star
    flux_Energy = (en_star+P_star) * momx_star/rho_star

    C_L = np.sqrt(gamma*P_L/rho_L) + np.abs(vx_L)
    C_R = np.sqrt(gamma*P_R/rho_R) + np.abs(vx_R)
    C = np.maximum( C_L, C_R )

    flux_Mass   -= C * 0.5 * (rho_L - rho_R)
    flux_Momx   -= C * 0.5 * (rho_L * vx_L - rho_R * vx_R)
    flux_Momy   -= C * 0.5 * (rho_L * vy_L - rho_R * vy_R)
    flux_Energy -= C * 0.5 * ( en_L - en_R )



