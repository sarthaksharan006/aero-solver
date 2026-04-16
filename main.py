import math
import numpy as np
import time
import matplotlib.pyplot as plt

# wing airfoil
def s1223():
    global R, x_max, alpha_stall, C_Lmax, alpha_zerolift, wing
    wing= "S1223"
    R = 0.121
    x_max = 0.198
    alpha_stall = 15
    C_Lmax = 2.2
    alpha_zerolift = -5

def clarky():
    global R, x_max, alpha_stall, C_Lmax, alpha_zerolift, wing
    wing = "clark-y"
    R = 0.117
    x_max = 0.28
    alpha_stall = 13
    C_Lmax = 1.3
    alpha_zerolift = -4

# vt airfoil
def naca0012_vt():
    global R_vt, x_maxvt
    R_vt = 0.122
    x_maxvt = 0.225

# ht airfoil
def naca0012_ht():
    global R_ht, x_maxht
    R_ht = 0.122
    x_maxht = 0.225

s1223()
naca0012_ht()
naca0012_vt()

# constants
u = 340
rho = 1.2
mu = 1.8*10**(-5)
RW = 0.8 #dry weight
PW = 1.35 #payload
W = (RW+PW)*9.81
A_max = 0.11*0.11
la = 1 #taper ratio
la_ht = 1 #lambda #graphically optimised
La = 0 # sweep
La_ht = 0
b = 1.2
T = 1.2*RW*9.81 

# to optimise, placeholders
r = 0.28
L_fuse = 1 #graphically optimised for balance b/w minimum weight of tail and fuselage and drag
r_ht = 0.2
r_vt = 0.2
x_p = 0.15*b

#### plane body
t = la*r
c = (t+r)/2
t_ht = la_ht*r_ht
c_ht = (t_ht+r_ht)/2

d_w = x_p+0.25*c

tail = "conventional"
if tail == "t":
    la_vt = 1
    t_vt = la_vt*r_vt
    c_vt= (t_vt+r_vt)/2
    l_ht = L_fuse-d_w-la_vt*r_vt+0.25*c_ht
    La_vt = math.radians(45)
else:
    la_vt = la_ht
    t_vt = la_vt*r_vt
    c_vt= (t_vt+r_vt)/2
    l_ht = L_fuse-d_w-0.75*c_vt
    La_vt = 0

l_vt = L_fuse-d_w-0.75*c_vt

MAC = (2/3)*( (1+la+la**2)/(1+la) )*r
MAC_ht = (2/3)*( (1+la_ht+la_ht**2)/(1+la_ht) )*r_ht
MAC_vt = (2/3)*( (1+la_vt+la_vt**2)/(1+la_vt) )*r_vt

S = b*c

Vol_ht = 0.65
Vol_vt = 0.06
S_ht = S*MAC*Vol_ht/l_ht
S_vt = S*b*Vol_vt/l_vt

b_ht = S_ht/c_ht
b_vt = S_vt/c_vt

AR = b**2/S
AR_ht = b_ht**2/S_ht
AR_vt = b_vt**2/S_vt

t = la*r
t_ht = la_ht*r_ht
t_vt = la_vt*r_vt



N_e = 1
k = 10**(-4)
x = 0.1
Q = 1.1

## Forces
def calcC_D0(L_fuse, V): # drag coefficient
    M = V/u
    ### Reynolds no.
    R_ewing = c*V*rho/mu
    R_efuse = L_fuse*V*rho/mu
    R_eht = c_ht*V*rho/mu
    R_evt = c_vt*V*rho/mu
    R_ecutoff = 38.21*L_fuse/k**1.053
    R_eeff = min(R_efuse, R_ecutoff)

    # wing
    FFw = (1 + 0.6*R/x_max + 100*(R**4)) * (1.34*M**0.18*(math.cos(La))**0.28)
    Cfw = x*(1.328/R_ewing**0.5) + (1-x)*( 0.455*(1+0.144*M**2)**(-0.65) / (math.log10(R_ewing))**2.58 )
    Swetw = 2*S*(1+0.25*R)
    # fuselage
    f = L_fuse/(4*A_max/math.pi)**0.5
    d = 2*(A_max/math.pi)**0.5
    FFf = 1+60/f**3+f/400
    Cff = ( 0.455*(1+0.144*M**2)**(-0.65) / (math.log10(R_eeff))**2.58 )
    Swetf = math.pi*d*L_fuse*(1-0.35*d/L_fuse)
    # ht
    FFht = (1 + 0.6*R_ht/x_maxht + 100*(R_ht**4)) * (1.34 * M**0.18 * (math.cos(La_ht))**0.28)
    Cfht = x*(1.328/R_eht**0.5) + (1-x)*( 0.455*(1+0.144*M**2)**(-0.65) / (math.log10(R_eht))**2.58 )
    Swetht = 2*S_ht*(1+0.25*R_ht)
    # vt
    FFvt = (1 + 0.6*R_vt/x_maxvt + 100*(R_vt**4)) * (1.34 * M**0.18 * (math.cos(La_vt))**0.28)
    Cfvt = x*(1.328/R_evt**0.5) + (1-x)*( 0.455*(1+0.144*M**2)**(-0.65) / (math.log10(R_evt))**2.58 )
    Swetvt = 2*S_vt*(1+0.25*R_vt)

    C_D0wing = FFw*Cfw*Swetw*Q/S
    C_D0fuselage = FFf*Cff*Swetf*Q/S
    C_D0ht = FFht*Cfht*Swetht*Q/S
    C_D0vt = FFvt*Cfvt*Swetvt*Q/S

    C_D0 = C_D0wing + C_D0fuselage + C_D0ht + C_D0vt
    return C_D0

def findV(V): #iteratively converging cruise velocity
    V_new=0
    while True:
        M = V/u
        e_0 = 1 / ( (1 + 0.12*M**6) * (1 + (0.142 + 0.005*(1 + 1.5*(la - 0.6)**2)*AR*(10*R)**0.33) / (math.cos(La)**2) + 0.1*(3*N_e + 1)/(4 + AR)**0.8) )
        C_D0 = calcC_D0(L_fuse, V)
        if C_D0>=0:
            V_new = ((2*W/(rho*S))**0.5) * (1/(C_D0 * math.pi * AR * e_0) )**0.25
        else:
            break
        if abs(V_new-V)<1e-3:
            break
        V = V_new
    return V
V = V_c = findV(10) # real V
M = V/u

e_0 = 1 / ( (1 + 0.12*M**6) * (1 + (0.142 + 0.005*(1 + 1.5*(la - 0.6)**2)*AR*(10*R)**0.33) / (math.cos(La)**2) + 0.1*(3*N_e + 1)/(4 + AR)**0.8) )
R_ewing = c*V*rho/mu
C_D0 = calcC_D0(L_fuse, V) #real CD0


C_Lcruise = W/(0.5*rho*(V**2)*S)
C_L = C_Lcruise
Lf = 0.5*rho*(V**2)*S*C_L

C_Dind = C_L**2/(math.pi*AR*e_0)
C_D = C_D0 + C_Dind
D = 0.5*rho*(V**2)*S*C_D
Vs = (2*W/(rho*S*C_Lmax))**0.5
WL = W/S
# Design after computation
epsilon = C_L/math.pi*AR*e_0
alpha_stall = alpha_stall+epsilon
alpha_zerolift = alpha_zerolift+epsilon
#theta = C_Lcruise*( (10+18*math.cos(La))/(AR*math.cos(La)) ) + alpha_zerolift
wash = 2
goal = Lf/D #maximise
print(f"ALL SI UNITS unless specified otherwise")
print("\n=== AIRFOIL PARAMETERS ===")
print(f"Wing ({wing}): R={R}, x_max={x_max}, downwash={epsilon:.4f}, effective alpha_stall={alpha_stall:.4f}, C_Lmax={C_Lmax}, effective alpha_zerolift={alpha_zerolift:.4f}")
print(f"HT (NACA0012): R_ht={R_ht}, x_maxht={x_maxht}")
print(f"VT (NACA0012): R_vt={R_vt}, x_maxvt={x_maxvt}")

print("\n=== GEOMETRIC PARAMETERS ===")
print(f"Wing: S={S:.4f}, b={b:.4f}, c={c:.4f}, MAC={MAC:.4f}, AR={AR:.4f}")
print(f"      Taper ratio (λ)={la:.4f} , Sweep (Λ)={La:.4f}°, Washout={wash}° ,Root chord={r:.4f}, Tip chord={t:.4f}")
print(f"HT: S_ht={S_ht:.4f}, b_ht={b_ht:.4f}, c_ht={c_ht:.4f}, MAC_ht={MAC_ht:.4f}, AR_ht={AR_ht:.4f}, Volume_ht={Vol_ht:.4f}")
print(f"    Taper ratio (λ)={la_ht:.4f}, Sweep (Λ)={La_ht:.4f}°, Root chord={r_ht:.4f}, Tip chord={t_ht:.4f}")
print(f"VT type(t-tail/conventional)={tail}")
print(f"VT: S_vt={S_vt:.4f}, b_vt={b_vt:.4f}, c_vt={c_vt:.4f}, MAC_vt={MAC_vt:.4f}, AR_vt={AR_vt:.4f},  Volume_vt={Vol_vt:.4f}")
print(f"    Taper ratio (λ)={la_vt:.4f}, Sweep (Λ)={math.degrees(La_vt):.4f}°, Root chord={r_vt:.4f}, Tip chord={t_vt:.4f}")
print(f"Tail moment arms: l_ht={l_ht:.4f} m, l_vt={l_vt:.4f} m for L_fuse={L_fuse:.4f} m")
print(f"Fuselage: L_fuse={L_fuse:.4f} m, Prop to LE dist: x_p={x_p:.4f} m")

print("\n=== FLIGHT CONDITIONS ===")
print(f"Cruise velocity: V={V:.4f} m/s")
print(f"Stall velocity: Vs={Vs:.4f} m/s")
print(f"Mach number: M={M:.4f}")
print(f"Raw Weight: RW={RW:.4f} kg")
print(f"Total Weight with payload: W={W/9.81:.4f} kg")
print(f"Reynolds number: R_ewing={R_ewing:.4f}")
print(f"Wing Loading: WL={WL/9.81:.4f} kg/m^2")


print("\n=== AERODYNAMIC COEFFICIENTS ===")
print(f"C_D0: {C_D0:.6f}")
print(f"C_L (cruise): {C_L:.6f}")
print(f"C_Dind: {C_Dind:.6f}")
print(f"C_D (total): {C_D:.6f}")
print(f"Oswald efficiency (e_0): {e_0:.6f}")
#print(f"Optimal incidence angle: {theta:.6f}")

print("\n=== FORCES ===")
print(f"Lift at cruise: Lf={Lf:.4f} N")
print(f"Drag at cruise: D={D:.4f} N")
print(f"Lift/Drag ratio: {goal:.4f}")
print(f"Max allowed Thrust: T={T:.4f} N")
print(f"Thrust/Drag ratio: {T/D:.4f}")
