import numpy as np
import pylab
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.spatial.transform import Rotation
from src.airProperties import AirProperties
from src.aeroDynamic import AeroDynamic

INITIAL_MASS = 26.409  #kg
FINALMASS = 20.4338 #kg
BURNINGOUTTIME = 15.0 #s
GROUND_GRAVITY = 9.81  #m/s^2
THRUST = 100 * 9.81 #N
HEIGHT = 100 #km
EARTH_RADIUS = 6371 #km
GRAVITY = GROUND_GRAVITY
MDOT = (INITIAL_MASS - FINALMASS)/BURNINGOUTTIME #kg/s

INITIALMASSCENTER = 0.2373 #m
FINALMASSCENTER = 0.2308 #m

MAC = 0.2573 

MASSCENTERMOVINGRATE = (INITIALMASSCENTER - FINALMASSCENTER)/BURNINGOUTTIME #m/s



INITIX = 0.1697
INITIY = 6.5455
INITIZ = 6.5455

FINALIX = 0.1354
FINALIY = 6.3605
FINALIZ = 6.3605

IXDOT = (INITIX - FINALIX)/BURNINGOUTTIME
IYDOT = (INITIY - FINALIY)/BURNINGOUTTIME
IZDOT = (INITIZ - FINALIZ)/BURNINGOUTTIME

TVCPOINT = 0.456+INITIALMASSCENTER
LENGTHARM = 0.456 #m
WIDTHARM = 0 #m

THETA_SP = np.deg2rad(10)
PSI_SP = np.deg2rad(5)
P_SP = 0

Q0_SP = np.cos(PSI_SP/2)*np.cos(THETA_SP/2)
Q1_SP = -np.sin(PSI_SP/2)*np.sin(THETA_SP/2)
Q2_SP = np.cos(PSI_SP/2)*np.sin(THETA_SP/2)
Q3_SP = np.sin(PSI_SP/2)*np.cos(THETA_SP/2)




KP1THETA = 0.02
KP2THETA = -5
KITHETA = -0.15

KP1PSI = 0.02
KP2PSI = -3
KIPSI = -0.04

KP2PHI = 50
KIPHI = 0.04


# KP1THETA = 0.2
# KP2THETA = -0.005
# KITHETA = -0.0000015
# KP1PSI = 0.2 
# KP2PSI = -0.5
# KIPSI = -0.005

# KP2PHI = 5
# KIPHI = 0.06

def quaternion_mutiply(p0,p1,p2,p3,q0,q1,q2,q3):
    return [p0*q0-(p1*q1+p2*q2+p3*q3), p0*q1+q0*p1+(p2*q3-p3*q2), p0*q2+q0*p2+(p3*q1-p1*q3), p0*q3+q0*p3+(p1*q2-p2*q1)]

def quaternion_rotation(q0,q1,q2,q3):
    theta = np.arccos(q0)
    sintheta = np.arsin(theta)
    return 2*np.arccos(q0), [q1/sintheta, q2/sintheta, q3/sintheta]

control=[]

def f(t,y):
    thrust = 100 * 9.81
######################################################system#############################################################################

    p, q, r, u, v, w, q0, q1, q2, q3,Xinertial,Yinertial,Zinertial,mass,masscenter,ix,iy,iz = y

    
    velocityBody = np.array([[u],[v],[w]])
    T = np.array([[q0**2 + q1**2 - q2**2 - q3**2, 2*(q1*q2 + q0*q3), 2*(q1*q3 - q0*q2)]
                  ,[2*(q1*q2 - q0*q3), q0**2 - q1**2 + q2**2 - q3**2, 2*(q2*q3 + q0*q1)]
                  ,[2*(q1*q3 + q0*q2), 2*(q2*q3 - q0*q1), q0**2 - q1**2 - q2**2 + q3**2]])
    eta = np.deg2rad(0)
    zeta =np.deg2rad(0)
    mtp = 0
    mtq = - (TVCPOINT - masscenter) * thrust * np.sin(eta)
    mtr = - (TVCPOINT - masscenter) * thrust * np.cos(eta) * np.sin(zeta)
    

    # TVC thrust
    ftx = thrust * np.cos(eta)*np.cos(zeta)
    fty = thrust * np.cos(eta)*np.sin(zeta)
    ftz = - thrust * np.sin(eta)

    # aerodynamiv force
    air = AirProperties(-Zinertial)
    aerodynamic = AeroDynamic()
    velocityMag = np.linalg.norm(velocityBody)
    AoA = np.arctan(w/u)
    beta= np.arcsin(v/velocityMag)

    dynamicPressure = 0.5*air.density*velocityMag**2

    # aerodynamic force from angle of attack
    Liftaoa = dynamicPressure*aerodynamic.CLalpha*aerodynamic.Sref
    Dragaoa = dynamicPressure*aerodynamic.CDalpha*aerodynamic.Sref
    sideSlipaoa = dynamicPressure*aerodynamic.CYalpha*aerodynamic.Sref

    fax_aoa = -(-Liftaoa*np.sin(AoA) + Dragaoa*np.sin(AoA))
    fay_aoa = sideSlipaoa
    faz_aoa = -(Liftaoa*np.cos(AoA) + Dragaoa*np.sin(AoA))

    # aerodynamic force from side slip angle

    Liftbeta = dynamicPressure*aerodynamic.CLbeta*aerodynamic.Sref
    Dragbeta = dynamicPressure*aerodynamic.CDbeta*aerodynamic.Sref
    sideSlipbeta = dynamicPressure*aerodynamic.CYbeta*aerodynamic.Sref

    fax_beta = sideSlipbeta*np.sin(beta) - Dragbeta*np.cos(beta)
    fay_beta = sideSlipbeta*np.cos(beta) + Dragbeta*np.sin(beta)
    faz_beta = -Liftbeta

    fax = fax_aoa + fax_beta
    fay = fay_aoa + fay_beta
    faz = faz_aoa + faz_beta

    fx = ftx + fax
    fy = fty + fay
    fz = ftz + faz
    dudt = (((r * v) - (q * w)) + 2*(q1*q3 - q0*q2) * GRAVITY 
            + fx / mass)
    
    dvdt = (((p * w) - (r * u)) + 2*(q2*q3 + q0*q1) * GRAVITY 
            + fy / mass )
    dwdt = (((q * u) - (p * v)) + (q0**2 - q1**2 - q2**2 + q3**2) * GRAVITY 
            + fz / mass )

    # aerodynamic moment
    ixdot = IXDOT
    iydot = IYDOT
    izdot = IZDOT

    maP = 0.0
    maQ = faz * (MAC-masscenter)
    maR = -fay * (MAC-masscenter)


    mp = mtp + maP
    mq = mtq + maQ
    mr = mtr + maR

    dpdt = ((iy-iz) * q * r 
            + mp )/ix  
    # dpdt=0
    dqdt = ((iz-ix) * p * r 
            + mq)/iy
    # dqdt = 0

    drdt = ((ix-iy) * p * q
            + mr)/iz
    
    massdot = MDOT
    masscenterdot = MASSCENTERMOVINGRATE
    Vinertial = np.matmul(T.transpose(),velocityBody)
    Vx=Vinertial[0].item()
    Vy=Vinertial[1].item()
    Vz=Vinertial[2].item()
    
    
    dq0dt = (-p * q1 - q * q2 - r * q3) / 2
    dq1dt = (p * q0 + r * q2 - q * q3) / 2
    dq2dt = (q * q0 - r * q1 + p * q3) / 2
    dq3dt = (r * q0 + q * q1 - p * q2) / 2


    return [dpdt, dqdt, drdt, dudt, dvdt, dwdt, dq0dt, dq1dt, dq2dt, dq3dt,Vx,Vy,Vz,massdot,masscenterdot,ixdot,iydot,izdot]


# y = [p, q, r, u, v, w, q0, q1, q2, q3,Xinertial,Yinertial,Zinertial,mass,masscenter,ix,iy,iz]
# k0 =  [0.0, 0.0, 0.0, 200, 0, 0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,1000.0, INITIAL_MASS, INITIALMASSCENTER, INITIX, INITIY, INITIZ]
y0 =  [0,0,0,200,0,0,0,0,0,1,0,0,1000,INITIAL_MASS,INITIALMASSCENTER,INITIX,INITIY,INITIZ]
# initial condition
t0 = 0
tfinal = 12
ts = np.linspace(t0, tfinal, 90000)

sol = solve_ivp(f, [t0, tfinal], y0, t_eval=ts,rtol=1e-13,atol=1e-14)

p, q, r, u, v, w, q0, q1, q2, q3,Xinertial,Yinertial,Zinertial,mass,masscenter,ix,iy,iz = sol.y           
phi = []
theta = []
psi = []
c_eta = []
c_zeta = []
for index in range(len(q0)):
   phi.append(np.arctan(2*(q2[index]*q3[index] + q0[index]*q1[index])/(q0[index]**2 - q1[index]**2 - q2[index]**2 + q3[index]**2)))
   theta.append(np.arcsin(-2 * (q1[index] * q3[index] - q0[index] * q2[index])))
   psi.append(np.arctan(2*(q1[index]*q2[index] + q0[index]*q3[index])/(q0[index]**2 + q1[index]**2 - q2[index]**2 - q3[index]**2)))
#    c_eta.append(KP2THETA * (r_theta_dot_sp[index]-q[index]) + KITHETA * theta_dot_err_sum[index])
#    c_zeta.append(np.arcsin((KP2PSI * (r_psi_dot_sp[index]-r[index]) + KIPSI * psi_dot_err_sum[index]) * IZ / 2 / THRUST / LENGTHARM / np.cos(c_eta[index])))




plt.figure(1,dpi=1200)
pylab.plot(sol.t, u, label='u(m/s)')
pylab.plot(sol.t, v, label='v(m/s)')
pylab.plot(sol.t, w, label='w(m/s)')

pylab.legend()
pylab.xlabel('t(s)')
pylab.savefig('RocketBodyV.png')


phi = np.rad2deg(phi)
theta = np.rad2deg(theta)
psi = np.rad2deg(psi)

plt.figure(2,dpi=1200)
# pylab.plot(sol.t, p, label='p')
pylab.plot(sol.t, phi, label='phi(deg)')
pylab.legend()
pylab.xlabel('t(s)')
pylab.savefig('Rocketeuler_phi.png')


plt.figure(3,dpi=1200)
# pylab.plot(sol.t, p, label='p')
pylab.plot(sol.t, theta, label='theta(deg)')
pylab.legend()
pylab.xlabel('t(s)')
pylab.savefig('Rocketeuler_theta.png')

plt.figure(4,dpi=1200)
pylab.plot(sol.t, psi, label='psi(deg)')
pylab.legend()
pylab.xlabel('t(s)')
pylab.savefig('Rocketeuler_psi.png')

plt.figure(5,dpi=1200)
# pylab.plot(sol.t, p, label='p')
pylab.plot(sol.t, p, label='p(rad/s)')
pylab.legend()
pylab.xlabel('t(s)')
pylab.savefig('Rocketp.png')


plt.figure(6,dpi=1200)
# pylab.plot(sol.t, p, label='p')
pylab.plot( control, label='r(rad/s)')
pylab.legend()
pylab.xlabel('t(s)')
# pylab.savefig('Rocketc.png')



plt.figure(7,dpi=1200)
# pylab.plot(sol.t, p, label='p')
pylab.plot(sol.t, p, label='p(rad/s)')
pylab.legend()
pylab.xlabel('t(s)')
pylab.savefig('Rocketp.png')


plt.figure(8,dpi=1200)
pylab.plot(sol.t, Xinertial, label='X(m/s)')
pylab.plot(sol.t, Yinertial, label='Y(m/s)')
pylab.plot(sol.t, Zinertial, label='Z(m/s)')

pylab.legend()
pylab.xlabel('t(s)')
pylab.savefig('RocketPotition.png')