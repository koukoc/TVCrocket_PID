import numpy as np
import pandas as pd
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


MASSCENTERMOVINGRATE = (INITIALMASSCENTER - FINALMASSCENTER)/BURNINGOUTTIME #m/s
WIND = np.array([[0],[0],[0]]) #m/s


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
    if t < 15:
        thrust = 100 * 9.81
        ixdot = IXDOT
        iydot = IYDOT
        izdot = IZDOT
        massdot = MDOT
        masscenterdot = MASSCENTERMOVINGRATE
    else:
        thrust=0
        ixdot = 0
        iydot = 0
        izdot = 0
        massdot = 0
        masscenterdot = 0
    
######################################################system#############################################################################

    p, q, r, u, v, w, q0, q1, q2, q3,Xinertial,Yinertial,Zinertial,mass,masscenter,ix,iy,iz = y

    # qmag = np.sqrt(q0**2+q1**2+q2**2+q3**2)
    # q0 = q0/qmag
    # q1 = q1/qmag
    # q2 = q2/qmag
    # q3 = q3/qmag
    
    velocityBody = np.array([[u],[v],[w]])
    T = np.array([[q0**2 + q1**2 - q2**2 - q3**2, 2*(q1*q2 + q0*q3), 2*(q1*q3 - q0*q2)]
                  ,[2*(q1*q2 - q0*q3), q0**2 - q1**2 + q2**2 - q3**2, 2*(q2*q3 + q0*q1)]
                  ,[2*(q1*q3 + q0*q2), 2*(q2*q3 - q0*q1), q0**2 - q1**2 - q2**2 + q3**2]])
    eta = np.deg2rad(0)
    zeta =np.deg2rad(3)
    mtp = 0
    mtq = - (TVCPOINT - masscenter) * thrust * np.sin(eta)
    mtr = - (TVCPOINT - masscenter) * thrust * np.cos(eta) * np.sin(zeta)
    

    # TVC thrust
    ftx = thrust * np.cos(eta)*np.cos(zeta)
    fty = thrust * np.cos(eta)*np.sin(zeta)
    ftz = - thrust * np.sin(eta)

    # aerodynamiv force
    air = AirProperties(-Zinertial)
    WIND_body = T @ WIND
    AirSpeed = velocityBody+WIND_body
    velocityMag = np.linalg.norm(velocityBody+WIND_body)
    AoA = np.arctan2((w+WIND_body[2].item()),(u+WIND_body[0].item()))
    beta= np.arcsin((v+WIND_body[1].item())/(velocityMag))
    # print('u',u,'v',v,'w',w,'AoA',AoA)
    phi=np.arctan(2*(q2*q3 + q0*q1)/(q0**2 - q1**2 - q2**2 + q3**2))


    aerodynamic = AeroDynamic(AoA,beta,p,AirSpeed,phi,air.SpeedofSound)

    dynamicPressure = 0.5*air.density*velocityMag**2

    # aerodynamic force from angle of attack
    Lift = dynamicPressure*aerodynamic.CL*aerodynamic.Sref
    Drag = dynamicPressure*aerodynamic.CD*aerodynamic.Sref
    sideSlip = dynamicPressure*aerodynamic.CY*aerodynamic.Sref

	# wind frame to bodyframe
    wind_to_body = np.array([
        [np.cos(beta)*np.cos(AoA), -np.sin(beta), np.cos(beta)*np.sin(AoA)],
        [np.sin(beta)*np.cos(AoA), np.cos(beta), np.sin(beta)*np.sin(AoA)],
        [-np.sin(AoA), 0, np.cos(AoA)]
    ])
    # aerodynamic force in bodyframe
    F_aero = wind_to_body @ np.array([-Drag,sideSlip,-Lift])   

    fx = ftx + F_aero[0].item()
    fy = fty + F_aero[1].item()
    fz = ftz + F_aero[2].item()
    
    dudt = (((r * v) - (q * w)) + 2*(q1*q3 - q0*q2) * GRAVITY 
            + fx / mass)
    
    dvdt = (((p * w) - (r * u)) + 2*(q2*q3 + q0*q1) * GRAVITY 
            + fy / mass )
    dwdt = (((q * u) - (p * v)) + (q0**2 - q1**2 - q2**2 + q3**2) * GRAVITY 
            + fz / mass )

    # aerodynamic moment
    M_roll = dynamicPressure*aerodynamic.CM_roll*aerodynamic.Sref*aerodynamic.L_ref
    M_pitch = dynamicPressure*aerodynamic.CM_pitch*aerodynamic.Sref*aerodynamic.L_ref
    M_yaw = dynamicPressure*aerodynamic.CM_yaw*aerodynamic.Sref*aerodynamic.L_ref
    
    # M_aero=np.array([M_roll, M_pitch, M_yaw])

    # print('beta:',beta,'mar:',maR)
    # print('p',p)
    mp = M_roll
    mq = mtq + M_pitch
    mr = mtr + M_yaw

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

    if Zinertial >= 0:
        Vx=0
        Vy=0
        Vz=0

    return [dpdt, dqdt, drdt, dudt, dvdt, dwdt, dq0dt, dq1dt, dq2dt, dq3dt,Vx,Vy,Vz,massdot,masscenterdot,ixdot,iydot,izdot]


# y = [p, q, r, u, v, w, q0, q1, q2, q3,Xinertial,Yinertial,Zinertial,mass,masscenter,ix,iy,iz]
# k0 =  [0.0, 0.0, 0.0, 200, 0, 0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,1000.0, INITIAL_MASS, INITIALMASSCENTER, INITIX, INITIY, INITIZ]
y0 =  [0.0,0.0,0.0
       ,200,1,0,
       0.7372773,0, 0.6755902, 0,
       0,0,-1000,
       INITIAL_MASS,
       INITIALMASSCENTER,
       INITIX,INITIY,INITIZ]
# initial condition
t0 = 0
tfinal = 90
ts = np.linspace(t0, tfinal, 9000)

sol = solve_ivp(f, [t0, tfinal], y0, t_eval=ts,rtol=1e-13,atol=1e-14,method='RK45')

p, q, r, u, v, w, q0, q1, q2, q3,Xinertial,Yinertial,Zinertial,mass,masscenter,ix,iy,iz = sol.y           
phi = []
theta = []
psi = []
AoA = []
Beta = []
for index in range(len(q0)):
   phi.append(np.arctan2(2*(q2[index]*q3[index] - q0[index]*q1[index]),(q0[index]**2 - q1[index]**2 - q2[index]**2 + q3[index]**2)))
   theta.append(np.arcsin(-2 * (q1[index] * q3[index] + q0[index] * q2[index])))
   psi.append(np.arctan(2*(q1[index]*q2[index] - q0[index]*q3[index])/(q0[index]**2 + q1[index]**2 - q2[index]**2 - q3[index]**2)))
   AoA.append(np.arctan(w[index]/u[index]))
   Beta.append(np.arcsin(v[index]/np.sqrt(u[index]**2+v[index]**2+w[index]**2)))
print('save data')


episodedata = pd.DataFrame({'time':sol.t,
                            'q0':q0,'q1':q1,'q2':q2,'q3':q3
                            ,'X':Xinertial,'Y':Yinertial,'Z':Zinertial}) 
episodedata.to_excel('rocketSim.xlsx', sheet_name='sheet1', index=False)


print('phi',phi[len(q0)-1],'theta',theta[len(q0)-1],'psi',psi[len(q0)-1])
print('saving picture')

plt.figure(1,dpi=400)
pylab.plot(sol.t, u, label='u(m/s)')
pylab.plot(sol.t, v, label='v(m/s)')
pylab.plot(sol.t, w, label='w(m/s)')

pylab.legend()
pylab.xlabel('t(s)')
pylab.savefig('RocketBodyV.png')


phi = np.rad2deg(phi)
theta = np.rad2deg(theta)
psi = np.rad2deg(psi)

plt.figure(2,dpi=400)
# pylab.plot(sol.t, p, label='p')
pylab.plot(sol.t, phi, label='phi(deg)')
pylab.legend()
pylab.xlabel('t(s)')
pylab.savefig('Rocketeuler_phi.png')


plt.figure(3,dpi=400)
# pylab.plot(sol.t, p, label='p')
pylab.plot(sol.t, theta, label='theta(deg)')
pylab.legend()
pylab.xlabel('t(s)')
pylab.savefig('Rocketeuler_theta.png')

plt.figure(4,dpi=400)
pylab.plot(sol.t, psi, label='psi(deg)')
pylab.legend()
pylab.xlabel('t(s)')
pylab.savefig('Rocketeuler_psi.png')

# plt.figure(5,dpi=1200)
# # pylab.plot(sol.t, p, label='p')
# pylab.plot(sol.t, p, label='p(rad/s)')
# pylab.legend()
# pylab.xlabel('t(s)')
# pylab.savefig('Rocketp.png')


# plt.figure(6,dpi=1200)
# # pylab.plot(sol.t, p, label='p')
# pylab.plot( control, label='r(rad/s)')
# pylab.legend()
# pylab.xlabel('t(s)')
# # pylab.savefig('Rocketc.png')



# plt.figure(7,dpi=1200)
# # pylab.plot(sol.t, p, label='p')
# pylab.plot(sol.t, p, label='p(rad/s)')
# pylab.legend()
# pylab.xlabel('t(s)')
# pylab.savefig('Rocketp.png')



# pylab.legend()
# pylab.xlabel('t(s)')
# pylab.savefig('RocketPotition.png')

plt.figure(9,dpi=400)
pylab.plot(sol.t, AoA, label='AoA(m/s)')
pylab.plot(sol.t, Beta, label='Beta(m/s)')

pylab.legend()
pylab.xlabel('t(s)')
pylab.savefig('RocketFlightPathAngle.png')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(Xinertial, Yinertial, -Zinertial)
ax.set_xlabel('X Position (m)')
ax.set_ylabel('Y Position (m)')
ax.set_zlabel('Z Position (m)')
ax.set_title('Rocket Trajectory')
ax.set_xlim([-2000, 3000])
ax.set_ylim([-2000, 3000])
pylab.savefig('Rocketpath.png')
plt.show()

