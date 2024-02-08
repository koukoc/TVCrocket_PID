import numpy as np
import pylab
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.spatial.transform import Rotation

INITIAL_MASS = 26.409  #kg
GROUND_GRAVITY = 9.81  #m/s^2
THRUST = 100 * 9.81 #N
HEIGHT = 100 #km
EARTH_RADIUS = 6371 #km
GRAVITY = GROUND_GRAVITY

IX = 0.183
IY = 4.479
IZ = 4.480
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


# y = [p  q  r  u  v  w  q0 q1 q2 q3 t1psi t1theta t2psi  t2theta t3psi  t3theta t4psi  t4theta theta_dot_err_sum psi_dot_err_sum]
def f(t,y):
    thrust = 100 * 9.81


######################################################system#############################################################################

    
    p, q, r, u, v, w, q0, q1, q2, q3, t1psi, t1theta, theta_dot_err_sum, psi_dot_err_sum, phi_dot_err_sum = y

    T = np.array([[q0**2 + q1**2 - q2**2 - q3**2, 2*(q1*q2 + q0*q3), 2*(q1*q3 - q0*q2)]
                  ,[2*(q1*q2 - q0*q3), q0**2 - q1**2 + q2**2 - q3**2, 2*(q2*q3 + q0*q1)]
                  ,[2*(q1*q3 + q0*q2), 2*(q2*q3 - q0*q1), q0**2 - q1**2 - q2**2 + q3**2]])
    
    
    dpdt = 0 
    # dpdt=0
    dqdt = ((IZ-IX) * p * r 
            - LENGTHARM * thrust * np.sin(t1theta) * np.cos(t1psi))/IY 
    # dqdt = 0

    
    drdt = ((IX-IY) * p * q
            - LENGTHARM * thrust * np.cos(t1theta) * np.sin(t1psi))/IZ


    dudt = (((r * v) - (q * w)) + 2*(q1*q3 - q0*q2) * GRAVITY 
            + thrust * np.cos(t1theta)*np.cos(t1psi) / INITIAL_MASS)
    
    dvdt = (((p * w) - (r * u)) + 2*(q2*q3 + q0*q1) * GRAVITY 
            + thrust * np.cos(t1theta)*np.sin(t1psi) / INITIAL_MASS )
    dwdt = (((q * u) - (p * v)) + (q0**2 - q1**2 - q2**2 + q3**2) * GRAVITY 
            - thrust * np.sin(t1theta)*np.cos(t1psi) / INITIAL_MASS )

    dq0dt = (-p * q1 - q * q2 - r * q3) / 2
    dq1dt = (p * q0 + r * q2 - q * q3) / 2
    dq2dt = (q * q0 - r * q1 + p * q3) / 2
    dq3dt = (r * q0 + q * q1 - p * q2) / 2

##########################################       controller    ###############################################################################



    Q0_SP = np.cos(PSI_SP/2)*np.cos(THETA_SP/2)
    Q1_SP = -np.sin(PSI_SP/2)*np.sin(THETA_SP/2)
    Q2_SP = np.cos(PSI_SP/2)*np.sin(THETA_SP/2)
    Q3_SP = np.sin(PSI_SP/2)*np.cos(THETA_SP/2)

    Q_diff=quaternion_mutiply(Q0_SP,Q1_SP,Q2_SP,Q3_SP,q0,-q1,-q2,-q3)
    # Q0_diff = Q0_SP * q0 + Q1_SP * q1 + Q2_SP * q2 + Q3_SP * q3
    # Q1_diff = Q0_SP * q1 - Q1_SP * q0 - Q2_SP * q3 + Q3_SP * q2
    # Q2_diff = Q0_SP * q2 - Q2_SP * q0 - Q3_SP * q1 + Q1_SP * q3
    # Q3_diff = Q0_SP * q3 - Q3_SP * q0 - Q1_SP * q2 + Q2_SP * q1


    # Q0_diff = Q0_SP * q0 + Q1_SP * q1 + Q2_SP * q2 + Q3_SP * q3
    # Q1_diff = Q0_SP * q1 - Q1_SP * q0 - Q2_SP * q3 + Q3_SP * q2
    # Q2_diff = Q0_SP * q2 - Q2_SP * q0 - Q3_SP * q1 + Q1_SP * q3
    # Q3_diff = Q0_SP * q3 - Q3_SP * q0 - Q1_SP * q2 + Q2_SP * q1
#     pitch controller
    theta = np.arcsin( -2 * (q1 * q3 - q0 * q2))
    # theta_err = THETA_SP - theta
    # theta_err = np.arcsin( -2 * (Q1_diff * Q3_diff - Q0_diff * Q2_diff))
    theta_err = -np.pi / 2 + 2 * np.arctan2( np.sqrt((1 + 2 * (Q_diff[0] * Q_diff[2] - Q_diff[1] * Q_diff[3]))),np.sqrt(1 - 2 * (Q_diff[0] * Q_diff[2] - Q_diff[1] * Q_diff[3])))
    theta_dot_sp = KP1THETA * theta_err
    theta_dot_err = theta_dot_sp - q
    c_Mtheta = (KP2THETA * theta_dot_err + KITHETA * theta_dot_err_sum)

#     yaw controller
    psi = np.arctan(2*(q1*q2 + q0*q3)/(q0**2 + q1**2 - q2**2 - q3**2))
    # psi_err = PSI_SP - psi
    psi_err = np.arctan(2*(Q_diff[1]*Q_diff[2] + Q_diff[0]*Q_diff[3])/(Q_diff[0]**2 + Q_diff[1]**2 - Q_diff[2]**2 - Q_diff[3]**2))
    psi_dot_sp = KP1PSI * psi_err
    psi_dot_err = psi_dot_sp - r
    c_Mpsi = (KP2PSI * psi_dot_err + KIPSI * psi_dot_err_sum)

#     roll controller
    
    phi_dot_err = P_SP - p
    
    control.append(psi_err)
    c_Mphi = (KP2PHI * phi_dot_err + KIPHI * phi_dot_err_sum)*IX
        
#########################################       Mixer       ################################################################################
#     roll Mixer
    # 要產生的roll moment
    roll_moment = c_Mphi
    c_t2theta =  - np.arcsin(roll_moment/WIDTHARM/thrust/2/IX)/2
    c_t3theta =  + np.arcsin(roll_moment/WIDTHARM/thrust/2/IX)/2
    c_t2psi =  + np.arcsin(roll_moment/WIDTHARM/thrust/2/IX)/2
    c_t3psi =  - np.arcsin(roll_moment/WIDTHARM/thrust/2/IX)/2



#     pitch Mixer
    c_theta = c_Mtheta
    
 

#     yaw Mixer
    # MZ = c_Mpsi *IZ
    # A = - LENGTHARM * thrust * np.cos(c_Mtheta) 
    # B = WIDTHARM * thrust * np.cos(c_Mtheta)
    # LENGTHARM * thrust * np.cos(t1theta) * np.sin(t1psi)
    c_psi = (c_Mpsi)
    # c_t1psi = np.arctan2(B, A) + np.arctan2(-MZ, np.sqrt(A**2 + B**2 - MZ**2))
    # print(np.sqrt(A**2 + B**2 - MZ**2))

  
    
    c_t1theta = c_theta



    c_t1psi = c_psi



###########################################     Mixer 2 Output      ####################################################################################


    dt1theta = -10*t1theta + 10 * c_t1theta
    dt1psi = -10*t1psi + 10 * c_t1psi

    control.append(c_t1psi)


    return [dpdt, dqdt, drdt, dudt, dvdt, dwdt, dq0dt, dq1dt, dq2dt, dq3dt, dt1psi, dt1theta, theta_dot_err, psi_dot_err, phi_dot_err]



# y = [p  q  r  u  v  w  q0 q1 q2 q3 t1psi t1theta theta_dot_err_sum psi_dot_err_sum phi_dot_err_sum]
y0 =  [0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0, 1.0, 0,    0,        0,           0, 0]      # initial condition
t0 = 0
tfinal = 600
ts = np.linspace(t0, tfinal, 90000)

sol = solve_ivp(f, [t0, tfinal], y0, t_eval=ts)

p, q, r, u, v, w, q0, q1, q2, q3, t1psi, t1theta, theta_dot_err_sum, psi_dot_err_sum,phi_dot_err = sol.y           
phi = []
theta = []
psi = []
c_t1theta = []
c_t1psi = []
r_theta_dot_sp = []
r_psi_dot_sp = []
for index in range(len(q0)):
   
   
   
   phi.append(np.arctan(2*(q2[index]*q3[index] + q0[index]*q1[index])/(q0[index]**2 - q1[index]**2 - q2[index]**2 + q3[index]**2)))
   theta.append(np.arcsin(-2 * (q1[index] * q3[index] - q0[index] * q2[index])))
   psi.append(np.arctan(2*(q1[index]*q2[index] + q0[index]*q3[index])/(q0[index]**2 + q1[index]**2 - q2[index]**2 - q3[index]**2)))
  
   r_psi_dot_sp.append(KP1PSI * (PSI_SP - psi[index]))
   r_theta_dot_sp.append(KP1THETA * (THETA_SP - theta[index]))
#    c_t1theta.append(KP2THETA * (r_theta_dot_sp[index]-q[index]) + KITHETA * theta_dot_err_sum[index])
#    c_t1psi.append(np.arcsin((KP2PSI * (r_psi_dot_sp[index]-r[index]) + KIPSI * psi_dot_err_sum[index]) * IZ / 2 / THRUST / LENGTHARM / np.cos(c_t1theta[index])))




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

plt.figure(11,dpi=1200)
# pylab.plot(sol.t, p, label='p')
pylab.plot(sol.t, phi, label='phi(deg)')
pylab.legend()
pylab.xlabel('t(s)')
pylab.savefig('Rocketeuler_phi.png')


plt.figure(8,dpi=1200)
# pylab.plot(sol.t, p, label='p')
pylab.plot(sol.t, theta, label='theta(deg)')
pylab.legend()
pylab.xlabel('t(s)')
pylab.savefig('Rocketeuler_theta.png')

plt.figure(2,dpi=1200)
pylab.plot(sol.t, psi, label='psi(deg)')
pylab.legend()
pylab.xlabel('t(s)')
pylab.savefig('Rocketeuler_psi.png')

plt.figure(9,dpi=1200)
# pylab.plot(sol.t, p, label='p')
pylab.plot(sol.t, p, label='p(rad/s)')
pylab.legend()
pylab.xlabel('t(s)')
pylab.savefig('Rocketp.png')


plt.figure(3,dpi=1200)
# pylab.plot(sol.t, p, label='p')
pylab.plot(sol.t, q, label='q(rad/s)')
pylab.plot(sol.t, r_theta_dot_sp, label='q_r(rad/s)')
pylab.legend()
pylab.xlabel('t(s)')
pylab.savefig('Rocketq.png')

plt.figure(6,dpi=1200)
# pylab.plot(sol.t, p, label='p')
pylab.plot(sol.t, r, label='r(rad/s)')
pylab.plot(sol.t, r_psi_dot_sp, label='r_r(rad/s)')
pylab.legend()
pylab.xlabel('t(s)')
pylab.savefig('Rocketr.png')



plt.figure(4,dpi=1200)
pylab.plot(sol.t, t1theta, label='theta_u(t)(rad)')
# pylab.plot(sol.t,c_t1theta, label='theta_r(t)(rad)')
pylab.legend()
pylab.xlabel('t(s)')
pylab.savefig('RocketTVC_theta.png')

plt.figure(5,dpi=1200)
pylab.plot(sol.t, t1psi, label='psi_u(t)(rad)')
# pylab.plot(sol.t,c_t1psi, label='psi_r(t)(rad)')
pylab.legend()
pylab.xlabel('t(s)')
pylab.savefig('RocketTVC_psi.png')

plt.figure(7,dpi=1200)
# pylab.plot(sol.t, p, label='p')
pylab.plot( control, label='r(rad/s)')
pylab.legend()
pylab.xlabel('t(s)')
# pylab.savefig('Rocketc.png')



plt.figure(10,dpi=1200)
# pylab.plot(sol.t, p, label='p')
pylab.plot(sol.t, p, label='p(rad/s)')
pylab.legend()
pylab.xlabel('t(s)')
pylab.savefig('Rocketp.png')
