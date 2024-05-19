import numpy as np
class AeroDynamic:
    # Example coefficients for moments (these should be determined experimentally or through detailed simulations)
    Cm_alpha = -30  # Pitching moment coefficient per unit alpha
    Cn_beta = 30   # Yawing moment coefficient per unit beta
    Cl_p = -0.012994     # Rolling moment coefficient per unit roll rate
    Cl_roll= -0.00  # Rolling moment coefficient per unit roll

    # rocket parameter
    L_ref = 0.205/2 # diameter of rocket
    Sref = L_ref**2*np.pi

    def __init__(self,alpha,beta,p,AirSpeed,rollAngle,SpeedofSound):
        self.CL = 11*alpha
        self.CD= 0.4 + 10.2 * np.arccos(AirSpeed[0].item()/np.linalg.norm(AirSpeed)) # angle between wind and rocket

        self.CY = -11 * beta
        self.CM_pitch = self.Cm_alpha * alpha 
        self.CM_yaw = self.Cn_beta * beta 
        self.CM_roll = self.Cl_p * p + self.Cl_roll * rollAngle

        pass