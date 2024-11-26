import numpy as np
import opensees.openseespy as ops

# TODO: Temporary patch pending version 0.1.11
ops.__all__.append("invoke")


# :::
# Degrading Bouc-Wen with updated stiffness
# :::
class deg_bw_opensees:
    def __init__(self, params):
        self.model = ops.Model(ndm=1, ndf=1)
        self.model.uniaxialMaterial("BoucWenMG", 1, #*params)
            eta    = params[0],
            k0     = params[1],  # true elastic stiffness
            sy0    = params[2],  # true yield stress
            sig    = params[3],
            lam    = params[4],
            mup    = params[5],
            sigp   = params[6],
            rsmax  = params[7],
            n      = params[8],
            alpha  = params[9],
            alpha1 = params[10],
            alpha2 = params[11],
            betam1 = params[12]
        )
    
    def set_trial_state(self, tstrain):
        self.tstrain = tstrain
        self.model.invoke("UniaxialMaterial", 1, f"{{strain {tstrain}}}")

    def commit_state(self):
        self.model.invoke("UniaxialMaterial", 1, "commit")
        self.cstate = self.tstate

    def __getattr__(self, name):
        if name == "tstate":
            return {
                "stress":    self.model.invoke("UniaxialMaterial", 1, "stress"),
                "strain":    self.tstrain,
                "stiffness": self.model.invoke("UniaxialMaterial", 1, "tangent")
            }
        else:
            raise AttributeError("No attribute named " + name)


class deg_bw_material_mod():
    '''
    This is a degrading bouc wen material class
    Modified to include reloading stiffness

    '''
    def __init__(self, params):
        self.params = {
            "eta1": params[0],
            "eta2": 1.0 - params[0],
            "k0": params[1],   # true elastic stiffness
            "sy0": params[2],  # true yield stress
            "sig": params[3],
            "lam": params[4],
            "mup": params[5],
            "sigp": params[6],
            "rsmax": params[7],
            "n": params[8],
            "alpha": params[9],
            "alpha1": params[10],
            "alpha2": params[11],
            "betam1": params[12]
        }
        
        self.cstate = {
            "strain": 0.0, # true strain
            "stress": 0.0, # true stress
            "stiffness": params[1], # true tangent stiffness
            "z": 0.0,
            "h": 0.0,
            "stress_el": 0.0,
            "stress_st": 0.0,
            "stress_y": 1.0,
            "rs": 0.0,
            "rk": 1.0,
            "rkmin": 1.0,
            "emax_pos": 0.0,
            "emax_neg": 0.0,
            "emax": 0.0,
            "regime": 0,
            "thetamax_pos": [0.0, 0.0],
            "thetamax_neg": [0.0, 0.0],
            "momenmax_pos": [0.0, 0.0],
            "momenmax_neg": [0.0, 0.0],
            "k2": 1.0,
        }
        
        self.tstate = {
            "strain": 0.0, # true strain
            "stress": 0.0, # true stress
            "stiffness": params[1], # true tangent stiffness
            "z": 0.0,
            "h": 0.0,
            "stress_el": 0.0,
            "stress_st": 0.0,
            "stress_y": 1.0,
            "rs": 0.0,
            "rk": 1.0,
            "rkmin": 1.0,
            "emax_pos": 0.0,
            "emax_neg": 0.0,
            "emax": 0.0,
            "regime": 0,
            "thetamax_pos": [0.0, 0.0],
            "thetamax_neg": [0.0, 0.0],
            "momenmax_pos": [0.0, 0.0],
            "momenmax_neg": [0.0, 0.0],
            "k2": 1.0,
        }
        self.conv_flag = True
        
    def __str__(self):
        print('Bouc_Wen material model')
        
        for vals in self.params:
            print(vals, ':=', self.params[vals])
            
        return ':::'
        
    def set_trial_state(self, tstrain):
        # (1) Model Parameters
        # Get elastic model parameters
        k_true = self.params["k0"]    # true elastic stiffness with units
        sy_true = self.params["sy0"]  # true yield stress with units
        ey_true = sy_true / k_true    # true yield strain
        
        k0 = 1.0   # normalized stiffness
        sy0 = 1.0  # normalized yield stress
        
        # Shape parameters
        eta1 = self.params["eta1"]
        eta2 = self.params["eta2"]
        n = self.params["n"]
        alpha = self.params["alpha"]
        
        # Pinching parameters
        sig = self.params["sig"]
        lam = self.params["lam"]
        mup = self.params["mup"]
        sigp = self.params["sigp"]
        rsmax = self.params["rsmax"]
        
        # Degradation parameters
        alpha1 = self.params["alpha1"]
        alpha2 = self.params["alpha2"]
        betam1 = self.params["betam1"]
        
        # (2) Load current state variables (all these are normalized)
        stress_el = self.cstate["stress_el"]
        stress_st = self.cstate["stress_st"]
        stress_y = self.cstate["stress_y"]
        
        stress_sig = sig * stress_y
        stress_bar = lam * stress_y
        
        # Elastic post-yield stiffness
        k_el = alpha * k0
        
        # Degradation state parameters
        rk = self.cstate["rk"]
        rkmin = self.cstate["rkmin"]
        h = self.cstate["h"]
        emax_pos = self.cstate["emax_pos"]
        emax_neg = self.cstate["emax_neg"]
        emax = self.cstate["emax"]
        rs = self.cstate["rs"]

        # Variables for unloading-reloading
        regime = self.cstate["regime"]
        thetamax_pos = self.cstate["thetamax_pos"] 
        thetamax_neg = self.cstate["thetamax_neg"]
        momenmax_pos = self.cstate["momenmax_pos"]
        momenmax_neg = self.cstate["momenmax_neg"]
        k2 = self.cstate["k2"]
        
        # Normalize responses
        cstress = self.cstate["stress"] / sy_true    # normalized commited stress
        strain_i = self.cstate["strain"] / ey_true  # normalized commited strain
        strain_ip1 = tstrain / ey_true              # normalized trial strain
        
        dstrain = strain_ip1 - strain_i   # dstrain in normalized space

        # Define regime (unload-reload or not)
        if cstress * dstrain < 0 and regime == 0: # if changing from loading to unloading
            regime = 1
            # if unloading happens, store the moment and rotation at the previous point (pivot)
            if cstress > 0:
                thetamax_pos = [strain_i, strain_i]
                momenmax_pos = [cstress, cstress]
            else:
                thetamax_neg = [strain_i, strain_i]
                momenmax_neg = [cstress, cstress]
        
        if regime == 1 and cstress > 0 and dstrain > 0 and strain_ip1 > 0:
            # reloading, store the values
            regime = 2
            thetamax_pos[0] = strain_i
            momenmax_pos[0] = cstress
            if thetamax_pos[1] - thetamax_pos[0] != 0:
                k2 = (momenmax_pos[1] - momenmax_pos[0]) / (thetamax_pos[1] - thetamax_pos[0])
            else:
                k2 = k0
                
        elif regime == 1 and cstress < 0 and dstrain < 0 and strain_ip1 < 0:
            regime = 2
            thetamax_neg[0] = strain_i
            momenmax_neg[0] = cstress
            if thetamax_neg[1] - thetamax_neg[0] != 0:
                k2 = (momenmax_neg[1] - momenmax_neg[0]) / (thetamax_neg[1] - thetamax_neg[0])
            else:
                k2 = k0
            
        if regime == 0 or regime == 1:

            # Option 1
            startPoint = stress_st
            Tz = startPoint
            Tzold = startPoint
            Tznew = stress_st + 0.01

            # Option 2
            #startPoint = 0.01
            #Tz = startPoint
            #Tzold = startPoint
            #Tznew = 1.0
            
            count = 0
            maxiter = 500
            
            while np.abs(Tzold - Tznew) > 1.0e-6 and count < maxiter:
                
                # Function f(Tz):
                A = 1 - (eta1 * np.sign(Tz * dstrain) + eta2) * np.abs(Tz / stress_y) ** n
                kh = (rk - alpha) * k0 * A
                
                s = max([rs * (emax_pos - emax_neg), 0.0001])
    
                B = np.exp(- 0.5 * ((Tz - stress_bar * np.sign(dstrain)) / (stress_sig)) ** 2)
    
                # Check B value for stability
                if B < 1.0e-20:
                    B = 1.0e-20
    
                ks = min(((1 / np.sqrt(2 * np.pi)) * (s / stress_sig) * B) ** (-1), 1000)  
                kr = (kh * ks) / (kh + ks)
                f = kr * dstrain - Tz + stress_st
                
                # And it's derivative f_z(Tz):
                A_z = - (eta1 * np.sign(Tz * dstrain) + eta2) * n * np.abs(Tz / stress_sig) ** (n - 1) * np.sign(Tz / stress_sig)
                B_z = - (Tz - stress_bar * np.sign(dstrain)) * B / (stress_sig ** 2)
                ks_z = - (np.sqrt(2 * np.pi) * stress_sig / s) * (B ** -2) * B_z
                kh_z = (rk - alpha) * k0 * A_z
                kr_z = (kh_z * ks ** 2 + ks_z * kh ** 2) / (kh + ks) ** 2
    
                f_z = - 1.0 + kr_z * dstrain
                    
                Tznew = Tz - f / f_z
                    
                # Replace Tz with new value, keep the old one for convergence check
                Tzold = Tz
                Tz = Tznew
                
                count += 1
                
                if Tz == np.nan:
                    Tz = 0.001
                
                if count == maxiter and self.conv_flag == True:
                    # print('MaxIter reached in PH iterations')
                    self.conv_flag = False
                    
            # Upon convergence of Tz...
            
            # Compute tangent stiffness
            try:
                k = k_el + (kh * ks) / (kh + ks)  # This one is the total stiffness of the parallel/series spring system
            except:
                k = k0
        
        else:
            k = k2  # use reloading stiffness
            Tz = self.cstate["z"]

        tstress = k * dstrain + cstress   # This is still a normalized stress value
        stress_el = stress_el + k_el * dstrain   # elastic normalized stress
        stress_st = tstress - stress_el           # elastic plastic stress

        # How to go back to regime 0?
        if regime == 2:
            if tstress > 0 and dstrain > 0 and tstress > momenmax_pos[1]: # if on Q1 and moment exceed pivot, then go back to regime 0
                regime = 0
                tstress = 0.1 * tstress + 0.9 * momenmax_pos[1]
                
            elif tstress > 0 and dstrain < 0 and tstress < momenmax_pos[0]:
                regime = 1
                momenmax_pos[0] = cstress
                
            elif tstress < 0 and dstrain < 0 and tstress < momenmax_neg[1]:
                regime = 0
                tstress = 0.1 * tstress + 0.9 * momenmax_neg[1]
                    
            elif tstress < 0 and dstrain > 0 and tstress > momenmax_neg[0]:
                regime = 1
                momenmax_neg[0] = cstress
        
        if tstress * strain_ip1 < 0:
            regime = 0
        
        # Stiffness Degradation
        rk_trial = (np.abs(tstress) + alpha1 * stress_y) / (k0 * np.abs(strain_i) + alpha1 * stress_y)
        
        rkmin = min([rk_trial, rkmin])
        rk = rk_trial + (1 - alpha2) * (rkmin - rk_trial)
        
        # Strength Degradation
        dh = tstress * dstrain
        h = h + dh
        stress_y = sy0 / (1 + betam1 * h)
        emax_pos = max([emax_pos, strain_i])
        emax_neg = min([emax_neg, strain_i])

        # Msig = sig * My
        # Mbar = lam * My

        if abs(strain_i) - emax > 0:
            drs = (1 / (np.sqrt(2 * np.pi) * sigp)) * np.exp(-0.5 * ((emax - mup) / (sigp)) ** 2) * (np.abs(strain_i) - emax)
            rs = rs + drs * rsmax
            emax = np.abs(strain_i)

        # De-normalize the state parameters:
        k = k * k_true
        tstress = tstress * sy_true
        tstrain = strain_ip1 * ey_true
            
        self.tstate = {
            "strain": tstrain,
            "stress": tstress,
            "stiffness": k,
            "z": Tz,
            "h": h,
            "stress_el": stress_el,
            "stress_st": stress_st,
            "stress_y": stress_y,
            "rs": rs,
            "rk": rk,
            "rkmin": rkmin,
            "emax_pos": emax_pos,
            "emax_neg": emax_neg,
            "emax": emax,
            "regime": regime,
            "thetamax_pos": thetamax_pos,
            "thetamax_neg": thetamax_neg,
            "momenmax_pos": momenmax_pos,
            "momenmax_neg": momenmax_neg,
            "k2": k2
        }
        #for vals in self.tstate:
        #    print(vals, '=', self.tstate[vals])
        #print('--')
        
    def commit_state(self):
        self.cstate = self.tstate

    def reset_state(self):
        self.cstate = {
            "strain": 0.0, # true strain
            "stress": 0.0, # true stress
            "stiffness": self.params["k0"], # true tangent stiffness
            "z": 0.0,
            "h": 0.0,
            "stress_el": 0.0,
            "stress_st": 0.0,
            "stress_y": self.params["sy0"],
            "rs": 0.0,
            "rk": 1.0,
            "rkmin": 1.0,
            "emax_pos": 0.0,
            "emax_neg": 0.0,
            "emax": 0.0,
            "regime": 0,
            "thetamax_pos": [0.0, 0.0],
            "thetamax_neg": [0.0, 0.0],
            "momenmax_pos": [0.0, 0.0],
            "momenmax_neg": [0.0, 0.0],
            "k2": 1.0
        }
        self.tstate = {
            "strain": 0.0, # true strain
            "stress": 0.0, # true stress
            "stiffness": self.params["k0"], # true tangent stiffness
            "z": 0.0,
            "h": 0.0,
            "stress_el": 0.0,
            "stress_st": 0.0,
            "stress_y": self.params["sy0"],
            "rs": 0.0,
            "rk": 1.0,
            "rkmin": 1.0,
            "emax_pos": 0.0,
            "emax_neg": 0.0,
            "emax": 0.0,
            "regime": 0,
            "thetamax_pos": [0.0, 0.0],
            "thetamax_neg": [0.0, 0.0],
            "momenmax_pos": [0.0, 0.0],
            "momenmax_neg": [0.0, 0.0],
            "k2": 1.0
        }

def material_test(material, strains):
    '''
    Function that tests a material model, using a history of displacements
    '''
    stress_vec = []
    strain_vec = []
    
    for tstrain in strains:
        # run state determination block and extract the
        material.set_trial_state(tstrain)
        material.commit_state()
        stress_vec.append(material.cstate["stress"])
        strain_vec.append(material.cstate["strain"])

    return strain_vec, stress_vec
    

def createCyclicStrain(maxStrains=[0.01, 0.02], dStrain=0.001):
    """
    Function to create cyclic test data

    Parameters
    ----------
    maxStrains : TYPE, optional
        DESCRIPTION. The default is [0.01, 0.02].
    dStrain : TYPE, optional
        DESCRIPTION. The default is 0.001.

    Returns
    -------
    cycles : TYPE
        DESCRIPTION.
    t : TYPE
        DESCRIPTION.

    """
    cycles = []
    ncycles = len(maxStrains)
    
    n = 0
    for maxStrain in maxStrains:
        upWards1 = np.arange(0, maxStrain, dStrain)
        dnWards1 = np.arange(maxStrain, -maxStrain, -dStrain)
        upWards2 = np.arange(-maxStrain, 0, dStrain)
        cycle = np.concatenate([upWards1, dnWards1, upWards2])
        cycles.append(cycle)
        
    cycles = np.concatenate(cycles)
    t = np.arange(0, len(cycles))
    
    return cycles, t




if __name__ == "__main__":
    import matplotlib.pyplot as plt
    params = [
        1.0, # eta
        10,  # k0
         1, # sy0
         1,
         0.1,
        *[1.0]*(13-5)
    ]

    a = deg_bw_material_mod(params)
    b = deg_bw_opensees(params)

    cycles, t = createCyclicStrain()

    fig, ax = plt.subplots()

    ax.plot(*material_test(a, cycles))
    ax.plot(*material_test(b, cycles), "x")

    plt.show()