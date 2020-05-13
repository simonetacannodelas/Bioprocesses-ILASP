from scipy.integrate import odeint
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class model:
    def __init__(self, dynamicQ = False):
        df = pd.read_excel("Kinetics_overview.xlsx")
        self.nuMaxGlu = float(df['Value2'][df.index[df['Parameters'] == "nuMaxGlu"]])  # s-1
        self.nuMaxXyl = float(df['Value2'][df.index[df['Parameters'] == "nuMaxXyl"]])  # s-1
        self.Ks_Glu = float(df['Value2'][df.index[df['Parameters'] == "Ks_Glu"]])  # kg Glu m-3
        self.Ks_Xyl = float(df['Value2'][df.index[df['Parameters'] == "Ks_Xyl"]])  # kg Xyl m-3
        self.Ki_Glu = float(df['Value2'][df.index[df['Parameters'] == "Ki_Glu"]])  # kg Glu m-3
        self.Ki_Xyl = float(df['Value2'][df.index[df['Parameters'] == "Ki_Xyl"]])  # kg Xyl m-3
        self.Ki_GluXyl = float(df['Value2'][df.index[df['Parameters'] == "Ki_GluXyl"]])  # kg Glu m-3
        self.Y_XGlu = float(df['Value2'][df.index[df['Parameters'] == "Y_XGlu"]])  # kg X/kg Glu
        self.Y_XXyl = float(df['Value2'][df.index[df['Parameters'] == "Y_XXyl"]])  # kg X/kg Xyl
        self.Ki_EtOHmaxGlu = float(df['Value2'][df.index[df['Parameters'] == "Ki_EtOHmaxGlu"]])  # kg Glu m-3
        self.Ki_EtOHmaxXyl = float(df['Value2'][df.index[df['Parameters'] == "Ki_EtOHmaxXyl"]])  # kg Xyl m-3
        self.Y_EtOHGlu = float(df['Value2'][df.index[df['Parameters'] == "Y_EtOHGlu"]])  # kg EtOH/kg Glu
        self.Y_EtOHXyl = float(df['Value2'][df.index[df['Parameters'] == "Y_EtOHXyl"]])  # kg EtOH/kg Xyl
        self.gammaG = float(df['Value2'][df.index[df['Parameters'] == "gammaG"]])  # no unit
        self.gammaX = float(df['Value2'][df.index[df['Parameters'] == "gammaX"]])  # no unit

        # Acetate-related parameters
        self.nuHAcMax = float(df['Value2'][df.index[df['Parameters'] == "nuHAcMax"]])  # s-1
        self.Ks_HAc = float(df['Value2'][df.index[df['Parameters'] == "Ks_HAc"]])  # kg HAc m-3
        self.Ki_HAcGlu = float(df['Value2'][df.index[df['Parameters'] == "Ki_HAcGlu"]])  # kg HAc m-3
        self.Ki_HAcXyl = float(df['Value2'][df.index[df['Parameters'] == "Ki_HAcXyl"]]) # kg HAc m-3
        self.Y_HAcHMF = float(df['Value2'][df.index[df['Parameters'] == "Y_HAcHMF"]])  # kg Ac/kg HMF

        # Furfural-related parameters
        self.nuFurMax = float(df['Value2'][df.index[df['Parameters'] == "nuFurMax"]])  # s-1
        self.Ks_Fur = float(df['Value2'][df.index[df['Parameters'] == "Ks_Fur"]])  # kg Furfural m-3
        self.Ki_FurGlu = float(df['Value2'][df.index[df['Parameters'] == "Ki_FurGlu"]])  # kg Furfural m-3
        self.Ki_FurXyl = float(df['Value2'][df.index[df['Parameters'] == "Ki_FurXyl"]])  # kg Furfural m-3
        self.Ki_FurHMF = float(df['Value2'][df.index[df['Parameters'] == "Ki_FurHMF"]])  # kg Furfural m-3
        self.Y_FurFA = float(df['Value2'][df.index[df['Parameters'] == "Y_FurFA"]])  # kg FA/kg Fur

        # Furfuryl alcohol-related parameters
        self.Ki_FAGlu = float(df['Value2'][df.index[df['Parameters'] == "Ki_FAGlu"]])  # kg FA m-3
        self.Ki_FAXyl = float(df['Value2'][df.index[df['Parameters'] == "Ki_FAXyl"]])  # kg FA m-3

        # HMF-related parameters
        self.nuHMFMax = float(df['Value2'][df.index[df['Parameters'] == "nuHMFMax"]])  # s-1
        self.Ks_HMF = float(df['Value2'][df.index[df['Parameters'] == "Ks_HMF"]])  # kg HMF m-3
        self.Ki_HMFGlu = float(df['Value2'][df.index[df['Parameters'] == "Ki_HMFGlu"]])  # kg HMF m-3
        self.Ki_HMFXyl = float(df['Value2'][df.index[df['Parameters'] == "Ki_HMFXyl"]])  # kg HMF m-3

        self.X_mass = 100  # g
        self.V0 = 100  # L
        self.X0 = self.X_mass / self.V0  # g/L
        self.Glu0 = 40  # g/L
        self.Xyl0 = 20  # g/L
        self.EtOH0 = 0  # g/L
        self.Fur0 = 1  # g/L
        self.HAc0 = 1  # g/L
        self.HMF0 = 0.5  # g/L
        self.FA0 = 0  # g/L
        self.Ferm_time = 50  # h
        self.n = 101  # steps for integration
        self.dynamicQ = dynamicQ
        self.Q_max = 150  # maximum operative flow
        self.n = 101  # steps for integration

        # Time vector to solve the ODE. It goes from 0 to the fermentation time and we are looking for one hour intervals to choose the optimum Q for every hour
        self.timepoints = np.linspace(0, self.Ferm_time, self.Ferm_time + 1)
        # This correction it is made to show the values with two decimals only to make it clearer
        self.timepoints = np.around(self.timepoints, decimals=2)

        self.t = np.linspace(0, self.Ferm_time, self.n)
        self.t = np.around(self.t, decimals=2)
        self.Q = 0
        self.Q_opt_Xylrates_data = None


    def rxn(self, C, t):

        # Initialization of the stoichiometric matrix

        s = np.zeros((5, 8))

        s[0, 0] = (self.Y_XGlu)
        s[1, 0] = (self.Y_XXyl)
        s[2, 0] = (0)
        s[3, 0] = (0)
        s[4, 0] = (0)

        s[0, 1] = (-1)
        s[1, 1] = (0)
        s[2, 1] = (0)
        s[3, 1] = (0)
        s[4, 1] = (0)

        s[0, 2] = (0)
        s[1, 2] = (-1)
        s[2, 2] = (0)
        s[3, 2] = (0)
        s[4, 2] = (0)

        s[0, 3] = (self.Y_EtOHGlu)
        s[1, 3] = (self.Y_EtOHXyl)
        s[2, 3] = (0)
        s[3, 3] = (0)
        s[4, 3] = (0)

        s[0, 4] = (0)
        s[1, 4] = (0)
        s[2, 4] = (-1)
        s[3, 4] = (0)
        s[4, 4] = (0)

        s[0, 5] = (0)
        s[1, 5] = (0)
        s[2, 5] = (0)
        s[3, 5] = (-1)
        s[4, 5] = (self.Y_HAcHMF)

        s[0, 6] = (0)
        s[1, 6] = (0)
        s[2, 6] = (0)
        s[3, 6] = (0)
        s[4, 6] = (-1)

        s[0, 7] = (0)
        s[1, 7] = (0)
        s[2, 7] = (self.Y_FurFA)
        s[3, 7] = (0)
        s[4, 7] = (0)

        # Initialization of the process vector used in the rxn function

        rho = np.zeros(5)

        # Initialization of the overall rates vector

        r = np.zeros(8)

        X, Glu, Xyl, EtOH, Fur, HAc, HMF, FA, V = C

#        Q = Flow(t)

        # Glucose uptake process
        rho[0] = self.nuMaxGlu * X * (Glu / (self.Ks_Glu + Glu + ((Glu ** 2) / self.Ki_Glu)) *
                                 (1 - (EtOH / self.Ki_EtOHmaxGlu) **self.gammaG) *
                                 (1 / (1 + (Fur / self.Ki_FurGlu))) *
                                 (1 / (1 + (HAc / self.Ki_HAcGlu))) *
                                 (1 / (1 + (HMF / self.Ki_HMFGlu))) *
                                 (1 / (1 + (FA / self.Ki_FAGlu))))
        # Xylose uptake process
        rho[1] = self.nuMaxXyl * X * (Xyl / (self.Ks_Xyl + Xyl + ((Xyl ** 2) / self.Ki_Xyl)) *
                                 (1 - (EtOH / self.Ki_EtOHmaxXyl) ** self.gammaX) *
                                 (1 / (1 + (Fur / self.Ki_FurXyl))) *
                                 (1 / (1 + (HAc / self.Ki_HAcXyl))) *
                                 (1 / (1 + (HMF / self.Ki_HMFXyl))) *
                                 (1 / (1 + (FA / self.Ki_FAXyl))) *
                                 (1 / (1 + (Glu / self.Ki_GluXyl))))
        # Fur uptake process
        rho[2] = self.nuFurMax * X * (Fur / (self.Ks_Fur + Fur))
        # HAc uptake process
        rho[3] = self.nuHAcMax * X * (HAc / (self.Ks_HAc + HAc))
        # HMF uptake process
        rho[4] = self.nuHMFMax * X * (HMF / (self.Ks_HMF + HMF)) * (1 / (1 + (Fur / self.Ki_FurGlu)))

        # Overall conversion rate (stoichiometric matrix * process vector)
        r[0] = s[0, 0] * rho[0] + s[1, 0] * rho[1] + s[2, 0] * rho[2] + s[3, 0] * \
               rho[3] + s[4, 0] * rho[4]
        r[1] = s[0, 1] * rho[0] + s[1, 1] * rho[1] + s[2, 1] * rho[2] + s[3, 1] * \
               rho[3] + s[4, 1] * rho[4]
        r[2] = s[0, 2] * rho[0] + s[1, 2] * rho[1] + s[2, 2] * rho[2] + s[3, 2] * \
               rho[3] + s[4, 2] * rho[4]
        r[3] = s[0, 3] * rho[0] + s[1, 3] * rho[1] + s[2, 3] * rho[2] + s[3, 3] * \
               rho[3] + s[4, 3] * rho[4]
        r[4] = s[0, 4] * rho[0] + s[1, 4] * rho[1] + s[2, 4] * rho[2] + s[3, 4] * \
               rho[3] + s[4, 4] * rho[4]
        r[5] = s[0, 5] * rho[0] + s[1, 5] * rho[1] + s[2, 5] * rho[2] + s[3, 5] * \
               rho[3] + s[4, 5] * rho[4]
        r[6] = s[0, 6] * rho[0] + s[1, 6] * rho[1] + s[2, 6] * rho[2] + s[3, 6] * \
               rho[3] + s[4, 6] * rho[4]
        r[7] = s[0, 7] * rho[0] + s[1, 7] * rho[1] + s[2, 7] * rho[2] + s[3, 7] * \
               rho[3] + s[4, 7] * rho[4]

        # if self.dynamicQ == True:
        #     if self.timepoints in self.t:
        #         self.timepoints = self.timepoints
        #     else:
        #         self.timepoints = self.t[(np.abs(self.t - self.timepoints)).argmin()]
        #         self.Q = self.Q_opt_Xylrates_data[0][round(self.timepoints, 2)]
        #     return self.Q # L/h
        # else:
        #     self.Q = 0#np.zeros(len(self.timepoints))
        dVdt = self.Q

        dXdt = r[0] - X * (self.Q / V)

        dGludt = r[1] + (self.Q / V) * (self.Glu0 - Glu)

        dXyldt = r[2] + (self.Q / V) * (self.Xyl0 - Xyl)

        dEtOHdt = r[3] - EtOH * (self.Q / V)

        dFurdt = r[4] + (self.Q / V) * (self.Fur0 - Fur)

        dHAcdt = r[5] + (self.Q / V) * (self.HAc0 - HAc)

        dHMFdt = r[6] + (self.Q / V) * (self.Fur0 - Fur)

        dFAdt = r[7] - FA * (self.Q / V)

        return [dXdt, dGludt, dXyldt, dEtOHdt, dFurdt, dHAcdt, dHMFdt, dFAdt,
                dVdt]

    def solve(self):
        t = np.linspace(0, self.Ferm_time, self.n)
        C0 = [self.X0, self.Glu0, self.Xyl0, self.EtOH0, self.Fur0, self.HAc0, self.HMF0, self.FA0, self.V0]
        C = odeint(self.rxn, C0, t)
        return t, C

M = model()
M.solve()
t= M.solve()[0]
C = M.solve()[1]

# plotting code
plt.plot(t, C[:, 0], 'b-', label="X")
plt.plot(t, C[:, 1], 'g-', label="Glu")
plt.plot(t, C[:, 2], 'r-', label="Xyl")
plt.plot(t, C[:, 3], 'c-', label="EtOH")
plt.plot(t, C[:, 4], 'm-', label="Fur")
plt.plot(t, C[:, 5], 'y-', label="HAc")
plt.plot(t, C[:, 6], 'k-', label="HMF")
plt.plot(t, C[:, 7], '0.75', label="FA")
plt.xlabel("Time (s)")
plt.ylabel("Conc. (g/L)")
plt.legend(loc="upper right")
plt.title('Batch simulation')
#plt.savefig('Sim_batch.png')
plt.show()



