from OOM_model_simoca import model
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# M = model()
# #M.dynamicQ=True
# M.solve()
# print(M.solve()[0])
# print(M.solve()[1])

class optimization:
    #We call the model to avoid the need to repeat all the equations etc.
    M = model()
    t = M.solve()[0]
    C0 = np.zeros((len(M.timepoints), 9))
    C0[0] = [M.X0, M.Glu0, M.Xyl0, M.EtOH0, M.Fur0, M.HAc0, M.HMF0, M.FA0, M.V0]
    # Vector containing all the possible values of Q to detect the optimal value from 0 to Qmax and n will define the amount of values to evaluate
    Q = np.linspace(0, M.Q_max, 5)
    Xylrates = np.zeros((len(M.timepoints), len(Q)))
    Q_opt_Xylrates = np.zeros(len(M.timepoints))
    for j in range(len(t)):
        for i in range(len(Q)):
            M.Q = Q [i]
            M.solve()
            C = M.solve()[1]
            Xylrates[j, i] = M.Y_EtOHXyl * (
                        M.nuMaxXyl * C[-1, 0] * (C[-1, 2] / (M.Ks_Xyl + C[-1, 2] + ((C[-1, 2] ** 2) / M.Ki_Xyl)) *
                                               (1 - (C[-1, 3] / M.Ki_EtOHmaxXyl) ** M.gammaX) *
                                               (1 / (1 + (C[-1, 4] / M.Ki_FurXyl))) *
                                               (1 / (1 + (C[-1, 5] / M.Ki_HAcXyl))) *
                                               (1 / (1 + (C[-1, 6] / M.Ki_HMFXyl))) *
                                               (1 / (1 + (C[-1, 7] / M.Ki_FAXyl))) *
                                               (1 / (1 + (C[-1, 1] / M.Ki_GluXyl))))) + \
                             M.Y_EtOHGlu * (M.nuMaxGlu * C[-1, 0] * (
                        C[-1, 1] / (M.Ks_Glu + C[-1, 1] + ((C[-1, 2] ** 1) / M.Ki_Glu)) *
                        (1 - (C[-1, 3] / M.Ki_EtOHmaxGlu) ** M.gammaG) *
                        (1 / (1 + (C[-1, 4] / M.Ki_FurGlu))) *
                        (1 / (1 + (C[-1, 5] / M.Ki_HAcGlu))) *
                        (1 / (1 + (C[-1, 6] / M.Ki_HMFGlu))) *
                        (1 / (1 + (C[-1, 7] / M.Ki_FAGlu)))))
            Xylrates_data = pd.DataFrame(Xylrates, columns=Q, index=M.timepoints)

            # Here comes the main point of the optimization. the initial concentrations of the next timepoint will be the values achieved with the simulation of the optimal Q
            if Xylrates[j, i] == Xylrates[j].max():
                C0[j + 1] = C[-1]

            Q_opt_Xylrates[j] = Q[np.argmax(Xylrates, axis=1)[j]]

            # Here a restriction of volume is included

            Q_opt_Xylrates_data = pd.DataFrame(Q_opt_Xylrates, index=M.timepoints)


    # Xylrates = np.zeros((len(M.timepoints), len(Q)))
    # Q_opt_Xylrates = np.zeros(len(M.timepoints))


    # M.dynamicQ=True

O = optimization()
print(O.Q_opt_Xylrates_data)


