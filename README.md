# Non-Commutative weak measurement

The simulation source code for arXiv:.

## Simulation

- ComupleteReadout.py

Simulation of the complete readout case Fig.(d).

- PartialReadout.py

Simulation of the partial readout case Fig.(c) and Fig.(d).

## Data folder

- EV_Crit_0: Complete readout case: $\braket{\psi|Z_iZ_j|\psi}^2$ (entanglement and SWSSB).

    "ZZX_EV+6+0.1+0.01.npy": 
    - ZZX_EV: $\braket{\psi|Z_iZ_j|\psi}^2$
    - 6: the number of qubits
    - 0.1: the strength of the ZZ weak measurement ($t_1=0.1$).
    - 0.01: the strength of the X weak measurement ($t_2=0.01$).

- EV_Crit_1: Partial readout case: $\text{Tr}(\rho Z_i Z_j)^2$ (entanglement) and $F(\rho,Z_iZ_j\rho Z_iZ_j)$ (SWSSB)

    "ZZX_EV_NoXM+6+0.1+0.01": discard the outcome of X measurements.
    - ZZX_EV_NoXM: $\text{Tr}(\rho Z_i Z_j)^2$
    - 6: the number of qubits
    - 0.1: the strength of the ZZ weak measurement ($t_1=0.1$).
    - 0.01: the strength of the X weak measurement ($t_2=0.01$).

    "ZZX_FL_NoXM+6+0.1+0.01": discard the outcome of X measurements.
    - ZZX_FL_NoXM: $F(\rho,Z_iZ_j\rho Z_iZ_j)$
    - 6: the number of qubits
    - 0.1: the strength of the ZZ weak measurement ($t_1=0.1$).
    - 0.01: the strength of the X weak measurement ($t_2=0.01$).

## Draw

- DrawFig2.ipynb
