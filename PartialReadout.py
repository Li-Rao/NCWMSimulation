from qiskit.quantum_info import Statevector, Pauli
from qiskit.circuit.library import HamiltonianGate, HGate
from qiskit.quantum_info import SparsePauliOp
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathos.multiprocessing import ProcessingPool
from functools import partial
import multiprocessing
import time
from qiskit.quantum_info import DensityMatrix, Clifford, Operator, partial_trace
from qiskit.quantum_info import state_fidelity

# Auxiliary operator: H gate and the transformation matrix UO that maps Z eigenstates to Y eigenstates
HO=1/np.sqrt(2)*np.array([[1.0, 1.0],[1.0,-1.0]])
UO=1/np.sqrt(2)*np.array([[1.0, 0-1.0j],[1,0+1.0j]])

# - init(N): Initialization function, returns the density matrix of the product state $\ket{+}^{\otimes N}$ for $N+1$ qubits, where qubit 0 is always used as the auxiliary qubit for measurement.
# - MeasureZ_reset(dm): Measurement function for the ZZ operator, performs measurement of the auxiliary qubit in the Y basis, resets the auxiliary qubit to $\ket{+}$, and returns the resulting density matrix.
# - NoMeasureX_reset(dm): Used for partial readout, i.e., X measurement without reading out, traces out the auxiliary qubit and introduces a new $\ket{+}$ state to replace the original auxiliary qubit, returning the resulting density matrix.
def init(N):
    dm = DensityMatrix.from_int(0,2**(N+1))
    for i in range(N+1):
        dm= dm.evolve(HO,[i])
    return dm

def MeasureZ_reset(dm):
    dm = dm.evolve(UO,[0])
    meas,dm = dm.measure([0])
    dm=dm.reset([0])
    dm=dm.evolve(HO,[0])
    return dm

def NoMeasureX_reset(dm):
    dm0=DensityMatrix.from_label('+')
    dm=dm0.expand(partial_trace(dm,[0]))
    return dm

# - evol_ZZ(dm,U1): Perform one layer of ZZ measurement and readout, return the resulting density matrix
# - evol_X(dm,U2): Perform one layer of X measurement without readout, return the resulting density matrix

def evol_ZZ(dm,U1):
    N = dm.num_qubits-1
    for i in range(N-1):
        dm=dm.evolve(U1,[0,i+1,i+2])
        dm = MeasureZ_reset(dm)
    return dm

def evol_X(dm,U2):
    N = dm.num_qubits-1
    for i in range(1,N+1):
        dm=dm.evolve(U2,[i,0])
        dm=NoMeasureX_reset(dm)
    return dm

# - Get_ExpectValue_ZZ(dm,N): Computes correlators <ZiZj>^2 (distance |i-j| from 1 to N-1), returns a list of length N-1
# - Get_fidelity(dm,N): Computes fidelity correlators (distance from 1 to N-1), returns a list of length N-1
def Get_ExpectValue_ZZ(dm,N):
    EV_list=[]
    for i in range(N-1):
        EV_list.append(np.real(dm.expectation_value(Pauli('ZZ'),[1,2+i]))**2)
    return EV_list

# def Get_fidelity(dm,N):
#    Fidelity_list=[]
#    for i in range(N-1):
#       dm2 = dm.evolve(Pauli('ZZ'),[1,2+i])
#       Fidelity_list.append(state_fidelity(dm,dm2))
#    return Fidelity_list

# Run_ZZX_NoXM(sa, Num_qubits, Depth, U1, U2): Evolution function. 'sa' represents the sample index used for parallelization, 'Num_qubits' is the length of the 1D chain, 'Depth' is the evolution depth, and 'U1', 'U2' are operators containing measurement strengths J and h used to simulate weak measurement evolution. Returns a 3-dimensional list.
def Run_ZZX_NoXM(sa,Num_qubits,Depth,U1,U2):
    N = Num_qubits
    dm=init(N)
    EV_Sam=[]
#    FL_Sam=[]

# evolve the system for a given number of Depths, and get the expectation values
    for j in range(Depth-1):
        dm=evol_ZZ(dm,U1)
        dm=evol_X(dm,U2)
        if j % 3 == 0:
            EV_Sam.append(Get_ExpectValue_ZZ(dm,N))
#            FL_Sam.append(Get_fidelity(dm,N))
    return EV_Sam
#    return FL_Sam

# Create a process pool based on the given sample size and number of processes, and return the results for different samples
def parallel_square_sum(Func, Sample, num_processes):
    numbers = list(range(1, Sample + 1))    
    # Create a process pool
    with multiprocessing.Pool(processes=num_processes) as pool:
        result = pool.map(Func, numbers)
    return result

if __name__ == "__main__":
    Sample = 4000  # Set the number of samples
    num_processes = multiprocessing.cpu_count()  # Use the number of CPU cores as the number of processes

    t1 = 0.1
    start_time = time.time()
    for i in range(7):
        t2 = np.round(0.01 + 0.03*i,2) # Set X weak measurement strength: 0.01, 0.04, 0.07, 0.10, 0.13, 0.16, 0.19
        U1 = HamiltonianGate(Pauli('ZZZ'),t1)
        U2 = HamiltonianGate(Pauli('ZX'),t2)
        Partial_Run = partial(Run_ZZX_NoXM, Num_qubits=6, Depth=50, U1=U1, U2=U2)  # Set the length of the 1D chain and evolution depth
        Fresult = parallel_square_sum(Partial_Run, Sample, num_processes) # Returns the sampling results, a three-dimensional list
        Final = np.sum(Fresult, axis=0) / Sample  # Average the results over all samples; the result is a 2D list, each row corresponds to the squared expectation value of ZZ at different distances for each evolution depth
# Set the file names for saving data
        name1 = "EV_Data/EV_Crit_1/" + "ZZX_EV_NoXM"+"+"+str(6)+ "+" + str(t1) + "+" + str(t2) + ".npy"
#        name2 = "EV_Data/EV_Crit_1/" + "ZZX_FL_NoXM"+"+"+str(6)+ "+" + str(t1) + "+" + str(t2) + ".npy" 
        np.save(name1,Final)
#        np.save(name2,Final)
    end_time = time.time()
    print(end_time - start_time)