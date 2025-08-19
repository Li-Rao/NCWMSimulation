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

# Auxiliary operator: H gate and the transformation matrix UO that maps Z eigenstates to Y eigenstates
HO=1/np.sqrt(2)*np.array([[1.0, 1.0],[1.0,-1.0]])
UO=1/np.sqrt(2)*np.array([[1.0, 0-1.0j],[1,0+1.0j]])


# N is the number of qubitsInitialize the state vector for N qubits in the |+> state.
def init(N): 
    stavec = Statevector.from_int(0,2**(N+1))
    for i in range(N+1):
        stavec= stavec.evolve(HO,[i])
    return stavec

# Reset the first qubit after measurement to the |+> state
def measure_reset(sv): 
    sv = sv.evolve(UO,[0])
    meas,sv = sv.measure([0])
    sv=sv.reset([0])
    sv=sv.evolve(HO,[0])
    return sv

# U1 is the Hamiltonian gate for ZZZ interaction. Evolves the state vector sv under ZZ WeakMeas for a round of time.
def evol_ZZ(sv,U1): 
    N = sv.num_qubits-1
    for i in range(N-1):
        sv=sv.evolve(U1,[0,i+1,i+2])
        sv = measure_reset(sv)
    sv = measure_reset(sv)
    return sv

# U2 is the Hamiltonian gate for ZX interaction. Evolves the state vector sv under X WeakMeas for a round of time.
def evol_X(sv,U2):
    N = sv.num_qubits-1
    for i in range(1,N+1):
        sv=sv.evolve(U2,[i,0])
        sv = measure_reset(sv)
    return sv 


# Returns a list of expectation values for each qubit pair (1,2), (1,3), ..., (1,N).
# The expectation value is calculated using the Pauli operator 'ZZ', which is squared to get the final result.
def Get_ExpectValue_ZZ(sv,N):
    EV_list=[]
    for i in range(N-1):
        EV_list.append(np.real(sv.expectation_value(Pauli('ZZ'),[1,2+i]))**2)
    return EV_list

# Run_ZZX function runs the evolution for a given number of qubits, depth, and Hamiltonian gates U1 and U2.
# It initializes the state vector, evolves it through the ZZ and X WeakMeas (complete readout) and collects expectation values at specified intervals.
# The results are returned as a order-2 list.
# sa is a sample number used for identification in the output for parallel processing.
def Run_ZZX(sa,Num_qubits,Depth,U1,U2):
    N = Num_qubits
    sv=init(N)
    EV_Sam=[]
    for j in range(Depth):
        sv=evol_ZZ(sv,U1)
        sv=evol_X(sv,U2)
        if j % 3 == 0:
            EV_Sam.append(Get_ExpectValue_ZZ(sv,N))
            print(j)
    print(10000+sa)
    return EV_Sam


# This function runs the Run_ZZX function in parallel for a range of sample numbers.
# It uses a multiprocessing pool to distribute the workload across multiple processes.
# Func is the function to be executed in parallel, Sample is the number of samples, and num_processes is the number of processes to use.
# returns a list of results from the parallel execution: order-3 list
def parallel_square_sum(Func, Sample, num_processes):
    numbers = list(range(1, Sample + 1))
    # Create a process pool
    with multiprocessing.Pool(processes=num_processes) as pool:
        result = pool.map(Func, numbers)  
    return result 

if __name__ == "__main__":
    Sample = 4000
    num_processes = multiprocessing.cpu_count()  # Use the number of CPU cores as the number of processes

    t1 = 0.1
    start_time = time.time()
    for i in range(7):
        t2 = np.round(0.01 + 0.03*i,2)  # Set X weak measurement strength: 0.01, 0.04, 0.07, 0.10, 0.13, 0.16, 0.19
        U1 = HamiltonianGate(Pauli('ZZZ'),t1)
        U2 = HamiltonianGate(Pauli('ZX'),t2)
        Partial_Run=partial(Run_ZZX,Num_qubits=6,Depth=50,U1=U1,U2=U2) # function for parallel execution: set the length of the 1D chain and evolution depth
        Fresult = parallel_square_sum(Partial_Run,Sample, num_processes) # parallel execution: returns the sampling results, a three-dimensional list
        Final = np.sum(Fresult,axis=0)/Sample # average over samples
        # Set the file names for saving data
        name1 = "EV_Data/EV_Crit_0/" + "ZZX_EV"+"+"+str(6)+ "+" + str(t1) + "+" + str(t2) + ".npy" 
        np.save(name1,Final)
    end_time = time.time()
    print(end_time - start_time)