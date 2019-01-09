import numpy as np



def read_oszicar(File):
    data = open(File, "r")
    Energy = np.array([])
    for line in data:
        x = line.split()
        if x[1] == "F=":
            Energy = np.append(Energy, x[4])
    Energy = Energy.astype(float)
    return Energy

def read_final_energy(data):
    last_line = data[-1]
    return last_line

def adsorption_energy(data, stoich, natoms, ads):
    adsorbtion_e = np.array([])
    for i in range(0, data.size):
        ae = (data[i] - ((natoms * ads) + stoich)) / natoms
        adsorbtion_e = np.append(adsorbtion_e, ae)
    return adsorbtion_e