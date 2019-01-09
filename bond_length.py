from pymatgen import Structure
import numpy as np
from pymatgen.core.bonds import CovalentBond

def get_bond_lengths(file, species_1, species_2, r):
    """Calcualte all bond lengths within threshold range
    between species A and B

    Parameters
    ----------
    file : str
        path to file
    species_1 : str
        Species A
    species_2 : str
        Species B
    r : float
        threshold bond length

    Returns
    -------
    nieghbours : array like
        array of bond lengths
    """
    stuc = read_structure(file)
    c1, c2 = get_index(stuc, species_1, species_2)
    dist = get_distances(stuc, c1, c2)
    nieghbours = get_nieghbour_distances(dist, r)
    return nieghbours

def read_structure(file):
    """Read a CONTCAR/POSCAR file and return
    a dictionary of sites

    Parameters
    ----------
    file : str
        path to file

    Returns
    -------
    struc : object
        pymatgen structure object
    """
    struc = Structure.from_file("CONTCAR")

    return struc


def get_index(struc, species_1, species_2):
    """Get the index of species A and species B.
    
    Parameters
    ----------
    struc : object
        pymatgen structure object
    species_1 : str
        Atom label
    species_2 : str
        Atom label

    Returns
    -------
    c1 : array like
        indexes of species 1 in pymatgen structure object
    c2 : array like
        indexes of species 2 in pymatgen structure object     
    """
    x = struc.as_dict()
    y = x['sites']
    c1 = np.array([])
    c2 = np.array([])
    for i in range(0, len(y)):
        if y[i]['label'] == species_1:
            c1 = np.append(c1, i)
        elif y[i]['label'] == species_2:
            c2 = np.append(c2, i)
    c1 = c1.astype(int)
    c2 = c2.astype(int)
    return c1, c2


def get_distances(struc, species_1, species_2):
    """Get the distances between species 1 and 2

    Parameters
    ----------
    struc : object
        pymatgen structure object
    species_1 : array like
        indexes of species 1
    species_2 : array like
        indexes of species 2

    Returns
    -------
    dist : array like
        array of bond distances
        
    """
    dist = np.array([])
    for i in range(species_1[0], species_1[-1]+1):
        for j in range(species_2[0], species_2[-1]+1):
            x = struc.get_distance(i, j)
            dist = np.append(dist, x)
    dist = np.reshape(dist, (species_1.size, species_2.size))
    return dist

def get_nieghbour_distances(dist, r):
    """Return the bond lengths less than threshold

    Parameters
    ----------
    dist : array like
        bond lengths
    r : float
        threshold bond length
    
    Returns
    -------
    array like
        bond legnths less than threshold
    """
    return dist[dist<r]
