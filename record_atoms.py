###Script for writing atoms in Pfin(N) with fixed maximum to local directory

from PfinN import * #Need functions for Pfin(N) subset arithmetic
from time import time
import multiprocessing as mp

atom_ref_path = 'atoms/atoms_max'

def ret_if_atom(X): #returns X if X is an atom, otherwise None
    if is_atom(X):
        return X
    else:
        pass

def atoms(n): #returns the list of atoms which have max n 
    p = mp.Pool()
    return [a for a in p.map(ret_if_atom,interval_subsets(n)) if a != None]

#Uses synchronous multiprocessing for the sake of keeping atom list "in order"
#Future improvements:
#**try async for better speed
#**avoid list implementation


def write_atoms(list_of_atoms,file): #for our purposes, list_of_atoms = atoms(n)
    for a in list_of_atoms:
        file.write(str(a)+'\n')
    file.close()


#### The following will record all atoms in Pfin(N) with max=n
#### to a file named 'atoms_maxn

import sys

n = int(sys.argv[1])

file = open('atoms/atoms_max'+str(n),'w')

write_atoms(atoms(n),file)

