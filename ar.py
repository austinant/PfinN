from atom_record import *
import sys

n = int(sys.argv[1])

file = open('atoms/atoms_max'+str(n),'w')
#file = open('testat','w')

write_atoms(atoms(n),file)
