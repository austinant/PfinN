import itertools as it
from time import time
import multiprocessing as mp


def powerset(X): #returns generator of all subsets of X.  Code borrowed from [1]
    s = list(X)
    return it.chain.from_iterable(
        map(set,it.combinations(s, r)) for r in range(len(s)+1))

def npowerset(X): #returns generator of all nonempty subsets of X.
    s = list(X)
    return it.chain.from_iterable(
        map(set,it.combinations(s, r)) for r in range(1,len(s)+1))

def subsets_cont(X,Y): #returns subsets of X containing Y
    X0 = X - Y
    for S in powerset(X0):
        yield S.union(Y)

def nonunits(X): #returns subsets of X containing 0 and some other element
    for S in npowerset(X-{0}):
        yield S.union({0})

def interval_subsets(n): #generator of subsets of [0,n] containing endpoints
    for S in powerset(set(range(1,n))):
        yield S.union({0,n})

####
####
####

def setsum(X,Y): #returns X+Y={x+y: (x,y) in (X,Y)}
    return set([x+y for x in X for y in Y])
    
def msum(*listofsets): #n-ary version of setsum
    S = listofsets[0]
    if len(listofsets) == 1:
        return S
    else:
        for i in range(1,len(listofsets)):
            S = setsum(S,listofsets[i])
        return S

def translate(X,c): #translates X by c
    return set([x+c for x in X])


def dilate(X,k): #dilates X by k
    return set([k*x for x in X])


def subsums(X): #returns set of subsums of a multiset or tuple X
    S = {0}
    for x in X:
        S = setsum(S,{0,x})
    return S


####
####
####

def cofactor(X,A): #returns the largest possible B with X = A + B
    B = X
    for a in A:
        B = B & translate(X,-a)
    return B

def divides(X,A): #decides whether A divides X
    return len(X) == len(setsum(A,cofactor(X,A)))


####
####
####

atom_ref_path = 'atoms/atoms_max'

def is_atom(X): #decides whether X is an atom
    for Y in nonunits(X-{max(X)}):
        if divides(X,Y):
            return False
    return True

def ret_if_atom(X): #returns X if X is an atom, otherwise None
    if is_atom(X):
        return X
    else:
        pass

def atoms(n): #returns the list of atoms which have max n 
    p = mp.Pool()
    return [a for a in p.map(ret_if_atom,interval_subsets(n)) if a != None]

def write_atoms(list_of_atoms,file):
    for a in list_of_atoms:
        file.write(str(a)+'\n')
    file.close()
