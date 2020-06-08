### Functions for investigating the natural power monoid, P_{fin,0}(\NN)
### Created by Austin Antoniou
###
### References:
### [1] https://docs.python.org/3/library/itertools.html
### [2] http://jeromekelleher.net/generating-integer-partitions.html
###

### Packages needed
import itertools as it
from time import time
#from numba import njit


################################################################################
### Functions for subsets ######################################################
################################################################################


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

def is_subset(X,Y): #decides whether X is a subset of Y
    for x in X:
        if not x in Y:
            return False
    return True


################################################################################
### Functions for setwise sums #################################################
################################################################################


def setsum(X,Y): #returns X+Y={x+y: (x,y) in (X,Y)}
    return set(x+y for x in X for y in Y)
    
def msum(*listofsets): #n-ary version of setsum
    S = listofsets[0]
    if len(listofsets) == 1:
        return S
    else:
        for i in range(1,len(listofsets)):
            S = setsum(S,listofsets[i])
        return S

def translate(X,c): #translates X by c
    return set(x+c for x in X)

def dilate(X,k): #dilates X by k
    return set(k*x for x in X)

def subsums(X): #returns set of subsums of a multiset or tuple X
    S = {0}
    for x in X:
        S = setsum(S,{0,x})
    return S

'''
### fast versions of these functions, using @njit decorator from numba

#@njit
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
#@njit
def translate(X,c): #translates X by c
    return set([x+c for x in X])

#@njit
def dilate(X,k): #dilates X by k
    return set([k*x for x in X])

#@njit
def subsums(X): #returns set of subsums of a multiset or tuple X
    S = {0}
    for x in X:
        S = setsum(S,{0,x})
    return S
'''

################################################################################
#Functions for testing sumset decompositions involving a given factor ##########
################################################################################


def cofactor(X,A): #returns the largest possible B with X = A + B
    B = X
    for a in A:
        B = B & translate(X,-a)
    return B

def cofactors(X,A): #returns generator of all B with X = A+B
    C = cofactor(X,A)
    for B in subsets_cont(C,{0,C}):
        if X == setsum(A,B):
            yield B

def vitals(X,A): #returns set of b contained in any B such that X = A+B
    vits = set()
    Xm = {}
    for x in X:
        Xm[x] = 0
    for a in A:
        for b in cofactor(X,A):
            if Xm[a+b] == 0:
                Xm[a+b] = b
            else:
                Xm[a+b] = -1
    return set(Xm[x] for x in X if Xm[x] >= 0)

def divides(X,A): #decides whether A divides X
    return len(X) == len(setsum(A,cofactor(X,A)))


################################################################################
#Functions for testing irreducibility and finding atoms ########################
################################################################################

atom_ref_path = 'atoms/atoms_max'

def is_atom(X): #decides whether X is an atom
    for Y in nonunits(X-{max(X)}):
        if divides(X,Y):
            return False
    return True

def atoms_in(X): #returns generator of atoms which are subsets of X
    for Y in nonunits(X):
        if is_atom(Y):
            yield Y
            
def atoms_dividing(X):
    for A in atoms_in(X):
        if divides(X,A):
            yield A

def record_atoms(n): #records a list of the atoms with max n
    cur_len = 1
    start = time()
    file = open(atom_ref_path+str(n),'w')
    for S in interval_subsets(n):
        if len(S) > cur_len:
            cur_len += 1
            print('Prev stage: '+str(round(time()-start))+'s. '+
                  'Now checking subsets of size '+str(cur_len))
            start = time()
        if is_atom(S):
            file.write(str(S)+'\n')
    file.close()

def atom_ref(n): #opens file containing strings of the atoms with max n
    return open(atom_ref_path+str(n),'r')


################################################################################
#Functions involving partitions ################################################
################################################################################

def partitions(n): #integer partitions of n, code borrowed from [2] above
    a = [0 for i in range(n + 1)]
    k = 1
    y = n - 1
    while k != 0:
        x = a[k - 1] + 1
        k -= 1
        while 2 * x <= y:
            a[k] = x
            y -= x
            k += 1
        l = k + 1
        while x <= y:
            a[k] = x
            a[l] = y
            yield tuple(sorted(a[:k + 2],reverse=True))
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        yield tuple(sorted(a[:k + 1],reverse=True))

'''def parts(n,k): #returns list of partitions of n into exactly k parts
    if k == 1:
        return [(n,)]
    P = []
    for l in range(1,int((n+k)/2)+1):
        for p in parts(l,k-1):
            P.append(tuple([n-l]+sorted(p,reverse=True)))
    return P'''


################################################################################
#Functions for finding factorizations, using partition type ####################
################################################################################


def facs_of_type(X,p): #returns all factorizations of X of p-type p
    F = set()
    atom_classes = []
    for m in p:
        f = open(atom_ref_path+str(m),'r')
        atom_classes.append([tuple(eval(a)) for a in f if is_subset(eval(a),X)])
        f.close()
    atom_tuples = it.product(*atom_classes)
    for aa in atom_tuples:
        if msum(*aa) == X:
            F.add(aa)
    return set(tuple(sorted(f,reverse=True)) for f in F)

def feasible_types(X): #generator with partitions which can occur as types of X
    for p in partitions(max(X)):
        if is_subset(subsums(p),X):
            yield p

def facs_by_type(X): #set of factorizations of X, assembled by p-types
    F = set()
    for p in feasible_types(X):
        F = F.union(facs_of_type(X,p))
    return F

def facs_by_force(X): #set of factorizations of X, not using p-types of X
    if is_atom(X):
        return {(tuple(sorted(X)),)}
    facs = set()
    for A in atoms_dividing(X):
        for B in nonunits(cofactor(X,A)):
            if setsum(A,B) == X:
                for f in facs_by_force(B):
                    facs.add(tuple(sorted(
                        (tuple(sorted(A)),)+f)))
    return facs


################################################################################

################################################################################
