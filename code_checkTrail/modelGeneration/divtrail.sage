
from sage.crypto.boolean_function import BooleanFunction
import sys

def SBOX_ANF(SBOX):
	"""
	Given a list SBOX of 2^n value
	return (P, y) where
	y is a list of n ANF for each output coordinate
	P is a common BooleanPolynomialRing for all element of y
	!!Has a memory leak, use custom_ANF if this is a problem
	"""

	#size of the sbox
	n = max([x.nbits() for x in SBOX])

	#Get the bit representation of each value
	SBOX_bits = [x.bits() + [0 for i in range(n-len(x.bits()))] for x in SBOX]

	#Create the boolean function corresponding to each bit of the output from its truth table
	#i.e. all bits of index i in SBOX-table for the i-th output bit
	B = [ BooleanFunction([x[i] for x in SBOX_bits]) for i in range(n)]

	#Set a common BooleanPolynomialRing for each output function
	y0 = B[0].algebraic_normal_form()
	P = y0.ring()

	#Compute the ANF of each output bit
	y = [P(b.algebraic_normal_form()) for b in B]

	return (P, y)

def vectorF2n(n):
	"""
	Create the list of all vector in {0,1}^n
	"""
	return [tuple(Integer(c).bits() + [0 for i in range(n-Integer(c).nbits())]) for c in range(2^n)]

def vecToInt(x):
	vx = 0
	for i in range(len(x)):
		if x[i] == 1:
			vx += (1 << i)

	return vx

def custom_ANF(B,SBOX):
	"""Compute the ANF of S in the ring B"""

	#size of the sbox
	n = max([x.nbits() for x in SBOX])
	#Get the bit representation of each value
	SBOX_bits = [x.bits() + [0 for i in range(n-len(x.bits()))] for x in SBOX]
	xvar = B.gens()

	ANF = [B(0) for _ in range(n)]
	F2n = vectorF2n(n)

	for i in range(n):
		#For each component of the ANF
		for u in F2n:
			#compute a_u for ANF[i]
			a_u = 0
			for x in F2n:
				if greater(u,x):
					vx = vecToInt(x)
					a_u ^^= SBOX_bits[vx][i]
			if a_u == 1:
				xu = prod(xvar[i] for i in range(n) if u[i] == 1)
				ANF[i] += xu

	return ANF

def greater(a,b):
	#return True if a[i] >= b[i] for all i
	#False otherwise
	for i in range(len(a)):
		if a[i] < b[i]:
			return False
	return True

def SboxDivTrail(y, k, reduceTable=True):
	"""
	input : 
	- y list of BooleanPolynomial representing the output ANF of the SBox
	- k the input division property
	output :
	 K the set of output division property
	"""

	n = len(k)
	len_out = len(y)
	P = y[0].ring()
	x = P.gens()

	S = set()
	for e in range(2^n):
		kbar = Integer(e).bits() + [0 for i in range(n-Integer(e).nbits())]
		if greater(kbar, k):
			S.add(tuple(kbar))

	F = set()
	for kbar in S:
		F.add(P(prod([x[i] for i in range(n) if kbar[i] == 1])))

	Kbar = set()

	for e in range(2^len_out):
		u = Integer(e).bits() + [0 for i in range(len_out-Integer(e).nbits())]
		puy = prod([y[i] for i in range(len_out) if u[i] == 1])
		puyMon = P(puy).monomials()
		contains = False
		for mon in F:
			if mon in puyMon:
				contains = True
				break

		if contains:
			Kbar.add(tuple(u))

	if reduceTable:
		K = []
		for kbar in Kbar:
			great = False
			for kbar2 in Kbar:
				if(kbar != kbar2 and greater(kbar, kbar2)):
					great = True
					break
			if not great:
				K.append(kbar)

		return K
	else:
		return list(Kbar)

def SboxDivTrailTable(y,reduceTable=True):
	"""
	Return a dict containing all possible division propagation of the SBOX, where y is a list containing the ANF of each output bits
	"""
	P = y[0].ring()
	n = P.n_variables()

	D = dict()
	for c in range(2^n):
		k = Integer(c).bits() + [0 for i in range(n-Integer(c).nbits())]
		k = tuple(k)
		D[k] = SboxDivTrail(y, k,reduceTable)

	return D

def cmpvec(x,y):
	if x == y:
		return int(0)
	elif sum(x) < sum(y):
		return int(-1)
	elif sum(x) == sum(y) and x > y:
		return int(-1)
	else:
		return int(1)

#Python 3 removed the use of a cmp function for sorting, wrapper to adapt it
def cmp_to_key(mycmp):
    '''Convert a cmp= function into a key= function'''
    class K(object):
        def __init__(self, obj, *args):
            self.obj = obj
        def __lt__(self, other):
            return mycmp(self.obj, other.obj) < 0
        def __gt__(self, other):
            return mycmp(self.obj, other.obj) > 0
        def __eq__(self, other):
            return mycmp(self.obj, other.obj) == 0
        def __le__(self, other):
            return mycmp(self.obj, other.obj) <= 0  
        def __ge__(self, other):
            return mycmp(self.obj, other.obj) >= 0
        def __ne__(self, other):
            return mycmp(self.obj, other.obj) != 0
    return K
		
def printDivTable(D):
	Lx = list(D.keys())
	Lx.sort(key=cmp_to_key(cmpvec))

	for dx in Lx:
		s = ""
		for ex in dx:
			s+=str(ex)
		s+= " : "
		Ly = D[dx]
		Ly.sort(key=cmp_to_key(cmpvec))
		for dy in Ly:
			for ey in dy:
				s+= str(ey)
			s += " "
		print(s)	

def printDivTableInt(D, nbit=4):
	T = [[] for _ in range(1 << nbit)]
	for u in D:
		uval = 0
		for i in range(nbit):
			if u[i] == 1:
				uval |= (1 << i)

		for v in D[u]:
			vval = 0
			for i in range(nbit):
				if v[i] == 1:
					vval |= (1 << i)
			T[uval].append(vval)

		T[uval].sort()

	for i in range(1 << nbit):
		print(str(i) + " : " + str(T[i]))

def expandMatrix(M, modulus=[1,1,0,1,1,0,0,0,1]):
	"""
	Expand the matrix M given on GF(2)[X]/modulus over GF(2)
	M is given as a matrix of integer, no need to coerce it to the actual field
	The modulus is given as a list of coefficient, in ascending degree order
	(i.e. what would be obtained with p.coefficients(sparse=False) where p is a Sage Polynomial)
	e.g : p = x^8 + x^4 + x^3 + x + 1 <=> modulus = [1, 1, 0, 1, 1, 0, 0, 0, 1]
	"""

	deg_mod = len(modulus)-1

	#Precomptuation of the matrices for multiplication by x^i
	#Matrix for multiplication by x
	Mx = Matrix(GF(2),deg_mod, deg_mod)
	#First subdiagonal set to 1
	for j in range(1,deg_mod):
		Mx[j,j-1] = 1
	#Last column set to the coefficients of the modulus
	for j in range(deg_mod):
		Mx[j,deg_mod-1] = modulus[j]

	#Compute each Mxi = matrix of the multiplication by x^i
	L_Mxi = [0 for i in range(deg_mod)]
	L_Mxi[0] = identity_matrix(GF(2),deg_mod)
	for i in range(1,deg_mod):
		L_Mxi[i] = Mx*L_Mxi[i-1]



	#Precomputation for the matrix conversion
	# multMatrix[a] will contains the matrix representing the multiplication of an element by a in GF(2^m) at a bit level
	multMatrix = dict()
	#Get the set of all coefficient appearing in M
	set_coeff = set([int(M[i][j]) for j in range(M.ncols()) for i in range(M.nrows())])
	#Compute each multiplication matrix we will need
	for a in set_coeff:
		Ma = matrix(GF(2),deg_mod,deg_mod)
		for b in range(deg_mod):
			if (a >> b) & 1:
				Ma += L_Mxi[b]
		multMatrix[a] = Ma

	#Now we can convert the diffusion matrix M on GF(2^m) to a matrix Mbit on GF2
	Mbit = matrix(GF(2), 0, M.nrows()*deg_mod)
	for i in range(M.nrows()):
		Mrow = matrix(GF(2), deg_mod, 0)
		for j in range(M.ncols()):
			Mrow = Mrow.augment(multMatrix[M[i][j]])
		Mbit = Mbit.stack(Mrow)

	return Mbit