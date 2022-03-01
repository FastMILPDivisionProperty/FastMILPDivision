#MILP related function, mostly to compute the inequalities for an Sbox and the search algorithm for division property

from sage.crypto.boolean_function import BooleanFunction
from gurobipy import *
import sys

def vectorF2n(n):
	"""
	Create the list of all vector in {0,1}^n
	"""
	return [tuple(Integer(c).bits() + [0 for i in range(n-Integer(c).nbits())]) for c in range(2^n)]

def ineqRepresentation(D):
	"""
	input : a dict D such that D[k] = [k1,...,kt] where k1,...,kt are the possible division propagation from k
	output : (Lineq, A) where Lineq is a list of linear inequalities and equalities representing Conv(A)
			  and A is the list of all vector k|k1,..., k|kt etc.
			  For a given E in Lineq, E[-1] == ">=" means that it's an inequality, and E[-1] == "==" means it's an equality
	Exemple :
	the vector [1,2,3,8,"=="] represent the inequality
	1*x[0] + 2*x[1] + 3*x[2] + 8 >= 0
	"""

	#create all the vector
	A = []
	for k in list(D.keys()):
		for k2 in D[k]:
			A.append(k+k2)

	#Create the corresponding polyhedron
	P = Polyhedron(vertices=A)

	#Get the inequalities
	#an inequality is represented as a list L of size 2n+1 (where n is the size of k)
	#the first item (i.e. L[0]) is the constant term, which is not really convenient, so we will put it at the end of the vector

	Ltmp = P.inequalities_list()
	Lineq = [l[1:] + [l[0]] for l in Ltmp]
	Lineq = [l + [">="] for l in Lineq]

	Ltmp = P.equations_list()
	Ltmp = [l[1:] + [l[0]] for l in Ltmp]
	Ltmp = [l + ["=="] for l in Ltmp]

	Lineq += Ltmp

	return(Lineq, A)

def ineqEvaluation(ineq, p):
	"""
	Return whether or not the inequality ineq is verified for point p
	ineq can also be an equation, ineq[-1] wil give the type
	"""

	if ineq[-1] == ">=":
		return sum([ineq[i]*p[i] for i in range(len(p))]) + ineq[-2] >= 0 #constant term of the inequality is at the index -2
	elif ineq[-1] == "==":
		return sum([ineq[i]*p[i] for i in range(len(p))]) + ineq[-2] == 0
	else:
		"Unknown relation character in inequality"


def ineqReduction(L, A):
	"""
	input : L is a list of linear inequalities representing Conv(A)
			 A is a list of points (tuples)		 
	output : a list of inequalities whose feasible solutions in {0,1}^n are A
	"""
	
	Ls = []
	n = len(A[0])

	#Create the list of all vector in {0,1}^n
	Btmp = vectorF2n(n)
	#Remove vector in A from B
	B = [b for b in Btmp if b not in A]
	Lbar = deepcopy(L)

	# print "Init"
	# print "B : " + str(B)
	# print "Lbar : " + str(Lbar)

	while(len(B) > 0):
		l = Lbar[0]
		Bs = [p for p in B if not ineqEvaluation(l, p)]

		#Get the inequality (l) in Lbar which maximize the number of points in B that do not satisfy this inequality (Bs)
		for ineq in Lbar:
			tmpBs = [p for p in B if not ineqEvaluation(ineq, p)]
			if len(tmpBs) > len(Bs):
				l = ineq
				Bs = tmpBs

		# print "Choosed to keep l = " + str(l)
		Ls.append(l)
		Lbar.remove(l)
		# s = ""
		for p in Bs:
			B.remove(p)
			# s += str(p) + ", "
		# print "Removed points " + s
		# print "Remaining points : " + str(B)
		# print ""

	return Ls

def sboxReducedInequalities(D):
	"""
	input : a dict D such that D[k] = [k1,...,kt] where k1,...,kt are the possible division propagation from k (obtained from SboxDivTrailTable)
	output : a reduced list of inequalities whose feasible solutions in {0,1}^n are exactly the possible division trail propagation
	"""

	(Lineq, A) = ineqRepresentation(D)
	return ineqReduction(Lineq, A)

def addSboxConstr(m, L, inputvar, outputvar):
	"""
	add to the model m the linear constraints obtained from L representing the propagation from inputvar to outputvar through one Sbox
	if n is the size of the input/output, then an inequality l in L can be seen with 3 part :
	 - l[:n] is the list of coefficient on the input variables
	 - l[n:2n] is the list of coefficient on the output variables
	 - l[2n] == l[-1] is the constant term of the inequality
	for example, with n = 2, l = [1,2,3,4,5], inputvar = [x0, x1] and outputvar = [y0, y2]
	then we have :
					 x0 + 2*x1 + 3*y0 + 4*y2 + 5 >= 0
	return a list containing the added Gurobi.Constr
	"""

	n = len(inputvar)
	c = []
	for l in L:
		if l[-1] == ">=":
			c.append(m.addConstr(quicksum([l[i]*inputvar[i] for i in range(n)] + [l[i]*outputvar[i-n] for i in range(n,2*n)]) + l[-2] >= 0))
		elif l[-1] == "==":
			c.append(m.addConstr(quicksum([l[i]*inputvar[i] for i in range(n)] + [l[i]*outputvar[i-n] for i in range(n,2*n)]) + l[-2] == 0))
		else:
			"Unknown relation character in inequality when adding Sbox constr"

	return c

def addXorConstr(m, x, y):
	"""
	Add an XOR constraint for y = XOR(x0,x1,...) in model m
	x is a list of variables
	y is a single variable
	"""
	m.addConstr(y == quicksum(x))

def addCopyConstr(m, x, y):
	"""
	Add a Copy constraint for (y0,y1,...) = Copy(x) in model m
	x is a single variable
	y is a list of variables
	"""
	m.addConstr(x == quicksum(y))

def addAndConstr(m, x, y, z):
	"""
	Add an And constraint for z = AND(x,y) in model m
	x,y,z are single variables
	"""
	m.addConstr(x <= z)
	m.addConstr(y <= z)
	m.addConstr(z <= x + y)

def addModAddConstr(m, x, y, z, cPrefix="c"):
	"""
	Add a modular addition constraint z = x + y mod 2^n (n = len(x)) in model m
	x,y,z are lists of variables, all of same size
	The carryPrefix is used on the naming of the carry variables
	"""
	assert (len(x) == len(y) and len(y) == len(z)), "Arguments for modular addition constraint not of the same size"

	n = len(x)

	#Algebraic computation of mod add is done as follow, ^ denotes xor
	# zi = x[i] ^ y[i] ^ c[i]
	# ci = (x[i-1] & y[i-1]) ^ (c[i-1] & (x[i-1] ^ y[i-1]))
	# c0 = 0

	#Modelization requires a bunch of copy operations and intermediate variables
	#In general, these are the required operations for 1 < i < n-1
	# - Copy of the base variables
	# (c0xi,c1xi,c2xi) = Copy(xi)
	# (c0yi,c1yi,c2yi) = Copy(yi)
	# (c0ci,c1ci) = Copy(ci)

	# - Intermediate variables
	# xyi = AND(c0xi, c0yi)
	# xXyi = XOR(c1xi, c1yi)
	# xXyci = AND(xXyi, c0ci)

	# - End computation
	# ci = XOR(xy[i-1], xXyc[i-1])
	# zi = XOR(c2xi, c2yi, c1ci)

	# However there are some small changes at the begining and at the end
	# - i = 0 (No carry on the first bit so much easier expression)
	# (c0x0, c1x0) = Copy(x0)
	# (c0y0, c1y0) = Copy(y0)
	# z0 = XOR(c1x0,c1y0)

	# - i = 1
	# (c0x1,c1x1,c2x1) = Copy(x1)
	# (c0y1,c1y1,c2y1) = Copy(y1)
	# (c0c1,c1c1) = Copy(c1)
	# xy1 = AND(c0x1, c0y1)
	# xXy1 = XOR(c1x1, c1y1)
	# xXyc1 = AND(xXy1, c0c1)
	# c1 = AND(c0x0,c0y0) (no carry c0 so easier expression, rest is the same)
	# z1 = XOR(c2x1,c2y1,c1c1)

	# - i = n-1 (No further carry computation so no need for copy and intermediates)
	# c[n-1] = XOR(xy[n-2], xXyc[n-2])
	# z[n-1] = XOR(x[n-1], y[n-1], c[n-1])

	#First, setup the intermediate variables
	#Carry variables
	c = [m.addVar(vtype=GRB.BINARY, name=cPrefix+"_"+str(i)) for i in range(1,n)]
	c = [None] + c #c0 = 0, not actually used but allows to have matching indexes
	m.update()

	#Copy variables
	cx = [[m.addVar(vtype=GRB.BINARY, name="copy"+str(j)+"_"+x[i].VarName) if not (j == 2 and i == 0) else None for i in range(n-1)] for j in range(3)]
	cy = [[m.addVar(vtype=GRB.BINARY, name="copy"+str(j)+"_"+y[i].VarName) if not (j == 2 and i == 0) else None for i in range(n-1)] for j in range(3)]
	cc = [[m.addVar(vtype=GRB.BINARY, name="copy"+str(j)+"_"+c[i].VarName) for i in range(1,n-1)] for j in range(2)]
	cc[0] = [None] + cc[0]
	cc[1] = [None] + cc[1] #no carry copy on i = 0
	m.update()

	#AND(x,y) variables
	xy = [m.addVar(vtype=GRB.BINARY, name=x[i].VarName + "_and_" + y[i].VarName) for i in range(1,n-1)]
	xy = [None] + xy #No xy0
	m.update()

	#XOR(x,y) variables
	xXy = [m.addVar(vtype=GRB.BINARY, name=x[i].VarName + "_xor_" + y[i].VarName) for i in range(1,n-1)] 
	xXy = [None] + xXy #No xXy0
	m.update()

	#AND(x+y,c) variables
	xXyc = [m.addVar(vtype=GRB.BINARY, name=xXy[i].VarName + "_and_" + c[i].VarName) for i in range(1,n-1)]
	xXyc = [None] + xXyc #No xXyc0
	m.update()


	#Now constraints
	#First i = 0

	addCopyConstr(m, x[0], [cx[0][0], cx[1][0]])
	addCopyConstr(m, y[0], [cy[0][0], cy[1][0]])
	addXorConstr(m, [cx[1][0], cy[1][0]], z[0])

	#i = 1
	addCopyConstr(m, x[1], [cx[0][1], cx[1][1], cx[2][1]])
	addCopyConstr(m, y[1], [cy[0][1], cy[1][1], cy[2][1]])
	addCopyConstr(m, c[1], [cc[0][1], cc[1][1]])
	addAndConstr(m, cx[0][1], cy[0][1], xy[1])
	addXorConstr(m, [cx[1][1], cy[1][1]], xXy[1])
	addAndConstr(m, xXy[1], cc[0][1], xXyc[1])
	addAndConstr(m, cx[0][0], cy[0][0], c[1])
	addXorConstr(m, [cx[2][1], cy[2][1], cc[1][1]], z[1])

	#i = 2 to n-2
	for i in range(2,n-1):
		addCopyConstr(m, x[i], [cx[0][i], cx[1][i], cx[2][i]])
		addCopyConstr(m, y[i], [cy[0][i], cy[1][i], cy[2][i]])
		addCopyConstr(m, c[i], [cc[0][i], cc[1][i]])
		addAndConstr(m, cx[0][i], cy[0][i], xy[i])
		addXorConstr(m, [cx[1][i], cy[1][i]], xXy[i])
		addAndConstr(m, xXy[i], cc[0][i], xXyc[i])
		addXorConstr(m, [xy[i-1], xXyc[i-1]], c[i])
		addXorConstr(m, [cx[2][i], cy[2][i], cc[1][i]], z[i])
	
	#i = n-1
	# c[n-1] = XOR(xy[n-2], xXyc[n-2])
	# z[n-1] = XOR(x[n-1], y[n-1], c[n-1])
	addXorConstr(m, [xy[n-2], xXyc[n-2]], c[n-1])
	addXorConstr(m, [x[n-1], y[n-1], c[n-1]], z[n-1])

def addModAddConstantConstr(m, x, z, cPrefix="ccst"):
	"""
	Add a modular addition constraint z = x + k mod 2^n (n = len(x)) in model m, where k is constant (e.g. key)
	x,z are lists of variables, all of same size
	The carryPrefix is used on the naming of the carry variables
	"""
	assert (len(x) == len(z)), "Arguments for modular addition constraint not of the same size"
	n = len(x)

	#Algebraic computation of mod add is done as follow, ^ denotes xor
	# zi = x[i] ^ y[i] ^ c[i]
	# ci = (x[i-1] & y[i-1]) ^ (c[i-1] & (x[i-1] ^ y[i-1]))
	# c0 = 0
	#However, here y is constant, and we know that an XOR or AND with a constant has no influence on bit-based DP
	#Hence, it can be simplified as (in term of constraints to set)
	# zi = x[i] ^ c[i]
	# ci = x[i-1] ^ (c[i-1] & x[i-1])
	# c0 = 0

	#Modelization requires a bunch of copy operations and intermediate variables
	#In general, these are the required operations for 1 < i < n-1
	# - Copy of the base variables
	# (c0xi,c1xi,c2xi) = Copy(xi)
	# (c0ci,c1ci) = Copy(ci)

	# - Intermediate variables
	# xci = AND(c0xi, c0ci)

	# - End computation
	# ci = XOR(c1x[i-1], xc[i-1])
	# zi = XOR(c2xi, c1ci)

	# However there are some small changes at the begining and at the end
	# - i = 0 (No carry on the first bit so much easier expression)
	# (c0x0, c1x0) = Copy(x0)
	# z0 = c1x0

	# - i = 1
	# (c0x1,c1x1,c2x1) = Copy(x1)
	# (c0c1,c1c1) = Copy(c1)
	# xc1 = AND(c0x1, c0c1)
	# c1 = c0x0 (no carry c0 so easier expression, rest is the same)
	# z1 = XOR(c2x1,c1c1)

	# - i = n-1 (No further carry computation so no need for copy and intermediates)
	# c[n-1] = XOR(c1x[n-2], xc[n-2])
	# z[n-1] = XOR(x[n-1], c[n-1])

	#First, setup the intermediate variables
	#Carry variables
	c = [m.addVar(vtype=GRB.BINARY, name=cPrefix+"_"+str(i)) for i in range(1,n)]
	c = [None] + c #c0 = 0, not actually used but allows to have matching indexes
	m.update()

	#Copy variables
	cx = [[m.addVar(vtype=GRB.BINARY, name="copy"+str(j)+"_"+x[i].VarName) if not (j == 2 and i == 0) else None for i in range(n-1)] for j in range(3)]
	cc = [[m.addVar(vtype=GRB.BINARY, name="copy"+str(j)+"_"+c[i].VarName) for i in range(1,n-1)] for j in range(2)]
	cc[0] = [None] + cc[0]
	cc[1] = [None] + cc[1] #no carry copy on i = 0
	m.update()

	#AND(x,c) variables
	xc = [m.addVar(vtype=GRB.BINARY, name=x[i].VarName + "_and_" + c[i].VarName) for i in range(1,n-1)]
	xc = [None] + xc #No xc0
	m.update()

	#Now constraints
	#First i = 0

	addCopyConstr(m, x[0], [cx[0][0], cx[1][0]])
	m.addConstr(z[0] == cx[1][0])

	#i = 1
	addCopyConstr(m, x[1], [cx[0][1], cx[1][1], cx[2][1]])
	addCopyConstr(m, c[1], [cc[0][1], cc[1][1]])
	addAndConstr(m, cx[0][1], cc[0][1], xc[1])
	m.addConstr(c[1] == cx[0][0])
	addXorConstr(m, [cx[2][1], cc[1][1]], z[1])

	#i = 2 to n-2
	for i in range(2,n-1):
		addCopyConstr(m, x[i], [cx[0][i], cx[1][i], cx[2][i]])
		addCopyConstr(m, c[i], [cc[0][i], cc[1][i]])
		addAndConstr(m, cx[0][i], cc[0][i], xc[i])

		addXorConstr(m, [cx[1][i-1], xc[i-1]], c[i])
		addXorConstr(m, [cx[2][i], cc[1][i]], z[i])
	
	#i = n-1
	addXorConstr(m, [cx[1][n-2], xc[n-2]], c[n-1])
	addXorConstr(m, [x[n-1], c[n-1]], z[n-1])



def addLinearConstrZR(m, x, y, M, n):
	"""
	Add constraints to modelize y = M(x) using the Zhang-Rijmen modelization
	#M must be given expanded on GF(2)
	#n is the size of the underlying field for M, i.e. the original matrix is over GF(2^n)
	--- the original matrix needs to be binary over GF(2^n) ---
	"""
	s = M.nrows()//n
	D = [[sig + i*n for i in range(s)] for sig in range(n)]

	for t in range(1,s):
		for sig in range(n):
			for seti in Combinations(D[sig], t):
				Lexpr = quicksum(-y[i] for i in seti)
				for k in range(n*s):
					a = GF(2)(0)
					for i in seti:
						a += M[i][k]
					Lexpr += int(a)*x[k]
				m.addConstr(Lexpr >= -(t-1))

	for i in range(n):
		m.addConstr(quicksum(x[i+k*n] for k in range(s)) == quicksum(y[i+k*n] for k in range(s)))

def addLinearConstrCX(m, x, y, M):
	"""
	Add constraints to modelize y = M(x) using the Copy+XOR modelization
	--- M needs to be given over GF(2) ---
	"""

	#First create the necessary copy variables
	nr = M.nrows()
	nc = M.ncols()
	c = [[0 for i in range(nc)] for j in range(nr)]
	for j in range(nc):
		for i in range(nr):
			if M[i][j] == 1:
				c[i][j] = m.addVar(vtype=GRB.BINARY, name="copy"+str(i)+"_"+x[j].VarName)

	m.update()

	#First step, copy constraints
	for j in range(nc):
		outvar = [c[i][j] for i in range(nr) if M[i][j] == 1]
		addCopyConstr(m, x[j], outvar)

	#Second step, xor constraints
	for i in range(nr):
		invar = [c[i][j] for j in range(nc) if M[i][j] == 1]
		addXorConstr(m, invar, y[i])



def addSimplifiedLinearConstr(m, x, y):
	"""
	Add the simplified linear constraint w(x) = w(y)
	"""
	m.addConstr(quicksum(x) == quicksum(y))

def addSimplifiedLinearConstr_v2(m,x,y,fullM,invFullM):
	n = fullM.nrows()
	m.addConstr(quicksum(x) == quicksum(y))
	for i in range(n):
		m.addConstr(y[i] <= quicksum(x[j] for j in range(n) if fullM[i][j] == 1))
		m.addConstr(x[i] <= quicksum(y[j] for j in range(n) if fullM[j][i] == 1))
		m.addConstr((1-y[i]) <= quicksum(1-x[j] for j in range(n) if invFullM[j][i] == 1))
		m.addConstr((1-x[i]) <= quicksum(1-y[j] for j in range(n) if invFullM[i][j] == 1))
	


def computeLBWeightSbox(divTable, sboxSize):
	"""
	Compute the lower bounds for the output weight of the sbox
	"""
	outLB = [sboxSize for _ in range(sboxSize+1)]
	for u in divTable:
		wu = sum(u)
		for v in divTable[u]:
			wv = sum(v)
			if wv < outLB[wu]:
				outLB[wu] = wv

	return outLB

def addSimplifiedSboxConstr(m, x, y, outLB, boundSuffix=""):

	xpts = [i for i in range(len(outLB))]

	bound0 = m.addVar(lb=0, ub=max(xpts), vtype=GRB.INTEGER, name="bound0_"+boundSuffix)
	bound1 = m.addVar(lb=0, ub=max(outLB), vtype=GRB.INTEGER, name="bound1_"+boundSuffix)

	m.addGenConstrPWL(bound0,bound1,xpts,outLB)
	m.addConstr(quicksum(x) == bound0)
	m.addConstr(quicksum(y) >= bound1)


#Printing fluff
def cmpvec(x,y):
	if x == y:
		return int(0)
	elif sum(x) < sum(y):
		return int(-1)
	elif sum(x) == sum(y) and x > y:
		return int(-1)
	else:
		return int(1)

def printDivSet(K):
	K.sort(cmpvec)
	for k in K:
		s = ""
		for i in range(len(k)):
			if i%4 == 0:
				s += " "
			s += str(k[i])
		print(s)