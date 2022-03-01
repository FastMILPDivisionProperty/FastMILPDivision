reset()
load("MILP_function.sage")
load("divtrail.sage")
load("genHIGHT.sage")
import sys
from itertools import combinations

def genTableF0(reduceTable=True):
	n = 8
	BPR = BooleanPolynomialRing(n, ["x"+str(i) for i in range(n)])
	(x0, x1, x2, x3, x4, x5, x6, x7) = BPR.gens()
	anf = [x1 + x6 + x7,
		   x0 + x2 + x7,
		   x0 + x1 + x3,
		   x1 + x2 + x4,
		   x2 + x3 + x5,
		   x3 + x4 + x6,
		   x4 + x5 + x7,
		   x0 + x5 + x6]
	D = SboxDivTrailTable(anf,reduceTable)
	if(reduceTable):
		save(D,"F0_table_reduced")
	else:
		save(D,"F0_table_notReduced")

def genTableF1(reduceTable=True):
	n = 8
	BPR = BooleanPolynomialRing(n, ["x"+str(i) for i in range(n)])
	(x0, x1, x2, x3, x4, x5, x6, x7) = BPR.gens()
	anf = [x2 + x4 + x5,
		   x3 + x5 + x6,
		   x4 + x6 + x7,
		   x0 + x5 + x7,
		   x0 + x1 + x6,
		   x1 + x2 + x7,
		   x0 + x2 + x3,
		   x1 + x3 + x4]

	D = SboxDivTrailTable(anf,reduceTable)
	if(reduceTable):
		save(D,"F1_table_reduced")
	else:
		save(D,"F1_table_notReduced")

# genTableF0(True)
# print "F0 true done"
# genTableF0(False)
# print "F0 false done"
# genTableF1(True)
# print "F1 true done"
# genTableF1(False)
# print "F1 false done"

n = 8
BPR = BooleanPolynomialRing(n, ["x"+str(i) for i in range(n)])
(x0, x1, x2, x3, x4, x5, x6, x7) = BPR.gens()
D = load("F1_table_notReduced.sobj")
ineqF0 = load("ineqF1.sobj")

F0 = matrix(GF(2), [[0,1,0,0,0,0,1,1],
					[1,0,1,0,0,0,0,1],
					[1,1,0,1,0,0,0,0],
					[0,1,1,0,1,0,0,0],
					[0,0,1,1,0,1,0,0],
					[0,0,0,1,1,0,1,0],
					[0,0,0,0,1,1,0,1],
					[1,0,0,0,0,1,1,0]])

F1 = matrix(GF(2), [[0,0,1,0,1,1,0,0],
					[0,0,0,1,0,1,1,0],
					[0,0,0,0,1,0,1,1],
					[1,0,0,0,0,1,0,1],
					[1,1,0,0,0,0,1,0],
					[0,1,1,0,0,0,0,1],
					[1,0,1,1,0,0,0,0],
					[0,1,0,1,1,0,0,0]])

vecF2n = [tuple([1 if ((u & (1 << i)) != 0) else 0 for i in range(n)]) for u in range(1 << n)]

ctr = 0
ctrFalseButOpti = 0
ctrTrueButInfea = 0
falsePositive = []
for u in vecF2n:
	for v in vecF2n:
		sys.stdout.write(str(vecF2n.index(u)) + "/" + str(len(vecF2n)) +" "+ str(vecF2n.index(v)) + "/" + str(len(vecF2n)) + "         \r")
		sys.stdout.flush()

		vInTableU = v in D[u]

		m = Model()
		m.Params.OutputFlag = 0
		xvar = [m.addVar(vtype=GRB.BINARY, name="x"+str(i)) for i in range(n)]
		yvar = [m.addVar(vtype=GRB.BINARY, name="y"+str(i)) for i in range(n)]

		m.update()
		# addF0Constr(m, xvar, yvar)
		addSboxConstr(m, ineqF0, xvar, yvar)
		# addLinearConstrZR(m, xvar, yvar, F0, 1)
		for i in range(n):
			m.addConstr(xvar[i] == u[i])
			m.addConstr(yvar[i] == v[i])

		m.update()
		m.optimize()

		
		if not ((m.Status == GRB.OPTIMAL and vInTableU) or (m.Status == GRB.INFEASIBLE and (not vInTableU))):
			print(("u = " + str(u)))
			print(("v = " + str(v)))
			print(("v in Table = " + str(vInTableU)))
			if m.Status == GRB.OPTIMAL:
				print("Status : OPTIMAL")
				ctrFalseButOpti += 1
				falsePositive.append(u+v)
			elif m.Status == GRB.INFEASIBLE:
				print("Status : INFEASIBLE")
				ctrTrueButInfea += 1
			else:
				print(("Status : " + str(m.Status)))
			print("")
			ctr += 1
print(ctr)
print(("Number of cases where w is in the table but infeasible model : " + str(ctrTrueButInfea)))
print(("Number of cases where w is NOT in the table, but optimal model : " + str(ctrFalseButOpti)))