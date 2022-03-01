import os
localdir = os.path.dirname(__file__)
if localdir == "":
	localdir = "."
os.chdir(localdir)

load("MILP_function.sage")
load("divtrail.sage")
import sys, getopt

def genJoltikModel(rMax, sboxModel="Hull", linModel="CX", lastLin=True,firstSbox=True):
	"""
	Generate the model for Joltik over #rMax rounds
	#lastLin defines if the linear layer is applied on the last round
	#firstSbox defines if the first Sbox layer is applied

	The #sboxModel parameter defines how the Sbox is modelized :
	- "Hull" use the Convex Hull technique
	- "Simple" use the simplified constraint with PWL

	The #linModel parameter defines how the linear layer is modelized :
	- "CX" modelize the linear layer using the classical copy+xor technique
	- "QM" use the Quin-McCluskey algorithm as in Abdelkhalek,Sasaki,Todo,Tolba,Youssef
	- "Simple" modelize the linear layer with the simplified constraint w(x) = w(y)
	"""

	if sboxModel not in ["Hull", "Simple"]:
		print("Unknown value for the sboxModel parameter, default to Hull")
		sboxModel = "Hull"

	if linModel not in ["CX", "Simple", "QM"]:
		print("Unknown value for the linModel parameter, default to CX")
		linModel = "CX"

	blockSize = 64
	nbSbox = 16
	sboxSize = 4

	S = [14,4,11,2,3,8,0,9,1,10,7,15,6,12,5,13]
	# P = [0,13,10,7,4,1,14,11,8,5,2,15,12,9,6,3]
	#Expand P at a bit level
	# Pbit = [0 for i in range(blockSize)]
	# for isbox in range(nbSbox):
	# 	for j in range(sboxSize):
	# 		Pbit[isbox*sboxSize+j] = P[isbox]*sboxSize+j
	Pbit = [0,1,2,3,52,53,54,55,40,41,42,43,28,29,30,31,16,17,18,19,4,5,6,7,56,57,58,59,44,45,46,47,32,33,34,35,20,21,22,23,8,9,10,11,60,61,62,63,48,49,50,51,36,37,38,39,24,25,26,27,12,13,14,15]
	fullM = matrix(GF(2), [[1,0,0,0,0,0,1,0,1,1,0,0,1,1,1,0],
						   [0,1,0,0,0,0,1,1,0,0,1,0,0,0,0,1],
						   [0,0,1,0,1,0,0,1,0,0,0,1,1,0,0,0],
						   [0,0,0,1,0,1,0,0,1,0,0,0,1,1,0,0],
						   [0,0,1,0,1,0,0,0,1,1,1,0,1,1,0,0],
						   [0,0,1,1,0,1,0,0,0,0,0,1,0,0,1,0],
						   [1,0,0,1,0,0,1,0,1,0,0,0,0,0,0,1],
						   [0,1,0,0,0,0,0,1,1,1,0,0,1,0,0,0],
						   [1,1,0,0,1,1,1,0,1,0,0,0,0,0,1,0],
						   [0,0,1,0,0,0,0,1,0,1,0,0,0,0,1,1],
						   [0,0,0,1,1,0,0,0,0,0,1,0,1,0,0,1],
						   [1,0,0,0,1,1,0,0,0,0,0,1,0,1,0,0],
						   [1,1,1,0,1,1,0,0,0,0,1,0,1,0,0,0],
						   [0,0,0,1,0,0,1,0,0,0,1,1,0,1,0,0],
						   [1,0,0,0,0,0,0,1,1,0,0,1,0,0,1,0],
						   [1,1,0,0,1,0,0,0,0,1,0,0,0,0,0,1]])


	#--- Precomputations ---

	#Compute the ANF of S
	# (BPR, anfS) = SBOX_ANF(S)
	BPR = BooleanPolynomialRing(4, ["x"+str(i) for i in range(4)])
	(x0, x1, x2, x3) = BPR.gens()
	anfS = [x0*x1*x2 + x0*x1 + x0*x2 + x0*x3 + x1*x2*x3 + x1*x3 + x1 + x2 + x3,
			x0*x1 + x0 + x1*x2*x3 + x1*x2 + x1*x3 + x2*x3 + x3 + 1,
			x1*x2 + x1 + x2 + x3 + 1,
			x0 + x2*x3 + x2 + x3 + 1]
	#Generate the inequalities for S
	divTableS = SboxDivTrailTable(anfS,False) #second argument indicates whether redundant propagations should be excluded
	ineqS = None
	outLB = None
	ineqMC = None

	if sboxModel == "Hull":
		ineqS = sboxReducedInequalities(divTableS)
	elif sboxModel == "Simple":
		#Generate the lower bounds for the output weight of S
		outLB = computeLBWeightSbox(divTableS,sboxSize)

	if linModel == "QM":
		ineqMC = load("JoltikMC_ineq.sobj")


	#--- Model generation ---

	#Create the model
	modelName = "Joltik_"+str(rMax)+"r_"+sboxModel+"_"+linModel
	if lastLin:
		modelName += "_lastLin"
	if not firstSbox:
		modelName += "_noFirstSbox"
	m = Model(modelName)

	#Create the variables
	x = []
	y = []
	z = []
	if lastLin:
		#x[r] --S--> y[r] --P--> z[r] --MC--> x[r+1]
		if firstSbox:
			x = [[m.addVar(vtype=GRB.BINARY, name="x"+str(r)+"_"+str(j)) for j in range(blockSize)] for r in range(rMax+1)]
		else: #empty list for first round
			x = [[]] + [[m.addVar(vtype=GRB.BINARY, name="x"+str(r)+"_"+str(j)) for j in range(blockSize)] for r in range(1,rMax+1)]

		y = [[m.addVar(vtype=GRB.BINARY, name="y"+str(r)+"_"+str(j)) for j in range(blockSize)] for r in range(rMax)]

		z = [[m.addVar(vtype=GRB.BINARY, name="z"+str(r)+"_"+str(j)) for j in range(blockSize)] for r in range(rMax)]
	else:
		#x[r] --S--> y[r] --P--> z[r] --MC--> x[r+1]
		#No MC on last round so stops at y[rMax-1]
		if firstSbox:
			x = [[m.addVar(vtype=GRB.BINARY, name="x"+str(r)+"_"+str(j)) for j in range(blockSize)] for r in range(rMax)]
		else: #empty list for first round
			x = [[]] + [[m.addVar(vtype=GRB.BINARY, name="x"+str(r)+"_"+str(j)) for j in range(blockSize)] for r in range(1,rMax)]

		y = [[m.addVar(vtype=GRB.BINARY, name="y"+str(r)+"_"+str(j)) for j in range(blockSize)] for r in range(rMax)]

		z = [[m.addVar(vtype=GRB.BINARY, name="z"+str(r)+"_"+str(j)) for j in range(blockSize)] for r in range(rMax-1)]
	m.update()

	#Check if we need the first Sbox layer
	roundFirstSbox = 0
	if not firstSbox:
		roundFirstSbox = 1

	#Modelize the Sbox layers
	#y[r] = S(x[r])
	if sboxModel == "Hull":
		for r in range(roundFirstSbox,rMax):
			for i in range(nbSbox):
				invar = [x[r][sboxSize*i + j] for j in range(sboxSize)]
				outvar = [y[r][sboxSize*i + j] for j in range(sboxSize)]
				addSboxConstr(m, ineqS, invar, outvar)

	elif sboxModel == "Simple":
		for r in range(roundFirstSbox,rMax):
			for i in range(nbSbox):
				invar = [x[r][sboxSize*i + j] for j in range(sboxSize)]
				outvar = [y[r][sboxSize*i + j] for j in range(sboxSize)]
				addSimplifiedSboxConstr(m,invar, outvar, outLB, boundSuffix="s"+str(r)+"_"+str(i))

	#Check if we need the last linear layer
	roundLastLin = rMax
	if not lastLin:
		roundLastLin = rMax-1

	#Modelize the Permutation layers
	#z[r] = P(y[r])
	for r in range(roundLastLin):
		for i in range(blockSize):
			m.addConstr(z[r][Pbit[i]] == y[r][i])

	#Modelize the linear layer
	#x[r+1] = MC(z[r])
	if linModel == "CX":
		#Modelize MC with the Copy+XOR technique
		for r in range(roundLastLin):
			for col in range(4):
				invar = [z[r][16*col + i] for i in range(16)]
				outvar = [x[r+1][16*col + i] for i in range(16)]
				addLinearConstrCX(m, invar, outvar, fullM)

	elif linModel == "QM":
		#Modelize MC with QM
		for r in range(roundLastLin):
			for col in range(4):
				invar = [z[r][16*col + i] for i in range(16)]
				outvar = [x[r+1][16*col + i] for i in range(16)]
				addSboxConstr(m, ineqMC, invar, outvar)

	elif linModel == "Simple":
		for r in range(roundLastLin):
			for col in range(4):
				invar = [z[r][16*col + i] for i in range(16)]
				outvar = [x[r+1][16*col + i] for i in range(16)]
				addSimplifiedLinearConstr(m, invar, outvar)

	m.update()
	m.write("./../models/"+modelName+".mps")



if __name__ == "__main__":
	rMax = 1
	sboxModel = "Hull"
	linModel = "CX"
	lastLin = True
	firstSbox = True
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hr:s:m:", ["lastLin=","firstSbox="])
	except getopt.GetoptError:
		print("genJoltik.sage -r <rMax> -s <sboxModel> -m <linModel> --lastLin <lastLin> --firstSbox <firstSbox>")
		print("<rMax> is the number of rounds")
		print("<sboxModel> is the technique used for the modelization of the Sbox")
		print("           Possible options : Hull Simple")
		print("<linModel> is the technique used for the modelization of MC")
		print("           Possible options : CX QM Simple")
		print("<lastLin> set if the last linear layer is applied")
		print("          Possible options : True False")
		print("<firstSbox> set if the first sbox layer is applied")
		print("          Possible options : True False")
		sys.exit(2)
	for opt, arg in opts:
		if opt == "-h":
			print("genJoltik.sage -r <rMax> -s <sboxModel> -m <linModel> --lastLin <lastLin> --firstSbox <firstSbox>")
			print("<rMax> is the number of rounds")
			print("<sboxModel> is the technique used for the modelization of the Sbox")
			print("           Possible options : Hull Simple")
			print("<linModel> is the technique used for the modelization of MC")
			print("           Possible options : CX QM Simple")
			print("<lastLin> set if the last linear layer is applied")
			print("          Possible options : True False")
			print("<firstSbox> set if the first sbox layer is applied")
			print("          Possible options : True False")
			sys.exit()
		elif opt == "-r":
			rMax = int(arg)
		elif opt == "-s":
			sboxModel = arg
		elif opt == "-m":
			linModel = arg
		elif opt == "--lastLin":
			if arg == "True":
				lastLin = True
			elif arg == "False":
				lastLin = False
		elif opt == "--firstSbox":
			if arg == "True":
				firstSbox = True
			elif arg == "False":
				firstSbox = False

	genJoltikModel(rMax,sboxModel,linModel,lastLin,firstSbox)