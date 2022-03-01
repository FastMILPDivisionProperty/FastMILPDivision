import os
localdir = os.path.dirname(__file__)
if localdir == "":
	localdir = "."
os.chdir(localdir)

load("MILP_function.sage")
load("divtrail.sage")
import sys, getopt

def genCRAFTModel(rMax, sboxModel="Hull", linModel="Lbox", lastLin=True, smartLin=True,firstSbox=True):
	"""
	Generate the model for CRAFT over #rMax rounds
	#lastLin defines if the linear layer is applied on the last round
	#smartLin only has effect if linModel=="Simple" and defines if MC is seen as independent Lboxes or as a single matrix
	#firstSbox defines if the first Sbox layer is applied

	The #sboxModel parameter defines how the Sbox is modelized :
	- "Hull" use the Convex Hull technique
	- "Simple" use the simplified constraint with PWL

	The #linModel parameter defines how the linear layer is modelized :
	- "Lbox" modelize the linear layer as the parallel application of the same linear Sbox (thanks to the specific form of the matrix)
	- "ZR" modelize the linear layer using the technique from Zhang and Rijmen
	- "CX" modelize the linear layer using the classical copy+xor technique
	- "Simple" modelize the linear layer with the simplified constraint w(x) = w(y)
	"""

	if sboxModel not in ["Hull", "Simple"]:
		print("Unknown value for the sboxModel parameter, default to Hull")
		sboxModel = "Hull"

	if linModel not in ["Lbox", "ZR", "CX", "Simple"]:
		print("Unknown value for the linModel parameter, default to Lbox")
		linModel = "Lbox"

	blockSize = 64
	nbSbox = 16
	sboxSize = 4

	S = [0xc, 0xa, 0xd, 0x3, 0xe, 0xb, 0xf, 0x7, 0x8, 0x9, 0x1, 0x5, 0x0, 0x2, 0x4, 0x6]
	# P = [15,10,9,4,3,6,5,8,7,2,1,12,11,14,13,0]
	#Expand P at a bit level
	# Pbit = [0 for i in range(blockSize)]
	# for isbox in range(nbSbox):
	# 	for j in range(sboxSize):
	# 		Pbit[isbox*sboxSize+j] = P[isbox]*sboxSize+j
	Pbit = [60,61,62,63,40,41,42,43,36,37,38,39,16,17,18,19,12,13,14,15,24,25,26,27,20,21,22,23,32,33,34,35,28,29,30,31,8,9,10,11,4,5,6,7,48,49,50,51,44,45,46,47,56,57,58,59,52,53,54,55,0,1,2,3]
	M = matrix(GF(2), [[1,0,1,1], [0,1,0,1], [0,0,1,0], [0,0,0,1]])
	
	#--- Precomputations ---

	#Compute the ANF of S
	# (BPR, anfS) = SBOX_ANF(S)
	BPR = BooleanPolynomialRing(4, ["x"+str(i) for i in range(4)])
	(x0, x1, x2, x3) = BPR.gens()
	anfS = [x0*x1*x2 + x0*x1*x3 + x0*x2 + x0*x3 + x1*x2*x3 + x1,
			x0*x2 + x0*x3 + x0 + x2*x3 + x2,
			x0*x1*x2 + x0*x1*x3 + x0*x3 + x0 + x1*x2*x3 + x3 + 1,
			x0*x1*x3 + x0*x1 + x1*x2*x3 + x1*x3 + x2*x3 + 1]
	#Generate the inequalities for S
	divTableS = SboxDivTrailTable(anfS,False) #second argument indicates whether redundant propagations should be excluded
	ineqS = None
	outLB = None

	if sboxModel == "Hull":
		ineqS = sboxReducedInequalities(divTableS)
	elif sboxModel == "Simple":
		#Generate the lower bounds for the output weight of S
		outLB = computeLBWeightSbox(divTableS,sboxSize)

	fullM = None
	if linModel == "ZR" or linModel == "CX":
		idmat = identity_matrix(GF(2),4,4)
		zmat = zero_matrix(GF(2),4,4)
		fullM = block_matrix([[idmat,zmat ,idmat,idmat], 
							  [zmat ,idmat,zmat ,idmat], 
							  [zmat ,zmat ,idmat,zmat ], 
							  [zmat ,zmat ,zmat ,idmat]], subdivide=False)

	elif linModel == "Lbox":
		#Compute the inequalities for the Lbox
		# V = vectorF2n(4)
		# L = [M*vector(x) for x in V]
		# MCbox = [None for _ in range(16)]
		# for i in range(16):
		# 	x = 0
		# 	for j in range(4):
		# 		if L[i][j] == 1:
		# 			x += 2^j
		# 	MCbox[i] = x

		# (BPR,anfMC) = SBOX_ANF(MCbox)
		anfMC = [x0 + x2 + x3, x1 + x3, x2, x3]
		divTableMC = SboxDivTrailTable(anfMC,False) #second argument indicates whether redundant propagations should be excluded
		ineqMC= sboxReducedInequalities(divTableMC)

	#--- Model generation ---

	#Create the model
	modelName = "CRAFT_"+str(rMax)+"r_"+sboxModel+"_"+linModel
	if lastLin:
		modelName += "_lastLin"
	if linModel=="Simple" and smartLin:
		modelName += "_smartLin"
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
	#The modelization depends on the method selected by the #linModel flag
	if linModel == "Lbox":
		#Modelize MC as linear sboxes
		for r in range(roundLastLin):
			for col in range(4): #For each column of the state
				for offset in range(4): #For each set of 4 bits
					invar = [z[r][16*col + offset + i*4] for i in range(4)]
					outvar = [x[r+1][16*col + offset + i*4] for i in range(4)]
					addSboxConstr(m,ineqMC,invar,outvar)

	elif linModel == "ZR":
		#Modelize MC with the Zhang-Rijmen technique
		for r in range(roundLastLin):
			for col in range(4):
				invar = [z[r][16*col + i] for i in range(16)]
				outvar = [x[r+1][16*col + i] for i in range(16)]
				addLinearConstrZR(m, invar, outvar, fullM, 4)

	elif linModel == "CX":
		#Modelize MC with the Copy+XOR technique
		for r in range(roundLastLin):
			for col in range(4):
				invar = [z[r][16*col + i] for i in range(16)]
				outvar = [x[r+1][16*col + i] for i in range(16)]
				addLinearConstrCX(m, invar, outvar, fullM)

	elif linModel == "Simple":
		#The modelization depends on the method selected by the #smartLin flag
		if smartLin:
			#Consider MC as independant Lboxes
			for r in range(roundLastLin):
				for col in range(4): #For each column of the state
					for offset in range(4): #For each set of 4 bits
						invar = [z[r][16*col + offset + i*4] for i in range(4)]
						outvar = [x[r+1][16*col + offset + i*4] for i in range(4)]
						addSimplifiedLinearConstr(m,invar,outvar)

		else:
			#Consider MC as a whole
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
	linModel = "Lbox"
	lastLin = True
	smartLin = True
	firstSbox = True
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hr:s:m:", ["lastLin=","smartLin=","firstSbox="])
	except getopt.GetoptError:
		print("genCRAFT.sage -r <rMax> -s <sboxModel> -m <linModel> --lastLin <lastLin> --smartLin <smartLin> --firstSbox <firstSbox>")
		print("<rMax> is the number of rounds")
		print("<sboxModel> is the technique used for the modelization of the Sbox")
		print("           Possible options : Hull Simple")
		print("<linModel> is the technique used for the modelization of MC")
		print("           Possible options : Lbox ZR CX Simple")
		print("<lastLin> set if the last linear layer is applied")
		print("          Possible options : True False")
		print("<smartLin> set if the linear layer is seen as independent Lboxes or as a single matrix")
		print("          Only has effect if the \"Simple\" option is selected for -m")
		print("          Possible options : True False")
		print("<firstSbox> set if the first sbox layer is applied")
		print("          Possible options : True False")
		sys.exit(2)
	for opt, arg in opts:
		if opt == "-h":
			print("genCRAFT.sage -r <rMax> -s <sboxModel> -m <linModel> --lastLin <lastLin> --smartLin <smartLin> --firstSbox <firstSbox>")
			print("<rMax> is the number of rounds")
			print("<sboxModel> is the technique used for the modelization of the Sbox")
			print("           Possible options : Hull Simple")
			print("<linModel> is the technique used for the modelization of MC")
			print("           Possible options : Lbox ZR CX Simple")
			print("<lastLin> set if the last linear layer is applied")
			print("          Possible options : True False")
			print("<smartLin> set if the linear layer is seen as independent Lboxes or as a single matrix")
			print("          Only has effect if the \"Simple\" option is selected for -m")
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
		elif opt == "--smartLin":
			if arg == "True":
				smartLin = True
			elif arg == "False":
				smartLin = False
		elif opt == "--firstSbox":
			if arg == "True":
				firstSbox = True
			elif arg == "False":
				firstSbox = False

	genCRAFTModel(rMax,sboxModel,linModel,lastLin,smartLin,firstSbox)