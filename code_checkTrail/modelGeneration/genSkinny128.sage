import os
localdir = os.path.dirname(__file__)
if localdir == "":
	localdir = "."
os.chdir(localdir)

load("MILP_function.sage")
load("divtrail.sage")
import sys, getopt

def genSkinny128Model(rMax, sboxModel="Simple", linModel="Lbox", lastLin=True, smartLin=True, firstSbox=True):
	"""
	Generate the model for Skinny over #rMax rounds
	#lastLin defines if the linear layer is applied on the last round
	#smartLin only has effect if linModel=="Simple" and defines if MC is seen as independent Lboxes or as a single matrix

	The #sboxModel parameter defines how the Sbox is modelized :
	- "QM" use the Quin-McCluskey algorithm as in Abdelkhalek,Sasaki,Todo,Tolba,Youssef
	- "Simple" use the simplified constraint with PWL

	The #linModel parameter defines how the linear layer is modelized :
	- "Lbox" modelize the linear layer as the parallel application of the same linear Sbox (thanks to the specific form of the matrix)
	- "ZR" modelize the linear layer using the technique from Zhang and Rijmen
	- "CX" modelize the linear layer using the classical copy+xor technique
	- "Simple" modelize the linear layer with the simplified constraint w(x) = w(y)
	"""

	if sboxModel not in ["Simple", "QM"]:
		print("Unknown value for the sboxModel parameter, default to Simple")
		sboxModel = "Simple"

	if linModel not in ["Lbox", "ZR", "CX", "Simple"]:
		print("Unknown value for the linModel parameter, default to Lbox")
		linModel = "Lbox"

	blockSize = 128
	nbSbox = 16
	sboxSize = 8

	# S = [0x65,0x4c,0x6a,0x42,0x4b,0x63,0x43,0x6b,0x55,0x75,0x5a,0x7a,0x53,0x73,0x5b,0x7b,
	# 	 0x35,0x8c,0x3a,0x81,0x89,0x33,0x80,0x3b,0x95,0x25,0x98,0x2a,0x90,0x23,0x99,0x2b,
	# 	 0xe5,0xcc,0xe8,0xc1,0xc9,0xe0,0xc0,0xe9,0xd5,0xf5,0xd8,0xf8,0xd0,0xf0,0xd9,0xf9,
	# 	 0xa5,0x1c,0xa8,0x12,0x1b,0xa0,0x13,0xa9,0x05,0xb5,0x0a,0xb8,0x03,0xb0,0x0b,0xb9,
	# 	 0x32,0x88,0x3c,0x85,0x8d,0x34,0x84,0x3d,0x91,0x22,0x9c,0x2c,0x94,0x24,0x9d,0x2d,
	# 	 0x62,0x4a,0x6c,0x45,0x4d,0x64,0x44,0x6d,0x52,0x72,0x5c,0x7c,0x54,0x74,0x5d,0x7d,
	# 	 0xa1,0x1a,0xac,0x15,0x1d,0xa4,0x14,0xad,0x02,0xb1,0x0c,0xbc,0x04,0xb4,0x0d,0xbd,
	# 	 0xe1,0xc8,0xec,0xc5,0xcd,0xe4,0xc4,0xed,0xd1,0xf1,0xdc,0xfc,0xd4,0xf4,0xdd,0xfd,
	# 	 0x36,0x8e,0x38,0x82,0x8b,0x30,0x83,0x39,0x96,0x26,0x9a,0x28,0x93,0x20,0x9b,0x29,
	# 	 0x66,0x4e,0x68,0x41,0x49,0x60,0x40,0x69,0x56,0x76,0x58,0x78,0x50,0x70,0x59,0x79,
	# 	 0xa6,0x1e,0xaa,0x11,0x19,0xa3,0x10,0xab,0x06,0xb6,0x08,0xba,0x00,0xb3,0x09,0xbb,
	# 	 0xe6,0xce,0xea,0xc2,0xcb,0xe3,0xc3,0xeb,0xd6,0xf6,0xda,0xfa,0xd3,0xf3,0xdb,0xfb,
	# 	 0x31,0x8a,0x3e,0x86,0x8f,0x37,0x87,0x3f,0x92,0x21,0x9e,0x2e,0x97,0x27,0x9f,0x2f,
	# 	 0x61,0x48,0x6e,0x46,0x4f,0x67,0x47,0x6f,0x51,0x71,0x5e,0x7e,0x57,0x77,0x5f,0x7f,
	# 	 0xa2,0x18,0xae,0x16,0x1f,0xa7,0x17,0xaf,0x01,0xb2,0x0e,0xbe,0x07,0xb7,0x0f,0xbf,
	# 	 0xe2,0xca,0xee,0xc6,0xcf,0xe7,0xc7,0xef,0xd2,0xf2,0xde,0xfe,0xd7,0xf7,0xdf,0xff]
	# P = [0,5,10,15,4,9,14,3,8,13,2,7,12,1,6,11]
	# #Expand P at a bit level
	# Pbit = [0 for i in range(blockSize)]
	# for isbox in range(nbSbox):
	# 	for j in range(sboxSize):
	# 		Pbit[isbox*sboxSize+j] = P[isbox]*sboxSize+j
	Pbit = [0,1,2,3,4,5,6,7,40,41,42,43,44,45,46,47,80,81,82,83,84,85,86,87,120,121,122,123,124,125,126,127,32,33,34,35,36,37,38,39,72,73,74,75,76,77,78,79,112,113,114,115,116,117,118,119,24,25,26,27,28,29,30,31,64,65,66,67,68,69,70,71,104,105,106,107,108,109,110,111,16,17,18,19,20,21,22,23,56,57,58,59,60,61,62,63,96,97,98,99,100,101,102,103,8,9,10,11,12,13,14,15,48,49,50,51,52,53,54,55,88,89,90,91,92,93,94,95]
	M = matrix(GF(2), [[1,0,1,1], [1,0,0,0], [0,1,1,0], [1,0,1,0]])

	#--- Precomputations ---

	#Compute the ANF of S
	# (BPR, anfS) = SBOX_ANF(S)
	#Generate the inequalities for S
	# divTableS = SboxDivTrailTable(anfS)
	divTableS = load("divTableSkinny128Sbox.sobj")
	ineqS = None
	outLB = None

	if sboxModel == "QM":
		ineqS = load("ineqSkinny128Sbox.sobj")

	elif sboxModel == "Simple":
		#Generate the lower bounds for the output weight of S
		outLB = computeLBWeightSbox(divTableS,sboxSize)

	fullM = None
	if linModel == "ZR" or linModel == "CX":
		idmat = identity_matrix(GF(2),8,8)
		zmat = zero_matrix(GF(2),8,8)
		fullM = block_matrix([[idmat,zmat ,idmat,idmat], 
							  [idmat,zmat ,zmat ,zmat], 
							  [zmat ,idmat,idmat,zmat], 
							  [idmat,zmat ,idmat,zmat]], subdivide=False)

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
		BPR = BooleanPolynomialRing(4, ["x"+str(i) for i in range(4)])
		(x0, x1, x2, x3) = BPR.gens()
		anfMC = [x0 + x2 + x3, x0, x1 + x2, x0 + x2]
		divTableMC = SboxDivTrailTable(anfMC,False) #second argument indicates whether redundant propagations should be excluded
		ineqMC= sboxReducedInequalities(divTableMC)

	#--- Model generation ---

	#Create the model
	modelName = "Skinny128_"+str(rMax)+"r_"+sboxModel+"_"+linModel
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
	if sboxModel == "QM":
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
				for offset in range(sboxSize): #For each set of 4 bits
					invar = [z[r][32*col + offset + i*sboxSize] for i in range(4)]
					outvar = [x[r+1][32*col + offset + i*sboxSize] for i in range(4)]
					addSboxConstr(m,ineqMC,invar,outvar)

	elif linModel == "ZR":
		#Modelize MC with the Zhang-Rijmen technique
		for r in range(roundLastLin):
			for col in range(4):
				invar = [z[r][32*col + i] for i in range(32)]
				outvar = [x[r+1][32*col + i] for i in range(32)]
				addLinearConstrZR(m, invar, outvar, fullM, sboxSize)

	elif linModel == "CX":
		#Modelize MC with the Copy+XOR technique
		for r in range(roundLastLin):
			for col in range(4):
				invar = [z[r][32*col + i] for i in range(32)]
				outvar = [x[r+1][32*col + i] for i in range(32)]
				addLinearConstrCX(m, invar, outvar, fullM)

	elif linModel == "Simple":
		#The modelization depends on the method selected by the #smartLin flag
		if smartLin:
			#Modelize MC as linear sboxes
			for r in range(roundLastLin):
				for col in range(4): #For each column of the state
					for offset in range(sboxSize): #For each set of 4 bits
						invar = [z[r][32*col + offset + i*sboxSize] for i in range(4)]
						outvar = [x[r+1][32*col + offset + i*sboxSize] for i in range(4)]
						addSimplifiedLinearConstr(m,invar,outvar)

		else:
			#Modelize MC with the Zhang-Rijmen technique
			for r in range(roundLastLin):
				for col in range(4):
					invar = [z[r][32*col + i] for i in range(32)]
					outvar = [x[r+1][32*col + i] for i in range(32)]
					addSimplifiedLinearConstr(m, invar, outvar)

	m.update()
	m.write("./../models/"+modelName+".mps")


if __name__ == "__main__":
	rMax = 1
	sboxModel = "Simple"
	linModel = "Lbox"
	lastLin = True
	smartLin = True
	firstSbox = True
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hr:s:m:", ["lastLin=","smartLin=","firstSbox="])
	except getopt.GetoptError:
		print("genSkinny64.sage -r <rMax> -s <sboxModel> -m <linModel> --lastLin <lastLin> --smartLin <smartLin> --firstSbox <firstSbox>")
		print("<rMax> is the number of rounds")
		print("<sboxModel> is the technique used for the modelization of the Sbox")
		print("           Possible options : QM Simple")
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
			print("genSkinny64.sage -r <rMax> -s <sboxModel> -m <linModel> --lastLin <lastLin> --firstSbox <firstSbox>")
			print("<rMax> is the number of rounds")
			print("<sboxModel> is the technique used for the modelization of the Sbox")
			print("           Possible options : QM Simple")
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

	genSkinny128Model(rMax,sboxModel,linModel,lastLin,smartLin,firstSbox)


