import os
localdir = os.path.dirname(__file__)
if localdir == "":
	localdir = "."
os.chdir(localdir)

load("MILP_function.sage")
load("divtrail.sage")
import sys, getopt

def genAESModel(rMax, sboxModel="Simple", linModel="Simple", lastLin=True,firstSbox=True):
	"""
	Generate the model for AES over #rMax rounds
	#lastLin defines if the linear layer is applied on the last round
	#firstSbox defines if the first Sbox layer is applied

	The #sboxModel parameter defines how the Sbox is modelized :
	- "QM" use the Quin-McCluskey algorithm as in Abdelkhalek,Sasaki,Todo,Tolba,Youssef
	- "Simple" use the simplified constraint with PWL

	The #linModel parameter defines how the linear layer is modelized :
	- "CX" modelize the linear layer using the classical copy+xor technique
	- "Simple" modelize the linear layer with the simplified constraint w(x) = w(y)
	"""

	if sboxModel not in ["Simple","QM"]:
		print("Unknown value for the sboxModel parameter, default to Simple")
		sboxModel = "Simple"

	if linModel not in ["CX", "Simple", "Simplev2"]:
		print("Unknown value for the linModel parameter, default to Simple")
		linModel = "Simple"

	blockSize = 128
	nbSbox = 16
	sboxSize = 8

	# S = [0x63,0x7c,0x77,0x7b,0xf2,0x6b,0x6f,0xc5,0x30,0x01,0x67,0x2b,0xfe,0xd7,0xab,0x76,
	# 	0xca,0x82,0xc9,0x7d,0xfa,0x59,0x47,0xf0,0xad,0xd4,0xa2,0xaf,0x9c,0xa4,0x72,0xc0,
	# 	0xb7,0xfd,0x93,0x26,0x36,0x3f,0xf7,0xcc,0x34,0xa5,0xe5,0xf1,0x71,0xd8,0x31,0x15,
	# 	0x04,0xc7,0x23,0xc3,0x18,0x96,0x05,0x9a,0x07,0x12,0x80,0xe2,0xeb,0x27,0xb2,0x75,
	# 	0x09,0x83,0x2c,0x1a,0x1b,0x6e,0x5a,0xa0,0x52,0x3b,0xd6,0xb3,0x29,0xe3,0x2f,0x84,
	# 	0x53,0xd1,0x00,0xed,0x20,0xfc,0xb1,0x5b,0x6a,0xcb,0xbe,0x39,0x4a,0x4c,0x58,0xcf,
	# 	0xd0,0xef,0xaa,0xfb,0x43,0x4d,0x33,0x85,0x45,0xf9,0x02,0x7f,0x50,0x3c,0x9f,0xa8,
	# 	0x51,0xa3,0x40,0x8f,0x92,0x9d,0x38,0xf5,0xbc,0xb6,0xda,0x21,0x10,0xff,0xf3,0xd2,
	# 	0xcd,0x0c,0x13,0xec,0x5f,0x97,0x44,0x17,0xc4,0xa7,0x7e,0x3d,0x64,0x5d,0x19,0x73,
	# 	0x60,0x81,0x4f,0xdc,0x22,0x2a,0x90,0x88,0x46,0xee,0xb8,0x14,0xde,0x5e,0x0b,0xdb,
	# 	0xe0,0x32,0x3a,0x0a,0x49,0x06,0x24,0x5c,0xc2,0xd3,0xac,0x62,0x91,0x95,0xe4,0x79,
	# 	0xe7,0xc8,0x37,0x6d,0x8d,0xd5,0x4e,0xa9,0x6c,0x56,0xf4,0xea,0x65,0x7a,0xae,0x08,
	# 	0xba,0x78,0x25,0x2e,0x1c,0xa6,0xb4,0xc6,0xe8,0xdd,0x74,0x1f,0x4b,0xbd,0x8b,0x8a,
	# 	0x70,0x3e,0xb5,0x66,0x48,0x03,0xf6,0x0e,0x61,0x35,0x57,0xb9,0x86,0xc1,0x1d,0x9e,
	# 	0xe1,0xf8,0x98,0x11,0x69,0xd9,0x8e,0x94,0x9b,0x1e,0x87,0xe9,0xce,0x55,0x28,0xdf,
	# 	0x8c,0xa1,0x89,0x0d,0xbf,0xe6,0x42,0x68,0x41,0x99,0x2d,0x0f,0xb0,0x54,0xbb,0x16]
	# P = [0,13,10,7,4,1,14,11,8,5,2,15,12,9,6,3]
	#Expand P at a bit level
	# Pbit = [0 for i in range(blockSize)]
	# for isbox in range(nbSbox):
	# 	for j in range(sboxSize):
	# 		Pbit[isbox*sboxSize+j] = P[isbox]*sboxSize+j
	Pbit = [0,1,2,3,4,5,6,7,104,105,106,107,108,109,110,111,80,81,82,83,84,85,86,87,56,57,58,59,60,61,62,63,32,33,34,35,36,37,38,39,8,9,10,11,12,13,14,15,112,113,114,115,116,117,118,119,88,89,90,91,92,93,94,95,64,65,66,67,68,69,70,71,40,41,42,43,44,45,46,47,16,17,18,19,20,21,22,23,120,121,122,123,124,125,126,127,96,97,98,99,100,101,102,103,72,73,74,75,76,77,78,79,48,49,50,51,52,53,54,55,24,25,26,27,28,29,30,31]

	fullM = matrix(GF(2), [[0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
						   [1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
						   [0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0],
						   [0,0,1,0,0,0,0,1,0,0,1,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0],
						   [0,0,0,1,0,0,0,1,0,0,0,1,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0],
						   [0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0],
						   [0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0],
						   [0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1],
						   [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0],
						   [0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0],
						   [0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0],
						   [0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,1,0,0,0,1,0,0,0,1,0,0,0,0],
						   [0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,1,0,0,1,0,0,0,0,1,0,0,0],
						   [0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0],
						   [0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0],
						   [0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1],
						   [1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1],
						   [0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,1],
						   [0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0],
						   [0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,1,0,0,0,1],
						   [0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,1,0,0,1],
						   [0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0],
						   [0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0],
						   [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1],
						   [1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
						   [1,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1],
						   [0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0],
						   [0,0,1,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1],
						   [0,0,0,1,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1],
						   [0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0],
						   [0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0],
						   [0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0]])
	invFullM = fullM.inverse()

	

	#--- Precomputations ---

	#Compute the ANF of S
	# (BPR, anfS) = SBOX_ANF(S)
	#Generate the inequalities for S
	# divTableS = SboxDivTrailTable(anfS)
	divTableS = load("divTableAESSbox.sobj")
	ineqS = None
	outLB = None

	if sboxModel == "QM":
		ineqS = load("ineqAESSbox.sobj")

	if sboxModel == "Simple":
		#Generate the lower bounds for the output weight of S
		outLB = computeLBWeightSbox(divTableS,sboxSize)


	#--- Model generation ---

	#Create the model
	modelName = "AES_"+str(rMax)+"r_"+sboxModel+"_"+linModel
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
	if linModel == "CX":
		#Modelize MC with the Copy+XOR technique
		for r in range(roundLastLin):
			for col in range(4):
				invar = [z[r][32*col + i] for i in range(32)]
				outvar = [x[r+1][32*col + i] for i in range(32)]
				addLinearConstrCX(m, invar, outvar, fullM)

	elif linModel == "Simple":
		for r in range(roundLastLin):
			for col in range(4):
				invar = [z[r][32*col + i] for i in range(32)]
				outvar = [x[r+1][32*col + i] for i in range(32)]
				addSimplifiedLinearConstr(m, invar, outvar)
	elif linModel == "Simplev2":
		for r in range(roundLastLin):
			for col in range(4):
				invar = [z[r][32*col + i] for i in range(32)]
				outvar = [x[r+1][32*col + i] for i in range(32)]
				addSimplifiedLinearConstr_v2(m, invar, outvar, fullM, invFullM)

	m.update()
	m.write("./../models/"+modelName+".mps")


if __name__ == "__main__":
	rMax = 1
	sboxModel = "Simple"
	linModel = "Simple"
	lastLin = True
	firstSbox = True
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hr:s:m:", ["lastLin=","firstSbox="])
	except getopt.GetoptError:
		print("genAES.sage -r <rMax> -s <sboxModel> -m <linModel> --lastLin <lastLin> --firstSbox <firstSbox>")
		print("<rMax> is the number of rounds")
		print("<sboxModel> is the technique used for the modelization of the Sbox")
		print("           Possible options : QM Simple")
		print("<linModel> is the technique used for the modelization of MC")
		print("           Possible options : CX Simple Simplev2")
		print("<lastLin> set if the last linear layer is applied")
		print("          Possible options : True False")
		print("<firstSbox> set if the first sbox layer is applied")
		print("          Possible options : True False")
		sys.exit(2)
	for opt, arg in opts:
		if opt == "-h":
			print("genAES.sage -r <rMax> -s <sboxModel> -m <linModel> --lastLin <lastLin> --firstSbox <firstSbox>")
			print("<rMax> is the number of rounds")
			print("<sboxModel> is the technique used for the modelization of the Sbox")
			print("           Possible options : QM Simple")
			print("<linModel> is the technique used for the modelization of MC")
			print("           Possible options : CX Simple Simplev2")
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

	genAESModel(rMax,sboxModel,linModel,lastLin,firstSbox)




