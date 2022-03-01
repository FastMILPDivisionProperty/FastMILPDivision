import os
localdir = os.path.dirname(__file__)
if localdir == "":
	localdir = "."
os.chdir(localdir)

load("MILP_function.sage")
load("divtrail.sage")
import sys, getopt

def genCLEFIAModel(rMax, sboxModel="Simple", linModel="Simple"):
	"""
	Generate the model for CLEFIA over #rMax rounds

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

	if linModel not in ["CX", "Simple"]:
		print("Unknown value for the linModel parameter, default to Simple")
		linModel = "Simple"

	blockSize = 128
	sboxSize = 8

	# S0 = [0x57,0x49,0xd1,0xc6,0x2f,0x33,0x74,0xfb,0x95,0x6d,0x82,0xea,0x0e,0xb0,0xa8,0x1c,
	# 	 0x28,0xd0,0x4b,0x92,0x5c,0xee,0x85,0xb1,0xc4,0x0a,0x76,0x3d,0x63,0xf9,0x17,0xaf,
	# 	 0xbf,0xa1,0x19,0x65,0xf7,0x7a,0x32,0x20,0x06,0xce,0xe4,0x83,0x9d,0x5b,0x4c,0xd8,
	# 	 0x42,0x5d,0x2e,0xe8,0xd4,0x9b,0x0f,0x13,0x3c,0x89,0x67,0xc0,0x71,0xaa,0xb6,0xf5,
	# 	 0xa4,0xbe,0xfd,0x8c,0x12,0x00,0x97,0xda,0x78,0xe1,0xcf,0x6b,0x39,0x43,0x55,0x26,
	# 	 0x30,0x98,0xcc,0xdd,0xeb,0x54,0xb3,0x8f,0x4e,0x16,0xfa,0x22,0xa5,0x77,0x09,0x61,
	# 	 0xd6,0x2a,0x53,0x37,0x45,0xc1,0x6c,0xae,0xef,0x70,0x08,0x99,0x8b,0x1d,0xf2,0xb4,
	# 	 0xe9,0xc7,0x9f,0x4a,0x31,0x25,0xfe,0x7c,0xd3,0xa2,0xbd,0x56,0x14,0x88,0x60,0x0b,
	# 	 0xcd,0xe2,0x34,0x50,0x9e,0xdc,0x11,0x05,0x2b,0xb7,0xa9,0x48,0xff,0x66,0x8a,0x73,
	# 	 0x03,0x75,0x86,0xf1,0x6a,0xa7,0x40,0xc2,0xb9,0x2c,0xdb,0x1f,0x58,0x94,0x3e,0xed,
	# 	 0xfc,0x1b,0xa0,0x04,0xb8,0x8d,0xe6,0x59,0x62,0x93,0x35,0x7e,0xca,0x21,0xdf,0x47,
	# 	 0x15,0xf3,0xba,0x7f,0xa6,0x69,0xc8,0x4d,0x87,0x3b,0x9c,0x01,0xe0,0xde,0x24,0x52,
	# 	 0x7b,0x0c,0x68,0x1e,0x80,0xb2,0x5a,0xe7,0xad,0xd5,0x23,0xf4,0x46,0x3f,0x91,0xc9,
	# 	 0x6e,0x84,0x72,0xbb,0x0d,0x18,0xd9,0x96,0xf0,0x5f,0x41,0xac,0x27,0xc5,0xe3,0x3a,
	# 	 0x81,0x6f,0x07,0xa3,0x79,0xf6,0x2d,0x38,0x1a,0x44,0x5e,0xb5,0xd2,0xec,0xcb,0x90,
	# 	 0x9a,0x36,0xe5,0x29,0xc3,0x4f,0xab,0x64,0x51,0xf8,0x10,0xd7,0xbc,0x02,0x7d,0x8e]

	# S1 = [0x6c,0xda,0xc3,0xe9,0x4e,0x9d,0x0a,0x3d,0xb8,0x36,0xb4,0x38,0x13,0x34,0x0c,0xd9,
	# 	 0xbf,0x74,0x94,0x8f,0xb7,0x9c,0xe5,0xdc,0x9e,0x07,0x49,0x4f,0x98,0x2c,0xb0,0x93,
	# 	 0x12,0xeb,0xcd,0xb3,0x92,0xe7,0x41,0x60,0xe3,0x21,0x27,0x3b,0xe6,0x19,0xd2,0x0e,
	# 	 0x91,0x11,0xc7,0x3f,0x2a,0x8e,0xa1,0xbc,0x2b,0xc8,0xc5,0x0f,0x5b,0xf3,0x87,0x8b,
	# 	 0xfb,0xf5,0xde,0x20,0xc6,0xa7,0x84,0xce,0xd8,0x65,0x51,0xc9,0xa4,0xef,0x43,0x53,
	# 	 0x25,0x5d,0x9b,0x31,0xe8,0x3e,0x0d,0xd7,0x80,0xff,0x69,0x8a,0xba,0x0b,0x73,0x5c,
	# 	 0x6e,0x54,0x15,0x62,0xf6,0x35,0x30,0x52,0xa3,0x16,0xd3,0x28,0x32,0xfa,0xaa,0x5e,
	# 	 0xcf,0xea,0xed,0x78,0x33,0x58,0x09,0x7b,0x63,0xc0,0xc1,0x46,0x1e,0xdf,0xa9,0x99,
	# 	 0x55,0x04,0xc4,0x86,0x39,0x77,0x82,0xec,0x40,0x18,0x90,0x97,0x59,0xdd,0x83,0x1f,
	# 	 0x9a,0x37,0x06,0x24,0x64,0x7c,0xa5,0x56,0x48,0x08,0x85,0xd0,0x61,0x26,0xca,0x6f,
	# 	 0x7e,0x6a,0xb6,0x71,0xa0,0x70,0x05,0xd1,0x45,0x8c,0x23,0x1c,0xf0,0xee,0x89,0xad,
	# 	 0x7a,0x4b,0xc2,0x2f,0xdb,0x5a,0x4d,0x76,0x67,0x17,0x2d,0xf4,0xcb,0xb1,0x4a,0xa8,
	# 	 0xb5,0x22,0x47,0x3a,0xd5,0x10,0x4c,0x72,0xcc,0x00,0xf9,0xe0,0xfd,0xe2,0xfe,0xae,
	# 	 0xf8,0x5f,0xab,0xf1,0x1b,0x42,0x81,0xd6,0xbe,0x44,0x29,0xa6,0x57,0xb9,0xaf,0xf2,
	# 	 0xd4,0x75,0x66,0xbb,0x68,0x9f,0x50,0x02,0x01,0x3c,0x7f,0x8d,0x1a,0x88,0xbd,0xac,
	# 	 0xf7,0xe4,0x79,0x96,0xa2,0xfc,0x6d,0xb2,0x6b,0x03,0xe1,0x2e,0x7d,0x14,0x95,0x1d]

	# M0_GF256 = Matrix([[1,2,4,6],[2,1,6,4],[4,6,1,2],[6,4,2,1]])
	# M0 = expandMatrix(M0_GF256, [1,0,1,1,1,0,0,0,1])
	M0 = matrix(GF(2), [[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1],
						[0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1],
						[0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,1,0,1,1,0,0,0,0,1,1],
						[0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,1,1,0,1,1,0,0,0,1,0],
						[0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,1,0,0,0,1,1,0,0,1,1,0,0,1,0],
						[0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,1,0,0,1],
						[0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0],
						[0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0],
						[0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0],
						[1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1],
						[0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,1,0,0,0,0,1,1,1,0,0,0,0,0,1,0],
						[0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,1,0,0,0,1,0,0,1,0,0,0,0,1,1],
						[0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,1,0,0,1,0,0,0,1,0,0,0,1,1],
						[0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,1,0,0,0,1,0,0,0,1],
						[0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0],
						[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0],
						[0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
						[0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
						[1,0,0,0,0,0,1,0,1,1,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1],
						[0,1,0,0,0,0,1,1,0,1,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1],
						[0,0,1,0,0,0,1,1,0,0,1,1,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1],
						[0,0,0,1,0,0,0,1,0,0,0,1,1,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0],
						[0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0],
						[0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0],
						[0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0],
						[1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
						[1,1,0,0,0,0,1,1,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0],
						[0,1,1,0,0,0,1,0,0,1,0,0,0,0,1,1,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,0],
						[0,0,1,1,0,0,1,0,0,0,1,0,0,0,1,1,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0],
						[0,0,0,1,1,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0],
						[0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0],
						[0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1]])

	# M1_GF256 = Matrix([[1,8,2,10],[8,1,10,2],[2,10,1,8],[10,2,8,1]])
	# M1 = expandMatrix(M1_GF256, [1,0,1,1,1,0,0,0,1])
	M1 = Matrix(GF(2), [[1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,1],
						[0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0],
						[0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,1,0,1,0,0,0,1,0,0],
						[0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,0,0,0,0,1,1,0,1,0,0,1,1,1],
						[0,0,0,0,1,0,0,0,0,1,0,0,0,1,1,1,0,0,0,1,0,0,0,1,0,1,0,1,0,1,1,0],
						[0,0,0,0,0,1,0,0,0,0,1,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,1,0,1,0,1,1],
						[0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,1,0,1],
						[0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,1,0],
						[0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1],
						[0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0],
						[0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,1],
						[1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,1,1,0,0,1,0,0,0,0,1],
						[0,1,0,0,0,1,1,1,0,0,0,0,1,0,0,0,0,1,0,1,0,1,1,0,0,0,0,1,0,0,0,1],
						[0,0,1,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,1,1,0,0,0,0,1,0,0,0],
						[0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,1,0,0,0,0,0,1,0,0],
						[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0],
						[0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
						[1,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
						[0,1,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,1],
						[0,0,1,0,0,0,0,1,1,0,1,0,0,1,1,1,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0],
						[0,0,0,1,0,0,0,1,0,1,0,1,0,1,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,1,1,1],
						[0,0,0,0,1,0,0,0,0,0,1,0,1,0,1,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,1,1],
						[0,0,0,0,0,1,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,1],
						[0,0,0,0,0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0],
						[0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0],
						[1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0],
						[0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0],
						[1,0,1,0,0,1,1,1,0,0,1,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0],
						[0,1,0,1,0,1,1,0,0,0,0,1,0,0,0,1,0,1,0,0,0,1,1,1,0,0,0,0,1,0,0,0],
						[0,0,1,0,1,0,1,1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,1,0,0],
						[0,0,0,1,0,1,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0],
						[0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1]])


	#Compute the ANF of S
	# (BPR, anfS) = SBOX_ANF(S)
	#Generate the inequalities for S
	# divTableS = SboxDivTrailTable(anfS)
	divTableS0 = load("divTableCLEFIAS0.sobj")
	divTableS1 = load("divTableCLEFIAS1.sobj")

	ineqS0 = None
	ineqS1 = None

	outLBS0 = None
	outLBS1 = None

	if sboxModel == "QM":
		ineqS0 = load("ineqCLEFIAS0Sbox.sobj")
		ineqS1 = load("ineqCLEFIAS1Sbox.sobj")

	elif sboxModel == "Simple":
		#Generate the lower bounds for the output weight of S
		outLBS0 = computeLBWeightSbox(divTableS0,sboxSize)
		outLBS1 = computeLBWeightSbox(divTableS1,sboxSize)

	modelName = "CLEFIA_"+str(rMax)+"r_"+sboxModel+"_"+linModel
	m = Model(modelName)
	x = [[[m.addVar(vtype=GRB.BINARY, name="x"+str(r)+"_"+str(i)+"_"+str(j))
		   for j in range(32)]
		   for i in range(4)]
		   for r in range(rMax+1)]

	#We also need copies of each even branch
	cx = [[[m.addVar(vtype=GRB.BINARY, name="cx"+str(r)+"_"+str(i)+"_"+str(j))
			 for j in range(32)] if i%2 == 0 else None
			 for i in range(4)]
			 for r in range(rMax)]

	#Variables for the output of Sboxes
	y = [[[m.addVar(vtype=GRB.BINARY, name="y"+str(r)+"_"+str(i)+"_"+str(j))
			 for j in range(32)] if i%2 == 0 else None
			 for i in range(4)]
			 for r in range(rMax)]

	#Variables for the output of Fi
	z = [[[m.addVar(vtype=GRB.BINARY, name="z"+str(r)+"_"+str(i)+"_"+str(j))
			 for j in range(32)] if i%2 == 0 else None
			 for i in range(4)]
			 for r in range(rMax)]

	m.update()
	#Constraints

	#Copy and XOR constraints
	#No alternatives for these ones
	for r in range(rMax):
		#(cx0,x3) = copy(x0)
		#x0 = x1 xor z0
		for i in range(32):
			addCopyConstr(m, x[r][0][i], [cx[r][0][i], x[r+1][3][i]])
			addXorConstr(m, [x[r][1][i], z[r][0][i]], x[r+1][0][i])

		#(cx2,x1) = copy(x2)
		#x2 = x3 xor z2
		for i in range(32):
			addCopyConstr(m, x[r][2][i], [cx[r][2][i], x[r+1][1][i]])
			addXorConstr(m, [x[r][3][i], z[r][2][i]], x[r+1][2][i])

	#Modelize the Sbox layers
	#y[r] = S(cx[r])
	if sboxModel == "QM":
		for r in range(rMax):
			#y[0] = (S0,S1,S0,S1)(cx[0])
			for i in range(4):
				invar  = [cx[r][0][sboxSize*i + j] for j in range(sboxSize)]
				outvar = [y[r][0][sboxSize*i + j] for j in range(sboxSize)]
				if(i%2 == 0):
					addSboxConstr(m, ineqS0, invar, outvar)
				else:
					addSboxConstr(m, ineqS1, invar, outvar)

			#y[2] = (S1,S0,S1,S0)(cx[2])
			for i in range(4):
				invar  = [cx[r][2][sboxSize*i + j] for j in range(sboxSize)]
				outvar = [y[r][2][sboxSize*i + j] for j in range(sboxSize)]
				if(i%2 == 0):
					addSboxConstr(m, ineqS1, invar, outvar)
				else:
					addSboxConstr(m, ineqS0, invar, outvar)

	elif sboxModel == "Simple":
		for r in range(rMax):
			#y[0] = (S0,S1,S0,S1)(cx[0])
			for i in range(4):
				invar  = [cx[r][0][sboxSize*i + j] for j in range(sboxSize)]
				outvar = [y[r][0][sboxSize*i + j] for j in range(sboxSize)]
				if(i%2 == 0):
					addSimplifiedSboxConstr(m,invar, outvar, outLBS0, boundSuffix="s0"+str(r)+"_"+str(i))
				else:
					addSimplifiedSboxConstr(m,invar, outvar, outLBS1, boundSuffix="s0"+str(r)+"_"+str(i))

			#y[2] = (S1,S0,S1,S0)(cx[2])
			for i in range(4):
				invar  = [cx[r][2][sboxSize*i + j] for j in range(sboxSize)]
				outvar = [y[r][2][sboxSize*i + j] for j in range(sboxSize)]
				if(i%2 == 0):
					addSimplifiedSboxConstr(m,invar, outvar, outLBS1, boundSuffix="s2"+str(r)+"_"+str(i))
				else:
					addSimplifiedSboxConstr(m,invar, outvar, outLBS0, boundSuffix="s2"+str(r)+"_"+str(i))

	#Linear constraints
	if linModel == "CX":
		#Modelize Mi with the Copy+XOR technique
		#z0 = M0*y0
		#z2 = M1*y2
		for r in range(rMax):
			addLinearConstrCX(m, y[r][0], z[r][0], M0)
			addLinearConstrCX(m, y[r][2], z[r][2], M1)

	elif linModel == "Simple":
		for r in range(rMax):
			addSimplifiedLinearConstr(m, y[r][0], z[r][0])
			addSimplifiedLinearConstr(m, y[r][2], z[r][2])

	m.update()
	m.write("./../models/"+modelName+".mps")

if __name__ == "__main__":
	rMax = 1
	sboxModel = "Simple"
	linModel = "Simple"
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hr:s:m:")
	except getopt.GetoptError:
		print("genARIA.sage -r <rMax> -s <sboxModel> -m <linModel>")
		print("<rMax> is the number of rounds")
		print("<sboxModel> is the technique used for the modelization of the Sbox")
		print("           Possible options : QM Simple")
		print("<linModel> is the technique used for the modelization of MC")
		print("           Possible options : CX Simple")
		sys.exit(2)
	for opt, arg in opts:
		if opt == "-h":
			print("genARIA.sage -r <rMax> -s <sboxModel> -m <linModel>")
			print("<rMax> is the number of rounds")
			print("<sboxModel> is the technique used for the modelization of the Sbox")
			print("           Possible options : QM Simple")
			print("<linModel> is the technique used for the modelization of MC")
			print("           Possible options : CX Simple")
			sys.exit()
		elif opt == "-r":
			rMax = int(arg)
		elif opt == "-s":
			sboxModel = arg
		elif opt == "-m":
			linModel = arg

	genCLEFIAModel(rMax,sboxModel,linModel)