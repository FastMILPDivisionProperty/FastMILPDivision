import os
localdir = os.path.dirname(__file__)
if localdir == "":
	localdir = "."
os.chdir(localdir)

load("MILP_function.sage")
load("divtrail.sage")
import sys, getopt

def genHIGHTmodel(rMax, arxModel="ARX", linModel="Lbox"):
	"""
	Generate the model for HIGHT over #rMax rounds

	The #arxModel parameter defines how the mod addition is modelized :
	- "ARX" use the technique from Sun,Wang,Liu,Wang

	The #linModel parameter defines how the linear layer is modelized :
	- "QM" use the Quin-McCluskey algorithm as in Abdelkhalek,Sasaki,Todo,Tolba,Youssef
	- "ZR" modelize the linear layer using the technique from Zhang and Rijmen
	- "CX" modelize the linear layer using the classical copy+xor technique
	- "Simple" modelize the linear layer with the simplified constraint w(x) = w(y)
	"""

	if arxModel not in ["ARX"]:
		print("Unknown value for the arxModel parameter, default to ARX")
		arxModel = "ARX"

	if linModel not in ["QM", "ZR", "CX", "Simple"]:
		print("Unknown value for the linModel parameter, default to QM")
		linModel = "QM"

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
	ineqF0 = None
	ineqF1 = None
	if linModel == "QM":
		ineqF0 = load("ineqF0.sobj")
		ineqF1 = load("ineqF1.sobj")

	modelName = "HIGHT_"+str(rMax)+"r_"+arxModel+"_"+linModel
	m = Model(modelName)
	x = [[[m.addVar(vtype=GRB.BINARY, name="x"+str(r)+"_"+str(i)+"_"+str(j))
		   for j in range(8)]
		   for i in range(8)]
		   for r in range(rMax+1)]

	#We also need copies of each even branch
	cx = [[[m.addVar(vtype=GRB.BINARY, name="cx"+str(r)+"_"+str(i)+"_"+str(j))
			 for j in range(8)] if i%2 == 0 else None
			 for i in range(8)]
			 for r in range(rMax)]

	#Variables for the output of Fi
	y = [[[m.addVar(vtype=GRB.BINARY, name="y"+str(r)+"_"+str(i)+"_"+str(j))
			 for j in range(8)] if i%2 == 0 else None
			 for i in range(8)]
			 for r in range(rMax)]
	#Variables for the output of mod add with key
	z = [[[m.addVar(vtype=GRB.BINARY, name="z"+str(r)+"_"+str(i)+"_"+str(j))
			 for j in range(8)] if i == 2 or i == 6 else None
			 for i in range(8)]
			 for r in range(rMax)]

	m.update()
	#Constraints
	#Alternate between F1+modAdd and F0+xor for each pair of block

	#Copy and XOR constraints
	#No alternatives for these ones
	for r in range(rMax):
		#(cx0,x1) = copy(x0)
		#y0 = F1(cx0)
		#x2 = x1 + y0 mod 256
		for i in range(8):
			addCopyConstr(m, x[r][0][i], [cx[r][0][i], x[r+1][1][i]])
		
		#(cx2,x3) = copy(x2)
		#y2 = F0(cx2)
		#x4 = x3 xor y2
		for i in range(8):
			addCopyConstr(m, x[r][2][i], [cx[r][2][i], x[r+1][3][i]])
			addXorConstr(m, [x[r][3][i], z[r][2][i]], x[r+1][4][i])

		#(cx4,x5) = copy(x4)
		#y4 = F1(cx4)
		#x6 = x5 + y4 mod 256
		for i in range(8):
			addCopyConstr(m, x[r][4][i], [cx[r][4][i], x[r+1][5][i]])

		#(cx6,x7) = copy(x6)
		#y6 = F0(cx6)
		#x0 = x7 xor y6
		for i in range(8):
			addCopyConstr(m, x[r][6][i], [cx[r][6][i], x[r+1][7][i]])
			addXorConstr(m, [x[r][7][i], z[r][6][i]], x[r+1][0][i])

	#ARX Constraints
	if arxModel == "ARX":
		for r in range(rMax):
			addModAddConstr(m, x[r][1], y[r][0], x[r+1][2], cPrefix="c"+str(r)+"_0")
			addModAddConstr(m, x[r][5], y[r][4], x[r+1][6], cPrefix="c"+str(r)+"_4")
			addModAddConstantConstr(m, y[r][2], z[r][2], cPrefix="cst"+str(r)+"_2")
			addModAddConstantConstr(m, y[r][6], z[r][6], cPrefix="cst"+str(r)+"_6")

	#Linear constraints
	if linModel == "CX":
		#Modelize Fi with the Copy+XOR technique
		for r in range(rMax):
			addLinearConstrCX(m, cx[r][0], y[r][0], F1)
			addLinearConstrCX(m, cx[r][2], y[r][2], F0)
			addLinearConstrCX(m, cx[r][4], y[r][4], F1)
			addLinearConstrCX(m, cx[r][6], y[r][6], F0)

	elif linModel == "ZR":
		#Modelize Fi with the Zhang-Rijmen technique
		for r in range(rMax):
			addLinearConstrZR(m, cx[r][0], y[r][0], F1, 1)
			addLinearConstrZR(m, cx[r][2], y[r][2], F0, 1)
			addLinearConstrZR(m, cx[r][4], y[r][4], F1, 1)
			addLinearConstrZR(m, cx[r][6], y[r][6], F0, 1)

	elif linModel == "QM":
		#Modelize Fi as linear sboxes
		for r in range(rMax):
			addSboxConstr(m, ineqF1, cx[r][0], y[r][0])
			addSboxConstr(m, ineqF0, cx[r][2], y[r][2])
			addSboxConstr(m, ineqF1, cx[r][4], y[r][4])
			addSboxConstr(m, ineqF0, cx[r][6], y[r][6])

	elif linModel == "Simple":
		for r in range(rMax):
			addSimplifiedLinearConstr(m, cx[r][0], y[r][0])
			addSimplifiedLinearConstr(m, cx[r][2], y[r][2])
			addSimplifiedLinearConstr(m, cx[r][4], y[r][4])
			addSimplifiedLinearConstr(m, cx[r][6], y[r][6])

	m.update()
	m.write("./../models/"+modelName+".mps")

	# return (m,x,cx,y)


if __name__ == "__main__":
	rMax = 1
	arxModel = "ARX"
	linModel = "QM"
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hr:a:m:")
	except getopt.GetoptError:
		print("genHIGHT.sage -r <rMax> -a <arxModel> -m <linModel>")
		print("<rMax> is the number of rounds")
		print("<arxModel> is the technique used for the modelization of the modular addition")
		print("           Possible options : ARX")
		print("<linModel> is the technique used for the modelization of MC")
		print("           Possible options : QM ZR CX Simple")
		sys.exit(2)
	for opt, arg in opts:
		if opt == "-h":
			print("genHIGHT.sage -r <rMax> -a <arxModel> -m <linModel>")
			print("<rMax> is the number of rounds")
			print("<arxModel> is the technique used for the modelization of the modular addition")
			print("           Possible options : ARX")
			print("<linModel> is the technique used for the modelization of MC")
			print("           Possible options : QM ZR CX Simple")
			sys.exit()
		elif opt == "-r":
			rMax = int(arg)
		elif opt == "-a":
			arxModel = arg
		elif opt == "-m":
			linModel = arg

	genHIGHTmodel(rMax,arxModel,linModel)
