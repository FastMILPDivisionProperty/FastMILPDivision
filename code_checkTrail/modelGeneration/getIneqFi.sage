reset()
for Fi in range(2):
	f = open("uv_F"+str(Fi)+".txt")
	print(("F" + str(Fi)))
	ineqF = []
	for l in f.readlines():
		L = l.split(" | ")
		L[1] = L[1].rstrip("\n")
		L = [list(L[0]), list(L[1])]
		u = L[0]
		v = L[1]
		for i in range(len(u)):
			if u[i] == '1':
				v[i] = '*'

		ineq = [0 for _ in range(len(v))]
		cst = -1
		for i in range(len(v)):
			if v[i] == '1':
				ineq[i] = -1
				cst += 1
			elif v[i] == '0':
				ineq[i] = 1
		ineq.append(cst)
		ineq.append(">=")
		
		print((l.rstrip("\n") + " -> " + "".join(v) + " -> " + " + ".join(map(str,ineq)) + " 0"))
		ineqF.append(ineq)

	save(ineqF, "ineqF"+str(Fi))


	f.close()