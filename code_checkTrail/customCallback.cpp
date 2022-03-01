#include "customCallback.hpp"

#define DISPLAY_CB_LOG false
using namespace std;
  using namespace std::chrono;

typedef unsigned int uint;

CustomCallback::CustomCallback(vector<BlockVarsTable> const & xvarsTable,
 							   vector<BlockVarsIneq> const & xvarsIneq,
 							   vector<BlockVarsMatrix> const & xvarsMatrix,
 							   uint64_t const xmaxNumberConstr) :
							   varsTable(xvarsTable),
							   varsIneq(xvarsIneq),
							   varsMatrix(xvarsMatrix),
							   maxNumberConstr(xmaxNumberConstr),
							   ctrGlobalConstr(0),
							   ctrSol(0),
							   // ctrNonInvMinor(0),
							   // histoNonInvMinors()
							   ctrConstrSbox(0),
							   ctrConstrLin(0),
							   ctrConstrSSB(0)
{}

void CustomCallback::callback(){
	try {
        if (where == GRB_CB_MIPSOL) { //If the solver found a solution
		    ctrSol++;
		    uint64_t ctrconstr = 0;

		    //Check the variables with ineq checks (if any)
	    	bool reachedMaxConstr = false;
	    	for(auto const & iot : varsIneq){
	    		//Extract the input values
	    		//vectors for fast access afterward
	    		uint64_t input_size = iot.input.size();
	    		vector<int> inputval(input_size);
	    		for(uint i = 0; i < input_size; i++)
	    			inputval[i] = int(round(getSolution(iot.input[i])));

	    		uint64_t output_size = iot.output.size();
	    		vector<int> outputval(output_size);
	    		for(uint i = 0; i < output_size; i++)
	    			outputval[i] = int(round(getSolution(iot.output[i])));

	    		//Evaluate each ineq
	    		for(auto const eq : iot.ineq){
	    			int lhs = 0;
	    			for(uint i = 0; i < input_size; i++)
	    				lhs += eq[i]*inputval[i];

	    			for(uint i = 0; i < output_size; i++)
	    				lhs += eq[input_size+i]*outputval[i];

	    			//Constant term is eq[-2]
	    			lhs += eq[eq.size()-2];

	    			//eq[-1] indicates the operator
	    			bool ineqCheck = true;
	    			if(eq[eq.size()-1] == 0)
	    				ineqCheck = (lhs == 0);
	    			else
	    				ineqCheck = (lhs >= 0);

	    			//If the ineq is not verified, add a constraint
	    			if(!ineqCheck){
	    				GRBLinExpr cutExpr = 0;
	    				for(uint i = 0; i < input_size; i++)
	    					cutExpr += eq[i]*iot.input[i];

	    				for(uint i = 0; i < output_size; i++)
	    					cutExpr += eq[input_size+i]*iot.output[i];

	    				//constant term
	    				cutExpr += eq[eq.size()-2];
	    				if(eq[eq.size()-1] == 0)
	    					addLazy(cutExpr == 0);
	    				else
	    					addLazy(cutExpr >= 0);

	    				ctrconstr++;
	    				ctrGlobalConstr++;
	    				if(iot.ctype == ConstrType::Sbox)
	    					ctrConstrSbox++;
	    				else if(iot.ctype == ConstrType::Lin)
	    					ctrConstrLin++;
	    				else if(iot.ctype == ConstrType::SSB)
	    					ctrConstrSSB++;
			    		if(ctrconstr >= maxNumberConstr){
			    			reachedMaxConstr = true;
			    			break;
			    		}
	    			}
	    		}
	    		if(reachedMaxConstr)
	    			break;
	    	}

		    if(ctrconstr < maxNumberConstr){
		    	//Check the variables with table checks (if any)
			    for(auto const & iot : varsTable){
			    	//Build the int value of the input
			    	uint32_t inputval = 0;
			    	uint32_t mask = 1;
			    	for(uint i = 0; i < iot.input.size(); i++){
			    		if(int(round(getSolution(iot.input[i]))))
			    			inputval |= mask;
			    		mask <<= 1;
			    	}

			    	//Build the int value of the output
			    	uint32_t outputval = 0;
			    	mask = 1;
			    	for(uint i = 0; i < iot.output.size(); i++){
			    		if(int(round(getSolution(iot.output[i]))))
			    			outputval |= mask;
			    		mask <<= 1;
			    	}

			    	//Check if it belongs to the table
			    	//If not, add a constraint to remove it
			    	// if(iot.table[inputval][outputval] == 0){
			    	auto const & table_input = iot.table[inputval];
			    	if(!binary_search(table_input.begin(), table_input.end(), outputval)){

			    		GRBLinExpr cutExpr(0);
			    		for(uint i = 0; i < iot.input.size(); i++){
			    			if(inputval & 1) cutExpr += (1 - iot.input[i]);
			    			else cutExpr += iot.input[i];
			    			inputval >>= 1;
			    		}
			    		for(uint i = 0; i < iot.output.size(); i++){
			    			if(outputval & 1) cutExpr += (1 - iot.output[i]);
			    			else cutExpr += iot.output[i]; 
			    			outputval >>= 1;
			    		}
			    		addLazy(cutExpr >= 1);
			    		ctrconstr++;
			    		ctrGlobalConstr++;
			    		if(iot.ctype == ConstrType::Sbox)
	    					ctrConstrSbox++;
	    				else if(iot.ctype == ConstrType::Lin)
	    					ctrConstrLin++;
	    				else if(iot.ctype == ConstrType::SSB)
	    					ctrConstrSSB++;
			    		if(ctrconstr >= maxNumberConstr)
			    			break;
			    	}
			    }
		    }

		    if(ctrconstr < maxNumberConstr){
		    	//Check the variables with matrix checks (if any)
		    	for(auto const & iot : varsMatrix){
		    		//Extract the input values
					//vectors for fast access afterward
					//Also build the row/column index vectors
					vector<uint> indexRow, indexCol;
					vector<uint> inputval(iot.input.size());
					for(uint i = 0; i < iot.input.size(); i++){
						inputval[i] = uint(round(getSolution(iot.input[i])));
						if(inputval[i] == 1)
							indexCol.emplace_back(i);
					}
					vector<uint> outputval(iot.output.size());
					for(uint i = 0; i < iot.output.size(); i++){
						outputval[i] = uint(round(getSolution(iot.output[i])));
						if(outputval[i] == 1)
							indexRow.emplace_back(i);
					}

					// //Check the corresponding minor
					// //If not invertible, remove this solution
					// if(!iot.mat.isInvertibleMinor(indexRow, indexCol)){
					// 	// ctrNonInvMinor++;
					// 	// histoNonInvMinors[indexRow.size()]++;

					// 	GRBLinExpr cutExpr(0);
					// 	for(uint i = 0; i < iot.input.size(); i++){
					// 		if(inputval[i] == 1) cutExpr += (1 - iot.input[i]);
					// 		else cutExpr += iot.input[i];
					// 	}
					// 	for(uint i = 0; i < iot.output.size(); i++){
					// 		if(outputval[i] == 1) cutExpr += (1 - iot.output[i]);
					// 		else cutExpr += iot.output[i];
					// 	}
					// 	addLazy(cutExpr >= 1);
					// 	ctrconstr++;
					// 	ctrGlobalConstr++;
					// 	if(ctrconstr >= maxNumberConstr)
					// 		break;
					// }
					// cout << "Minor of size " << indexRow.size() << " x " << indexCol.size() << endl;
					// cout << "indexRow : ";
					// for(auto const & tmpir : indexRow) cout << tmpir << ",";
					// cout << endl << "indexCol : ";
					// for(auto const & tmpir : indexCol) cout << tmpir << ",";

					auto checkIndex = iot.mat.isInvertibleMinorWithAugment(indexRow, indexCol);
					if(!checkIndex.first){ //Minor not invertible, invalid trail
						if(checkIndex.second.size() == 0){
							//With a properly made model, this shouldn't happen
							//Still implement out of safety, it's basically when the minor is not a square matrix

							GRBLinExpr cutExpr(0);
							for(uint i = 0; i < iot.input.size(); i++){
								if(inputval[i] == 1) cutExpr += (1 - iot.input[i]);
								else cutExpr += iot.input[i];
							}
							for(uint i = 0; i < iot.output.size(); i++){
								if(outputval[i] == 1) cutExpr += (1 - iot.output[i]);
								else cutExpr += iot.output[i];
							}
							addLazy(cutExpr >= 1);
							ctrconstr++;
							ctrGlobalConstr++;
							ctrConstrLin++;
							if(ctrconstr >= maxNumberConstr)
								break;
						}
						else{
							//Else we add the constraints resulting from the linear combinations leading to a null row in the minor
							// v1 + ... + vk  - (k - 1) <=  u1 + ... + un

							for(auto const & rc : checkIndex.second){
								GRBLinExpr cutExpr(0);
								for(auto const i : rc.first){
									//Indexes for rows, i.e. output vars
									cutExpr += iot.output[i];
								}
								for(auto const i : rc.second){
									//Indexes for columns, i.e. input vars
									cutExpr -= iot.input[i];
								}
								addLazy(cutExpr <= (rc.first.size() - 1));
								ctrconstr++;
								ctrGlobalConstr++;
								ctrConstrLin++;
								if(ctrconstr >= maxNumberConstr)
									break;
							}
							if(ctrconstr >= maxNumberConstr)
								break;
						}
					}
		    	}
		    }

		    if(DISPLAY_CB_LOG){
		    	cout << "Solution number : " << ctrSol << ", " << ctrGlobalConstr << " constraints added" << endl;
		    	// cout << ctrNonInvMinor << " non-invertible minors so far, distribution : " << endl;
		    	// for(auto const & sn : histoNonInvMinors)
		    	// 	cout << sn.first << " : " << sn.second << endl;
		    	// cout << endl;
		    }
		}
	} catch (GRBException e) {
    cout << "Error number: " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
	} catch (...) {
	cout << "Error during callback" << endl;
	}
}

tuple<bool,uint64_t,uint64_t> solveWithoutCallback(GRBModel & m,
						  vector<BlockVarsTable> const & varsTable,
						  vector<BlockVarsIneq> const & varsIneq,
						  vector<BlockVarsMatrix> const & varsMatrix,
						  uint64_t const maxNumberConstr){
	/*
		Apply the same strategy as in the callback, but without the callback
	*/

	uint64_t ctrGlobalConstr = 0;
	uint64_t ctrSol = 0;

	while(true){
		m.optimize();
		if(m.get(GRB_IntAttr_Status) == GRB_OPTIMAL){
			//Verify if the solution pass all the checks
			ctrSol++;
			uint64_t ctrconstr = 0;

		    //Check the variables with ineq checks (if any)
	    	bool reachedMaxConstr = false;
	    	for(auto const & iot : varsIneq){
	    		//Extract the input values
	    		//vectors for fast access afterward
	    		vector<int> inputval(iot.input.size());
	    		for(uint i = 0; i < iot.input.size(); i++)
	    			inputval[i] = int(round(iot.input[i].get(GRB_DoubleAttr_X)));

	    		vector<uint> outputval(iot.output.size());
	    		for(uint i = 0; i < iot.output.size(); i++)
	    			outputval[i] = int(round(iot.output[i].get(GRB_DoubleAttr_X)));

	    		//Evaluate each ineq
	    		for(auto const eq : iot.ineq){
	    			int lhs = 0;
	    			for(uint i = 0; i < inputval.size(); i++)
	    				lhs += eq[i]*inputval[i];

	    			for(uint i = 0; i < outputval.size(); i++)
	    				lhs += eq[inputval.size()+i]*outputval[i];

	    			//Constant term is eq[-2]
	    			lhs += eq[eq.size()-2];

	    			//eq[-1] indicates the operator
	    			bool ineqCheck = true;
	    			if(eq[eq.size()-1] == 0)
	    				ineqCheck = (lhs == 0);
	    			else
	    				ineqCheck = (lhs >= 0);

	    			//If the ineq is not verified, add a constraint
	    			if(!ineqCheck){

	    				GRBLinExpr cutExpr = 0;
	    				for(uint i = 0; i < iot.input.size(); i++)
	    					cutExpr += eq[i]*iot.input[i];

	    				for(uint i =0; i < iot.output.size(); i++)
	    					cutExpr += eq[iot.input.size()+i]*iot.output[i];

	    				//constant term
	    				cutExpr += eq[eq.size()-2];
	    				if(eq[eq.size()-1] == 0)
	    					m.addConstr(cutExpr == 0);
	    				else
	    					m.addConstr(cutExpr >= 0);

	    				ctrconstr++;
	    				ctrGlobalConstr++;
			    		if(ctrconstr >= maxNumberConstr){
			    			reachedMaxConstr = true;
			    			break;
			    		}
	    			}
	    		}
	    		if(reachedMaxConstr)
	    			break;
	    	}

		    if(ctrconstr < maxNumberConstr){
		    	//Check the variables with table checks (if any)
			    for(auto const & iot : varsTable){
			    	//Build the int value of the input
			    	uint64_t inputval = 0;
			    	uint64_t mask = 1;
			    	for(uint i = 0; i < iot.input.size(); i++){
			    		if(int(round(iot.input[i].get(GRB_DoubleAttr_X))))
			    			inputval |= mask;
			    		mask <<= 1;
			    	}

			    	//Build the int value of the output
			    	uint64_t outputval = 0;
			    	mask = 1;
			    	for(uint i = 0; i < iot.output.size(); i++){
			    		if(int(round(iot.output[i].get(GRB_DoubleAttr_X))))
			    			outputval |= mask;
			    		mask <<= 1;
			    	}

			    	//Check if it belongs to the table
			    	//If not, add a constraint to remove it
			    	// if(iot.table[inputval][outputval] == 0){
			    	auto const & table_input = iot.table[inputval];
			    	if(!binary_search(table_input.begin(), table_input.end(), outputval)){

			    		GRBLinExpr cutExpr(0);
			    		for(uint i = 0; i < iot.input.size(); i++){
			    			if(inputval & 1) cutExpr += (1 - iot.input[i]);
			    			else cutExpr += iot.input[i];
			    			inputval >>= 1;
			    		}
			    		for(uint i = 0; i < iot.output.size(); i++){
			    			if(outputval & 1) cutExpr += (1 - iot.output[i]);
			    			else cutExpr += iot.output[i];
			    			outputval >>= 1;
			    		}
			    		m.addConstr(cutExpr >= 1);
			    		ctrconstr++;
			    		ctrGlobalConstr++;
			    		if(ctrconstr >= maxNumberConstr)
			    			break;
			    	}
			    }
		    }

		    if(ctrconstr < maxNumberConstr){
		    	//Check the variables with matrix checks (if any)
		    	for(auto const & iot : varsMatrix){
		    		//Extract the input values
					//vectors for fast access afterward
					//Also build the row/column index vectors
					vector<uint> indexRow, indexCol;
					vector<uint> inputval(iot.input.size());
					for(uint i = 0; i < iot.input.size(); i++){
						inputval[i] = uint(round(iot.input[i].get(GRB_DoubleAttr_X)));
						if(inputval[i] == 1)
							indexCol.emplace_back(i);
					}
					vector<uint> outputval(iot.output.size());
					for(uint i = 0; i < iot.output.size(); i++){
						outputval[i] = uint(round(iot.output[i].get(GRB_DoubleAttr_X)));
						if(outputval[i] == 1)
							indexRow.emplace_back(i);
					}

					//Check the corresponding minor
					//If not invertible, remove this solution

					// clock_t begin = clock();
					// auto begin = high_resolution_clock::now();
					auto isInvert = iot.mat.isInvertibleMinor(indexRow, indexCol);
					// auto end = high_resolution_clock::now();
					// duration<double> time_span = duration_cast<duration<double>>(end - begin);
					// cout << "Time spent : " << time_span.count() << " - ";
					// totalSec += (double)(time_span.count());
					// cout << "Total sec : " << totalSec << " sec" << endl;

					if(!isInvert){
						GRBLinExpr cutExpr(0);
			    		for(uint i = 0; i < iot.input.size(); i++){
			    			if(inputval[i] == 1) cutExpr += (1 - iot.input[i]);
			    			else cutExpr += iot.input[i];
			    		}
			    		for(uint i = 0; i < iot.output.size(); i++){
			    			if(outputval[i] == 1) cutExpr += (1 - iot.output[i]);
			    			else cutExpr += iot.output[i];
			    		}
			    		m.addConstr(cutExpr >= 1);
			    		ctrconstr++;
			    		ctrGlobalConstr++;
			    		if(ctrconstr >= maxNumberConstr)
			    			break;
					}
		    	}
		    }

		    if(ctrconstr == 0){
		    	//We did not add any constraints, so the solution passed the checks
		    	return make_tuple(true,ctrGlobalConstr,ctrSol);
		    }
		    else{
		    	//Solution was not valid update the model and solve again
		    	if(DISPLAY_CB_LOG)
		    		cout << "Solution number : " << ctrSol << ", " << ctrGlobalConstr << " constraints added" << endl;
		    	m.update();
		    }
		}


		else if(m.get(GRB_IntAttr_Status) == GRB_INFEASIBLE){
			return make_tuple(false,ctrGlobalConstr,ctrSol);
		}

		else{
			cout << "Solving stopped with status : " << m.get(GRB_IntAttr_Status) << endl;
			return make_tuple(false,ctrGlobalConstr,ctrSol);
		}
	}
}