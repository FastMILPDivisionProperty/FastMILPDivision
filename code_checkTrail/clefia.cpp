#include "clefia.hpp"

using namespace std;

#define DISPLAY_GEN_OUTPUT false
#define DISPLAY_SOLVER_OUTPUT false
#define DISPLAY_NUMBER_ADDED_CONSTRAINT false
#define QUICK_LAZY_COUNT_CLEFIA false

GRBModel getModelCLEFIA(modelInfo & MI){
/*
	Get the model for CLEFIA over #rMax round

	The #sboxModel parameter defines how the Sbox is modelized :
	- "QM" use the Quin-McCluskey algorithm as in Abdelkhalek,Sasaki,Todo,Tolba,Youssef
	- "Simple" use the simplified constraint with PWL

	The #linModel parameter defines how the linear layer is modelized :
	- "CX" modelize the linear layer using the classical copy+xor technique
	- "Simple" modelize the linear layer with the simplified constraint w(x) = w(y)
*/
    string modelName = "CLEFIA_"+to_string(MI.rMax)+"r_"+MI.sboxModel+"_"+MI.linModel;
    modelName += ".mps";

    if(!fileExist("./models/"+modelName)){
        string cmd = "sage ./modelGeneration/genCLEFIA.sage";
        cmd += " -r " + to_string(MI.rMax);
        cmd += " -s " + MI.sboxModel;
        cmd += " -m " + MI.linModel;
        auto sysstdout = exec(cmd.c_str());
        if(DISPLAY_GEN_OUTPUT) cout << sysstdout << endl;
    }

    if(!DISPLAY_SOLVER_OUTPUT)
        MI.gurobiEnv.set(GRB_IntParam_OutputFlag,0);

    GRBModel m(MI.gurobiEnv, "./models/"+modelName);
    return m;
}

modelData getVariablesCLEFIA(GRBModel & m,
                          modelInfo const & MI){
     /*
    Read the model and return the corresponding blocks of variables
    - #MI defines various parameters for how the trails are checked (see aux_function.hpp)
    */

    modelData MD;
    MD.rMax = MI.rMax;

    MD.allTable = vector<vector<vector<uint32_t>>>(6);
    MD.allIneq = vector<vector<vector<int>>>(6);
    MD.allMatrix = vector<Matrix>(2);

    if(MI.sboxCheckType == CheckType::Table){ //Sbox table
        MD.allTable[0] = readDivTableFromFile("./checkData/CLEFIAS0Sbox_table.bin");
    	MD.allTable[1] = readDivTableFromFile("./checkData/CLEFIAS1Sbox_table.bin");
    }
    else if(MI.sboxCheckType == CheckType::Ineq){
        MD.allIneq[0] = readIneqFromFile("./checkData/CLEFIAS0Sbox_ineq.bin");
        MD.allIneq[1] = readIneqFromFile("./checkData/CLEFIAS1Sbox_ineq.bin");
    }

    if(MI.MCCheckType == CheckType::Matrix){
    	MD.allMatrix[0] = Matrix("./checkData/matrixCLEFIAM0.bin");
    	MD.allMatrix[1] = Matrix("./checkData/matrixCLEFIAM1.bin");
    }
    else if(MI.MCCheckType == CheckType::Table){
    	cerr << "Error : Table Check for CLEFIA MC not implemented" << endl;
    	exit(1);
        // MD.allTable[2] = //Table for M0
        // MD.allTable[3] = //Table for M1
    }
    else if(MI.MCCheckType == CheckType::Ineq){
    	cerr << "Error : Ineq Check for CLEFIA MC not implemented" << endl;
    	exit(1);
        // MD.allIneq[2]   = //Ineq for M0
        // MD.allIneq[3]   = //Ineq for M1
    }

    if(MI.SSBCheckType != CheckType::None){
    	cerr << "Error : SSB Check for CLEFIA not implemented" << endl;
    	exit(1);
    }
    // if(SSBCheckType == CheckType::Table)
    //     MD.allTable[4] = //Table for SSB0
    //     MD.allTable[5] = //Table for SSB1
    // else if(SSBCheckType == CheckType::Ineq)
    //     MD.allIneq[4]   = //Ineq for SSB0
    //     MD.allIneq[5]   = //Ineq for SSB1

    auto & tableSbox0 = MD.allTable[0];
    auto & tableSbox1 = MD.allTable[1];

    auto & tableMC0   = MD.allTable[2];
    auto & tableMC1   = MD.allTable[3];

    auto & tableSSB0  = MD.allTable[4];
    auto & tableSSB1  = MD.allTable[5];

    auto & ineqSbox0  = MD.allIneq[0];
    auto & ineqSbox1  = MD.allIneq[1];

    auto & ineqMC0    = MD.allIneq[2];
    auto & ineqMC1    = MD.allIneq[3];

    auto & ineqSSB0   = MD.allIneq[4];
    auto & ineqSSB1   = MD.allIneq[5];

    auto & MCmatrix0 = MD.allMatrix[0];
    auto & MCmatrix1 = MD.allMatrix[1];

    //Sbox variables
    if(MI.sboxCheckType != CheckType::None){
        for(auto const & r : MI.roundOrder){
        	for(uint i = 0; i < 4; i++){
	        	vector<GRBVar> in0(8);
	        	vector<GRBVar> out0(8);
	        	vector<GRBVar> in1(8);
	        	vector<GRBVar> out1(8);
	        	for(uint j = 0; j < 8; j++){
	        		in0[j] = m.getVarByName("cx"+to_string(r)+"_0_"+to_string(8*i+j));
	        		out0[j] = m.getVarByName("y"+to_string(r)+"_0_"+to_string(8*i+j));

					in1[j] = m.getVarByName("cx"+to_string(r)+"_2_"+to_string(8*i+j));
	        		out1[j] = m.getVarByName("y"+to_string(r)+"_2_"+to_string(8*i+j));
	        	}

	        	if(MI.sboxCheckType == CheckType::Table){
	        		if(i%2 == 0){
		        		MD.allVarsTable.emplace_back(move(in0), move(out0), tableSbox0,ConstrType::Sbox);
		                MD.allVarsTable.emplace_back(move(in1), move(out1), tableSbox1,ConstrType::Sbox);
		            }
		            else{
		            	MD.allVarsTable.emplace_back(move(in0), move(out0), tableSbox1,ConstrType::Sbox);
		                MD.allVarsTable.emplace_back(move(in1), move(out1), tableSbox0,ConstrType::Sbox);
		            }
	            }
	            else if(MI.sboxCheckType == CheckType::Ineq){
	                if(i%2 == 0){
		        		MD.allVarsIneq.emplace_back(move(in0), move(out0), ineqSbox0,ConstrType::Sbox);
		                MD.allVarsIneq.emplace_back(move(in1), move(out1), ineqSbox1,ConstrType::Sbox);
		            }
		            else{
		            	MD.allVarsIneq.emplace_back(move(in0), move(out0), ineqSbox1,ConstrType::Sbox);
		                MD.allVarsIneq.emplace_back(move(in1), move(out1), ineqSbox0,ConstrType::Sbox);
		            }
	            }
        	}
        }
    }

    //MC variables

    if(MI.MCCheckType != CheckType::None){
        for(auto const & r : MI.roundOrder){
            vector<GRBVar> in0(32);
            vector<GRBVar> out0(32);
            vector<GRBVar> in1(32);
            vector<GRBVar> out1(32);

            for(uint j = 0; j < 32; j++){
                in0[j] = m.getVarByName("y"+to_string(r)+"_0_"+to_string(j));
                out0[j] = m.getVarByName("z"+to_string(r)+"_0_"+to_string(j));
                in1[j] = m.getVarByName("y"+to_string(r)+"_2_"+to_string(j));
                out1[j] = m.getVarByName("z"+to_string(r)+"_2_"+to_string(j));
            }
            if(MI.MCCheckType == CheckType::Table){
                MD.allVarsTable.emplace_back(move(in0), move(out0), tableMC0,ConstrType::Lin);
                MD.allVarsTable.emplace_back(move(in1), move(out1), tableMC1,ConstrType::Lin);
            }
            else if(MI.MCCheckType == CheckType::Ineq){
                MD.allVarsIneq.emplace_back(move(in0), move(out0), ineqMC0,ConstrType::Lin);
                MD.allVarsIneq.emplace_back(move(in1), move(out1), ineqMC1,ConstrType::Lin);
            }
            else if(MI.MCCheckType == CheckType::Matrix){
            	MD.allVarsMatrix.emplace_back(move(in0), move(out0), MCmatrix0);
            	MD.allVarsMatrix.emplace_back(move(in1), move(out1), MCmatrix1);
            }
        }
    }

    //SSB variables
    if(MI.SSBCheckType != CheckType::None){
        for(auto const & r : MI.roundOrder){
    		vector<GRBVar> in0(32);
            vector<GRBVar> out0(32);
            vector<GRBVar> in1(32);
            vector<GRBVar> out1(32);

            for(uint j = 0; j < 32; j++){
                in0[j] = m.getVarByName("cx"+to_string(r)+"_0_"+to_string(j));
                out0[j] = m.getVarByName("z"+to_string(r)+"_0_"+to_string(j));
                in1[j] = m.getVarByName("cx"+to_string(r)+"_2_"+to_string(j));
                out1[j] = m.getVarByName("z"+to_string(r)+"_2_"+to_string(j));
            }
            if(MI.SSBCheckType == CheckType::Table){
                MD.allVarsTable.emplace_back(move(in0), move(out0), tableSSB0,ConstrType::SSB);
                MD.allVarsTable.emplace_back(move(in1), move(out1), tableSSB1,ConstrType::SSB);
            }
            else if(MI.SSBCheckType == CheckType::Ineq){
                MD.allVarsIneq.emplace_back(move(in0), move(out0), ineqSSB0,ConstrType::SSB);
                MD.allVarsIneq.emplace_back(move(in1), move(out1), ineqSSB1,ConstrType::SSB);
            }
        }
    }

    return MD;
}

uint64_t globalCtrCheckCLEFIA = 0;
uint64_t globalCtrSolCLEFIA = 0;
uint64_t globalCtrConstrSboxCLEFIA = 0;
uint64_t globalCtrConstrLinCLEFIA = 0;
uint64_t globalCtrConstrSSBCLEFIA = 0;
bool existTrailCLEFIA(GRBModel & m,
                   modelData const & MD,
                   vector<uint8_t> const & input,
                   vector<uint8_t> const & output,
                   bool const useCallback,
                   uint const maxNumberConstr){
/*
    Return true if there is a trail from #input to #output using #model for CLEFIA
    - #MD contains the variables to be checked
    - if #useCallback is true then the solving using CB + LC, otherwise repeated solving
    - #maxNumberConstr is the max number of contraints added for each solution found
*/

    uint ctr = 0;
    for(uint i = 0; i < 4; i++){
    	for(uint j = 0; j < 32; j++){
	        m.addConstr(m.getVarByName("x0_"+to_string(i)+"_"+to_string(j)) == input[ctr]);
	        m.addConstr(m.getVarByName("x"+to_string(MD.rMax)+"_"+to_string(i)+"_"+to_string(j)) == output[ctr]);
	        ctr++;
	    }
    }
    m.update();

    if(useCallback){
        m.set(GRB_IntParam_LazyConstraints, 1);
        CustomCallback cb(MD.allVarsTable, MD.allVarsIneq, MD.allVarsMatrix, maxNumberConstr);
        m.setCallback(&cb);

        m.update();
        m.optimize();
        if(DISPLAY_NUMBER_ADDED_CONSTRAINT){
            cout << cb.ctrSol << " solutions examined, ";
            cout << cb.ctrGlobalConstr << " constraints added" << endl;
        }
        globalCtrCheckCLEFIA++;
        globalCtrSolCLEFIA += cb.ctrSol;
        globalCtrConstrSboxCLEFIA += cb.ctrConstrSbox;
		globalCtrConstrLinCLEFIA += cb.ctrConstrLin;
		globalCtrConstrSSBCLEFIA += cb.ctrConstrSSB;

        if(m.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
            return true;
        else if(m.get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
            return false;
        else{
            cout << "Model terminated with status : " << m.get(GRB_IntAttr_Status) << endl;
            return false;
        }
    }
    else{
        auto tmp = solveWithoutCallback(m, MD.allVarsTable, MD.allVarsIneq, MD.allVarsMatrix, maxNumberConstr);
        if(DISPLAY_NUMBER_ADDED_CONSTRAINT){
            cout << get<2>(tmp) << " solutions examined, ";
            cout << get<1>(tmp) << " constraints added" << endl;
        }
        return get<0>(tmp);
    }
}

bool existTrailCLEFIA(modelInfo & MI,
                   vector<uint8_t> const & input,
                   vector<uint8_t> const & output){
/*
    Return true if there is a trail from #input to #output for CLEFIA
    #MI should contain the necessary information to generate the model and parametrize the solving
*/

    GRBModel m = getModelCLEFIA(MI);
    auto MD = getVariablesCLEFIA(m, MI);

    return existTrailCLEFIA(m,MD,input,output,MI.useCallback,MI.maxNumberConstr);
}

void createMatrixFileCLEFIA(){
	//Create the matrix file ./checkData/matrixCLEFIAMi.bin

	vector<vector<uint8_t>> val0({{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1},
								  {0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1},
								  {0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,1,0,1,1,0,0,0,0,1,1},
								  {0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,1,1,0,1,1,0,0,0,1,0},
								  {0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,1,0,0,0,1,1,0,0,1,1,0,0,1,0},
								  {0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,1,0,0,1},
								  {0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0},
								  {0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0},
								  {0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0},
								  {1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1},
								  {0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,1,0,0,0,0,1,1,1,0,0,0,0,0,1,0},
								  {0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,1,0,0,0,1,0,0,1,0,0,0,0,1,1},
								  {0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,1,0,0,1,0,0,0,1,0,0,0,1,1},
								  {0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,1,0,0,0,1,0,0,0,1},
								  {0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0},
								  {0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0},
								  {0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
								  {0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
								  {1,0,0,0,0,0,1,0,1,1,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1},
								  {0,1,0,0,0,0,1,1,0,1,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1},
								  {0,0,1,0,0,0,1,1,0,0,1,1,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1},
								  {0,0,0,1,0,0,0,1,0,0,0,1,1,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0},
								  {0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0},
								  {0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0},
								  {0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0},
								  {1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0},
								  {1,1,0,0,0,0,1,1,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0},
								  {0,1,1,0,0,0,1,0,0,1,0,0,0,0,1,1,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,0},
								  {0,0,1,1,0,0,1,0,0,0,1,0,0,0,1,1,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0},
								  {0,0,0,1,1,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0},
								  {0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0},
								  {0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1}});
	Matrix M0(32,32);
	for(uint i = 0; i < 32; i++){
		for(uint j = 0; j < 32; j++)
			M0.set(i,j,val0[i][j]);
	}

	M0.saveToFile("./checkData/matrixCLEFIAM0.bin");

	vector<vector<uint8_t>> val1({{1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,1},
								  {0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0},
								  {0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,1,0,1,0,0,0,1,0,0},
								  {0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,0,0,0,0,1,1,0,1,0,0,1,1,1},
								  {0,0,0,0,1,0,0,0,0,1,0,0,0,1,1,1,0,0,0,1,0,0,0,1,0,1,0,1,0,1,1,0},
								  {0,0,0,0,0,1,0,0,0,0,1,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,1,0,1,0,1,1},
								  {0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,1,0,1},
								  {0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,1,0},
								  {0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1},
								  {0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0},
								  {0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,1},
								  {1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,1,1,0,0,1,0,0,0,0,1},
								  {0,1,0,0,0,1,1,1,0,0,0,0,1,0,0,0,0,1,0,1,0,1,1,0,0,0,0,1,0,0,0,1},
								  {0,0,1,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,1,1,0,0,0,0,1,0,0,0},
								  {0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,1,0,0,0,0,0,1,0,0},
								  {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0},
								  {0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0},
								  {1,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0},
								  {0,1,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,1},
								  {0,0,1,0,0,0,0,1,1,0,1,0,0,1,1,1,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0},
								  {0,0,0,1,0,0,0,1,0,1,0,1,0,1,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,1,1,1},
								  {0,0,0,0,1,0,0,0,0,0,1,0,1,0,1,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,1,1},
								  {0,0,0,0,0,1,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,1},
								  {0,0,0,0,0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0},
								  {0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0},
								  {1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0},
								  {0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0},
								  {1,0,1,0,0,1,1,1,0,0,1,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0},
								  {0,1,0,1,0,1,1,0,0,0,0,1,0,0,0,1,0,1,0,0,0,1,1,1,0,0,0,0,1,0,0,0},
								  {0,0,1,0,1,0,1,1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,1,0,0},
								  {0,0,0,1,0,1,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0},
								  {0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1}});
	Matrix M1(32,32);
	for(uint i = 0; i < 32; i++){
		for(uint j = 0; j < 32; j++)
			M1.set(i,j,val1[i][j]);
	}

	M1.saveToFile("./checkData/matrixCLEFIAM1.bin");
}

void genSboxFileCLEFIA(){
	//Generate the tables and ineq files for CLEFIA

	vector<uint32_t> S0(
		{0x57,0x49,0xd1,0xc6,0x2f,0x33,0x74,0xfb,0x95,0x6d,0x82,0xea,0x0e,0xb0,0xa8,0x1c,
		 0x28,0xd0,0x4b,0x92,0x5c,0xee,0x85,0xb1,0xc4,0x0a,0x76,0x3d,0x63,0xf9,0x17,0xaf,
		 0xbf,0xa1,0x19,0x65,0xf7,0x7a,0x32,0x20,0x06,0xce,0xe4,0x83,0x9d,0x5b,0x4c,0xd8,
		 0x42,0x5d,0x2e,0xe8,0xd4,0x9b,0x0f,0x13,0x3c,0x89,0x67,0xc0,0x71,0xaa,0xb6,0xf5,
		 0xa4,0xbe,0xfd,0x8c,0x12,0x00,0x97,0xda,0x78,0xe1,0xcf,0x6b,0x39,0x43,0x55,0x26,
		 0x30,0x98,0xcc,0xdd,0xeb,0x54,0xb3,0x8f,0x4e,0x16,0xfa,0x22,0xa5,0x77,0x09,0x61,
		 0xd6,0x2a,0x53,0x37,0x45,0xc1,0x6c,0xae,0xef,0x70,0x08,0x99,0x8b,0x1d,0xf2,0xb4,
		 0xe9,0xc7,0x9f,0x4a,0x31,0x25,0xfe,0x7c,0xd3,0xa2,0xbd,0x56,0x14,0x88,0x60,0x0b,
		 0xcd,0xe2,0x34,0x50,0x9e,0xdc,0x11,0x05,0x2b,0xb7,0xa9,0x48,0xff,0x66,0x8a,0x73,
		 0x03,0x75,0x86,0xf1,0x6a,0xa7,0x40,0xc2,0xb9,0x2c,0xdb,0x1f,0x58,0x94,0x3e,0xed,
		 0xfc,0x1b,0xa0,0x04,0xb8,0x8d,0xe6,0x59,0x62,0x93,0x35,0x7e,0xca,0x21,0xdf,0x47,
		 0x15,0xf3,0xba,0x7f,0xa6,0x69,0xc8,0x4d,0x87,0x3b,0x9c,0x01,0xe0,0xde,0x24,0x52,
		 0x7b,0x0c,0x68,0x1e,0x80,0xb2,0x5a,0xe7,0xad,0xd5,0x23,0xf4,0x46,0x3f,0x91,0xc9,
		 0x6e,0x84,0x72,0xbb,0x0d,0x18,0xd9,0x96,0xf0,0x5f,0x41,0xac,0x27,0xc5,0xe3,0x3a,
		 0x81,0x6f,0x07,0xa3,0x79,0xf6,0x2d,0x38,0x1a,0x44,0x5e,0xb5,0xd2,0xec,0xcb,0x90,
		 0x9a,0x36,0xe5,0x29,0xc3,0x4f,0xab,0x64,0x51,0xf8,0x10,0xd7,0xbc,0x02,0x7d,0x8e});

	cout << "S0" << endl;
	genAndSaveTableIneq(vector<vector<uint32_t>>({S0}),8,8,"./checkData/CLEFIAS0Sbox",true);

	vector<uint32_t> S1(
		{0x6c,0xda,0xc3,0xe9,0x4e,0x9d,0x0a,0x3d,0xb8,0x36,0xb4,0x38,0x13,0x34,0x0c,0xd9,
		 0xbf,0x74,0x94,0x8f,0xb7,0x9c,0xe5,0xdc,0x9e,0x07,0x49,0x4f,0x98,0x2c,0xb0,0x93,
		 0x12,0xeb,0xcd,0xb3,0x92,0xe7,0x41,0x60,0xe3,0x21,0x27,0x3b,0xe6,0x19,0xd2,0x0e,
		 0x91,0x11,0xc7,0x3f,0x2a,0x8e,0xa1,0xbc,0x2b,0xc8,0xc5,0x0f,0x5b,0xf3,0x87,0x8b,
		 0xfb,0xf5,0xde,0x20,0xc6,0xa7,0x84,0xce,0xd8,0x65,0x51,0xc9,0xa4,0xef,0x43,0x53,
		 0x25,0x5d,0x9b,0x31,0xe8,0x3e,0x0d,0xd7,0x80,0xff,0x69,0x8a,0xba,0x0b,0x73,0x5c,
		 0x6e,0x54,0x15,0x62,0xf6,0x35,0x30,0x52,0xa3,0x16,0xd3,0x28,0x32,0xfa,0xaa,0x5e,
		 0xcf,0xea,0xed,0x78,0x33,0x58,0x09,0x7b,0x63,0xc0,0xc1,0x46,0x1e,0xdf,0xa9,0x99,
		 0x55,0x04,0xc4,0x86,0x39,0x77,0x82,0xec,0x40,0x18,0x90,0x97,0x59,0xdd,0x83,0x1f,
		 0x9a,0x37,0x06,0x24,0x64,0x7c,0xa5,0x56,0x48,0x08,0x85,0xd0,0x61,0x26,0xca,0x6f,
		 0x7e,0x6a,0xb6,0x71,0xa0,0x70,0x05,0xd1,0x45,0x8c,0x23,0x1c,0xf0,0xee,0x89,0xad,
		 0x7a,0x4b,0xc2,0x2f,0xdb,0x5a,0x4d,0x76,0x67,0x17,0x2d,0xf4,0xcb,0xb1,0x4a,0xa8,
		 0xb5,0x22,0x47,0x3a,0xd5,0x10,0x4c,0x72,0xcc,0x00,0xf9,0xe0,0xfd,0xe2,0xfe,0xae,
		 0xf8,0x5f,0xab,0xf1,0x1b,0x42,0x81,0xd6,0xbe,0x44,0x29,0xa6,0x57,0xb9,0xaf,0xf2,
		 0xd4,0x75,0x66,0xbb,0x68,0x9f,0x50,0x02,0x01,0x3c,0x7f,0x8d,0x1a,0x88,0xbd,0xac,
		 0xf7,0xe4,0x79,0x96,0xa2,0xfc,0x6d,0xb2,0x6b,0x03,0xe1,0x2e,0x7d,0x14,0x95,0x1d});

	cout << "S1" << endl;
	genAndSaveTableIneq(vector<vector<uint32_t>>({S1}),8,8,"./checkData/CLEFIAS1Sbox",true);

}

void checkNormalDistinguisherCLEFIA(){
	/*
		Search again for a known (no linear combination) distinguisher over CLEFIA to compare execution time
	*/

	vector<uint8_t> input({0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1});

	//These parameters have no effect for CLEFIA but are required for the modelInfo object
	string arxModel = "ARX";
	CheckType ARXCheckType = CheckType::None;
	bool smartLin = true;
	unsigned int startingRound = 1;

	//Fixed parameters
	uint rMax = 10;
	bool firstSbox = true;
	bool lastLin = false;
	CheckType sboxCheckType = CheckType::None;
	CheckType MCCheckType = CheckType::None;
	CheckType SSBCheckType = CheckType::None;
	bool useCallback = true;
	uint maxNumberConstr = 1;
	vector<uint> roundOrder = roundOrderAscend(rMax);

	//Variable parameters
	//Even with the simple model, no need to use the checks, so limited number of variable parameters
	vector<string> listSboxModel({"Simple", "QM"});
	vector<string> listLinModel({"Simple", "CX"});
	// vector<string> listSboxModel({"QM"});
	// vector<string> listLinModel({"CX"});

	//The searchs
	for(auto const & sboxModel : listSboxModel){
		for(auto const & linModel : listLinModel){

			globalCtrCheckCLEFIA = 0;
			globalCtrSolCLEFIA = 0;
			globalCtrConstrSboxCLEFIA = 0;
			globalCtrConstrLinCLEFIA = 0;
			globalCtrConstrSSBCLEFIA = 0;

			modelInfo MI(rMax,sboxModel,linModel,arxModel,firstSbox,lastLin,smartLin,startingRound,sboxCheckType,MCCheckType,SSBCheckType,ARXCheckType,useCallback,maxNumberConstr,roundOrder);

			auto start = chrono::high_resolution_clock::now();
			//Check each output bit (all balanced)
			uint ctrBalanced = 0;
			vector<bool> balancedBits(128);
			auto durBalanced = start - start;

			for(uint i = 0; i < 128; i++){
				vector<uint8_t> output(128,0);
				output[i] = 1;
				auto startBit = chrono::high_resolution_clock::now();
				balancedBits[i] = !existTrailCLEFIA(MI,input,output);
				auto endBit = chrono::high_resolution_clock::now();
				if(balancedBits[i]){
					ctrBalanced++;
					durBalanced += endBit - startBit;
				}
				cout << i << " " << flush;
			}
			cout << endl;
			auto end = chrono::high_resolution_clock::now();
			cout << "sboxModel : " << sboxModel << " - linModel : " << linModel << endl;
			cout << ctrBalanced << " balanced bits" << endl;
			for(uint i = 0; i < 128; i++){
				if(balancedBits[i]) cout << i << " ";
			}
			cout << endl;
			cout << chrono::duration<double>(end - start).count() << " seconds" << endl;
			cout << "Time for only balanced bits : " << chrono::duration<double>(durBalanced).count() << " second" << endl;

			cout << globalCtrCheckCLEFIA << " input/output pairs checked" << endl;
			cout << globalCtrSolCLEFIA << " solutions examined" << endl;
			cout << globalCtrConstrSboxCLEFIA << " constraints added for the sbox" << endl;
			cout << globalCtrConstrLinCLEFIA  << " constraints added for the lin layer" << endl;
			cout << globalCtrConstrSSBCLEFIA << " constraints added for the SSB" << endl;
		}
	}
}

void timingCLEFIA11r(){
	/*
		Search again for a known (no linear combination) distinguisher over CLEFIA to compare execution time
	*/

	vector<uint8_t> input({0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1});

	//These parameters have no effect for CLEFIA but are required for the modelInfo object
	string arxModel = "ARX";
	CheckType ARXCheckType = CheckType::None;
	bool smartLin = true;
	unsigned int startingRound = 1;

	//Fixed parameters
	uint rMax = 11;
	bool firstSbox = true;
	bool lastLin = false;
	// CheckType sboxCheckType = CheckType::None;
	CheckType MCCheckType = CheckType::Matrix;
	CheckType SSBCheckType = CheckType::None;
	bool useCallback = true;
	// uint maxNumberConstr = 1;
	// vector<uint> roundOrder = roundOrderAscend(rMax);

	//Variable parameters
	//Even with the simple model, no need to use the checks, so limited number of variable parameters
	vector<string> listSboxModel({"Simple","QM"});
	vector<string> listLinModel({"Simple","CX"});
	vector<uint> listMaxNbConstr({100,10,1});
	vector<vector<uint>> listRoundOrder({roundOrderOutsideIn(rMax), roundOrderAscend(rMax), roundOrderDescend(rMax)});

	//The searchs
	for(auto const & maxNumberConstr : listMaxNbConstr){
		for(auto const & roundOrder : listRoundOrder){
			for(auto const & sboxModel : listSboxModel){
				CheckType sboxCheckType = CheckType::None;
				if(sboxModel == "Simple")
					sboxCheckType = CheckType::Ineq;
				for(auto const & linModel : listLinModel){

					globalCtrCheckCLEFIA = 0;
				 	globalCtrSolCLEFIA = 0;
	      	 	 	globalCtrConstrSboxCLEFIA = 0;
				 	globalCtrConstrLinCLEFIA = 0;
				 	globalCtrConstrSSBCLEFIA = 0;

					cout << "Starting round : " << startingRound << " - ";
					cout << "Round Order : ";
				    for(auto const & tmpr : roundOrder) cout << tmpr << " ";
				    cout << endl;
				    cout << " sboxModel : " << sboxModel << " - linModel : " << linModel << endl;
					cout << "sboxCheckType : " << sboxCheckType << " - MCCheckType : " << MCCheckType << endl;
				    cout << "maxNumberConstr : " << maxNumberConstr << endl;

					modelInfo MI(rMax,sboxModel,linModel,arxModel,firstSbox,lastLin,smartLin,startingRound,sboxCheckType,MCCheckType,SSBCheckType,ARXCheckType,useCallback,maxNumberConstr,roundOrder);

					auto start = chrono::high_resolution_clock::now();
					auto durBalanced = start - start;
					//Check each output bit (half balanced)
					vector<bool> balancedBits(128);
					uint ctrBalanced = 0;
					for(uint i = 0; i < 128; i++){
						vector<uint8_t> output(128,0);
						output[i] = 1;
						auto startBit = chrono::high_resolution_clock::now();
						balancedBits[i] = !existTrailCLEFIA(MI,input,output);
						auto endBit = chrono::high_resolution_clock::now();
						if(balancedBits[i]){
							ctrBalanced++;
							durBalanced += endBit - startBit;
						}
					}
					auto end = chrono::high_resolution_clock::now();
					cout << ctrBalanced << " balanced bits" << endl;
					cout << chrono::duration<double>(end - start).count() << " seconds" << endl;
					cout << "Time for only balanced bits : " << chrono::duration<double>(durBalanced).count() << " second" << endl;
					cout << globalCtrCheckCLEFIA << " input/output pairs checked" << endl;
				 	cout << globalCtrSolCLEFIA << " solutions examined" << endl;
				 	cout << globalCtrConstrSboxCLEFIA << " constraints added for the sbox" << endl;
				 	cout << globalCtrConstrLinCLEFIA  << " constraints added for the lin layer" << endl;
				 	cout << globalCtrConstrSSBCLEFIA << " constraints added for the SSB" << endl;
				}
			}
		}
	}
	cout << endl;
}


void searchCLEFIA10r(){


	vector<uint8_t> input({0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1});

	//These parameters have no effect for CLEFIA but are required for the modelInfo object
	string arxModel = "ARX";
	CheckType ARXCheckType = CheckType::None;
	bool smartLin = true;
	unsigned int startingRound = 1;

	//Fixed parameters
	uint rMax = 10;
	bool firstSbox = true;
	bool lastLin = false;
	CheckType sboxCheckType = CheckType::None;
	CheckType MCCheckType = CheckType::None;
	CheckType SSBCheckType = CheckType::None;
	bool useCallback = true;
	uint maxNumberConstr = 1;
	vector<uint> roundOrder = roundOrderAscend(rMax);

	string sboxModel = "Simple";
	string linModel = "Simple";

	globalCtrCheckCLEFIA = 0;
	globalCtrSolCLEFIA = 0;
	globalCtrConstrSboxCLEFIA = 0;
	globalCtrConstrLinCLEFIA = 0;
	globalCtrConstrSSBCLEFIA = 0;

	modelInfo MI(rMax,sboxModel,linModel,arxModel,firstSbox,lastLin,smartLin,startingRound,sboxCheckType,MCCheckType,SSBCheckType,ARXCheckType,useCallback,maxNumberConstr,roundOrder);

	auto start = chrono::high_resolution_clock::now();
	//Check each output bit (all balanced)
	uint ctrBalanced = 0;
	vector<bool> balancedBits(128);
	auto durBalanced = start - start;

	for(uint i = 0; i < 128; i++){
		vector<uint8_t> output(128,0);
		output[i] = 1;
		auto startBit = chrono::high_resolution_clock::now();
		balancedBits[i] = !existTrailCLEFIA(MI,input,output);
		auto endBit = chrono::high_resolution_clock::now();
		if(balancedBits[i]){
			ctrBalanced++;
			durBalanced += endBit - startBit;
		}
		cout << i << " " << flush;
	}
	cout << endl;
	auto end = chrono::high_resolution_clock::now();
	cout << "sboxModel : " << sboxModel << " - linModel : " << linModel << endl;
	cout << ctrBalanced << " balanced bits" << endl;
	for(uint i = 0; i < 128; i++){
		if(balancedBits[i]) cout << i << " ";
	}
	cout << endl;
	cout << chrono::duration<double>(end - start).count() << " seconds" << endl;
	cout << "Time for only balanced bits : " << chrono::duration<double>(durBalanced).count() << " second" << endl;

	cout << globalCtrCheckCLEFIA << " input/output pairs checked" << endl;
	cout << globalCtrSolCLEFIA << " solutions examined" << endl;
	cout << globalCtrConstrSboxCLEFIA << " constraints added for the sbox" << endl;
	cout << globalCtrConstrLinCLEFIA  << " constraints added for the lin layer" << endl;
	cout << globalCtrConstrSSBCLEFIA << " constraints added for the SSB" << endl;

}

void searchCLEFIA11r(){

	uint rMax = 11;
	string sboxModel = "QM";
	string linModel = "CX";
	string arxModel = "ARX";
	bool firstSbox = true;
	bool lastLin = false;
	bool smartLin = true;
	CheckType sboxCheckType = CheckType::None;
	CheckType MCCheckType = CheckType::Matrix;
	CheckType SSBCheckType = CheckType::None;
	CheckType ARXCheckType = CheckType::None;
	bool useCallback = true;
	uint startingRound = 1;
	uint maxNumberConstr = 100;
	vector<uint> roundOrder = roundOrderOutsideIn(rMax);

	globalCtrCheckCLEFIA = 0;
	globalCtrSolCLEFIA = 0;
	globalCtrConstrSboxCLEFIA = 0;
	globalCtrConstrLinCLEFIA = 0;
	globalCtrConstrSSBCLEFIA = 0;

	for(uint indexInput = 0; indexInput < 128; indexInput++){
		vector<uint8_t> input(128,1);
		input[indexInput] = 0;
		modelInfo MI(rMax,sboxModel,linModel,arxModel,firstSbox,lastLin,smartLin,startingRound,sboxCheckType,MCCheckType,SSBCheckType,ARXCheckType,useCallback,maxNumberConstr,roundOrder);

		auto start = chrono::high_resolution_clock::now();
		auto durBalanced = start - start;
		//Check each output bit (half balanced)
		vector<bool> balancedBits(128);
		uint ctrBalanced = 0;
		for(uint i = 0; i < 128; i++){
			vector<uint8_t> output(128,0);
			output[i] = 1;
			auto startBit = chrono::high_resolution_clock::now();
			balancedBits[i] = !existTrailCLEFIA(MI,input,output);
			auto endBit = chrono::high_resolution_clock::now();
			if(balancedBits[i]){
				ctrBalanced++;
				durBalanced += endBit - startBit;
			}
		}
		auto end = chrono::high_resolution_clock::now();
		cout << "Constant input bit : " << indexInput << endl;
		cout << ctrBalanced << " balanced bits" << endl;
		cout << chrono::duration<double>(end - start).count() << " seconds" << endl;
		cout << "Time for only balanced bits : " << chrono::duration<double>(durBalanced).count() << " second" << endl;
		if(QUICK_LAZY_COUNT_CLEFIA)
			break;
	}

	cout << globalCtrCheckCLEFIA << " input/output pairs checked" << endl;
	cout << globalCtrSolCLEFIA << " solutions examined" << endl;
	cout << globalCtrConstrSboxCLEFIA << " constraints added for the sbox" << endl;
	cout << globalCtrConstrLinCLEFIA  << " constraints added for the lin layer" << endl;
	cout << globalCtrConstrSSBCLEFIA << " constraints added for the SSB" << endl;


}