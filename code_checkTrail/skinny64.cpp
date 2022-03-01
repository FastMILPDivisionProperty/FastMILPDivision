#include "skinny64.hpp"
#include "DivTable.hpp"

using namespace std;

#define DISPLAY_GEN_OUTPUT false
#define DISPLAY_SOLVER_OUTPUT false
#define DISPLAY_NUMBER_ADDED_CONSTRAINT false

vector<vector<pair<vector<uint8_t>, vector<uint8_t>>>> allInOutPairsSKINNY64(unsigned pos_in, unsigned pos_out, unsigned R) { // 0 <= i < 4  - 0 <= j < 4
  vector<uint16_t> sbox = {0xc, 6, 9, 0,   1, 0xa, 2, 0xb,
                            3,   8, 5, 0xd, 4, 0xe, 7, 0xf}; // SKINNY
  vector<uint16_t> linear = {
    0b0001'0001'0000'0001,
    0b0010'0010'0000'0010,
    0b0100'0100'0000'0100,
    0b1000'1000'0000'1000,
    0b0000'0000'0000'0001,
    0b0000'0000'0000'0010,
    0b0000'0000'0000'0100,
    0b0000'0000'0000'1000,
    0b0000'0001'0001'0000,
    0b0000'0010'0010'0000,
    0b0000'0100'0100'0000,
    0b0000'1000'1000'0000,
    0b0000'0001'0000'0001,
    0b0000'0010'0000'0010,
    0b0000'0100'0000'0100,
    0b0000'1000'0000'1000
  };

  vector<uint8_t> perm_inv = {0, 13, 10, 7, 4, 1, 14, 11, 8, 5, 2, 15, 12, 9, 6, 3};
  vector<uint8_t> perm (16);
  for (unsigned i = 0; i < 16; ++i) perm[perm_inv[i]] = i;
  vector<uint8_t> permBit({0,1,2,3,20,21,22,23,40,41,42,43,60,61,62,63,16,17,18,19,36,37,38,39,56,57,58,59,12,13,14,15,32,33,34,35,52,53,54,55,8,9,10,11,28,29,30,31,48,49,50,51,4,5,6,7,24,25,26,27,44,45,46,47});
  vector<uint8_t> permBit_inv(64);
  for(uint i = 0; i < 64; i++)
  	permBit_inv[permBit[i]] = i;

  auto anf_s = getANF_Parallel_S(4, sbox);

  auto anf_kos = anf_s;
  for (unsigned k = 0; k < 8; ++k) {
    anf_kos[k] += Monome(uint32_t(1) << (16 + k));
    //anf_kos[4*k+1] += Monome(uint32_t(1) << (16 + k + 4));
  }
  auto anf_l = convertLtoPolynome(linear);
  auto anf_lokos = composition(anf_l, anf_kos);
  auto anf_solokos = composition(anf_s, anf_lokos);

  auto anf_sol = composition(anf_s, anf_l);

  auto anf_cokos = anf_kos;
  anf_cokos[9] += Monome(0); //add 1
  auto anf_locokos = composition(anf_l, anf_cokos);
  auto anf_solocokos = composition(anf_s, anf_locokos);

  auto v_in = (pos_in == 0) ? loadTransitionsIN(anf_solocokos) : loadTransitionsIN(anf_solokos);
  auto v_out = (R%2 == 0) ? ((pos_out == 0) ? loadTransitionsOUT(anf_solocokos) : loadTransitionsOUT(anf_solokos) ): loadTransitionsOUT(anf_sol);

  vector<vector<pair<vector<uint8_t>, vector<uint8_t>>>> res;
  for (auto pin : v_in) {
    for (auto pout : v_out) {
      vector<pair<vector<uint8_t>, vector<uint8_t>>> tmp;
      for (auto xin : pin.first) {
        for (auto xout : pout.first) {
          vector<uint8_t> in_8t (64, 1);
          vector<uint8_t> out_8t (64, 0);
          for (unsigned b = 0; b < 16; ++b) in_8t[16*pos_in + b] = (xin >> b) & 1;
          for (unsigned b = 0; b < 16; ++b) out_8t[16*pos_out + b] = (xout >> b) & 1;
          //Permute out to fit the modeling
          vector<uint8_t> pout_8t(64);
          for(uint i = 0; i < 64; i++)
          	pout_8t[permBit_inv[i]] = out_8t[i];
          tmp.emplace_back(move(in_8t), move(pout_8t));
        }
      }
      res.emplace_back(move(tmp));
    }
  }

  return res;

}

GRBModel getModelSkinny64(modelInfo & MI){
/*
	Get the model for Skinny64 over #rMax round
	#lastLin defines if the linear layer is applied on the last round
	#smartLin only has effect if linModel=="Simple" and defines if MC is seen as independent Lboxes or as a single matrix

	The #sboxModel parameter defines how the Sbox is modelized :
	- "Hull" use the Convex Hull technique
	- "Simple" use the simplified constraint with PWL

	The #linModel parameter defines how the linear layer is modelized :
	- "Lbox" modelize the linear layer as the parallel application of the same linear Sbox (thanks to the specific form of the matrix)
	- "ZR" modelize the linear layer using the technique from Zhang and Rijmen
	- "CX" modelize the linear layer using the classical copy+xor technique
	- "Simple" modelize the linear layer with the simplified constraint w(x) = w(y)
*/
    string modelName = "Skinny64_"+to_string(MI.rMax)+"r_"+MI.sboxModel+"_"+MI.linModel;
    if(MI.lastLin) modelName += "_lastLin";
    if(MI.linModel=="Simple" && MI.smartLin) modelName += "_smartLin";
    if(!MI.firstSbox) modelName += "_noFirstSbox";
    modelName += ".mps";

    if(!fileExist("./models/"+modelName)){
        string cmd = "sage ./modelGeneration/genSkinny64.sage";
        cmd += " -r " + to_string(MI.rMax);
        cmd += " -s " + MI.sboxModel;
        cmd += " -m " + MI.linModel;
        if(MI.lastLin) cmd += " --lastLin True";
        else cmd += " --lastLin False";
        if(MI.smartLin) cmd += " --smartLin True";
        else cmd += " --smartLin False";
        if(MI.firstSbox) cmd += " --firstSbox True";
        else cmd += " --firstSbox False";
        auto sysstdout = exec(cmd.c_str());
        if(DISPLAY_GEN_OUTPUT) cout << sysstdout << endl;
    }

    if(!DISPLAY_SOLVER_OUTPUT)
        MI.gurobiEnv.set(GRB_IntParam_OutputFlag,0);

    GRBModel m(MI.gurobiEnv, "./models/"+modelName);
    return m;
}

modelData getVariablesSkinny64(GRBModel & m,
                        	   modelInfo const & MI){
/*
	Read the model and return the corresponding blocks of variables
	- #MI defines various parameters for how the trails are checked (see aux_function.hpp)
*/

    modelData MD;
    MD.rMax = MI.rMax;
    MD.firstSbox = MI.firstSbox;
    MD.lastLin = MI.lastLin;

    MD.allTable = vector<vector<vector<uint32_t>>>(4);
    MD.allIneq = vector<vector<vector<int>>>(4);
    MD.allMatrix = vector<Matrix>(1);

    if(MI.sboxCheckType == CheckType::Table) //Sbox table
        MD.allTable[0] = vector<vector<uint32_t>>({{0},{1,2,8},{1,2,4},{1,2,12},{1,2,4,8},{1,6,10,12},{1,2,4},{1,6,12},{1,2,4,8},{1,6,12},{1,2,12},{3,6,13},{1,2,8},{3,14},{1,2,12},{15}});
    else if(MI.sboxCheckType == CheckType::Ineq)
        MD.allIneq[0] = vector<vector<int>>({{1, 1, 1, 1, -1, -1, -1, -1, 0, 1},
											 {-2, 0, -1, -1, 3, 2, 1, 2, 0, 1},
											 {0, -1, 0, 0, 0, 0, 1, -1, 1, 1},
											 {1, 0, 0, 0, -1, -1, -1, 1, 1, 1},
											 {-3, -1, -1, -2, 1, 1, 2, -1, 4, 1},
											 {0, 0, 0, 3, -2, -1, -1, -1, 2, 1},
											 {-1, 0, -1, -1, 0, 1, -1, 1, 2, 1},
											 {2, 1, 0, 0, -1, -1, -1, -1, 1, 1},
											 {-1, -1, 0, -1, 2, 2, 1, 1, 0, 1},
											 {0, 0, 1, 1, 0, -2, -1, -1, 2, 1},
											 {0, 1, 1, 0, -1, -1, 0, -1, 1, 1},
											 {1, 0, 0, 1, -1, 0, -1, -1, 1, 1}});

    if(MI.smartLin){
    	if(MI.MCCheckType == CheckType::Table)
		    MD.allTable[1] = vector<vector<uint32_t>>({{0},{1,2,8},{4},{5,6,12},{1,4,8},{3,5,6,10,12},{5,12},{7,14},{1},{3,9},{5},{7,13},{5,9},{7,11,13},{13},{15}});
		else if(MI.MCCheckType == CheckType::Ineq)
		    MD.allIneq[1] = vector<vector<int>>({{1, 1, 1, 1, -1, -1, -1, -1, 0, 0},
												 {1, 0, 1, 1, -1, -1, 0, -1, 0, 1},
												 {0, 0, 0, -1, 1, 0, 0, 0, 0, 1},
												 {1, 0, 0, 0, 0, -1, 0, 0, 0, 1},
												 {0, 0, 0, 1, -1, 0, 0, -1, 1, 1},
												 {-1, 0, 0, -1, 1, 1, 0, 1, 0, 1}});
		else if(MI.MCCheckType == CheckType::Matrix)
    		MD.allMatrix[0] = Matrix("./checkData/matrixSkinny64_small.bin");
    }
     else{
    	if(MI.MCCheckType == CheckType::Matrix)
    		MD.allMatrix[0] = Matrix("./checkData/matrixSkinny64.bin");
	    else if(MI.MCCheckType == CheckType::Table){
	    	cerr << "Error : Table Check for Skinny64 MC (non smartLin) not implemented, probably not worth it compared to using smartLin=true" << endl;
	    	exit(1);
	        // MD.allTable[1] = //Table for MC
	    }
	    else if(MI.MCCheckType == CheckType::Ineq){
	    	cerr << "Error : Ineq Check for Skinny64 MC (non smartLin) not implemented, probably not worth it compared to using smartLin=true" << endl;
	    	exit(1);
	        // MD.allIneq[1]   = //Ineq for MC
	    }
    }

    if(MI.SSBCheckType == CheckType::Table){
        MD.allTable[2] = readDivTableFromFile("./checkData/Skinny64SSB0_table.bin");
        MD.allTable[3] = readDivTableFromFile("./checkData/Skinny64SSB1_table.bin");
    }
    else if(MI.SSBCheckType == CheckType::Ineq){
        MD.allIneq[2]   = readIneqFromFile("./checkData/Skinny64SSB0_ineq.bin");
        MD.allIneq[3]   = readIneqFromFile("./checkData/Skinny64SSB1_ineq.bin");
    }

    auto & tableSbox = MD.allTable[0];
    auto & tableMC   = MD.allTable[1];

    auto & ineqSbox  = MD.allIneq[0];
    auto & ineqMC    = MD.allIneq[1];

    auto & MCmatrix = MD.allMatrix[0];

    auto & tableSSB0  = MD.allTable[2];
    auto & tableSSB1  = MD.allTable[3];

    auto & ineqSSB0   = MD.allIneq[2];
    auto & ineqSSB1   = MD.allIneq[3];

    //Sbox variables
    uint roundFirstSbox = 0;
    if(!MI.firstSbox) roundFirstSbox = 1;

    if(MI.sboxCheckType != CheckType::None){
        for(auto const & r : MI.roundOrder){
        	if(r >= roundFirstSbox){
	            for(uint i = 0; i < 16; i++){
	                vector<GRBVar> in(4);
	                vector<GRBVar> out(4);
	                for(uint j = 0; j < 4; j++){
	                    in[j] = m.getVarByName("x"+to_string(r)+"_"+to_string(4*i+j));
	                    out[j] = m.getVarByName("y"+to_string(r)+"_"+to_string(4*i+j));
	                }
	                if(MI.sboxCheckType == CheckType::Table)
	                    MD.allVarsTable.emplace_back(move(in), move(out), tableSbox,ConstrType::Sbox);
	                else if(MI.sboxCheckType == CheckType::Ineq)
	                    MD.allVarsIneq.emplace_back(move(in), move(out), ineqSbox,ConstrType::Sbox);
	            }
	        }
        }
    }

    //MC variables
    uint roundLastLin = MI.rMax;
    if(!MI.lastLin) roundLastLin--;

    if(MI.MCCheckType != CheckType::None){
    	if(MI.smartLin){
    		//Consider MC as independent Lboxes
	        for(auto const & r : MI.roundOrder){
	        	if(r < roundLastLin){
		            for(uint col = 0; col < 4; col++){ //For each column
		            	for(uint offset = 0; offset < 4; offset++){ //For each set of 4 bits
			                vector<GRBVar> in(4);
			                vector<GRBVar> out(4);
			                for(uint j = 0; j < 4; j++){
			                    in[j] = m.getVarByName("z"+to_string(r)+"_"+to_string(16*col + offset + j*4));
			                    out[j] = m.getVarByName("x"+to_string(r+1)+"_"+to_string(16*col + offset + j*4));
			                }
			                if(MI.MCCheckType == CheckType::Table)
			                    MD.allVarsTable.emplace_back(move(in), move(out), tableMC,ConstrType::Lin);
			                else if(MI.MCCheckType == CheckType::Ineq)
			                    MD.allVarsIneq.emplace_back(move(in), move(out), ineqMC,ConstrType::Lin);
			                else if(MI.MCCheckType == CheckType::Matrix)
	                			MD.allVarsMatrix.emplace_back(move(in), move(out), MCmatrix);
			            }
		            }
		        }
	        }
    	}
    	else{
    		//Consider the matrix as a whole
	        for(auto const & r : MI.roundOrder){
	        	if(r < roundLastLin){
		            for(uint i = 0; i < 4; i++){
		                vector<GRBVar> in(16);
		                vector<GRBVar> out(16);
		                for(uint j = 0; j < 16; j++){
		                    in[j] = m.getVarByName("z"+to_string(r)+"_"+to_string(16*i+j));
		                    out[j] = m.getVarByName("x"+to_string(r+1)+"_"+to_string(16*i+j));
		                }
		                if(MI.MCCheckType == CheckType::Table)
		                    MD.allVarsTable.emplace_back(move(in), move(out), tableMC,ConstrType::Lin);
		                else if(MI.MCCheckType == CheckType::Ineq)
		                    MD.allVarsIneq.emplace_back(move(in), move(out), ineqMC,ConstrType::Lin);
		                else if(MI.MCCheckType == CheckType::Matrix)
	                		MD.allVarsMatrix.emplace_back(move(in), move(out), MCmatrix);
		            }
		        }
	        }
	    }
    }

    //SSB variables
    //Input/Output nibbles for each SSB are as follow
    //{ 0, 13, 10,  7} -> {0,1,2,3}
    //{ 4,  1, 14, 11} -> {4,5,6,7}
    //{ 8,  5,  2, 15} -> {8,9,10,11}
    //{12,  9,  6,  3} -> {12,13,14,15}
    static const uint8_t indexIn[4][4] = {{0,13,10,7},{ 4,1,14,11},{8,5,2,15},{12,9,6,3}};
    static const uint8_t indexOut[4][4] = {{0,1,2,3},{4,5,6,7},{8,9,10,11},{12,13,14,15}};
    roundLastLin--; //Last SSB starts one round before the last linear layer

    if(MI.SSBCheckType != CheckType::None){
        for(auto const & r : MI.roundOrder){
        	if(r < roundLastLin&& r >= roundFirstSbox){
	            for(uint i = 0; i < 4; i++){
	                vector<GRBVar> in(16);
	                vector<GRBVar> out(16);
	                uint ctr = 0;
	                for(auto const j : indexIn[i]){
	                    for(uint k = 0; k < 4; k++){
	                        in[ctr] = m.getVarByName("x"+to_string(r)+"_"+to_string(4*j+k));
	                        ctr++;
	                    }
	                }
	                ctr = 0;
	                for(auto const j : indexOut[i]){
	                    for(uint k = 0; k < 4; k++){
	                        out[ctr] = m.getVarByName("y"+to_string(r+1)+"_"+to_string(4*j+k));
	                        ctr++;
	                    }
	                }
	                if(i == 0){
		                if(MI.SSBCheckType == CheckType::Table)
		                    MD.allVarsTable.emplace_back(move(in), move(out), tableSSB0,ConstrType::SSB);
		                else if(MI.SSBCheckType == CheckType::Ineq)
		                    MD.allVarsIneq.emplace_back(move(in), move(out), ineqSSB0,ConstrType::SSB);
		            }
		            else{
		            	if(MI.SSBCheckType == CheckType::Table)
		                    MD.allVarsTable.emplace_back(move(in), move(out), tableSSB1,ConstrType::SSB);
		                else if(MI.SSBCheckType == CheckType::Ineq)
		                    MD.allVarsIneq.emplace_back(move(in), move(out), ineqSSB1,ConstrType::SSB);
		            }
	            }
	        }
        }
    }

    return MD;
}

uint64_t globalCtrCheckSkinny = 0;
uint64_t globalCtrSolSkinny = 0;
uint64_t globalCtrConstrSboxSkinny = 0;
uint64_t globalCtrConstrLinSkinny = 0;
uint64_t globalCtrConstrSSBSkinny = 0;
bool existTrailSkinny64(GRBModel & m,
                   modelData const & MD,
                   vector<uint8_t> const & input,
                   vector<uint8_t> const & output,
                   bool const useCallback,
                   uint const maxNumberConstr){
/*
    Return true if there is a trail from #input to #output using #model for Skinny64
    - #MD contains the variables to be checked
    - if #useCallback is true then the solving using CB + LC, otherwise repeated solving
    - #maxNumberConstr is the max number of contraints added for each solution found
*/

    string outputvarPrefix = "x"+to_string(MD.rMax)+"_";
    if(!MD.lastLin)
        outputvarPrefix = "y"+to_string(MD.rMax-1)+"_";
    string inputvarPrefix = "x0_";
    if(!MD.firstSbox)
    	inputvarPrefix = "y0_";

    for(uint i = 0; i < 64; i++){
        m.addConstr(m.getVarByName(inputvarPrefix+to_string(i)) == input[i]);
        m.addConstr(m.getVarByName(outputvarPrefix+to_string(i)) == output[i]);
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
        globalCtrCheckSkinny++;
        globalCtrSolSkinny += cb.ctrSol;
        globalCtrConstrSboxSkinny += cb.ctrConstrSbox;
		globalCtrConstrLinSkinny += cb.ctrConstrLin;
		globalCtrConstrSSBSkinny += cb.ctrConstrSSB;


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

bool existTrailSkinny64(modelInfo & MI,
                   vector<uint8_t> const & input,
                   vector<uint8_t> const & output){
/*
    Return true if there is a trail from #input to #output for Skinny64
    #MI should contain the necessary information to generate the model and parametrize the solving
*/

    GRBModel m = getModelSkinny64(MI);
    auto MD = getVariablesSkinny64(m, MI);

    return existTrailSkinny64(m,MD,input,output,MI.useCallback,MI.maxNumberConstr);
}

void createMatrixFileSkinny64(){
	//Create the matrix file ./checkData/matrixSkinny64.bin

	vector<vector<uint8_t>> val({{1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0},
								 {0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,0},
								 {0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0},
								 {0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1},
								 {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
								 {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
								 {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
								 {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0},
								 {0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0},
								 {0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0},
								 {0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0},
								 {0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0},
								 {1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
								 {0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0},
								 {0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0},
								 {0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0}});
	Matrix M(16,16);
	for(uint i = 0; i < 16; i++){
		for(uint j = 0; j < 16; j++)
			M.set(i,j,val[i][j]);
	}

	M.saveToFile("./checkData/matrixSkinny64.bin");

	val = vector<vector<uint8_t>>({{1,0,1,1},
								   {1,0,0,0},
								   {0,1,1,0},
								   {1,0,1,0}});

	M = Matrix(4,4);
	for(uint i = 0; i < 4; i++){
		for(uint j = 0; j < 4; j++)
			M.set(i,j,val[i][j]);
	}
	M.saveToFile("./checkData/matrixSkinny64_small.bin");
}

void checkNormalDistinguisherSkinny64(){
	/*
		Search for a normal (no linear combination) distinguisher over 10 rounds of Skinny64 
	*/

	vector<uint8_t> input({1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1});

	//These parameters have no effect for skinny64 but are required for the modelInfo object
	string arxModel = "ARX";
	unsigned int startingRound = 1;
	CheckType ARXCheckType = CheckType::None;

	//Fixed parameters
	CheckType SSBCheckType = CheckType::None; //no need to check SSB for this one
	uint rMax = 10;
	bool smartLin = true; //probably not useful to not have it honestly
	bool firstSbox = true;
	bool lastLin = true;

	//Variable parameters
	// string sboxModel = "QM";
	// string linModel = "CX";
	// CheckType sboxCheckType = CheckType::None;
	// CheckType MCCheckType = CheckType::Matrix;
	// bool useCallback = false;
	// uint maxNumberConstr = 4;
	// vector<uint> roundOrder = roundOrderOutsideIn(rMax);

	vector<string> listSboxModel({"Hull", "Simple"});
	vector<string> listLinModel({"Lbox", "ZR", "CX", "Simple"});
	vector<bool> listUseCallback({true});
	vector<vector<uint>> listroundOrderIfNeeded({roundOrderOutsideIn(rMax), roundOrderAscend(rMax), roundOrderDescend(rMax)});
	// vector<uint> listMaxNbConstr({500,100,10,1});

	uint ctrParameterSet = 0;
	uint firstParameterSet = 60;

	for(auto const & useCallback : listUseCallback){
	 for(auto const & linModel : listLinModel){
	  vector<CheckType> listMCCheckType({CheckType::Ineq});
	  if(linModel == "Lbox" || linModel == "ZR") //No need for checks on MC
	   listMCCheckType = vector<CheckType>({CheckType::None});
	  for(auto const & MCCheckType : listMCCheckType){
	  	for(auto const & sboxModel : listSboxModel){
	    vector<CheckType> listSboxCheckType({CheckType::Ineq});
	    if(sboxModel == "Hull" || sboxModel == "QM") //No need for checks on Sbox
	     listSboxCheckType = vector<CheckType>({CheckType::None});
	    for(auto const & sboxCheckType : listSboxCheckType){
	      vector<uint> listMaxNbConstr({1});
	      if(!(sboxCheckType == CheckType::None && MCCheckType == CheckType::None && SSBCheckType == CheckType::None))
	        listMaxNbConstr = vector<uint>({500,100,10,1});
	      for(auto const & maxNumberConstr : listMaxNbConstr){
	      	vector<vector<uint>> listroundOrder(1);
	      	if(!(sboxCheckType == CheckType::None && MCCheckType == CheckType::None && SSBCheckType == CheckType::None))
	      		listroundOrder = vector<vector<uint>>(listroundOrderIfNeeded);
	      	for(auto const & roundOrder : listroundOrder){

	       	if(ctrParameterSet >= firstParameterSet){

	       		globalCtrCheckSkinny = 0;
				globalCtrSolSkinny = 0;
				globalCtrConstrSboxSkinny = 0;
				globalCtrConstrLinSkinny = 0;
				globalCtrConstrSSBSkinny = 0;

		       	modelInfo MI(rMax,sboxModel,linModel,arxModel,firstSbox,lastLin,smartLin,startingRound,sboxCheckType,MCCheckType,SSBCheckType,ARXCheckType,useCallback,maxNumberConstr,roundOrder);
		       	cout << "Parameter set " << ctrParameterSet << endl;
		       	if(useCallback) cout << "Using callback" << endl;
				else cout << "Using iterative model" << endl;
				cout << "roundOrder : ";
				for(auto const & tmpr : roundOrder) cout << tmpr << " ";
				cout << endl;
				cout << "linModel : " << linModel << " - MCCheckType : " << MCCheckType << endl;
				cout << "sboxModel : " << sboxModel << " - sboxCheckType : " << sboxCheckType << endl;
				cout << "maxNumberConstr : " << maxNumberConstr << endl;

				auto start = chrono::high_resolution_clock::now();
				auto durBalanced = start - start;
				//Check each output bit 
				vector<bool> balancedBits(64);
				uint ctrBalanced = 0;
				for(uint i = 0; i < 64; i++){
					vector<uint8_t> output(64,0);
					output[i] = 1;
					auto startBit = chrono::high_resolution_clock::now();
					balancedBits[i] = !existTrailSkinny64(MI,input,output);
					auto endBit = chrono::high_resolution_clock::now();
					if(balancedBits[i]){
						ctrBalanced++;
						durBalanced += endBit - startBit;
					}
				}
				auto end = chrono::high_resolution_clock::now();
				cout << ctrBalanced << " balanced bits" << endl;
				cout << chrono::duration<double>(end - start).count() << " seconds, ";
				cout << "Time for only balanced bits : " << chrono::duration<double>(durBalanced).count() << " second" << endl;
				cout << endl;
				cout << globalCtrCheckSkinny << " input/output pairs checked" << endl;
				cout << globalCtrSolSkinny << " solutions examined" << endl;
				cout << globalCtrConstrSboxSkinny << " constraints added for the sbox" << endl;
				cout << globalCtrConstrLinSkinny  << " constraints added for the lin layer" << endl;
				cout << globalCtrConstrSSBSkinny << " constraints added for the SSB" << endl;
			}
			ctrParameterSet++;
	}}}}}}}
}


void timingSearchSkinny64(){
	/*
	*/

	//These parameters have no effect for skinny64 but are required for the modelInfo object
	string arxModel = "ARX";
	unsigned int startingRound = 1;
	CheckType ARXCheckType = CheckType::None;

	//Fixed parameters
	CheckType SSBCheckType = CheckType::Ineq;
	uint rMax = 8;
	bool smartLin = true; //probably not useful to not have it honestly
	bool firstSbox = false;
	bool lastLin = true;

	//Variable parameters
	// string sboxModel = "QM";
	// string linModel = "CX";
	// CheckType sboxCheckType = CheckType::None;
	// CheckType MCCheckType = CheckType::Matrix;
	// bool useCallback = false;
	// uint maxNumberConstr = 4;
	// vector<uint> roundOrder = roundOrderOutsideIn(rMax);

	vector<string> listSboxModel({"Hull", "Simple"});
	vector<string> listLinModel({"ZR", "CX", "Simple", "Lbox"});
	vector<bool> listUseCallback({true});
	vector<vector<uint>> listroundOrderIfNeeded({roundOrderOutsideIn(rMax), roundOrderAscend(rMax), roundOrderDescend(rMax)});
	auto const allTest = allInOutPairsSKINNY64(0,0,0);

	uint ctrParameterSet = 0;
	uint firstParameterSet = 0;

	for(auto const & useCallback : listUseCallback){
	 for(auto const & linModel : listLinModel){
	  vector<CheckType> listMCCheckType({CheckType::Ineq});
	  if(linModel == "Lbox" || linModel == "ZR") //No need for checks on MC
	   listMCCheckType = vector<CheckType>({CheckType::None});
	  for(auto const & MCCheckType : listMCCheckType){
	  	for(auto const & sboxModel : listSboxModel){
	    vector<CheckType> listSboxCheckType({CheckType::Ineq});
	    if(sboxModel == "Hull" || sboxModel == "QM") //No need for checks on Sbox
	     listSboxCheckType = vector<CheckType>({CheckType::None});
	    for(auto const & sboxCheckType : listSboxCheckType){
	      vector<uint> listMaxNbConstr({1});
	      if(!(sboxCheckType == CheckType::None && MCCheckType == CheckType::None && SSBCheckType == CheckType::None))
	        listMaxNbConstr = vector<uint>({500,100,10,1});
	      for(auto const & maxNumberConstr : listMaxNbConstr){
	      	vector<vector<uint>> listroundOrder(1);
	      	if(!(sboxCheckType == CheckType::None && MCCheckType == CheckType::None && SSBCheckType == CheckType::None))
	      		listroundOrder = vector<vector<uint>>(listroundOrderIfNeeded);
	      	for(auto const & roundOrder : listroundOrder){
	      	 if(ctrParameterSet >= firstParameterSet){

	      	 	 globalCtrCheckSkinny = 0;
				 globalCtrSolSkinny = 0;
	      	 	 globalCtrConstrSboxSkinny = 0;
				 globalCtrConstrLinSkinny = 0;
				 globalCtrConstrSSBSkinny = 0;

		       	 modelInfo MI(rMax,sboxModel,linModel,arxModel,firstSbox,lastLin,smartLin,startingRound,sboxCheckType,MCCheckType,SSBCheckType,ARXCheckType,useCallback,maxNumberConstr,roundOrder);
		       	 cout << "Parameter set " << ctrParameterSet << endl;
		       	 if(useCallback) cout << "Using callback" << endl;
				 else cout << "Using iterative model" << endl;
				 cout << "roundOrder : ";
				 for(auto const & tmpr : roundOrder) cout << tmpr << " ";
				 cout << endl;
				 cout << "linModel : " << linModel << " - MCCheckType : " << MCCheckType << endl;
				 cout << "sboxModel : " << sboxModel << " - sboxCheckType : " << sboxCheckType << endl;
				 cout << "SSBCheckType : " << SSBCheckType << endl;
				 cout << "maxNumberConstr : " << maxNumberConstr << endl;

				 auto start = chrono::high_resolution_clock::now();
				 auto haveDist = checkForDistinguisher(allTest, existTrailSkinny64, MI, false);
				 auto end = chrono::high_resolution_clock::now();
	 
				 cout << "Time : " << chrono::duration<double>(end - start).count() << endl;
				 if(haveDist)cout << "Distinguisher !!!" << endl << endl;
				 else cout << "No distinguisher..." << endl << endl;
				 cout << globalCtrCheckSkinny << " input/output pairs checked" << endl;
				 cout << globalCtrSolSkinny << " solutions examined" << endl;
				 cout << globalCtrConstrSboxSkinny << " constraints added for the sbox" << endl;
				 cout << globalCtrConstrLinSkinny  << " constraints added for the lin layer" << endl;
				 cout << globalCtrConstrSSBSkinny << " constraints added for the SSB" << endl;
			 }
			 ctrParameterSet++;
	}}}}}}}
}


void searchDistinguisherSkinny64_11r(){

// SB |
// MC | Precalc
// SB |

// MC |
// SB |
// MC |
// SB |
// MC |
// SB |
// MC | MILP
// SB |
// MC |
// SB |
// MC |
// SB |
// MC |
// SB |
// MC |

// SB |
// MC | precalc
// SB |

	string arxModel = "ARX";
	unsigned int startingRound = 1;
	CheckType ARXCheckType = CheckType::None;

	CheckType SSBCheckType = CheckType::Ineq;
	uint rMax = 8; //1 partial round + 7 rounds -> 11 rounds with precalc
	bool smartLin = true;
	bool firstSbox = false;
	bool lastLin = true;

	string sboxModel = "Hull";
	string linModel = "Lbox";
	CheckType sboxCheckType = CheckType::None;
	CheckType MCCheckType = CheckType::None;
	bool useCallback = true;
	uint maxNumberConstr = 100;
	vector<uint> roundOrder = roundOrderAscend(rMax);

	globalCtrCheckSkinny = 0;
	globalCtrSolSkinny = 0;
	globalCtrConstrSboxSkinny = 0;
	globalCtrConstrLinSkinny = 0;
	globalCtrConstrSSBSkinny = 0;

	auto globalStart = chrono::high_resolution_clock::now();
	for(uint i = 0; i < 1; i++){
		for(uint j = 0; j < 4; j++){
			cout << "Skinny64 Search i = " << i << " j = " << j << endl;
			//Get the list of input/output to test
			auto const allTest = allInOutPairsSKINNY64(i,j,0);

			//Generate the modelInfo object
			modelInfo MI(rMax,sboxModel,linModel,arxModel,firstSbox,lastLin,smartLin,startingRound,sboxCheckType,MCCheckType,SSBCheckType,ARXCheckType,useCallback,maxNumberConstr,roundOrder);

			//Check if we can find a distinguisher
			auto start = chrono::high_resolution_clock::now();
			auto haveDist = checkForDistinguisher(allTest, existTrailSkinny64, MI,false);
			auto end = chrono::high_resolution_clock::now();
			cout << "Time : " << chrono::duration<double>(end - start).count() << endl;
			if(haveDist){
				cout << "Distinguisher !!!" << endl;
			}
			else
				cout << "No distinguisher..." << endl;
			cout << endl;
		}
	}
	auto globalEnd = chrono::high_resolution_clock::now();
	cout << "Total time : " << chrono::duration<double>(globalEnd - globalStart).count() << endl;
	cout << globalCtrCheckSkinny << " input/output pairs checked" << endl;
	cout << globalCtrSolSkinny << " solutions examined" << endl;
	cout << globalCtrConstrSboxSkinny << " constraints added for the sbox" << endl;
	cout << globalCtrConstrLinSkinny  << " constraints added for the lin layer" << endl;
	cout << globalCtrConstrSSBSkinny << " constraints added for the SSB" << endl;
}


void searchDistinguisherSkinny64_12r(){

	string arxModel = "ARX";
	unsigned int startingRound = 1;
	CheckType ARXCheckType = CheckType::None;

	CheckType SSBCheckType = CheckType::Ineq;
	uint rMax = 9; //1 partial round + 8 rounds -> 12 rounds with precalc
	bool smartLin = true;
	bool firstSbox = false;
	bool lastLin = true;

	string sboxModel = "Hull";
	string linModel = "Lbox";
	CheckType sboxCheckType = CheckType::None;
	CheckType MCCheckType = CheckType::None;
	bool useCallback = true;
	uint maxNumberConstr = 500;
	vector<uint> roundOrder = roundOrderOutsideIn(rMax);

	globalCtrCheckSkinny = 0;
	globalCtrSolSkinny = 0;
	globalCtrConstrSboxSkinny = 0;
	globalCtrConstrLinSkinny = 0;
	globalCtrConstrSSBSkinny = 0;
	auto globalStart = chrono::high_resolution_clock::now();
	for(uint i = 0; i < 4; i++){
		for(uint j = 0; j < 4; j++){
			cout << "Skinny64 Search i = " << i << " j = " << j << endl;
			//Get the list of input/output to test
			auto const allTest = allInOutPairsSKINNY64(i,j,0);

			//Generate the modelInfo object
			modelInfo MI(rMax,sboxModel,linModel,arxModel,firstSbox,lastLin,smartLin,startingRound,sboxCheckType,MCCheckType,SSBCheckType,ARXCheckType,useCallback,maxNumberConstr,roundOrder);

			//Check if we can find a distinguisher
			auto start = chrono::high_resolution_clock::now();
			auto haveDist = checkForDistinguisher(allTest, existTrailSkinny64, MI,false);
			auto end = chrono::high_resolution_clock::now();
			cout << "Time : " << chrono::duration<double>(end - start).count() << endl;
			if(haveDist){
				cout << "Distinguisher !!!" << endl;
			}
			else
				cout << "No distinguisher..." << endl;
			cout << endl;
		}
	}
	auto globalEnd = chrono::high_resolution_clock::now();
	cout << "Total time : " << chrono::duration<double>(globalEnd - globalStart).count() << endl;
	cout << globalCtrCheckSkinny << " input/output pairs checked" << endl;
	cout << globalCtrSolSkinny << " solutions examined" << endl;
	cout << globalCtrConstrSboxSkinny << " constraints added for the sbox" << endl;
	cout << globalCtrConstrLinSkinny  << " constraints added for the lin layer" << endl;
	cout << globalCtrConstrSSBSkinny << " constraints added for the SSB" << endl;
}
