#include "led.hpp"
#include "DivTable.hpp"

using namespace std;

#define DISPLAY_GEN_OUTPUT false
#define DISPLAY_SOLVER_OUTPUT false	
#define DISPLAY_NUMBER_ADDED_CONSTRAINT false
#define TIME_LIMIT_LED 30
#define QUICK_LAZY_COUNT_LED true

vector<vector<pair<vector<uint8_t>, vector<uint8_t>>>> allInOutPairsLED(unsigned pos_in, unsigned pos_out, unsigned R) { // 0 <= i < 4  - 0 <= j < 4
  vector<uint16_t> sbox = {12, 5, 6, 11, 9, 0, 10, 13, 3, 14, 15, 8, 4, 7, 1, 2};
  vector<uint16_t> linear = {
    0b1000100000010100,
    0b1001100100101100,
    0b0010001001001001,
    0b0100010010000010,
    0b1100010111000010,
    0b0101111001010110,
    0b1011110110111100,
    0b0110101001101001,
    0b0011101011101011,
    0b0100111100111101,
    0b1000111001111010,
    0b0001110111110101,
    0b1011111110001000,
    0b1101000110011001,
    0b1010001100100010,
    0b0101011101000100
  };

  vector<uint8_t> perm = {0, 13, 10, 7, 4, 1, 14, 11, 8, 5, 2, 15, 12, 9, 6, 3};
  vector<uint8_t> perm_inv (16);
  for (unsigned i = 0; i < 16; ++i) perm_inv[perm[i]] = i;

  vector<uint8_t> permBit({0,1,2,3,52,53,54,55,40,41,42,43,28,29,30,31,16,17,18,19,4,5,6,7,56,57,58,59,44,45,46,47,32,33,34,35,20,21,22,23,8,9,10,11,60,61,62,63,48,49,50,51,36,37,38,39,24,25,26,27,12,13,14,15});
  vector<uint8_t> permBit_inv(64);
    for(uint i = 0; i < 64; i++)
  	  permBit_inv[permBit[i]] = i;
  

  auto anf_s = getANF_Parallel_S(4, sbox);
  auto anf_l = convertLtoPolynome(linear);

  auto anf_los = composition(anf_l, anf_s);
  auto anf_solos = composition(anf_s, anf_los);

  vector<vector<Polynome>> anf (4);
  anf[2] = anf_solos;
  anf[3] = anf_solos;

  uint8_t cst0 [4] = {(64>>4)&0xf, ((64>>4)&0xf)^1, (64 & 0xf)^2, (64 & 0xf)^3};
  anf[0] = anf_los;
  for (unsigned i = 0; i < 4; ++i) {
    for (unsigned j = 0; j < 4; ++j) if (((cst0[i] >> j) & 1) != 0) anf[0][4*i + j] += Monome(0);
  }
  anf[0] = composition(anf_s, anf[0]);
  anf[1] = anf_los;
  for (unsigned i = 0; i < 2; ++i) {
    for (unsigned j = 0; j < 3; ++j) {
      anf[1][4*i + j] += Monome(1u << (16 + (3*i + j)));
      anf[1][4*(i+2) + j] += Monome(1u << (16 + (3*i + j)));
    }
  }
  anf[1] = composition(anf_s, anf[1]);

  auto anf_sol = composition(anf_s, anf_l);

  auto v_in = loadTransitionsIN(anf[pos_in]);
  auto v_out = (R%2 == 0) ? loadTransitionsOUT(anf[pos_out]) : loadTransitionsOUT(anf_sol);

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

GRBModel getModelLED(modelInfo & MI){
	/*
	Get the model for LED over #rMax round
	#lastLin defines if the linear layer is applied on the last round
	#firstSbox defines if the first Sbox layer is applied

	The #sboxModel parameter defines how the Sbox is modelized :
	- "Hull" use the Convex Hull technique
	- "Simple" use the simplified constraint with PWL

	The #linModel parameter defines how the linear layer is modelized :
	- "CX" modelize the linear layer using the classical copy+xor technique
	- "QM" use the Quin-McCluskey algorithm as in Abdelkhalek,Sasaki,Todo,Tolba,Youssef
	- "Simple" modelize the linear layer with the simplified constraint w(x) = w(y)
*/
    string modelName = "LED_"+to_string(MI.rMax)+"r_"+MI.sboxModel+"_"+MI.linModel;
    if(MI.lastLin) modelName += "_lastLin";
    if(!MI.firstSbox) modelName += "_noFirstSbox";
    modelName += ".mps";

    if(!fileExist("./models/"+modelName)){
    	cout << "Generating the model..." << endl;
        string cmd = "sage ./modelGeneration/genLED.sage";
        cmd += " -r " + to_string(MI.rMax);
        cmd += " -s " + MI.sboxModel;
        cmd += " -m " + MI.linModel;
        if(MI.lastLin) cmd += " --lastLin True";
        else cmd += " --lastLin False";
        if(MI.firstSbox) cmd += " --firstSbox True";
        else cmd += " --firstSbox False";
        auto sysstdout = exec(cmd.c_str());
        if(DISPLAY_GEN_OUTPUT) cout << sysstdout << endl;
        cout << "Model generated" << endl;
    }

    if(!DISPLAY_SOLVER_OUTPUT)
        MI.gurobiEnv.set(GRB_IntParam_OutputFlag,0);

    GRBModel m(MI.gurobiEnv, "./models/"+modelName);
    return m;
}

modelData getVariablesLED(GRBModel & m,
                          modelInfo const & MI){
     /*
    Read the model and return the corresponding blocks of variables
    - #MI defines various parameters for how the trails are checked (see aux_function.hpp)
    */

    modelData MD;
    MD.rMax = MI.rMax;
    MD.firstSbox = MI.firstSbox;
    MD.lastLin = MI.lastLin;

    MD.allTable = vector<vector<vector<uint32_t>>>(5);
    MD.allIneq = vector<vector<vector<int>>>(5);
    MD.allMatrix = vector<Matrix>(1);

    if(MI.sboxCheckType == CheckType::Table) //Sbox table
        MD.allTable[0] = vector<vector<uint32_t>>({{0},{1,2,4,8},{1,2,4,8},{2,4,8},{1,2,4,8},{2,4,8},{1,2,8},{2,8},{1,2,4,8},{2,4,8},{2,4,8},{2,4,8},{2,4,8},{2,4,8},{5,11,14},{15}});
    else if(MI.sboxCheckType == CheckType::Ineq)
        MD.allIneq[0] = vector<vector<int>>({{1, 1, 1, 1, -1, -1, -1, -1, 0, 1},
                                              {-4, -2, -2, -2, -3, 1, 4, 1, 7, 1},
                                              {-2, 0, 0, 0, 2, -1, -1, -1, 3, 1},
                                              {0, -1, -1, -2, 2, 3, 3, 3, 0, 1},
                                              {1, 1, 1, 1, -2, 1, -2, -2, 1, 1},
                                              {1, 0, 0, 0, -1, -2, -1, 1, 2, 1},
                                              {-2, -1, -1, 0, -1, 1, 0, 1, 3, 1},
                                              {0, 0, 0, 0, 1, -1, 1, -1, 1, 1},
                                              {0, -2, -2, 0, 2, 1, -1, 1, 3, 1},
                                              {-1, 0, -1, 0, 1, 2, 2, 2, 0, 1}});

    if(MI.MCCheckType == CheckType::Matrix)
    	MD.allMatrix[0] = Matrix("./checkData/matrixLED.bin");
	else if(MI.MCCheckType == CheckType::Table){
    	cerr << "Error : Table Check for LED MC not implemented" << endl;
    	exit(1);
        // MD.allTable[1] = //Table for MC
    }
    else if(MI.MCCheckType == CheckType::Ineq){
        MD.allIneq[1]   = readIneqFromFile("./checkData/LEDMC_ineq.bin");
    }

    if(MI.SSBCheckType == CheckType::Table){
        MD.allTable[2] = readDivTableFromFile("./checkData/LEDSSB0_table.bin");
        MD.allTable[3] = readDivTableFromFile("./checkData/LEDSSB1_table.bin");
        MD.allTable[4] = readDivTableFromFile("./checkData/LEDSSB2_table.bin");
    }
    else if(MI.SSBCheckType == CheckType::Ineq){
        MD.allIneq[2]   = readIneqFromFile("./checkData/LEDSSB0_ineq.bin");
        MD.allIneq[3]   = readIneqFromFile("./checkData/LEDSSB1_ineq.bin");
        MD.allIneq[4]   = readIneqFromFile("./checkData/LEDSSB2_ineq.bin");
    }

    auto & tableSbox = MD.allTable[0];
    auto & tableMC   = MD.allTable[1];

    auto & ineqSbox  = MD.allIneq[0];
    auto & ineqMC    = MD.allIneq[1];

    auto & MCmatrix = MD.allMatrix[0];

    auto & tableSSB0  = MD.allTable[2];
    auto & tableSSB1  = MD.allTable[3];
    auto & tableSSB2  = MD.allTable[4];

    auto & ineqSSB0   = MD.allIneq[2];
    auto & ineqSSB1   = MD.allIneq[3];
    auto & ineqSSB2   = MD.allIneq[4];

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

    //SSB variables
    //Input/Output nibbles for each SSB are as follow
    //{0 ,5 ,10,15} -> {0,1,2,3}
    //{4 ,9 ,14,3} -> {4,5,6,7}
    //{8 ,13,2 ,7} -> {8,9,10,11}
    //{12,1 ,6 ,11} -> {12,13,14,15}
    static const uint8_t indexIn[4][4] = {{0,5,10,15},{4,9,14,3},{8,13,2,7},{12,1,6,11}};
    static const uint8_t indexOut[4][4] = {{0,1,2,3},{4,5,6,7},{8,9,10,11},{12,13,14,15}};
    roundLastLin--; //Last SSB starts one round before the last linear layer

    if(MI.SSBCheckType != CheckType::None){
        for(auto const & r : MI.roundOrder){
        	if(r < roundLastLin && r >= roundFirstSbox){
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
		            else if(i == 1){
		                if(MI.SSBCheckType == CheckType::Table)
		                    MD.allVarsTable.emplace_back(move(in), move(out), tableSSB1,ConstrType::SSB);
		                else if(MI.SSBCheckType == CheckType::Ineq)
		                    MD.allVarsIneq.emplace_back(move(in), move(out), ineqSSB1,ConstrType::SSB);
		            }
		            else{
		            	if(MI.SSBCheckType == CheckType::Table)
		                    MD.allVarsTable.emplace_back(move(in), move(out), tableSSB2,ConstrType::SSB);
		                else if(MI.SSBCheckType == CheckType::Ineq)
		                    MD.allVarsIneq.emplace_back(move(in), move(out), ineqSSB2,ConstrType::SSB);
		            }
	            }
	        }
        }
    }

    return MD;
}


uint64_t globalCtrCheckLED = 0;
uint64_t globalCtrSolLED = 0;
uint64_t globalCtrConstrSboxLED = 0;
uint64_t globalCtrConstrLinLED = 0;
uint64_t globalCtrConstrSSBLED = 0;
bool existTrailLED(GRBModel & m,
                   modelData const & MD,
                   vector<uint8_t> const & input,
                   vector<uint8_t> const & output,
                   bool const useCallback,
                   uint const maxNumberConstr){
/*
    Return true if there is a trail from #input to #output using #model for LED
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

    // GRBLinExpr outWeight = 0;
    // for(uint i = 0; i < 64; i++)
    // 	outWeight += m.getVarByName(outputvarPrefix+to_string(i));
    // m.addConstr(outWeight == 1);

    m.update();

    if(useCallback){
        m.set(GRB_IntParam_LazyConstraints, 1);
        m.set(GRB_DoubleParam_TimeLimit, TIME_LIMIT_LED);
        CustomCallback cb(MD.allVarsTable, MD.allVarsIneq, MD.allVarsMatrix, maxNumberConstr);
        m.setCallback(&cb);

        m.update();
        m.optimize();
        if(DISPLAY_NUMBER_ADDED_CONSTRAINT){
            cout << cb.ctrSol << " solutions examined, ";
            cout << cb.ctrGlobalConstr << " constraints added" << endl;
        }
        globalCtrCheckLED++;
        globalCtrSolLED += cb.ctrSol;
        globalCtrConstrSboxLED += cb.ctrConstrSbox;
		globalCtrConstrLinLED += cb.ctrConstrLin;
		globalCtrConstrSSBLED += cb.ctrConstrSSB;

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

bool existTrailLED(modelInfo & MI,
                   vector<uint8_t> const & input,
                   vector<uint8_t> const & output){
/*
    Return true if there is a trail from #input to #output for LED
    #MI should contain the necessary information to generate the model and parametrize the solving
*/

    GRBModel m = getModelLED(MI);
    auto MD = getVariablesLED(m, MI);
    return existTrailLED(m,MD,input,output,MI.useCallback,MI.maxNumberConstr);
}

bool existTrailLEDAllInOne(modelInfo & MI, 
						   vector<uint8_t> const & input){
	//Return true if there is a trail from #input to any vector of weight 1
	GRBModel m = getModelLED(MI);
    auto MD = getVariablesLED(m, MI);
    return existTrailLEDAllInOne(m,MD,input,MI.useCallback,MI.maxNumberConstr);
}

bool existTrailLEDAllInOne(GRBModel & m,
                   modelData const & MD,
                   vector<uint8_t> const & input,
                   bool const useCallback,
                   uint const maxNumberConstr){
/*
    //Return true if there is a trail from #input to any vector of weight 1 using #model for LED
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
    }

    GRBLinExpr outWeight = 0;
    for(uint i = 0; i < 64; i++)
    	outWeight += m.getVarByName(outputvarPrefix+to_string(i));
    m.addConstr(outWeight == 1);

    m.update();

    if(useCallback){
        m.set(GRB_IntParam_LazyConstraints, 1);
        m.set(GRB_DoubleParam_TimeLimit, TIME_LIMIT_LED);
        CustomCallback cb(MD.allVarsTable, MD.allVarsIneq, MD.allVarsMatrix, maxNumberConstr);
        m.setCallback(&cb);

        m.update();
        m.optimize();
        if(DISPLAY_NUMBER_ADDED_CONSTRAINT){
            cout << cb.ctrSol << " solutions examined, ";
            cout << cb.ctrGlobalConstr << " constraints added" << endl;
        }
        globalCtrCheckLED++;
        globalCtrSolLED += cb.ctrSol;
        globalCtrConstrSboxLED += cb.ctrConstrSbox;
		globalCtrConstrLinLED += cb.ctrConstrLin;
		globalCtrConstrSSBLED += cb.ctrConstrSSB;

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

void createMatrixFileLED(){
	//Create the matrix file ./checkData/matrixLED.bin

	vector<vector<uint8_t>> val({{0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1},
								 {0,0,1,1,0,1,0,0,1,0,0,1,1,0,0,1},
								 {1,0,0,1,0,0,1,0,0,1,0,0,0,1,0,0},
								 {0,1,0,0,0,0,0,1,0,0,1,0,0,0,1,0},
								 {0,1,0,0,0,0,1,1,1,0,1,0,0,0,1,1},
								 {0,1,1,0,1,0,1,0,0,1,1,1,1,0,1,0},
								 {0,0,1,1,1,1,0,1,1,0,1,1,1,1,0,1},
								 {1,0,0,1,0,1,1,0,0,1,0,1,0,1,1,0},
								 {1,1,0,1,0,1,1,1,0,1,0,1,1,0,0,1},
								 {1,0,1,1,1,1,0,0,1,1,1,1,1,1,0,1},
								 {0,1,0,1,1,1,1,0,0,1,1,1,0,1,1,0},
								 {1,0,1,0,1,1,1,1,1,0,1,1,0,0,1,1},
								 {0,0,0,1,0,0,0,1,1,1,1,1,1,1,0,1},
								 {1,0,0,1,1,0,0,1,1,0,0,0,1,0,1,1},
								 {0,1,0,0,0,1,0,0,1,1,0,0,0,1,0,1},
								 {0,0,1,0,0,0,1,0,1,1,1,0,1,0,1,0}});
	Matrix M(16,16);
	for(uint i = 0; i < 16; i++){
		for(uint j = 0; j < 16; j++)
			M.set(i,j,val[i][j]);
	}

	M.saveToFile("./checkData/matrixLED.bin");
}

void computeSboxMatrixLED(){
	// vector<vector<uint8_t>> Mb({{0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1},
	// 							{0,0,1,1,0,1,0,0,1,0,0,1,1,0,0,1},
	// 							{1,0,0,1,0,0,1,0,0,1,0,0,0,1,0,0},
	// 							{0,1,0,0,0,0,0,1,0,0,1,0,0,0,1,0},
	// 							{0,1,0,0,0,0,1,1,1,0,1,0,0,0,1,1},
	// 							{0,1,1,0,1,0,1,0,0,1,1,1,1,0,1,0},
	// 							{0,0,1,1,1,1,0,1,1,0,1,1,1,1,0,1},
	// 							{1,0,0,1,0,1,1,0,0,1,0,1,0,1,1,0},
	// 							{1,1,0,1,0,1,1,1,0,1,0,1,1,0,0,1},
	// 							{1,0,1,1,1,1,0,0,1,1,1,1,1,1,0,1},
	// 							{0,1,0,1,1,1,1,0,0,1,1,1,0,1,1,0},
	// 							{1,0,1,0,1,1,1,1,1,0,1,1,0,0,1,1},
	// 							{0,0,0,1,0,0,0,1,1,1,1,1,1,1,0,1},
	// 							{1,0,0,1,1,0,0,1,1,0,0,0,1,0,1,1},
	// 							{0,1,0,0,0,1,0,0,1,1,0,0,0,1,0,1},
	// 							{0,0,1,0,0,0,1,0,1,1,1,0,1,0,1,0}});

	// vector<uint16_t> M(16);
	// for(uint i = 0; i < 16; i++){
	// 	uint16_t x = 0;
	// 	for(uint j = 0; j < 16; j++){
	// 		if(Mb[i][j] == 1)
	// 			x += (1 << j);
	// 	}
	// 	M[i] = x;
	// }
	vector<uint16_t> M({34836, 39212, 8777, 17538, 50626, 24150, 48572, 27241, 39659, 48957, 28282, 52725, 49032, 53657, 41762, 22340});

	uint32_t bound = (1ul << 16);
	vector<uint32_t> S(bound,0);
	for(uint32_t x = 0; x < bound; x++)
		S[x] = apply16bitMC16bitState(M,x);

	genAndSaveTableIneq(vector<vector<uint32_t>>({S}),16,16,"./LEDMC",true);
}

void checkNormalDistinguisherLED(){
	/*
		Search for a normal (no linear combination) distinguisher over 7 rounds of LED 
	*/

	// uint rMax = 7;
	// vector<uint8_t> input({1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1});

	uint rMax = 6;
	vector<uint8_t> input({1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1});

	//These parameters have no effect for LED but are required for the modelInfo object
	string arxModel = "ARX";
	unsigned int startingRound = 1;
	CheckType ARXCheckType = CheckType::None;

	//Fixed parameters
	CheckType SSBCheckType = CheckType::None; //no need to check SSB for this one
	bool smartLin = true; //probably not useful to not have it honestly
	bool firstSbox = true;
	bool lastLin = false;

	//Variable parameters
	// string sboxModel = "QM";
	// string linModel = "CX";
	// CheckType sboxCheckType = CheckType::None;
	// CheckType MCCheckType = CheckType::Matrix;
	// bool useCallback = false;
	// uint maxNumberConstr = 4;
	// vector<uint> roundOrder = roundOrderOutsideIn(rMax);


	vector<string> listSboxModel({"Hull", "Simple"});
	vector<string> listLinModel({"QM", "Simple", "CX"});
	vector<bool> listUseCallback({true});
	vector<vector<uint>> listroundOrderIfNeeded({roundOrderOutsideIn(rMax), roundOrderAscend(rMax), roundOrderDescend(rMax)});
	// vector<uint> listMaxNbConstr({500,100,10,1});

	uint ctrParameterSet = 0;
	uint firstParameterSet = 0;

	for(auto const & useCallback : listUseCallback){
	 for(auto const & linModel : listLinModel){
	  vector<CheckType> listMCCheckType({CheckType::Ineq});
	  if(linModel == "Lbox" || linModel == "ZR" || linModel == "QM") //No need for checks on MC
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

	       		globalCtrCheckLED = 0;
				globalCtrSolLED = 0;
				globalCtrConstrSboxLED = 0;
				globalCtrConstrLinLED = 0;
				globalCtrConstrSSBLED = 0;

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
					balancedBits[i] = !existTrailLED(MI,input,output);
					auto endBit = chrono::high_resolution_clock::now();
					if(balancedBits[i]){
						cout << i << " balanced" << endl;
						ctrBalanced++;
						durBalanced += endBit - startBit;
					}
				}

				//sum(output) == 1
				// bool trailToWeight1 = existTrailLEDAllInOne(MI,input);
				// int ctrBalanced = -1;
				// if(!trailToWeight1)
				// 	ctrBalanced = 64;

				auto end = chrono::high_resolution_clock::now();
				cout << ctrBalanced << " balanced bits" << endl;
				cout << chrono::duration<double>(end - start).count() << " seconds, ";
				cout << "Time for only balanced bits : " << chrono::duration<double>(durBalanced).count() << " second" << endl;
				cout << endl;
				cout << globalCtrCheckLED << " input/output pairs checked" << endl;
				cout << globalCtrSolLED << " solutions examined" << endl;
				cout << globalCtrConstrSboxLED << " constraints added for the sbox" << endl;
				cout << globalCtrConstrLinLED  << " constraints added for the lin layer" << endl;
				cout << globalCtrConstrSSBLED << " constraints added for the SSB" << endl;
			}
			ctrParameterSet++;
	}}}}}}}
}

void timingTestSearchNewLED(){
	/*
	*/

//8r
// SB | 
// MC | precalc
// SB | 

// MC | 
// SB | 
// MC | 
// SB | 
// MC | MILP (1 partial + 4 full)
// SB | 
// MC | 
// SB | 
// MC | 

// SB | 
// MC | precalc
// SB | 

	//These parameters have no effect for Midori64 but are required for the modelInfo object
	string arxModel = "ARX";
	unsigned int startingRound = 1;
	CheckType ARXCheckType = CheckType::None;
	bool smartLin = true;

	//Fixed parameters
	CheckType SSBCheckType = CheckType::Ineq; 
	uint rMax = 5;
	bool firstSbox = false;
	bool lastLin = true;
	// auto const allTest = allInOutPairsLED(0,0,0);


	//Variable parameters
	// string sboxModel = "QM";
	// string linModel = "CX";
	// CheckType sboxCheckType = CheckType::None;
	// CheckType MCCheckType = CheckType::Matrix;
	// bool useCallback = false;
	// uint maxNumberConstr = 4;
	// vector<uint> roundOrder = roundOrderOutsideIn(rMax);

	vector<bool> listUseCallback({true});
	vector<uint> listMaxNbConstr({500,100,10,1});
	vector<vector<uint>> listroundOrder({roundOrderOutsideIn(rMax), roundOrderAscend(rMax), roundOrderDescend(rMax)});
	vector<string> listSboxModel({"Hull", "Simple"});
	// vector<string> listLinModel({"QM", "CX", "Simple"});
	vector<string> listLinModel({"QM"});
	
	for(auto const & useCallback : listUseCallback){
	 for(auto const & maxNumberConstr : listMaxNbConstr){
	   for(auto const & roundOrder : listroundOrder){
	    for(auto const & linModel : listLinModel){
	     vector<CheckType> listMCCheckType({CheckType::Ineq});
	     if(linModel == "QM") //No need for checks on MC
	      listMCCheckType = vector<CheckType>({CheckType::None});
	     for(auto const & MCCheckType : listMCCheckType){
	     	for(auto const & sboxModel : listSboxModel){
	       vector<CheckType> listSboxCheckType({CheckType::Ineq});
	       if(sboxModel == "Hull") //No need for checks on Sbox
	        listSboxCheckType = vector<CheckType>({CheckType::None});
	       for(auto const & sboxCheckType : listSboxCheckType){

	       	 globalCtrCheckLED = 0;
			 globalCtrSolLED = 0;
	      	 globalCtrConstrSboxLED = 0;
			 globalCtrConstrLinLED = 0;
			 globalCtrConstrSSBLED = 0;
	        
	       	 modelInfo MI(rMax,sboxModel,linModel,arxModel,firstSbox,lastLin,smartLin,startingRound,sboxCheckType,MCCheckType,SSBCheckType,ARXCheckType,useCallback,maxNumberConstr,roundOrder);
	       	 if(useCallback) cout << "Using callback" << endl;
			 else cout << "Using iterative model" << endl;
			 cout << "roundOrder : ";
			 for(auto const & tmpr : roundOrder) cout << tmpr << " ";
			 cout << endl;
			 cout << "linModel : " << linModel << " - MCCheckType : " << MCCheckType << endl;
			 cout << "sboxModel : " << sboxModel << " - sboxCheckType : " << sboxCheckType << endl;
			 cout << "maxNumberConstr : " << maxNumberConstr << endl;

			 auto start = chrono::high_resolution_clock::now();
			 // auto haveDist = checkForDistinguisher(allTest, existTrailLED, MI);
			 auto haveDist = iterativeSearchSetPair(0,0,0,existTrailLED, MI, true, 100);
			 auto end = chrono::high_resolution_clock::now();

			 cout << "Time : " << chrono::duration<double>(end - start).count() << endl;
			 if(haveDist)cout << "Distinguisher !!!" << endl;
			 else cout << "No distinguisher..." << endl;
			 cout << endl;
			 cout << globalCtrCheckLED << " input/output pairs checked" << endl;
			 cout << globalCtrSolLED << " solutions examined" << endl;
			 cout << globalCtrConstrSboxLED << " constraints added for the sbox" << endl;
			 cout << globalCtrConstrLinLED  << " constraints added for the lin layer" << endl;
			 cout << globalCtrConstrSSBLED << " constraints added for the SSB" << endl;
	}}}}}}}
}

void searchNewDistLED(){

	string arxModel = "ARX";
	unsigned int startingRound = 1;
	CheckType ARXCheckType = CheckType::None;
	bool smartLin = true;

	//Fixed parameters
	CheckType SSBCheckType = CheckType::Ineq; 
	uint rMax = 5;
	bool firstSbox = false;
	bool lastLin = true;

	string sboxModel = "Hull";
	string linModel = "Simple";
	CheckType sboxCheckType = CheckType::None;
	CheckType MCCheckType = CheckType::Ineq;
	bool useCallback = true;
	uint maxNumberConstr = 10;
	vector<uint> roundOrder = roundOrderDescend(rMax);

	globalCtrCheckLED = 0;
	globalCtrSolLED = 0;
	globalCtrConstrSboxLED = 0;
	globalCtrConstrLinLED = 0;
	globalCtrConstrSSBLED = 0;

	auto globalStart = chrono::high_resolution_clock::now();
	for(uint i = 0; i < 4; i++){
		for(uint j = 0; j < 4; j++){
			cout << "LED Search i = " << i << " j = " << j << endl;

			//Generate the modelInfo object
			modelInfo MI(rMax,sboxModel,linModel,arxModel,firstSbox,lastLin,smartLin,startingRound,sboxCheckType,MCCheckType,SSBCheckType,ARXCheckType,useCallback,maxNumberConstr,roundOrder);
			//Check if we can find a distinguisher
			auto start = chrono::high_resolution_clock::now();
			auto haveDist = iterativeSearchSetPair(i,j,0,existTrailLED, MI, false);
			auto end = chrono::high_resolution_clock::now();

			cout << "Time : " << chrono::duration<double>(end - start).count() << endl;
			if(haveDist)cout << "Distinguisher !!!" << endl;
			else cout << "No distinguisher..." << endl;
			cout << endl;
			if(QUICK_LAZY_COUNT_LED)
				break;
		}
		if(QUICK_LAZY_COUNT_LED)
			break;
	}
	auto globalEnd = chrono::high_resolution_clock::now();
	cout << "Total time : " << chrono::duration<double>(globalEnd - globalStart).count() << endl;

	cout << globalCtrCheckLED << " input/output pairs checked" << endl;
	cout << globalCtrSolLED << " solutions examined" << endl;
	cout << globalCtrConstrSboxLED << " constraints added for the sbox" << endl;
	cout << globalCtrConstrLinLED  << " constraints added for the lin layer" << endl;
	cout << globalCtrConstrSSBLED << " constraints added for the SSB" << endl;

}

//Because there are a lot of input/output to go through, we cannot store them all, so go through them iteratively
bool iterativeSearchSetPair(unsigned pos_in, 
							unsigned pos_out, 
							unsigned R,
							t_existTrail existTrail,
							modelInfo & MI,
							bool const fastSearch,
							uint const boundNbSet) {
// boundNbSet = 0 : test all sets
// otherwise test only boundNbSet

	if(boundNbSet > 0)
		cout << "Warning, starting a search with a bounded number of set tested !!" << endl;

	// 0 <= i < 4  - 0 <= j < 4
  vector<uint16_t> sbox = {12, 5, 6, 11, 9, 0, 10, 13, 3, 14, 15, 8, 4, 7, 1, 2};
  vector<uint16_t> linear = {
    0b1000100000010100,
    0b1001100100101100,
    0b0010001001001001,
    0b0100010010000010,
    0b1100010111000010,
    0b0101111001010110,
    0b1011110110111100,
    0b0110101001101001,
    0b0011101011101011,
    0b0100111100111101,
    0b1000111001111010,
    0b0001110111110101,
    0b1011111110001000,
    0b1101000110011001,
    0b1010001100100010,
    0b0101011101000100
  };

  vector<uint8_t> perm = {0, 13, 10, 7, 4, 1, 14, 11, 8, 5, 2, 15, 12, 9, 6, 3};
  vector<uint8_t> perm_inv (16);
  for (unsigned i = 0; i < 16; ++i) perm_inv[perm[i]] = i;

  vector<uint8_t> permBit({0,1,2,3,52,53,54,55,40,41,42,43,28,29,30,31,16,17,18,19,4,5,6,7,56,57,58,59,44,45,46,47,32,33,34,35,20,21,22,23,8,9,10,11,60,61,62,63,48,49,50,51,36,37,38,39,24,25,26,27,12,13,14,15});
  vector<uint8_t> permBit_inv(64);
    for(uint i = 0; i < 64; i++)
  	  permBit_inv[permBit[i]] = i;

  auto anf_s = getANF_Parallel_S(4, sbox);
  auto anf_l = convertLtoPolynome(linear);

  auto anf_los = composition(anf_l, anf_s);
  auto anf_solos = composition(anf_s, anf_los);

  vector<vector<Polynome>> anf (4);
  anf[2] = anf_solos;
  anf[3] = anf_solos;

  uint8_t cst0 [4] = {(64>>4)&0xf, ((64>>4)&0xf)^1, (64 & 0xf)^2, (64 & 0xf)^3};
  anf[0] = anf_los;
  for (unsigned i = 0; i < 4; ++i) {
    for (unsigned j = 0; j < 4; ++j) if (((cst0[i] >> j) & 1) != 0) anf[0][4*i + j] += Monome(0);
  }
  anf[0] = composition(anf_s, anf[0]);
  anf[1] = anf_los;
  for (unsigned i = 0; i < 2; ++i) {
    for (unsigned j = 0; j < 3; ++j) {
      anf[1][4*i + j] += Monome(1u << (16 + (3*i + j)));
      anf[1][4*(i+2) + j] += Monome(1u << (16 + (3*i + j)));
    }
  }
  anf[1] = composition(anf_s, anf[1]);

  auto anf_sol = composition(anf_s, anf_l);

  auto v_in = loadTransitionsIN(anf[pos_in]);
  cout << "v_in : " << v_in.size() << endl;
  auto v_out = (R%2 == 0) ? loadTransitionsOUT(anf[pos_out]) : loadTransitionsOUT(anf_sol);
  cout << "v_out : " << v_out.size() << endl;

  // vector<vector<pair<vector<uint8_t>, vector<uint8_t>>>> res;
  uint ctrSet = 0;
  uint ctrDist = 0;
  bool ret = false;
  map<pair<vector<uint8_t>, vector<uint8_t>>, bool> checkedPairs;
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
      vector<vector<pair<vector<uint8_t>, vector<uint8_t>>>> res;
      res.emplace_back(move(tmp));

      cout << "Set " << ctrSet << ", " << res[0].size() << " pairs to check" << endl;
      ctrSet++;
      auto start = chrono::high_resolution_clock::now();
      auto haveDist = checkForDistinguisher(res, existTrail, MI, checkedPairs);
      auto end = chrono::high_resolution_clock::now();
      cout << "Time for this check : " << chrono::duration<double>(end - start).count() << endl;
      if(haveDist){
      	ctrDist++;
      	ret = true;
      	if(fastSearch) return true;
      }
      if(boundNbSet != 0 && ctrSet == boundNbSet){
      	cout << "Reached set limit bound" << endl;
      	return ret;
      }
    }
  }
  cout << "Found " << ctrDist << " distinguishers in total" << endl;

  return ret;

}