#include "aria.hpp"
#include "DivTable.hpp"

using namespace std;

#define DISPLAY_GEN_OUTPUT false
#define DISPLAY_SOLVER_OUTPUT false
#define DISPLAY_NUMBER_ADDED_CONSTRAINT false
#define TIME_LIMIT_ARIA 7200
#define QUICK_LAZY_COUNT_ARIA true


vector<vector<pair<vector<uint8_t>, vector<uint8_t>>>> allInOutPairsARIA(unsigned i, unsigned j, unsigned R) {
  vector<uint16_t> SB1 = {
	0x63,0x7c,0x77,0x7b,0xf2,0x6b,0x6f,0xc5,0x30,0x01,0x67,0x2b,0xfe,0xd7,0xab,0x76,
	0xca,0x82,0xc9,0x7d,0xfa,0x59,0x47,0xf0,0xad,0xd4,0xa2,0xaf,0x9c,0xa4,0x72,0xc0,
	0xb7,0xfd,0x93,0x26,0x36,0x3f,0xf7,0xcc,0x34,0xa5,0xe5,0xf1,0x71,0xd8,0x31,0x15,
	0x04,0xc7,0x23,0xc3,0x18,0x96,0x05,0x9a,0x07,0x12,0x80,0xe2,0xeb,0x27,0xb2,0x75,
	0x09,0x83,0x2c,0x1a,0x1b,0x6e,0x5a,0xa0,0x52,0x3b,0xd6,0xb3,0x29,0xe3,0x2f,0x84,
	0x53,0xd1,0x00,0xed,0x20,0xfc,0xb1,0x5b,0x6a,0xcb,0xbe,0x39,0x4a,0x4c,0x58,0xcf,
	0xd0,0xef,0xaa,0xfb,0x43,0x4d,0x33,0x85,0x45,0xf9,0x02,0x7f,0x50,0x3c,0x9f,0xa8,
	0x51,0xa3,0x40,0x8f,0x92,0x9d,0x38,0xf5,0xbc,0xb6,0xda,0x21,0x10,0xff,0xf3,0xd2,
	0xcd,0x0c,0x13,0xec,0x5f,0x97,0x44,0x17,0xc4,0xa7,0x7e,0x3d,0x64,0x5d,0x19,0x73,
	0x60,0x81,0x4f,0xdc,0x22,0x2a,0x90,0x88,0x46,0xee,0xb8,0x14,0xde,0x5e,0x0b,0xdb,
	0xe0,0x32,0x3a,0x0a,0x49,0x06,0x24,0x5c,0xc2,0xd3,0xac,0x62,0x91,0x95,0xe4,0x79,
	0xe7,0xc8,0x37,0x6d,0x8d,0xd5,0x4e,0xa9,0x6c,0x56,0xf4,0xea,0x65,0x7a,0xae,0x08,
	0xba,0x78,0x25,0x2e,0x1c,0xa6,0xb4,0xc6,0xe8,0xdd,0x74,0x1f,0x4b,0xbd,0x8b,0x8a,
	0x70,0x3e,0xb5,0x66,0x48,0x03,0xf6,0x0e,0x61,0x35,0x57,0xb9,0x86,0xc1,0x1d,0x9e,
	0xe1,0xf8,0x98,0x11,0x69,0xd9,0x8e,0x94,0x9b,0x1e,0x87,0xe9,0xce,0x55,0x28,0xdf,
	0x8c,0xa1,0x89,0x0d,0xbf,0xe6,0x42,0x68,0x41,0x99,0x2d,0x0f,0xb0,0x54,0xbb,0x16};

  vector<uint16_t> SB2 = {
	0xe2,0x4e,0x54,0xfc,0x94,0xc2,0x4a,0xcc,0x62,0x0d,0x6a,0x46,0x3c,0x4d,0x8b,0xd1,
	0x5e,0xfa,0x64,0xcb,0xb4,0x97,0xbe,0x2b,0xbc,0x77,0x2e,0x03,0xd3,0x19,0x59,0xc1,
	0x1d,0x06,0x41,0x6b,0x55,0xf0,0x99,0x69,0xea,0x9c,0x18,0xae,0x63,0xdf,0xe7,0xbb,
	0x00,0x73,0x66,0xfb,0x96,0x4c,0x85,0xe4,0x3a,0x09,0x45,0xaa,0x0f,0xee,0x10,0xeb,
	0x2d,0x7f,0xf4,0x29,0xac,0xcf,0xad,0x91,0x8d,0x78,0xc8,0x95,0xf9,0x2f,0xce,0xcd,
	0x08,0x7a,0x88,0x38,0x5c,0x83,0x2a,0x28,0x47,0xdb,0xb8,0xc7,0x93,0xa4,0x12,0x53,
	0xff,0x87,0x0e,0x31,0x36,0x21,0x58,0x48,0x01,0x8e,0x37,0x74,0x32,0xca,0xe9,0xb1,
	0xb7,0xab,0x0c,0xd7,0xc4,0x56,0x42,0x26,0x07,0x98,0x60,0xd9,0xb6,0xb9,0x11,0x40,
	0xec,0x20,0x8c,0xbd,0xa0,0xc9,0x84,0x04,0x49,0x23,0xf1,0x4f,0x50,0x1f,0x13,0xdc,
	0xd8,0xc0,0x9e,0x57,0xe3,0xc3,0x7b,0x65,0x3b,0x02,0x8f,0x3e,0xe8,0x25,0x92,0xe5,
	0x15,0xdd,0xfd,0x17,0xa9,0xbf,0xd4,0x9a,0x7e,0xc5,0x39,0x67,0xfe,0x76,0x9d,0x43,
	0xa7,0xe1,0xd0,0xf5,0x68,0xf2,0x1b,0x34,0x70,0x05,0xa3,0x8a,0xd5,0x79,0x86,0xa8,
	0x30,0xc6,0x51,0x4b,0x1e,0xa6,0x27,0xf6,0x35,0xd2,0x6e,0x24,0x16,0x82,0x5f,0xda,
	0xe6,0x75,0xa2,0xef,0x2c,0xb2,0x1c,0x9f,0x5d,0x6f,0x80,0x0a,0x72,0x44,0x9b,0x6c,
	0x90,0x0b,0x5b,0x33,0x7d,0x5a,0x52,0xf3,0x61,0xa1,0xf7,0xb0,0xd6,0x3f,0x7c,0x6d,
	0xed,0x14,0xe0,0xa5,0x3d,0x22,0xb3,0xf8,0x89,0xde,0x71,0x1a,0xaf,0xba,0xb5,0x81};

  vector<uint16_t> SB3 (256);
  for (unsigned x = 0; x < 256; ++x) SB3[SB1[x]] = x;

  vector<uint16_t> SB4 (256);
  for (unsigned x = 0; x < 256; ++x) SB4[SB2[x]] = x;

  vector<vector<uint16_t>> sboxes ({SB1, SB2, SB3, SB4});
  
  //First round is always (S1,S2,invS1,invS2)
  auto in = loadTransitionsIN(getANF_S(sboxes[(i%4)]));

  auto out = (R%2 == 0) ? loadTransitionsOUT(getANF_S(sboxes[((j+2)%4)])) : loadTransitionsOUT(getANF_S(sboxes[(j%4)]));

  cout << "in.size() : " << in.size() << endl;
  cout << "out.size() : " << out.size() << endl;

  vector<vector<pair<vector<uint8_t>, vector<uint8_t>>>> res;
  for (auto pin : in) {
    for (auto pout : out) {
      vector<pair<vector<uint8_t>, vector<uint8_t>>> tmp;
      for (auto xin : pin.first) {
        for (auto xout : pout.first) {
          vector<uint8_t> in_8t (128, 1);
          vector<uint8_t> out_8t (128, 0);
          for (unsigned b = 0; b < 8; ++b) in_8t[8*i + b] = (xin >> b) & 1;
          for (unsigned b = 0; b < 8; ++b) out_8t[8*j + b] = (xout >> b) & 1;
          tmp.emplace_back(move(in_8t), move(out_8t));
        }
      }
      res.emplace_back(move(tmp));
    }
  }

  return res;

}

GRBModel getModelARIA(modelInfo & MI){
/*
	Get the model for ARIA over #rMax round
	#lastLin defines if the linear layer is applied on the last round
	#startingRound defines the index of the first round
	#firstSbox defines if the first Sbox layer is applied

	The #sboxModel parameter defines how the Sbox is modelized :
	- "Simple" use the simplified constraint with PWL

	The #linModel parameter defines how the linear layer is modelized :
	- "CX" modelize the linear layer using the classical copy+xor technique
	- "Simple" modelize the linear layer with the simplified constraint w(x) = w(y)
*/
    string modelName = "ARIA_"+to_string(MI.rMax)+"r_"+MI.sboxModel+"_"+MI.linModel;
    if(MI.lastLin) modelName += "_lastLin";
    if(!MI.firstSbox) modelName += "_noFirstSbox";
    modelName += "_start"+to_string(MI.startingRound);
    modelName += ".mps";

    if(!fileExist("./models/"+modelName)){
        string cmd = "sage ./modelGeneration/genARIA.sage";
        cmd += " -r " + to_string(MI.rMax);
        cmd += " -s " + MI.sboxModel;
        cmd += " -m " + MI.linModel;
        if(MI.lastLin) cmd += " --lastLin True";
        else cmd += " --lastLin False";
        cmd += " --startingRound " + to_string(MI.startingRound);
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

modelData getVariablesARIA(GRBModel & m,
                           modelInfo const & MI){
     /*
    Read the model and return the corresponding blocks of variables
    - #MI defines various parameters for how the trails are checked (see aux_function.hpp)
    */

    modelData MD;
    MD.rMax = MI.rMax;
    MD.lastLin = MI.lastLin;
    MD.firstSbox = MI.firstSbox;


    //ARIA uses 4 different sboxes but no SSB
    MD.allTable = vector<vector<vector<uint32_t>>>(5);
    MD.allIneq = vector<vector<vector<int>>>(5);
    MD.allMatrix = vector<Matrix>(1);

    if(MI.sboxCheckType == CheckType::Table){ //Sbox table
        MD.allTable[0] = readDivTableFromFile("./checkData/ARIAS1Sbox_table.bin"); //S1
        MD.allTable[1] = readDivTableFromFile("./checkData/ARIAS2Sbox_table.bin"); //S2
        MD.allTable[2] = readDivTableFromFile("./checkData/ARIAinvS1Sbox_table.bin"); //invS1
        MD.allTable[3] = readDivTableFromFile("./checkData/ARIAinvS2Sbox_table.bin"); //invS2
    }
    else if(MI.sboxCheckType == CheckType::Ineq){
        MD.allIneq[0] = readIneqFromFile("./checkData/ARIAS1Sbox_ineq.bin"); //S1
        MD.allIneq[1] = readIneqFromFile("./checkData/ARIAS2Sbox_ineq.bin"); //S2
        MD.allIneq[2] = readIneqFromFile("./checkData/ARIAinvS1Sbox_ineq.bin"); //invS1
        MD.allIneq[3] = readIneqFromFile("./checkData/ARIAinvS2Sbox_ineq.bin"); //invS2
    }

    if(MI.MCCheckType == CheckType::Matrix)
    	MD.allMatrix[0] = Matrix("./checkData/matrixARIA.bin");
    else if(MI.MCCheckType == CheckType::Table){
        MD.allTable[4] = readDivTableFromFile("./checkData/ARIALbox_table.bin");
    }
    else if(MI.MCCheckType == CheckType::Ineq){
        MD.allIneq[4]  = readIneqFromFile("./checkData/ARIALbox_ineq.bin");
    }


    auto & tableS1    = MD.allTable[0];
    auto & tableS2    = MD.allTable[1];
    auto & tableinvS1 = MD.allTable[2];
    auto & tableinvS2 = MD.allTable[3];

    auto & ineqS1     = MD.allIneq[0];
    auto & ineqS2     = MD.allIneq[1];
    auto & ineqinvS1  = MD.allIneq[2];
    auto & ineqinvS2  = MD.allIneq[3];

    auto & tableMC   = MD.allTable[4];
    auto & ineqMC    = MD.allIneq[4];
    auto & MCmatrix  = MD.allMatrix[0];

    //Sbox variables
    uint roundFirstSbox = 0;
    if(!MI.firstSbox) roundFirstSbox = 1;

    if(MI.sboxCheckType != CheckType::None){
        for(auto const & r : MI.roundOrder){
        	if(r >= roundFirstSbox){
	        	if((r+MI.startingRound)%2 == 0){
	        		//Sbox layer invS1 invS2 S1 S2 ...
	        		for(uint i = 0; i < 16; i++){
	        			vector<GRBVar> in(8);
						vector<GRBVar> out(8);
						for(uint j = 0; j < 8; j++){
						    in[j] = m.getVarByName("x"+to_string(r)+"_"+to_string(8*i+j));
						    out[j] = m.getVarByName("y"+to_string(r)+"_"+to_string(8*i+j));
						}
						uint imod4 = i%4;
						if(imod4 == 0){ //invS1
							if(MI.sboxCheckType == CheckType::Table)
			                    MD.allVarsTable.emplace_back(move(in), move(out), tableinvS1,ConstrType::Sbox);
			                else if(MI.sboxCheckType == CheckType::Ineq)
			                    MD.allVarsIneq.emplace_back(move(in), move(out), ineqinvS1,ConstrType::Sbox);
						}
						else if(imod4 == 1){//invS2
							if(MI.sboxCheckType == CheckType::Table)
			                    MD.allVarsTable.emplace_back(move(in), move(out), tableinvS2,ConstrType::Sbox);
			                else if(MI.sboxCheckType == CheckType::Ineq)
			                    MD.allVarsIneq.emplace_back(move(in), move(out), ineqinvS2,ConstrType::Sbox);
						}
						else if(imod4 == 2){
							if(MI.sboxCheckType == CheckType::Table)
			                    MD.allVarsTable.emplace_back(move(in), move(out), tableS1,ConstrType::Sbox);
			                else if(MI.sboxCheckType == CheckType::Ineq)
			                    MD.allVarsIneq.emplace_back(move(in), move(out), ineqS1,ConstrType::Sbox);
						}
						else{
							if(MI.sboxCheckType == CheckType::Table)
			                    MD.allVarsTable.emplace_back(move(in), move(out), tableS2,ConstrType::Sbox);
			                else if(MI.sboxCheckType == CheckType::Ineq)
			                    MD.allVarsIneq.emplace_back(move(in), move(out), ineqS2,ConstrType::Sbox);
						}
	        		}
	        	}
	        	else{
	        		//Sbox layer S1 S2 invS1 invS2...
	        		for(uint i = 0; i < 16; i++){
	        			vector<GRBVar> in(8);
						vector<GRBVar> out(8);
						for(uint j = 0; j < 8; j++){
						    in[j] = m.getVarByName("x"+to_string(r)+"_"+to_string(8*i+j));
						    out[j] = m.getVarByName("y"+to_string(r)+"_"+to_string(8*i+j));
						}
						uint imod4 = i%4;
						if(imod4 == 0){ //invS1
							if(MI.sboxCheckType == CheckType::Table)
			                    MD.allVarsTable.emplace_back(move(in), move(out), tableS1,ConstrType::Sbox);
			                else if(MI.sboxCheckType == CheckType::Ineq)
			                    MD.allVarsIneq.emplace_back(move(in), move(out), ineqS1,ConstrType::Sbox);
						}
						else if(imod4 == 1){//invS2
							if(MI.sboxCheckType == CheckType::Table)
			                    MD.allVarsTable.emplace_back(move(in), move(out), tableS2,ConstrType::Sbox);
			                else if(MI.sboxCheckType == CheckType::Ineq)
			                    MD.allVarsIneq.emplace_back(move(in), move(out), ineqS2,ConstrType::Sbox);
						}
						else if(imod4 == 2){
							if(MI.sboxCheckType == CheckType::Table)
			                    MD.allVarsTable.emplace_back(move(in), move(out), tableinvS1,ConstrType::Sbox);
			                else if(MI.sboxCheckType == CheckType::Ineq)
			                    MD.allVarsIneq.emplace_back(move(in), move(out), ineqinvS1,ConstrType::Sbox);
						}
						else{
							if(MI.sboxCheckType == CheckType::Table)
			                    MD.allVarsTable.emplace_back(move(in), move(out), tableinvS2,ConstrType::Sbox);
			                else if(MI.sboxCheckType == CheckType::Ineq)
			                    MD.allVarsIneq.emplace_back(move(in), move(out), ineqinvS2,ConstrType::Sbox);
						}
	        		}
	        	}
	        }
        }
    }

    //MC variables
    uint roundLastLin = MI.rMax;
    if(!MI.lastLin) roundLastLin--;

    if(MI.MCCheckType != CheckType::None){
    	if(MI.MCCheckType == CheckType::Matrix){
	        for(auto const & r : MI.roundOrder){
	        	if(r < roundLastLin){
		            vector<GRBVar> in(128);
		            vector<GRBVar> out(128);
		            for(uint i = 0; i < 128; i++){
		                in[i] = m.getVarByName("y"+to_string(r)+"_"+to_string(i));
		                out[i] = m.getVarByName("x"+to_string(r+1)+"_"+to_string(i));
		            }
		            // if(MI.MCCheckType == CheckType::Table)
		            //     MD.allVarsTable.emplace_back(move(in), move(out), tableMC);
		            // else if(MI.MCCheckType == CheckType::Ineq)
		            //     MD.allVarsIneq.emplace_back(move(in), move(out), ineqMC);
		            if(MI.MCCheckType == CheckType::Matrix)
		                MD.allVarsMatrix.emplace_back(move(in), move(out), MCmatrix);
		        }
	        }
	    }
	    else if(MI.MCCheckType == CheckType::Ineq || MI.MCCheckType == CheckType::Table){
	    	//Consider MC as independent Lboxes
	    	for(auto const & r : MI.roundOrder){
	        	if(r < roundLastLin){
	        		for(uint offset = 0; offset < 8; offset++){ //For each set of 16 bits
			            vector<GRBVar> in(16);
			            vector<GRBVar> out(16);
			            for(uint i = 0; i < 16; i++){
			                in[i] = m.getVarByName("y"+to_string(r)+"_"+to_string(offset + i*8));
			                out[i] = m.getVarByName("x"+to_string(r+1)+"_"+to_string(offset + i*8));
			            }
			            if(MI.MCCheckType == CheckType::Table)
			                MD.allVarsTable.emplace_back(move(in), move(out), tableMC,ConstrType::Lin);
			            else if(MI.MCCheckType == CheckType::Ineq)
			                MD.allVarsIneq.emplace_back(move(in), move(out), ineqMC,ConstrType::Lin);
			        }
		        }
	        }
	    }
    }

    return MD;
}

uint64_t globalCtrCheckARIA = 0;
uint64_t globalCtrSolARIA = 0;
uint64_t globalCtrConstrSboxARIA = 0;
uint64_t globalCtrConstrLinARIA = 0;
uint64_t globalCtrConstrSSBARIA = 0;
bool existTrailARIA(GRBModel & m,
                   modelData const & MD,
                   vector<uint8_t> const & input,
                   vector<uint8_t> const & output,
                   bool const useCallback,
                   uint const maxNumberConstr){
/*
    Return true if there is a trail from #input to #output using #model for ARIA
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

    for(uint i = 0; i < 128; i++){
        m.addConstr(m.getVarByName(inputvarPrefix+to_string(i)) == input[i]);
        m.addConstr(m.getVarByName(outputvarPrefix+to_string(i)) == output[i]);
    }
    m.update();

    if(useCallback){
        m.set(GRB_IntParam_LazyConstraints, 1);
        m.set(GRB_DoubleParam_TimeLimit, TIME_LIMIT_ARIA);

        CustomCallback cb(MD.allVarsTable, MD.allVarsIneq, MD.allVarsMatrix, maxNumberConstr);
        m.setCallback(&cb);

        m.update();
        m.optimize();
        if(DISPLAY_NUMBER_ADDED_CONSTRAINT){
            cout << cb.ctrSol << " solutions examined, ";
            cout << cb.ctrGlobalConstr << " constraints added" << endl;
        }
        globalCtrCheckARIA++;
        globalCtrSolARIA += cb.ctrSol;
        globalCtrConstrSboxARIA += cb.ctrConstrSbox;
		globalCtrConstrLinARIA += cb.ctrConstrLin;
		globalCtrConstrSSBARIA += cb.ctrConstrSSB;

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

bool existTrailARIA(modelInfo & MI,
                   vector<uint8_t> const & input,
                   vector<uint8_t> const & output){
/*
    Return true if there is a trail from #input to #output for ARIA
    #MI should contain the necessary information to generate the model and parametrize the solving
*/

    GRBModel m = getModelARIA(MI);
    auto MD = getVariablesARIA(m, MI);

    return existTrailARIA(m,MD,input,output,MI.useCallback,MI.maxNumberConstr);
}

void genSboxFileARIA(){
	//Generate the tables and ineq files for ARIA

	vector<uint32_t> S1(
		{0x63,0x7C,0x77,0x7B,0xF2,0x6B,0x6F,0xC5,0x30,0x01,0x67,0x2B,0xFE,0xD7,0xAB,0x76,
		 0xCA,0x82,0xC9,0x7D,0xFA,0x59,0x47,0xF0,0xAD,0xD4,0xA2,0xAF,0x9C,0xA4,0x72,0xC0,
		 0xB7,0xFD,0x93,0x26,0x36,0x3F,0xF7,0xCC,0x34,0xA5,0xE5,0xF1,0x71,0xD8,0x31,0x15,
		 0x04,0xC7,0x23,0xC3,0x18,0x96,0x05,0x9A,0x07,0x12,0x80,0xE2,0xEB,0x27,0xB2,0x75,
		 0x09,0x83,0x2C,0x1A,0x1B,0x6E,0x5A,0xA0,0x52,0x3B,0xD6,0xB3,0x29,0xE3,0x2F,0x84,
		 0x53,0xD1,0x00,0xED,0x20,0xFC,0xB1,0x5B,0x6A,0xCB,0xBE,0x39,0x4A,0x4C,0x58,0xCF,
		 0xD0,0xEF,0xAA,0xFB,0x43,0x4D,0x33,0x85,0x45,0xF9,0x02,0x7F,0x50,0x3C,0x9F,0xA8,
		 0x51,0xA3,0x40,0x8F,0x92,0x9D,0x38,0xF5,0xBC,0xB6,0xDA,0x21,0x10,0xFF,0xF3,0xD2,
		 0xCD,0x0C,0x13,0xEC,0x5F,0x97,0x44,0x17,0xC4,0xA7,0x7E,0x3D,0x64,0x5D,0x19,0x73,
		 0x60,0x81,0x4F,0xDC,0x22,0x2A,0x90,0x88,0x46,0xEE,0xB8,0x14,0xDE,0x5E,0x0B,0xDB,
		 0xE0,0x32,0x3A,0x0A,0x49,0x06,0x24,0x5C,0xC2,0xD3,0xAC,0x62,0x91,0x95,0xE4,0x79,
		 0xE7,0xC8,0x37,0x6D,0x8D,0xD5,0x4E,0xA9,0x6C,0x56,0xF4,0xEA,0x65,0x7A,0xAE,0x08,
		 0xBA,0x78,0x25,0x2E,0x1C,0xA6,0xB4,0xC6,0xE8,0xDD,0x74,0x1F,0x4B,0xBD,0x8B,0x8A,
		 0x70,0x3E,0xB5,0x66,0x48,0x03,0xF6,0x0E,0x61,0x35,0x57,0xB9,0x86,0xC1,0x1D,0x9E,
		 0xE1,0xF8,0x98,0x11,0x69,0xD9,0x8E,0x94,0x9B,0x1E,0x87,0xE9,0xCE,0x55,0x28,0xDF,
		 0x8C,0xA1,0x89,0x0D,0xBF,0xE6,0x42,0x68,0x41,0x99,0x2D,0x0F,0xB0,0x54,0xBB,0x16});

	cout << "S1" << endl;
	genAndSaveTableIneq(vector<vector<uint32_t>>({S1}),8,8,"./checkData/ARIAS1Sbox",true);

	vector<uint32_t> invS1(
		{0x52,0x09,0x6A,0xD5,0x30,0x36,0xA5,0x38,0xBF,0x40,0xA3,0x9E,0x81,0xF3,0xD7,0xFB,
		 0x7C,0xE3,0x39,0x82,0x9B,0x2F,0xFF,0x87,0x34,0x8E,0x43,0x44,0xC4,0xDE,0xE9,0xCB,
		 0x54,0x7B,0x94,0x32,0xA6,0xC2,0x23,0x3D,0xEE,0x4C,0x95,0x0B,0x42,0xFA,0xC3,0x4E,
		 0x08,0x2E,0xA1,0x66,0x28,0xD9,0x24,0xB2,0x76,0x5B,0xA2,0x49,0x6D,0x8B,0xD1,0x25,
		 0x72,0xF8,0xF6,0x64,0x86,0x68,0x98,0x16,0xD4,0xA4,0x5C,0xCC,0x5D,0x65,0xB6,0x92,
		 0x6C,0x70,0x48,0x50,0xFD,0xED,0xB9,0xDA,0x5E,0x15,0x46,0x57,0xA7,0x8D,0x9D,0x84,
		 0x90,0xD8,0xAB,0x00,0x8C,0xBC,0xD3,0x0A,0xF7,0xE4,0x58,0x05,0xB8,0xB3,0x45,0x06,
		 0xD0,0x2C,0x1E,0x8F,0xCA,0x3F,0x0F,0x02,0xC1,0xAF,0xBD,0x03,0x01,0x13,0x8A,0x6B,
		 0x3A,0x91,0x11,0x41,0x4F,0x67,0xDC,0xEA,0x97,0xF2,0xCF,0xCE,0xF0,0xB4,0xE6,0x73,
		 0x96,0xAC,0x74,0x22,0xE7,0xAD,0x35,0x85,0xE2,0xF9,0x37,0xE8,0x1C,0x75,0xDF,0x6E,
		 0x47,0xF1,0x1A,0x71,0x1D,0x29,0xC5,0x89,0x6F,0xB7,0x62,0x0E,0xAA,0x18,0xBE,0x1B,
		 0xFC,0x56,0x3E,0x4B,0xC6,0xD2,0x79,0x20,0x9A,0xDB,0xC0,0xFE,0x78,0xCD,0x5A,0xF4,
		 0x1F,0xDD,0xA8,0x33,0x88,0x07,0xC7,0x31,0xB1,0x12,0x10,0x59,0x27,0x80,0xEC,0x5F,
		 0x60,0x51,0x7F,0xA9,0x19,0xB5,0x4A,0x0D,0x2D,0xE5,0x7A,0x9F,0x93,0xC9,0x9C,0xEF,
		 0xA0,0xE0,0x3B,0x4D,0xAE,0x2A,0xF5,0xB0,0xC8,0xEB,0xBB,0x3C,0x83,0x53,0x99,0x61,
		 0x17,0x2B,0x04,0x7E,0xBA,0x77,0xD6,0x26,0xE1,0x69,0x14,0x63,0x55,0x21,0x0C,0x7D});

	cout << "invS1" << endl;
	genAndSaveTableIneq(vector<vector<uint32_t>>({invS1}),8,8,"./checkData/ARIAinvS1Sbox",true);

	vector<uint32_t> S2(
		{0xE2,0x4E,0x54,0xFC,0x94,0xC2,0x4A,0xCC,0x62,0x0D,0x6A,0x46,0x3C,0x4D,0x8B,0xD1,
		 0x5E,0xFA,0x64,0xCB,0xB4,0x97,0xBE,0x2B,0xBC,0x77,0x2E,0x03,0xD3,0x19,0x59,0xC1,
		 0x1D,0x06,0x41,0x6B,0x55,0xF0,0x99,0x69,0xEA,0x9C,0x18,0xAE,0x63,0xDF,0xE7,0xBB,
		 0x00,0x73,0x66,0xFB,0x96,0x4C,0x85,0xE4,0x3A,0x09,0x45,0xAA,0x0F,0xEE,0x10,0xEB,
		 0x2D,0x7F,0xF4,0x29,0xAC,0xCF,0xAD,0x91,0x8D,0x78,0xC8,0x95,0xF9,0x2F,0xCE,0xCD,
		 0x08,0x7A,0x88,0x38,0x5C,0x83,0x2A,0x28,0x47,0xDB,0xB8,0xC7,0x93,0xA4,0x12,0x53,
		 0xFF,0x87,0x0E,0x31,0x36,0x21,0x58,0x48,0x01,0x8E,0x37,0x74,0x32,0xCA,0xE9,0xB1,
		 0xB7,0xAB,0x0C,0xD7,0xC4,0x56,0x42,0x26,0x07,0x98,0x60,0xD9,0xB6,0xB9,0x11,0x40,
		 0xEC,0x20,0x8C,0xBD,0xA0,0xC9,0x84,0x04,0x49,0x23,0xF1,0x4F,0x50,0x1F,0x13,0xDC,
		 0xD8,0xC0,0x9E,0x57,0xE3,0xC3,0x7B,0x65,0x3B,0x02,0x8F,0x3E,0xE8,0x25,0x92,0xE5,
		 0x15,0xDD,0xFD,0x17,0xA9,0xBF,0xD4,0x9A,0x7E,0xC5,0x39,0x67,0xFE,0x76,0x9D,0x43,
		 0xA7,0xE1,0xD0,0xF5,0x68,0xF2,0x1B,0x34,0x70,0x05,0xA3,0x8A,0xD5,0x79,0x86,0xA8,
		 0x30,0xC6,0x51,0x4B,0x1E,0xA6,0x27,0xF6,0x35,0xD2,0x6E,0x24,0x16,0x82,0x5F,0xDA,
		 0xE6,0x75,0xA2,0xEF,0x2C,0xB2,0x1C,0x9F,0x5D,0x6F,0x80,0x0A,0x72,0x44,0x9B,0x6C,
		 0x90,0x0B,0x5B,0x33,0x7D,0x5A,0x52,0xF3,0x61,0xA1,0xF7,0xB0,0xD6,0x3F,0x7C,0x6D,
		 0xED,0x14,0xE0,0xA5,0x3D,0x22,0xB3,0xF8,0x89,0xDE,0x71,0x1A,0xAF,0xBA,0xB5,0x81});

	cout << "S2" << endl;
	genAndSaveTableIneq(vector<vector<uint32_t>>({S2}),8,8,"./checkData/ARIAS2Sbox",true);

	vector<uint32_t> invS2(
		{0x30,0x68,0x99,0x1B,0x87,0xB9,0x21,0x78,0x50,0x39,0xDB,0xE1,0x72,0x09,0x62,0x3C,
		 0x3E,0x7E,0x5E,0x8E,0xF1,0xA0,0xCC,0xA3,0x2A,0x1D,0xFB,0xB6,0xD6,0x20,0xC4,0x8D,
		 0x81,0x65,0xF5,0x89,0xCB,0x9D,0x77,0xC6,0x57,0x43,0x56,0x17,0xD4,0x40,0x1A,0x4D,
		 0xC0,0x63,0x6C,0xE3,0xB7,0xC8,0x64,0x6A,0x53,0xAA,0x38,0x98,0x0C,0xF4,0x9B,0xED,
		 0x7F,0x22,0x76,0xAF,0xDD,0x3A,0x0B,0x58,0x67,0x88,0x06,0xC3,0x35,0x0D,0x01,0x8B,
		 0x8C,0xC2,0xE6,0x5F,0x02,0x24,0x75,0x93,0x66,0x1E,0xE5,0xE2,0x54,0xD8,0x10,0xCE,
		 0x7A,0xE8,0x08,0x2C,0x12,0x97,0x32,0xAB,0xB4,0x27,0x0A,0x23,0xDF,0xEF,0xCA,0xD9,
		 0xB8,0xFA,0xDC,0x31,0x6B,0xD1,0xAD,0x19,0x49,0xBD,0x51,0x96,0xEE,0xE4,0xA8,0x41,
		 0xDA,0xFF,0xCD,0x55,0x86,0x36,0xBE,0x61,0x52,0xF8,0xBB,0x0E,0x82,0x48,0x69,0x9A,
		 0xE0,0x47,0x9E,0x5C,0x04,0x4B,0x34,0x15,0x79,0x26,0xA7,0xDE,0x29,0xAE,0x92,0xD7,
		 0x84,0xE9,0xD2,0xBA,0x5D,0xF3,0xC5,0xB0,0xBF,0xA4,0x3B,0x71,0x44,0x46,0x2B,0xFC,
		 0xEB,0x6F,0xD5,0xF6,0x14,0xFE,0x7C,0x70,0x5A,0x7D,0xFD,0x2F,0x18,0x83,0x16,0xA5,
		 0x91,0x1F,0x05,0x95,0x74,0xA9,0xC1,0x5B,0x4A,0x85,0x6D,0x13,0x07,0x4F,0x4E,0x45,
		 0xB2,0x0F,0xC9,0x1C,0xA6,0xBC,0xEC,0x73,0x90,0x7B,0xCF,0x59,0x8F,0xA1,0xF9,0x2D,
		 0xF2,0xB1,0x00,0x94,0x37,0x9F,0xD0,0x2E,0x9C,0x6E,0x28,0x3F,0x80,0xF0,0x3D,0xD3,
		 0x25,0x8A,0xB5,0xE7,0x42,0xB3,0xC7,0xEA,0xF7,0x4C,0x11,0x33,0x03,0xA2,0xAC,0x60});

	cout << "invS2" << endl;
	genAndSaveTableIneq(vector<vector<uint32_t>>({invS2}),8,8,"./checkData/ARIAinvS2Sbox",true);
}

void createMatrixFileARIA(){
	//Create the matrix file ./checkData/matrixARIA.bin

	//Slightly compressed in the code, the binary matrix is 128 x 128...
	//The values were obtained by hardcoding the binary matrix and priting it with M.printRawHex()
	//The hex values printed are actually the value of the chunks in the internal representation, so more compact code (still ugly and large)
	vector<vector<uint64_t>> val({{0x0001000101000000,0x0001010000000101},
								  {0x0002000202000000,0x0002020000000202},
								  {0x0004000404000000,0x0004040000000404},
								  {0x0008000808000000,0x0008080000000808},
								  {0x0010001010000000,0x0010100000001010},
								  {0x0020002020000000,0x0020200000002020},
								  {0x0040004040000000,0x0040400000004040},
								  {0x0080008080000000,0x0080800000008080},
								  {0x0100010000010000,0x0100000100000101},
								  {0x0200020000020000,0x0200000200000202},
								  {0x0400040000040000,0x0400000400000404},
								  {0x0800080000080000,0x0800000800000808},
								  {0x1000100000100000,0x1000001000001010},
								  {0x2000200000200000,0x2000002000002020},
								  {0x4000400000400000,0x4000004000004040},
								  {0x8000800000800000,0x8000008000008080},
								  {0x0001000100000100,0x0100000101010000},
								  {0x0002000200000200,0x0200000202020000},
								  {0x0004000400000400,0x0400000404040000},
								  {0x0008000800000800,0x0800000808080000},
								  {0x0010001000001000,0x1000001010100000},
								  {0x0020002000002000,0x2000002020200000},
								  {0x0040004000004000,0x4000004040400000},
								  {0x0080008000008000,0x8000008080800000},
								  {0x0100010000000001,0x0001010001010000},
								  {0x0200020000000002,0x0002020002020000},
								  {0x0400040000000004,0x0004040004040000},
								  {0x0800080000000008,0x0008080008080000},
								  {0x1000100000000010,0x0010100010100000},
								  {0x2000200000000020,0x0020200020200000},
								  {0x4000400000000040,0x0040400040400000},
								  {0x8000800000000080,0x0080800080800000},
								  {0x0000010000010001,0x0101000001000001},
								  {0x0000020000020002,0x0202000002000002},
								  {0x0000040000040004,0x0404000004000004},
								  {0x0000080000080008,0x0808000008000008},
								  {0x0000100000100010,0x1010000010000010},
								  {0x0000200000200020,0x2020000020000020},
								  {0x0000400000400040,0x4040000040000040},
								  {0x0000800000800080,0x8080000080000080},
								  {0x0000000101000100,0x0101000000010100},
								  {0x0000000202000200,0x0202000000020200},
								  {0x0000000404000400,0x0404000000040400},
								  {0x0000000808000800,0x0808000000080800},
								  {0x0000001010001000,0x1010000000101000},
								  {0x0000002020002000,0x2020000000202000},
								  {0x0000004040004000,0x4040000000404000},
								  {0x0000008080008000,0x8080000000808000},
								  {0x0100000000010001,0x0000010100010100},
								  {0x0200000000020002,0x0000020200020200},
								  {0x0400000000040004,0x0000040400040400},
								  {0x0800000000080008,0x0000080800080800},
								  {0x1000000000100010,0x0000101000101000},
								  {0x2000000000200020,0x0000202000202000},
								  {0x4000000000400040,0x0000404000404000},
								  {0x8000000000800080,0x0000808000808000},
								  {0x0001000001000100,0x0000010101000001},
								  {0x0002000002000200,0x0000020202000002},
								  {0x0004000004000400,0x0000040404000004},
								  {0x0008000008000800,0x0000080808000008},
								  {0x0010000010001000,0x0000101010000010},
								  {0x0020000020002000,0x0000202020000020},
								  {0x0040000040004000,0x0000404040000040},
								  {0x0080000080008000,0x0000808080000080},
								  {0x0100000100000101,0x0100010000010000},
								  {0x0200000200000202,0x0200020000020000},
								  {0x0400000400000404,0x0400040000040000},
								  {0x0800000800000808,0x0800080000080000},
								  {0x1000001000001010,0x1000100000100000},
								  {0x2000002000002020,0x2000200000200000},
								  {0x4000004000004040,0x4000400000400000},
								  {0x8000008000008080,0x8000800000800000},
								  {0x0001010000000101,0x0001000101000000},
								  {0x0002020000000202,0x0002000202000000},
								  {0x0004040000000404,0x0004000404000000},
								  {0x0008080000000808,0x0008000808000000},
								  {0x0010100000001010,0x0010001010000000},
								  {0x0020200000002020,0x0020002020000000},
								  {0x0040400000004040,0x0040004040000000},
								  {0x0080800000008080,0x0080008080000000},
								  {0x0001010001010000,0x0100010000000001},
								  {0x0002020002020000,0x0200020000000002},
								  {0x0004040004040000,0x0400040000000004},
								  {0x0008080008080000,0x0800080000000008},
								  {0x0010100010100000,0x1000100000000010},
								  {0x0020200020200000,0x2000200000000020},
								  {0x0040400040400000,0x4000400000000040},
								  {0x0080800080800000,0x8000800000000080},
								  {0x0100000101010000,0x0001000100000100},
								  {0x0200000202020000,0x0002000200000200},
								  {0x0400000404040000,0x0004000400000400},
								  {0x0800000808080000,0x0008000800000800},
								  {0x1000001010100000,0x0010001000001000},
								  {0x2000002020200000,0x0020002000002000},
								  {0x4000004040400000,0x0040004000004000},
								  {0x8000008080800000,0x0080008000008000},
								  {0x0101000000010100,0x0000000101000100},
								  {0x0202000000020200,0x0000000202000200},
								  {0x0404000000040400,0x0000000404000400},
								  {0x0808000000080800,0x0000000808000800},
								  {0x1010000000101000,0x0000001010001000},
								  {0x2020000000202000,0x0000002020002000},
								  {0x4040000000404000,0x0000004040004000},
								  {0x8080000000808000,0x0000008080008000},
								  {0x0101000001000001,0x0000010000010001},
								  {0x0202000002000002,0x0000020000020002},
								  {0x0404000004000004,0x0000040000040004},
								  {0x0808000008000008,0x0000080000080008},
								  {0x1010000010000010,0x0000100000100010},
								  {0x2020000020000020,0x0000200000200020},
								  {0x4040000040000040,0x0000400000400040},
								  {0x8080000080000080,0x0000800000800080},
								  {0x0000010101000001,0x0001000001000100},
								  {0x0000020202000002,0x0002000002000200},
								  {0x0000040404000004,0x0004000004000400},
								  {0x0000080808000008,0x0008000008000800},
								  {0x0000101010000010,0x0010000010001000},
								  {0x0000202020000020,0x0020000020002000},
								  {0x0000404040000040,0x0040000040004000},
								  {0x0000808080000080,0x0080000080008000},
								  {0x0000010100010100,0x0100000000010001},
								  {0x0000020200020200,0x0200000000020002},
								  {0x0000040400040400,0x0400000000040004},
								  {0x0000080800080800,0x0800000000080008},
								  {0x0000101000101000,0x1000000000100010},
								  {0x0000202000202000,0x2000000000200020},
								  {0x0000404000404000,0x4000000000400040},
								  {0x0000808000808000,0x8000000000800080}});

	Matrix M(128,128);
	for(uint i = 0; i < 128; i++){
		for(uint j = 0; j < 2; j++)
			M.get(i)[j] = val[i][j];
	}

	M.saveToFile("./checkData/matrixARIA.bin");
}

void createMatrixQMAria(){
	vector<uint16_t> M({
		0b0110001101011000,
		0b1001001110100100,
		0b1001110001010010,
		0b0110110010100001,
		0b1100100100100101,
		0b1100011000011010,
		0b0011011010000101,
		0b0011100101001010,
		0b1010010010010011,
		0b0101100001100011,
		0b1010000101101100,
		0b0101001010011100,
		0b0001101011000110,
		0b0010010111001001,
		0b0100101000111001,
		0b1000010100110110});

	vector<uint32_t> S(1 << 16);
	for(uint32_t x = 0; x < (1ULL << 16); x++){
		uint32_t y = 0;
		for(uint i = 0; i < 16; i++){
			if(__builtin_parity(x & M[i]) != 0)
				y |= (1 << i);
		}
		S[x] = y;
	}

	genAndSaveTableIneq(vector<vector<uint32_t>>({S}),16,16,"./checkData/ARIALbox",true);

}


void checkNormalDistinguisherARIA(){
	/*
		Search for a normal (no linear combination) distinguisher over 4 rounds of ARIA (new)
	*/

	//These parameters have no effect for ARIA but are required for the modelInfo object
	string arxModel = "ARX";
	CheckType ARXCheckType = CheckType::None;
	bool smartLin = true;

	//Fixed parameters
	uint rMax = 4;
	bool firstSbox = true;
	bool lastLin = false;
	CheckType sboxCheckType = CheckType::None;
	CheckType MCCheckType = CheckType::None;
	CheckType SSBCheckType = CheckType::None;
	bool useCallback = true;
	uint maxNumberConstr = 100;
	vector<uint> roundOrder = roundOrderAscend(rMax);

	//Variable parameters
	//Even with the simple model, no need to use the checks, so limited number of variable parameters
	vector<string> listSboxModel({"QM", "Simple"});
	vector<string> listLinModel({"CX", "Simple"});
	vector<uint> listStartinground({0,1}); //only the parity of the round matter
	//127 bits of data or 120 bits of data
	vector<vector<uint8_t>> listInput({{0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},
		{0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}});


	//The searchs
	for(auto const & input : listInput){
		cout << "Input : ";
		for(auto const & tmpv : input)
			cout << uint(tmpv);
		cout << endl;
		for(auto const & startingRound : listStartinground){
			for(auto const & sboxModel : listSboxModel){
				for(auto const & linModel : listLinModel){

					globalCtrCheckARIA = 0;
					globalCtrSolARIA = 0;
					globalCtrConstrSboxARIA = 0;
					globalCtrConstrLinARIA = 0;
					globalCtrConstrSSBARIA = 0;

					modelInfo MI(rMax,sboxModel,linModel,arxModel,firstSbox,lastLin,smartLin,startingRound,sboxCheckType,MCCheckType,SSBCheckType,ARXCheckType,useCallback,maxNumberConstr,roundOrder);

					auto start = chrono::high_resolution_clock::now();
					//Check each output bit (all balanced)
					uint ctrBalanced = 0;
					for(uint i = 0; i < 128; i++){
						vector<uint8_t> output(128,0);
						output[i] = 1;
						if(!existTrailARIA(MI,input,output)) ctrBalanced++;
					}
					auto end = chrono::high_resolution_clock::now();
					cout << "Starting round : " << startingRound << " - ";
					cout << " sboxModel : " << sboxModel << " - linModel : " << linModel << endl;
					cout << ctrBalanced << " balanced bits" << endl;
					cout << chrono::duration<double>(end - start).count() << " seconds" << endl;
					cout << globalCtrCheckARIA << " input/output pairs checked" << endl;
					cout << globalCtrSolARIA << " solutions examined" << endl;
					cout << globalCtrConstrSboxARIA << " constraints added for the sbox" << endl;
					cout << globalCtrConstrLinARIA  << " constraints added for the lin layer" << endl;
					cout << globalCtrConstrSSBARIA << " constraints added for the SSB" << endl;
				}
			}
		}
	}
}

void testSearchARIA(){


// 1 SB | precalc, Sbox Type 1

//   MC |
// 2 SB |
//   MC |
// 3 SB | MILP
//   MC |
// 4 SB |
//   MC |

// 5 SB | precalc

	//These parameters have no effect for ARIA but are required for the modelInfo object
	string arxModel = "ARX";
	unsigned int startingRound = 1;
	CheckType ARXCheckType = CheckType::None;
	CheckType SSBCheckType = CheckType::None; //No SSB for ARIA

	//Fixed parameters
	uint rMax = 4; //1 partial round at the start (only MC) + 3 full rounds
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

	vector<string> listSboxModel({"QM", "Simple"});
	vector<string> listLinModel({"QM", "CX", "Simple"});
	vector<bool> listUseCallback({true});
	vector<vector<uint>> listroundOrderIfNeeded({roundOrderOutsideIn(rMax), roundOrderAscend(rMax), roundOrderDescend(rMax)});
	auto const allTest = allInOutPairsARIA(0,0,5);

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

	      	 	 globalCtrCheckARIA = 0;
				 globalCtrSolARIA = 0;
	      	 	 globalCtrConstrSboxARIA = 0;
				 globalCtrConstrLinARIA = 0;
				 globalCtrConstrSSBARIA = 0;

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
				 auto haveDist = checkForDistinguisher(allTest, existTrailARIA, MI, false);
				 auto end = chrono::high_resolution_clock::now();
	 
				 cout << "Time : " << chrono::duration<double>(end - start).count() << endl;
				 if(haveDist)cout << "Distinguisher !!!" << endl << endl;
				 else cout << "No distinguisher..." << endl << endl;
				 cout << globalCtrCheckARIA << " input/output pairs checked" << endl;
				 cout << globalCtrSolARIA << " solutions examined" << endl;
				 cout << globalCtrConstrSboxARIA << " constraints added for the sbox" << endl;
				 cout << globalCtrConstrLinARIA  << " constraints added for the lin layer" << endl;
				 cout << globalCtrConstrSSBARIA << " constraints added for the SSB" << endl;
			 }
			 ctrParameterSet++;
	}}}}}}}
}

void searchDistinguisher4RoundsARIA(){

	//These parameters have no effect for ARIA but are required for the modelInfo object
	string arxModel = "ARX";
	CheckType ARXCheckType = CheckType::None;
	bool smartLin = true;

	//Fixed parameters
	uint rMax = 4;
	bool firstSbox = true;
	bool lastLin = false;
	CheckType sboxCheckType = CheckType::None;
	CheckType MCCheckType = CheckType::None;
	CheckType SSBCheckType = CheckType::None;
	bool useCallback = true;
	uint maxNumberConstr = 100;
	vector<uint> roundOrder = roundOrderAscend(rMax);

	vector<uint8_t> input({0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1});
	unsigned int startingRound = 0;
	string sboxModel = "Simple";
	string linModel = "Simple";

	globalCtrCheckARIA = 0;
	globalCtrSolARIA = 0;
	globalCtrConstrSboxARIA = 0;
	globalCtrConstrLinARIA = 0;
	globalCtrConstrSSBARIA = 0;

	modelInfo MI(rMax,sboxModel,linModel,arxModel,firstSbox,lastLin,smartLin,startingRound,sboxCheckType,MCCheckType,SSBCheckType,ARXCheckType,useCallback,maxNumberConstr,roundOrder);

	auto start = chrono::high_resolution_clock::now();
	//Check each output bit (all balanced)
	uint ctrBalanced = 0;
	for(uint i = 0; i < 128; i++){
		vector<uint8_t> output(128,0);
		output[i] = 1;
		if(!existTrailARIA(MI,input,output)) ctrBalanced++;
	}
	auto end = chrono::high_resolution_clock::now();
	cout << "Starting round : " << startingRound << " - ";
	cout << " sboxModel : " << sboxModel << " - linModel : " << linModel << endl;
	cout << ctrBalanced << " balanced bits" << endl;
	cout << chrono::duration<double>(end - start).count() << " seconds" << endl;
	cout << globalCtrCheckARIA << " input/output pairs checked" << endl;
	cout << globalCtrSolARIA << " solutions examined" << endl;
	cout << globalCtrConstrSboxARIA << " constraints added for the sbox" << endl;
	cout << globalCtrConstrLinARIA  << " constraints added for the lin layer" << endl;
	cout << globalCtrConstrSSBARIA << " constraints added for the SSB" << endl;

}

void searchDistinguisher5RoundsARIA(){
	/*
		Search for a distinguisher over 5 rounds ARIA
		Parameters are fixed and were determined with experiments on a subset of cases
	*/

	uint rMax = 4; //1 partial round at the start (only MC) + 3 full rounds
	string sboxModel = "QM";
	string linModel = "CX";
	string arxModel = "ARX";
	bool firstSbox = false;
	bool lastLin = true;
	bool smartLin = true;
	unsigned int startingRound = 1;
	CheckType sboxCheckType = CheckType::None;
	CheckType MCCheckType = CheckType::Ineq;
	CheckType SSBCheckType = CheckType::None;
	CheckType ARXCheckType = CheckType::None;
	bool useCallback = true;
	uint maxNumberConstr = 10;
	vector<uint> roundOrder = roundOrderAscend(rMax);

	globalCtrCheckARIA = 0;
	globalCtrSolARIA = 0;
	globalCtrConstrSboxARIA = 0;
	globalCtrConstrLinARIA = 0;
	globalCtrConstrSSBARIA = 0;
	auto globalStart = chrono::high_resolution_clock::now();
	for(uint i = 0; i < 16; i++){
		for(uint j = 0; j < 16; j++){
			cout << "ARIA Search i = " << i << " j = " << j << endl;
			//Get the list of input/output to test
			auto const allTest = allInOutPairsARIA(i,j,5);

			//Generate the modelInfo object
			modelInfo MI(rMax,sboxModel,linModel,arxModel,firstSbox,lastLin,smartLin,startingRound,sboxCheckType,MCCheckType,SSBCheckType,ARXCheckType,useCallback,maxNumberConstr,roundOrder);

			//Check if we can find a distinguisher
			auto start = chrono::high_resolution_clock::now();
			auto haveDist = checkForDistinguisher(allTest, existTrailARIA, MI, false);
			auto end = chrono::high_resolution_clock::now();
			cout << "Time : " << chrono::duration<double>(end - start).count() << endl;
			if(haveDist){
				cout << "Distinguisher !!!" << endl;
			}
			else
				cout << "No distinguisher..." << endl;

			if(QUICK_LAZY_COUNT_ARIA)
				break;
		}
		if(QUICK_LAZY_COUNT_ARIA)
			break;
	}
	auto globalEnd = chrono::high_resolution_clock::now();
	cout << "Total time : " << chrono::duration<double>(globalEnd - globalStart).count() << endl;
	cout << globalCtrCheckARIA << " input/output pairs checked" << endl;
	cout << globalCtrSolARIA << " solutions examined" << endl;
	cout << globalCtrConstrSboxARIA << " constraints added for the sbox" << endl;
	cout << globalCtrConstrLinARIA  << " constraints added for the lin layer" << endl;
	cout << globalCtrConstrSSBARIA << " constraints added for the SSB" << endl;
}