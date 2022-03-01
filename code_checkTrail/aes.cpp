#include "aes.hpp"
#include "DivTable.hpp"

using namespace std;

#define DISPLAY_GEN_OUTPUT false
#define DISPLAY_SOLVER_OUTPUT false
#define DISPLAY_NUMBER_ADDED_CONSTRAINT false

// map<uint64_t, uint64_t> globalHisto;
// uint64_t globalCtrSolAES = 0;
// uint64_t globalCtrConstr = 0;

vector<vector<pair<vector<uint8_t>, vector<uint8_t>>>> allInOutPairsAES(unsigned i, unsigned j) {
  vector<uint16_t> sbox = {
  //0     1    2      3     4    5     6     7      8    9     A      B    C     D     E     F
  0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
  0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
  0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
  0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
  0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
  0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
  0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
  0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
  0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
  0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
  0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
  0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
  0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
  0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
  0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
  0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16 };

  auto anf_s = getANF_S(sbox);

  auto in = loadTransitionsIN(anf_s);
  auto out = loadTransitionsOUT(anf_s);

  vector<uint8_t> permBit({0,1,2,3,4,5,6,7,104,105,106,107,108,109,110,111,80,81,82,83,84,85,86,87,56,57,58,59,60,61,62,63,32,33,34,35,36,37,38,39,8,9,10,11,12,13,14,15,112,113,114,115,116,117,118,119,88,89,90,91,92,93,94,95,64,65,66,67,68,69,70,71,40,41,42,43,44,45,46,47,16,17,18,19,20,21,22,23,120,121,122,123,124,125,126,127,96,97,98,99,100,101,102,103,72,73,74,75,76,77,78,79,48,49,50,51,52,53,54,55,24,25,26,27,28,29,30,31});
  vector<uint8_t> permBit_inv(128);
  for(uint i = 0; i < 128; i++)
  	permBit_inv[permBit[i]] = i;

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
          //Permute out to fit the modeling
          vector<uint8_t> pout_8t(128);
          for(uint i = 0; i < 128; i++)
          	pout_8t[permBit_inv[i]] = out_8t[i];
          tmp.emplace_back(move(in_8t), move(pout_8t));
        }
      }
      res.emplace_back(move(tmp));
    }
  }

  return res;

}


GRBModel getModelAES(modelInfo & MI){
/*
	Get the model for AES over #rMax round
	#lastLin defines if the linear layer is applied on the last round
	#firstSbox defines if the first Sbox layer is applied

	The #sboxModel parameter defines how the Sbox is modelized :
	- "QM" use the Quin-McCluskey algorithm as in Abdelkhalek,Sasaki,Todo,Tolba,Youssef
	- "Simple" use the simplified constraint with PWL

	The #linModel parameter defines how the linear layer is modelized :
	- "CX" modelize the linear layer using the classical copy+xor technique
	- "Simple" modelize the linear layer with the simplified constraint w(x) = w(y)
*/
    string modelName = "AES_"+to_string(MI.rMax)+"r_"+MI.sboxModel+"_"+MI.linModel;
    if(MI.lastLin) modelName += "_lastLin";
    if(!MI.firstSbox) modelName += "_noFirstSbox";
    modelName += ".mps";

    if(!fileExist("./models/"+modelName)){
        string cmd = "sage ./modelGeneration/genAES.sage";
        cmd += " -r " + to_string(MI.rMax);
        cmd += " -s " + MI.sboxModel;
        cmd += " -m " + MI.linModel;
        if(MI.lastLin) cmd += " --lastLin True";
        else cmd += " --lastLin False";
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

modelData getVariablesAES(GRBModel & m,
                          modelInfo const & MI){
     /*
    Read the model and return the corresponding blocks of variables
    - #MI defines various parameters for how the trails are checked (see aux_function.hpp)
    */

    modelData MD;
    MD.rMax = MI.rMax;
    MD.firstSbox = MI.firstSbox;
    MD.lastLin = MI.lastLin;

    MD.allTable = vector<vector<vector<uint32_t>>>(3);
    MD.allIneq = vector<vector<vector<int>>>(3);
    MD.allMatrix = vector<Matrix>(1);

    if(MI.sboxCheckType == CheckType::Table) //Sbox table
        MD.allTable[0] = readDivTableFromFile("./checkData/AESSbox_table.bin");
    else if(MI.sboxCheckType == CheckType::Ineq)
        MD.allIneq[0] = readIneqFromFile("./checkData/AESSbox_ineq.bin");

    if(MI.MCCheckType == CheckType::Matrix)
    	MD.allMatrix[0] = Matrix("./checkData/matrixAES.bin");
    else if(MI.MCCheckType == CheckType::Table){
    	cerr << "Error : Table Check for AES MC not implemented" << endl;
    	exit(1);
        // MD.allTable[1] = //Table for MC
    }
    else if(MI.MCCheckType == CheckType::Ineq){
    	cerr << "Error : Ineq Check for AES MC not implemented" << endl;
    	exit(1);
        // MD.allIneq[1]   = //Ineq for MC
    }

    if(MI.SSBCheckType != CheckType::None){
    	cerr << "Error : SSB Check for AES not implemented" << endl;
    	exit(1);
    }
    // if(SSBCheckType == CheckType::Table)
    //     MD.allTable[2] = //Table for SSB
    // else if(SSBCheckType == CheckType::Ineq)
    //     MD.allIneq[2]   = //Ineq for SSB

    auto & tableSbox = MD.allTable[0];
    auto & tableMC   = MD.allTable[1];
    auto & tableSSB  = MD.allTable[2];
    auto & ineqSbox  = MD.allIneq[0];
    auto & ineqMC    = MD.allIneq[1];
    auto & ineqSSB   = MD.allIneq[2];

    auto & MCmatrix = MD.allMatrix[0];

    //Sbox variables
    uint roundFirstSbox = 0;
    if(!MI.firstSbox) roundFirstSbox = 1;

    if(MI.sboxCheckType != CheckType::None){
        for(auto const & r : MI.roundOrder){
        	if(r >= roundFirstSbox){
	            for(uint i = 0; i < 16; i++){
	                vector<GRBVar> in(8);
	                vector<GRBVar> out(8);
	                for(uint j = 0; j < 8; j++){
	                    in[j] = m.getVarByName("x"+to_string(r)+"_"+to_string(8*i+j));
	                    out[j] = m.getVarByName("y"+to_string(r)+"_"+to_string(8*i+j));
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
	                vector<GRBVar> in(32);
	                vector<GRBVar> out(32);
	                for(uint j = 0; j < 32; j++){
	                    in[j] = m.getVarByName("z"+to_string(r)+"_"+to_string(32*i+j));
	                    out[j] = m.getVarByName("x"+to_string(r+1)+"_"+to_string(32*i+j));
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
	                vector<GRBVar> in(32);
	                vector<GRBVar> out(32);
	                uint ctr = 0;
	                for(auto const j : indexIn[i]){
	                    for(uint k = 0; k < 8; k++){
	                        in[ctr] = m.getVarByName("x"+to_string(r)+"_"+to_string(8*j+k));
	                        ctr++;
	                    }
	                }
	                ctr = 0;
	                for(auto const j : indexOut[i]){
	                    for(uint k = 0; k < 8; k++){
	                        out[ctr] = m.getVarByName("y"+to_string(r+1)+"_"+to_string(8*j+k));
	                        ctr++;
	                    }
	                }
	                if(MI.SSBCheckType == CheckType::Table)
	                    MD.allVarsTable.emplace_back(move(in), move(out), tableSSB,ConstrType::SSB);
	                else if(MI.SSBCheckType == CheckType::Ineq)
	                    MD.allVarsIneq.emplace_back(move(in), move(out), ineqSSB,ConstrType::SSB);
	            }
	        }
        }
    }

    return MD;
}

uint64_t globalCtrCheckAES = 0;
uint64_t globalCtrSolAES = 0;
uint64_t globalCtrConstrSboxAES = 0;
uint64_t globalCtrConstrLinAES = 0;
uint64_t globalCtrConstrSSBAES = 0;
bool existTrailAES(GRBModel & m,
                   modelData const & MD,
                   vector<uint8_t> const & input,
                   vector<uint8_t> const & output,
                   bool const useCallback,
                   uint const maxNumberConstr){
/*
    Return true if there is a trail from #input to #output using #model for AES
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
        CustomCallback cb(MD.allVarsTable, MD.allVarsIneq, MD.allVarsMatrix, maxNumberConstr);
        m.setCallback(&cb);

        m.update();
        m.optimize();
        if(DISPLAY_NUMBER_ADDED_CONSTRAINT){
            cout << cb.ctrSol << " solutions examined, ";
            cout << cb.ctrGlobalConstr << " constraints added" << endl;
            // cout << cb.ctrNonInvMinor << " non-invertible minors" << endl;
            // cout << "Distribution of non-invertible minors size :" << endl;
            // for(auto const & sn : cb.histoNonInvMinors)
            // 	cout << sn.first << " : " << sn.second << endl;
            
            // globalCtrSolAES += cb.ctrSol;
            // globalCtrConstr += cb.ctrGlobalConstr;
            // for(auto const & sn : cb.histoNonInvMinors)
            // 	globalHisto[sn.first] += sn.second;
        }
        globalCtrCheckAES++;
        globalCtrSolAES += cb.ctrSol;
        globalCtrConstrSboxAES += cb.ctrConstrSbox;
		globalCtrConstrLinAES += cb.ctrConstrLin;
		globalCtrConstrSSBAES += cb.ctrConstrSSB;

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

bool existTrailAES(modelInfo & MI,
                   vector<uint8_t> const & input,
                   vector<uint8_t> const & output){
/*
    Return true if there is a trail from #input to #output for AES
    #MI should contain the necessary information to generate the model and parametrize the solving
*/

    GRBModel m = getModelAES(MI);
    auto MD = getVariablesAES(m, MI);

    return existTrailAES(m,MD,input,output,MI.useCallback,MI.maxNumberConstr);
}

void genSboxFileAES(){
	//Generate the tables and ineq files for AES

	vector<uint32_t> S(
		{0x63,0x7c,0x77,0x7b,0xf2,0x6b,0x6f,0xc5,0x30,0x01,0x67,0x2b,0xfe,0xd7,0xab,0x76,
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
		 0x8c,0xa1,0x89,0x0d,0xbf,0xe6,0x42,0x68,0x41,0x99,0x2d,0x0f,0xb0,0x54,0xbb,0x16});

	genAndSaveTableIneq(vector<vector<uint32_t>>({S}),8,8,"./checkData/AESSbox",true);
}

void createMatrixFileAES(){
	//Create the matrix file ./checkData/matrixAES.bin

	vector<vector<uint8_t>> val({{0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
								 {1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0},
								 {0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0},
								 {0,0,1,0,0,0,0,1,0,0,1,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0},
								 {0,0,0,1,0,0,0,1,0,0,0,1,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0},
								 {0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0},
								 {0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0},
								 {0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1},
								 {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0},
								 {0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0},
								 {0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0},
								 {0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,1,0,0,0,1,0,0,0,1,0,0,0,0},
								 {0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,1,0,0,1,0,0,0,0,1,0,0,0},
								 {0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0},
								 {0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0},
								 {0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1},
								 {1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1},
								 {0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,1},
								 {0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0},
								 {0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,1,0,0,0,1},
								 {0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,1,0,0,1},
								 {0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0},
								 {0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0},
								 {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1},
								 {1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
								 {1,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1},
								 {0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0},
								 {0,0,1,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1},
								 {0,0,0,1,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1},
								 {0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0},
								 {0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0},
								 {0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0}});
	Matrix M(32,32);
	for(uint i = 0; i < 32; i++){
		for(uint j = 0; j < 32; j++)
			M.set(i,j,val[i][j]);
	}

	M.saveToFile("./checkData/matrixAES.bin");
}

void checkNormalDistinguisherAES(){
	/*
		Search again for a known (no linear combination) distinguisher over AES to compare execution time
	*/

	vector<uint8_t> input({1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1});

	//These parameters have no effect for AES but are required for the modelInfo object
	string arxModel = "ARX";
	CheckType ARXCheckType = CheckType::None;
	bool smartLin = true;
	unsigned int startingRound = 1;

	//Fixed parameters
	uint rMax = 4;
	bool firstSbox = true;
	bool lastLin = false;
	CheckType sboxCheckType = CheckType::None;
	CheckType MCCheckType = CheckType::None;
	CheckType SSBCheckType = CheckType::None;
	bool useCallback = false;
	uint maxNumberConstr = 1;
	vector<uint> roundOrder = roundOrderAscend(rMax);

	//Variable parameters
	//Even with the simple model, no need to use the checks, so limited number of variable parameters
	vector<string> listSboxModel({"QM", "Simple"});
	vector<string> listLinModel({"CX", "Simple"});

	//The searchs
	for(auto const & sboxModel : listSboxModel){
		for(auto const & linModel : listLinModel){

			globalCtrCheckAES = 0;
			globalCtrSolAES = 0;
			globalCtrConstrSboxAES = 0;
			globalCtrConstrLinAES = 0;
			globalCtrConstrSSBAES = 0;

			modelInfo MI(rMax,sboxModel,linModel,arxModel,firstSbox,lastLin,smartLin,startingRound,sboxCheckType,MCCheckType,SSBCheckType,ARXCheckType,useCallback,maxNumberConstr,roundOrder);

			auto start = chrono::high_resolution_clock::now();
			//Check each output bit (all balanced)
			uint ctrBalanced = 0;
			for(uint i = 0; i < 128; i++){
				vector<uint8_t> output(128,0);
				output[i] = 1;
				if(!existTrailAES(MI,input,output)) ctrBalanced++;
			}
			auto end = chrono::high_resolution_clock::now();
			cout << "sboxModel : " << sboxModel << " - linModel : " << linModel << endl;
			cout << ctrBalanced << " balanced bits" << endl;
			cout << chrono::duration<double>(end - start).count() << " seconds" << endl;
			cout << globalCtrCheckAES << " input/output pairs checked" << endl;
			cout << globalCtrSolAES << " solutions examined" << endl;
			cout << globalCtrConstrSboxAES << " constraints added for the sbox" << endl;
			cout << globalCtrConstrLinAES  << " constraints added for the lin layer" << endl;
			cout << globalCtrConstrSSBAES << " constraints added for the SSB" << endl;
		}
	}
}

void timingSearchAES(){
	//These parameters have no effect for AES but are required for the modelInfo object
	string arxModel = "ARX";
	unsigned int startingRound = 1;
	CheckType ARXCheckType = CheckType::None;
	CheckType SSBCheckType = CheckType::None; //No SSB for AES

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
	vector<string> listLinModel({"CX", "Simple"});
	vector<bool> listUseCallback({true});
	vector<vector<uint>> listroundOrderIfNeeded({roundOrderOutsideIn(rMax), roundOrderAscend(rMax), roundOrderDescend(rMax)});
	auto const allTest = allInOutPairsAES(0,0);

	uint ctrParameterSet = 0;
	uint firstParameterSet = 0;

	for(auto const & useCallback : listUseCallback){
	 for(auto const & linModel : listLinModel){
	  vector<CheckType> listMCCheckType({CheckType::Matrix});
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

	      	 	 globalCtrCheckAES = 0;
				 globalCtrSolAES = 0;
	      	 	 globalCtrConstrSboxAES = 0;
				 globalCtrConstrLinAES = 0;
				 globalCtrConstrSSBAES = 0;

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
				 auto haveDist = checkForDistinguisher(allTest, existTrailAES, MI, false);
				 auto end = chrono::high_resolution_clock::now();
	 
				 cout << "Time : " << chrono::duration<double>(end - start).count() << endl;
				 if(haveDist)cout << "Distinguisher !!!" << endl << endl;
				 else cout << "No distinguisher..." << endl << endl;
				 cout << globalCtrCheckAES << " input/output pairs checked" << endl;
				 cout << globalCtrSolAES << " solutions examined" << endl;
				 cout << globalCtrConstrSboxAES << " constraints added for the sbox" << endl;
				 cout << globalCtrConstrLinAES  << " constraints added for the lin layer" << endl;
				 cout << globalCtrConstrSSBAES << " constraints added for the SSB" << endl;
			 }
			 ctrParameterSet++;
	}}}}}}}
}

void searchDistinguisher4RoundsAES(){

	vector<uint8_t> input({1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1});

	//These parameters have no effect for AES but are required for the modelInfo object
	string arxModel = "ARX";
	CheckType ARXCheckType = CheckType::None;
	bool smartLin = true;
	unsigned int startingRound = 1;

	//Fixed parameters
	uint rMax = 4;
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

	globalCtrCheckAES = 0;
	globalCtrSolAES = 0;
	globalCtrConstrSboxAES = 0;
	globalCtrConstrLinAES = 0;
	globalCtrConstrSSBAES = 0;

	modelInfo MI(rMax,sboxModel,linModel,arxModel,firstSbox,lastLin,smartLin,startingRound,sboxCheckType,MCCheckType,SSBCheckType,ARXCheckType,useCallback,maxNumberConstr,roundOrder);

	auto start = chrono::high_resolution_clock::now();
	//Check each output bit (all balanced)
	uint ctrBalanced = 0;
	for(uint i = 0; i < 128; i++){
		vector<uint8_t> output(128,0);
		output[i] = 1;
		if(!existTrailAES(MI,input,output)) ctrBalanced++;
	}
	auto end = chrono::high_resolution_clock::now();
	cout << "sboxModel : " << sboxModel << " - linModel : " << linModel << endl;
	cout << ctrBalanced << " balanced bits" << endl;
	cout << chrono::duration<double>(end - start).count() << " seconds" << endl;
	cout << globalCtrCheckAES << " input/output pairs checked" << endl;
	cout << globalCtrSolAES << " solutions examined" << endl;
	cout << globalCtrConstrSboxAES << " constraints added for the sbox" << endl;
	cout << globalCtrConstrLinAES  << " constraints added for the lin layer" << endl;
	cout << globalCtrConstrSSBAES << " constraints added for the SSB" << endl;

}

void searchDistinguisher5RoundsAES(){
	/*
		Search for a distinguisher over 5 rounds AES
		Parameters are fixed and were determined with experiments on a subset of cases
	*/

	uint rMax = 4; //1 partial round at the start (only MC) + 3 full rounds
	string sboxModel = "QM";
	string linModel = "Simple";
	string arxModel = "ARX";
	bool firstSbox = false;
	bool lastLin = true;
	bool smartLin = true;
	unsigned int startingRound = 1;
	CheckType sboxCheckType = CheckType::None;
	CheckType MCCheckType = CheckType::Matrix;
	CheckType SSBCheckType = CheckType::None;
	CheckType ARXCheckType = CheckType::None;
	bool useCallback = true;
	uint maxNumberConstr = 10;
	vector<uint> roundOrder = roundOrderDescend(rMax);

	globalCtrCheckAES = 0;
	globalCtrSolAES = 0;
	globalCtrConstrSboxAES = 0;
	globalCtrConstrLinAES = 0;
	globalCtrConstrSSBAES = 0;
	auto globalStart = chrono::high_resolution_clock::now();
	for(uint i = 0; i < 16; i++){
		for(uint j = 0; j < 4; j++){
			cout << "AES Search i = " << i << " j = " << j << endl;
			//Get the list of input/output to test
			auto const allTest = allInOutPairsAES(i,j);

			//Generate the modelInfo object
			modelInfo MI(rMax,sboxModel,linModel,arxModel,firstSbox,lastLin,smartLin,startingRound,sboxCheckType,MCCheckType,SSBCheckType,ARXCheckType,useCallback,maxNumberConstr,roundOrder);

			//Check if we can find a distinguisher
			auto start = chrono::high_resolution_clock::now();
			auto haveDist = checkForDistinguisher(allTest, existTrailAES, MI, false);
			auto end = chrono::high_resolution_clock::now();
			cout << "Time : " << chrono::duration<double>(end - start).count() << endl;
			if(haveDist){
				cout << "Distinguisher !!!" << endl;
			}
			else
				cout << "No distinguisher..." << endl;
		}
	}
	auto globalEnd = chrono::high_resolution_clock::now();
	cout << "Total time : " << chrono::duration<double>(globalEnd - globalStart).count() << endl;
	cout << globalCtrCheckAES << " input/output pairs checked" << endl;
	cout << globalCtrSolAES << " solutions examined" << endl;
	cout << globalCtrConstrSboxAES << " constraints added for the sbox" << endl;
	cout << globalCtrConstrLinAES  << " constraints added for the lin layer" << endl;
	cout << globalCtrConstrSSBAES << " constraints added for the SSB" << endl;
}