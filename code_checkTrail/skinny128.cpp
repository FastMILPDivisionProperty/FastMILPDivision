#include "skinny128.hpp"
#include "DivTable.hpp"

using namespace std;

#define DISPLAY_GEN_OUTPUT false
#define DISPLAY_SOLVER_OUTPUT false
#define DISPLAY_NUMBER_ADDED_CONSTRAINT true

vector<vector<pair<vector<uint8_t>, vector<uint8_t>>>> allInOutPairsSKINNY128(unsigned i, unsigned j) {
  vector<uint16_t> sbox = 
    {0x65 , 0x4c , 0x6a , 0x42 , 0x4b , 0x63 , 0x43 , 0x6b , 0x55 , 0x75 , 0x5a , 0x7a , 0x53 , 0x73 , 0x5b , 0x7b ,0x35 , 0x8c , 0x3a , 0x81 , 0x89 , 0x33 , 0x80 , 0x3b , 0x95 , 0x25 , 0x98 , 0x2a , 0x90 , 0x23 , 0x99 , 0x2b ,0xe5 , 0xcc , 0xe8 , 0xc1 , 0xc9 , 0xe0 , 0xc0 , 0xe9 , 0xd5 , 0xf5 , 0xd8 , 0xf8 , 0xd0 , 0xf0 , 0xd9 , 0xf9 ,0xa5 , 0x1c , 0xa8 , 0x12 , 0x1b , 0xa0 , 0x13 , 0xa9 , 0x05 , 0xb5 , 0x0a , 0xb8 , 0x03 , 0xb0 , 0x0b , 0xb9 ,0x32 , 0x88 , 0x3c , 0x85 , 0x8d , 0x34 , 0x84 , 0x3d , 0x91 , 0x22 , 0x9c , 0x2c , 0x94 , 0x24 , 0x9d , 0x2d ,0x62 , 0x4a , 0x6c , 0x45 , 0x4d , 0x64 , 0x44 , 0x6d , 0x52 , 0x72 , 0x5c , 0x7c , 0x54 , 0x74 , 0x5d , 0x7d ,0xa1 , 0x1a , 0xac , 0x15 , 0x1d , 0xa4 , 0x14 , 0xad , 0x02 , 0xb1 , 0x0c , 0xbc , 0x04 , 0xb4 , 0x0d , 0xbd ,0xe1 , 0xc8 , 0xec , 0xc5 , 0xcd , 0xe4 , 0xc4 , 0xed , 0xd1 , 0xf1 , 0xdc , 0xfc , 0xd4 , 0xf4 , 0xdd , 0xfd ,0x36 , 0x8e , 0x38 , 0x82 , 0x8b , 0x30 , 0x83 , 0x39 , 0x96 , 0x26 , 0x9a , 0x28 , 0x93 , 0x20 , 0x9b , 0x29 ,0x66 , 0x4e , 0x68 , 0x41 , 0x49 , 0x60 , 0x40 , 0x69 , 0x56 , 0x76 , 0x58 , 0x78 , 0x50 , 0x70 , 0x59 , 0x79 ,0xa6 , 0x1e , 0xaa , 0x11 , 0x19 , 0xa3 , 0x10 , 0xab , 0x06 , 0xb6 , 0x08 , 0xba , 0x00 , 0xb3 , 0x09 , 0xbb ,0xe6 , 0xce , 0xea , 0xc2 , 0xcb , 0xe3 , 0xc3 , 0xeb , 0xd6 , 0xf6 , 0xda , 0xfa , 0xd3 , 0xf3 , 0xdb , 0xfb ,0x31 , 0x8a , 0x3e , 0x86 , 0x8f , 0x37 , 0x87 , 0x3f , 0x92 , 0x21 , 0x9e , 0x2e , 0x97 , 0x27 , 0x9f , 0x2f ,0x61 , 0x48 , 0x6e , 0x46 , 0x4f , 0x67 , 0x47 , 0x6f , 0x51 , 0x71 , 0x5e , 0x7e , 0x57 , 0x77 , 0x5f , 0x7f ,0xa2 , 0x18 , 0xae , 0x16 , 0x1f , 0xa7 , 0x17 , 0xaf , 0x01 , 0xb2 , 0x0e , 0xbe , 0x07 , 0xb7 , 0x0f , 0xbf ,0xe2 , 0xca , 0xee , 0xc6 , 0xcf , 0xe7 , 0xc7 , 0xef , 0xd2 , 0xf2 , 0xde , 0xfe , 0xd7 , 0xf7 , 0xdf , 0xff};

  vector<uint8_t> permBit({0,1,2,3,4,5,6,7,40,41,42,43,44,45,46,47,80,81,82,83,84,85,86,87,120,121,122,123,124,125,126,127,32,33,34,35,36,37,38,39,72,73,74,75,76,77,78,79,112,113,114,115,116,117,118,119,24,25,26,27,28,29,30,31,64,65,66,67,68,69,70,71,104,105,106,107,108,109,110,111,16,17,18,19,20,21,22,23,56,57,58,59,60,61,62,63,96,97,98,99,100,101,102,103,8,9,10,11,12,13,14,15,48,49,50,51,52,53,54,55,88,89,90,91,92,93,94,95});
  vector<uint8_t> permBit_inv(128);
  for(uint i = 0; i < 128; i++)
  	permBit_inv[permBit[i]] = i;


  auto anf_s = getANF_S(sbox);

  auto in = loadTransitionsIN(anf_s);
  auto out = loadTransitionsOUT(anf_s);

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

GRBModel getModelSkinny128(modelInfo & MI){
/*
	Get the model for Skinny128 over #rMax round
	#lastLin defines if the linear layer is applied on the last round
	#smartLin only has effect if linModel=="Simple" and defines if MC is seen as independent Lboxes or as a single matrix
	#firstSbox defines if the first Sbox layer is applied

	The #sboxModel parameter defines how the Sbox is modelized :
	- "QM" use the Quin-McCluskey algorithm as in Abdelkhalek,Sasaki,Todo,Tolba,Youssef
	- "Simple" use the simplified constraint with PWL

	The #linModel parameter defines how the linear layer is modelized :
	- "Lbox" modelize the linear layer as the parallel application of the same linear Sbox (thanks to the specific form of the matrix)
	- "ZR" modelize the linear layer using the technique from Zhang and Rijmen
	- "CX" modelize the linear layer using the classical copy+xor technique
	- "Simple" modelize the linear layer with the simplified constraint w(x) = w(y)
*/
    string modelName = "Skinny128_"+to_string(MI.rMax)+"r_"+MI.sboxModel+"_"+MI.linModel;
    if(MI.lastLin) modelName += "_lastLin";
    if(MI.linModel=="Simple" && MI.smartLin) modelName += "_smartLin";
    if(!MI.firstSbox) modelName += "_noFirstSbox";
    modelName += ".mps";

    if(!fileExist("./models/"+modelName)){
        string cmd = "sage ./modelGeneration/genSkinny128.sage";
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

modelData getVariablesSkinny128(GRBModel & m,
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
        MD.allTable[0] = readDivTableFromFile("./checkData/Skinny128Sbox_table.bin");
    else if(MI.sboxCheckType == CheckType::Ineq)
        MD.allIneq[0] = readIneqFromFile("./checkData/Skinny128Sbox_ineq.bin");

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
			MD.allMatrix[0] = Matrix("./checkData/matrixSkinny128_small.bin");
    }
    else{
    	if(MI.MCCheckType == CheckType::Matrix)
    		MD.allMatrix[0] = Matrix("./checkData/matrixSkinny128.bin");
	    else if(MI.MCCheckType == CheckType::Table){
	    	cerr << "Error : Table Check for Skinny128 MC (non smartLin) not implemented, probably not worth it compared to using smartLin=true" << endl;
	    	exit(1);
	        // MD.allTable[1] = //Table for MC
	    }
	    else if(MI.MCCheckType == CheckType::Ineq){
	    	cerr << "Error : Ineq Check for Skinny128 MC (non smartLin) not implemented, probably not worth it compared to using smartLin=true" << endl;
	    	exit(1);
	        // MD.allIneq[1]   = //Ineq for MC
	    }
    }

    if(MI.SSBCheckType != CheckType::None){
    	cerr << "Error : SSB Check for Skinny128 not implemented" << endl;
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
	                    MD.allVarsTable.emplace_back(move(in), move(out), tableSbox);
	                else if(MI.sboxCheckType == CheckType::Ineq)
	                    MD.allVarsIneq.emplace_back(move(in), move(out), ineqSbox);
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
		            	for(uint offset = 0; offset < 8; offset++){ //For each set of 4 bits
			                vector<GRBVar> in(4);
			                vector<GRBVar> out(4);
			                for(uint j = 0; j < 4; j++){
			                    in[j] = m.getVarByName("z"+to_string(r)+"_"+to_string(32*col + offset + j*8));
			                    out[j] = m.getVarByName("x"+to_string(r+1)+"_"+to_string(32*col + offset + j*8));
			                }
			                if(MI.MCCheckType == CheckType::Table)
			                    MD.allVarsTable.emplace_back(move(in), move(out), tableMC);
			                else if(MI.MCCheckType == CheckType::Ineq)
			                    MD.allVarsIneq.emplace_back(move(in), move(out), ineqMC);
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
		                vector<GRBVar> in(32);
		                vector<GRBVar> out(32);
		                for(uint j = 0; j < 32; j++){
		                    in[j] = m.getVarByName("z"+to_string(r)+"_"+to_string(32*i+j));
		                    out[j] = m.getVarByName("x"+to_string(r+1)+"_"+to_string(32*i+j));
		                }
		                if(MI.MCCheckType == CheckType::Table)
		                    MD.allVarsTable.emplace_back(move(in), move(out), tableMC);
		                else if(MI.MCCheckType == CheckType::Ineq)
		                    MD.allVarsIneq.emplace_back(move(in), move(out), ineqMC);
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
	                    MD.allVarsTable.emplace_back(move(in), move(out), tableSSB);
	                else if(MI.SSBCheckType == CheckType::Ineq)
	                    MD.allVarsIneq.emplace_back(move(in), move(out), ineqSSB);
	            }
	        }
        }
    }

    return MD;
}

bool existTrailSkinny128(GRBModel & m,
                   modelData const & MD,
                   vector<uint8_t> const & input,
                   vector<uint8_t> const & output,
                   bool const useCallback,
                   uint const maxNumberConstr){
/*
    Return true if there is a trail from #input to #output using #model for Skinny128
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
        }

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

bool existTrailSkinny128(modelInfo & MI,
                   vector<uint8_t> const & input,
                   vector<uint8_t> const & output){
/*
    Return true if there is a trail from #input to #output for Skinny128
    #MI should contain the necessary information to generate the model and parametrize the solving
*/

    GRBModel m = getModelSkinny128(MI);
    auto MD = getVariablesSkinny128(m, MI);

    return existTrailSkinny128(m,MD,input,output,MI.useCallback,MI.maxNumberConstr);
}

void genSboxFileSkinny128(){
	//Generate the tables and ineq files for Skinny 128

	vector<uint32_t> S(
		{0x65,0x4c,0x6a,0x42,0x4b,0x63,0x43,0x6b,0x55,0x75,0x5a,0x7a,0x53,0x73,0x5b,0x7b,
		 0x35,0x8c,0x3a,0x81,0x89,0x33,0x80,0x3b,0x95,0x25,0x98,0x2a,0x90,0x23,0x99,0x2b,
		 0xe5,0xcc,0xe8,0xc1,0xc9,0xe0,0xc0,0xe9,0xd5,0xf5,0xd8,0xf8,0xd0,0xf0,0xd9,0xf9,
		 0xa5,0x1c,0xa8,0x12,0x1b,0xa0,0x13,0xa9,0x05,0xb5,0x0a,0xb8,0x03,0xb0,0x0b,0xb9,
		 0x32,0x88,0x3c,0x85,0x8d,0x34,0x84,0x3d,0x91,0x22,0x9c,0x2c,0x94,0x24,0x9d,0x2d,
		 0x62,0x4a,0x6c,0x45,0x4d,0x64,0x44,0x6d,0x52,0x72,0x5c,0x7c,0x54,0x74,0x5d,0x7d,
		 0xa1,0x1a,0xac,0x15,0x1d,0xa4,0x14,0xad,0x02,0xb1,0x0c,0xbc,0x04,0xb4,0x0d,0xbd,
		 0xe1,0xc8,0xec,0xc5,0xcd,0xe4,0xc4,0xed,0xd1,0xf1,0xdc,0xfc,0xd4,0xf4,0xdd,0xfd,
		 0x36,0x8e,0x38,0x82,0x8b,0x30,0x83,0x39,0x96,0x26,0x9a,0x28,0x93,0x20,0x9b,0x29,
		 0x66,0x4e,0x68,0x41,0x49,0x60,0x40,0x69,0x56,0x76,0x58,0x78,0x50,0x70,0x59,0x79,
		 0xa6,0x1e,0xaa,0x11,0x19,0xa3,0x10,0xab,0x06,0xb6,0x08,0xba,0x00,0xb3,0x09,0xbb,
		 0xe6,0xce,0xea,0xc2,0xcb,0xe3,0xc3,0xeb,0xd6,0xf6,0xda,0xfa,0xd3,0xf3,0xdb,0xfb,
		 0x31,0x8a,0x3e,0x86,0x8f,0x37,0x87,0x3f,0x92,0x21,0x9e,0x2e,0x97,0x27,0x9f,0x2f,
		 0x61,0x48,0x6e,0x46,0x4f,0x67,0x47,0x6f,0x51,0x71,0x5e,0x7e,0x57,0x77,0x5f,0x7f,
		 0xa2,0x18,0xae,0x16,0x1f,0xa7,0x17,0xaf,0x01,0xb2,0x0e,0xbe,0x07,0xb7,0x0f,0xbf,
		 0xe2,0xca,0xee,0xc6,0xcf,0xe7,0xc7,0xef,0xd2,0xf2,0xde,0xfe,0xd7,0xf7,0xdf,0xff});

	genAndSaveTableIneq(vector<vector<uint32_t>>({S}),8,8,"./checkData/Skinny128Sbox",true);
}

void createMatrixFileSkinny128(){
	//Create the matrix file ./checkData/matrixSkinny128.bin

	vector<vector<uint8_t>> val({{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
								 {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0},
								 {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0},
								 {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0},
								 {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0},
								 {0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0},
								 {0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0},
								 {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1},
								 {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
								 {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
								 {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
								 {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
								 {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
								 {0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
								 {0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
								 {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
								 {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
								 {0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
								 {0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
								 {0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0},
								 {0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0},
								 {0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0},
								 {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0},
								 {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0},
								 {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
								 {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
								 {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
								 {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0},
								 {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0},
								 {0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0},
								 {0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0},
								 {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0}});
	Matrix M(32,32);
	for(uint i = 0; i < 32; i++){
		for(uint j = 0; j < 32; j++)
			M.set(i,j,val[i][j]);
	}

	M.saveToFile("./checkData/matrixSkinny128.bin");

	val = vector<vector<uint8_t>>({{1,0,1,1},
								   {1,0,0,0},
								   {0,1,1,0},
								   {1,0,1,0}});

	M = Matrix(4,4);
	for(uint i = 0; i < 4; i++){
		for(uint j = 0; j < 4; j++)
			M.set(i,j,val[i][j]);
	}
	M.saveToFile("./checkData/matrixSkinny128_small.bin");
}

void checkNormalDistinguisherSkinny128(){
	/*
		Search for a normal (no linear combination) distinguisher over 10 rounds of Skinny128
	*/

	vector<uint8_t> input({1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1});

	//These parameters have no effect for Skinny128 but are required for the modelInfo object
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

	vector<string> listSboxModel({"QM", "Simple"});
	vector<string> listLinModel({"Lbox", "ZR", "CX", "Simple"});
	vector<bool> listUseCallback({false,true});
	vector<vector<uint>> listroundOrder({roundOrderOutsideIn(rMax), roundOrderAscend(rMax), roundOrderDescend(rMax)});
	vector<uint> listMaxNbConstr({100,50,25,10,1});

	for(auto const & useCallback : listUseCallback){
	 for(auto const & roundOrder : listroundOrder){
	  for(auto const & linModel : listLinModel){
	   vector<CheckType> listMCCheckType({CheckType::Ineq, CheckType::Table, CheckType::Matrix});
	   if(linModel == "Lbox" || linModel == "ZR") //No need for checks on MC
	    listMCCheckType = vector<CheckType>({CheckType::None});
	   for(auto const & MCCheckType : listMCCheckType){
	   	for(auto const & sboxModel : listSboxModel){
	     vector<CheckType> listSboxCheckType({CheckType::Ineq, CheckType::Table});
	     if(sboxModel == "Hull") //No need for checks on Sbox
	      listSboxCheckType = vector<CheckType>({CheckType::None});
	     for(auto const & sboxCheckType : listSboxCheckType){
	       for(auto const & maxNumberConstr : listMaxNbConstr){

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
			auto durBalanced = start - start;
			//Check each output bit (half balanced)
			vector<bool> balancedBits(128);
			uint ctrBalanced = 0;
			for(uint i = 0; i < 128; i++){
				vector<uint8_t> output(128,0);
				output[i] = 1;
				auto startBit = chrono::high_resolution_clock::now();
				balancedBits[i] = !existTrailSkinny128(MI,input,output);
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
	}}}}}}}
}

void timingTestSearchNewSkinny128(){
	/*
	*/

	//These parameters have no effect for Skinny128 but are required for the modelInfo object
	string arxModel = "ARX";
	unsigned int startingRound = 1;
	CheckType ARXCheckType = CheckType::None;

	//Fixed parameters
	CheckType SSBCheckType = CheckType::None; //No SSB check, 32-bit SSB :'(
	uint rMax = 9;
	bool smartLin = true; //probably not useful to not have it honestly
	bool firstSbox = false;
	bool lastLin = true;

	//Variable parameters
	// string sboxModel = "QM";
	string linModel = "Lbox";
	CheckType sboxCheckType = CheckType::None;
	CheckType MCCheckType = CheckType::None;
	// bool useCallback = false;
	// uint maxNumberConstr = 4;
	// vector<uint> roundOrder = roundOrderOutsideIn(rMax);

	vector<string> listSboxModel({"QM", "Simple"});
	vector<string> listLinModel({"Lbox"}); //Simple and CX are just really really bad here. ZR seems worse than Lbox
	vector<bool> listUseCallback({true,false});
	vector<vector<uint>> listroundOrder({roundOrderOutsideIn(rMax), roundOrderAscend(rMax), roundOrderDescend(rMax)});
	vector<uint> listMaxNbConstr({500,100,25,10,1});
	// vector<CheckType> listSSBCheckType({CheckType::None}); //No SSB check, 32-bit SSB :'(
	auto const allTest = allInOutPairsSKINNY128(1,0);

	for(auto const & useCallback : listUseCallback){
	 for(auto const & sboxModel : listSboxModel){
	  if(sboxModel == "QM"){ //No need for checks on Sbox
		sboxCheckType = CheckType::None;
		listMaxNbConstr = vector<uint>({1});
		listroundOrder = vector<vector<uint>>({roundOrderOutsideIn(rMax)});
	  }
	  else if(sboxModel == "Simple"){
		sboxCheckType = CheckType::Ineq;
		listMaxNbConstr = vector<uint>({100});
		listroundOrder = vector<vector<uint>>({roundOrderOutsideIn(rMax), roundOrderAscend(rMax), roundOrderDescend(rMax)});
	  }
	  for(auto const & roundOrder : listroundOrder){
	   for(auto const & maxNumberConstr : listMaxNbConstr){
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
			auto haveDist = checkForDistinguisher(allTest, existTrailSkinny128, MI, true);
			auto end = chrono::high_resolution_clock::now();

			cout << "Time : " << chrono::duration<double>(end - start).count() << endl;
			if(haveDist)cout << "Distinguisher !!!" << endl;
			else cout << "No distinguisher..." << endl;
	}}}}
}


void searchDistinguisherSkinny128(){


	string arxModel = "ARX";
	unsigned int startingRound = 1;
	CheckType ARXCheckType = CheckType::None;

	CheckType SSBCheckType = CheckType::None;
	uint rMax = 9; 
	bool smartLin = true;
	bool firstSbox = false;
	bool lastLin = true;

	string sboxModel = "QM";
	string linModel = "Lbox";
	CheckType sboxCheckType = CheckType::None;
	CheckType MCCheckType = CheckType::None;
	bool useCallback = false;
	uint maxNumberConstr = 100;
	vector<uint> roundOrder = roundOrderAscend(rMax);

	auto globalStart = chrono::high_resolution_clock::now();
	for(uint i = 0; i < 16; i++){
		for(uint j = 0; j < 16; j++){
			cout << "Skinny128 Search i = " << i << " j = " << j << endl;
			//Get the list of input/output to test
			auto const allTest = allInOutPairsSKINNY128(i,j);

			//Generate the modelInfo object
			modelInfo MI(rMax,sboxModel,linModel,arxModel,firstSbox,lastLin,smartLin,startingRound,sboxCheckType,MCCheckType,SSBCheckType,ARXCheckType,useCallback,maxNumberConstr,roundOrder);

			//Check if we can find a distinguisher
			auto start = chrono::high_resolution_clock::now();
			auto haveDist = checkForDistinguisher(allTest, existTrailSkinny128, MI);
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

}