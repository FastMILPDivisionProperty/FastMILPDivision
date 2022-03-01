#include "camellia.hpp"

using namespace std;

#define DISPLAY_GEN_OUTPUT false
#define DISPLAY_SOLVER_OUTPUT false
#define DISPLAY_NUMBER_ADDED_CONSTRAINT false
#define QUICK_LAZY_COUNT_CAMELLIA false

GRBModel getModelCamellia(modelInfo & MI){
/*
	Get the model for Camellia over #rMax round
	#startingRound defines the index of the first round

	The #sboxModel parameter defines how the Sbox is modelized :
	- "QM" use the Quin-McCluskey algorithm as in Abdelkhalek,Sasaki,Todo,Tolba,Youssef
	- "Simple" use the simplified constraint with PWL

	The #linModel parameter defines how the linear layer is modelized :
	- "CX" modelize the linear layer using the classical copy+xor technique
	- "Simple" modelize the linear layer with the simplified constraint w(x) = w(y)
*/
    string modelName = "Camellia_"+to_string(MI.rMax)+"r_"+MI.sboxModel+"_"+MI.linModel;
    modelName += "_start"+to_string(MI.startingRound);
    modelName += ".mps";

    if(!fileExist("./models/"+modelName)){
        string cmd = "sage ./modelGeneration/genCamellia.sage";
        cmd += " -r " + to_string(MI.rMax);
        cmd += " -s " + MI.sboxModel;
        cmd += " -m " + MI.linModel;
        cmd += " --startingRound " + to_string(MI.startingRound);
        auto sysstdout = exec(cmd.c_str());
        if(DISPLAY_GEN_OUTPUT) cout << sysstdout << endl;
    }

    if(!DISPLAY_SOLVER_OUTPUT)
        MI.gurobiEnv.set(GRB_IntParam_OutputFlag,0);

    GRBModel m(MI.gurobiEnv, "./models/"+modelName);
    return m;
}

modelData getVariablesCamellia(GRBModel & m,
                           modelInfo const & MI){
     /*
    Read the model and return the corresponding blocks of variables
    - #MI defines various parameters for how the trails are checked (see aux_function.hpp)
    */

    modelData MD;
    MD.rMax = MI.rMax;

    //Camellia uses 4 different sboxes but no SSB (well, too big)
    MD.allTable = vector<vector<vector<uint32_t>>>(7);
    MD.allIneq = vector<vector<vector<int>>>(7);
    MD.allMatrix = vector<Matrix>(3);

    if(MI.sboxCheckType == CheckType::Table){ //Sbox table
        MD.allTable[0] = readDivTableFromFile("./checkData/CamelliaS1Sbox_table.bin"); //S1
        MD.allTable[1] = readDivTableFromFile("./checkData/CamelliaS2Sbox_table.bin"); //S2
        MD.allTable[2] = readDivTableFromFile("./checkData/CamelliaS3Sbox_table.bin"); //S3
        MD.allTable[3] = readDivTableFromFile("./checkData/CamelliaS4Sbox_table.bin"); //S4
    }
    else if(MI.sboxCheckType == CheckType::Ineq){
        MD.allIneq[0] = readIneqFromFile("./checkData/CamelliaS1Sbox_ineq.bin"); //S1
        MD.allIneq[1] = readIneqFromFile("./checkData/CamelliaS2Sbox_ineq.bin"); //S2
        MD.allIneq[2] = readIneqFromFile("./checkData/CamelliaS3Sbox_ineq.bin"); //S3
        MD.allIneq[3] = readIneqFromFile("./checkData/CamelliaS4Sbox_ineq.bin"); //S4
    }

    if(MI.MCCheckType == CheckType::Matrix){
    	MD.allMatrix[0] = Matrix("./checkData/matrixCamelliaP.bin");
    	MD.allMatrix[1] = Matrix("./checkData/matrixCamelliaFL.bin");
    	MD.allMatrix[2] = Matrix("./checkData/matrixCamelliainvFL.bin");
    }
    else if(MI.MCCheckType == CheckType::Table){
    	cerr << "Error : Table Check for Camellia MC not implemented, unlikely to be" << endl;
    	exit(1);
        // MD.allTable[4] = //Table for P
        // MD.allTable[5] = //Table for FL
        // MD.allTable[6] = //Table for invFL
    }
    else if(MI.MCCheckType == CheckType::Ineq){
    	cerr << "Error : Ineq Check for Camellia MC not implemented, unlikely to be" << endl;
    	exit(1);
        // MD.allIneq[4]   = //Ineq for P
        // MD.allIneq[5]   = //Ineq for FL
        // MD.allIneq[6]   = //Ineq for invFL
    }


    auto & tableS1 = MD.allTable[0];
    auto & tableS2 = MD.allTable[1];
    auto & tableS3 = MD.allTable[2];
    auto & tableS4 = MD.allTable[3];

    auto & ineqS1  = MD.allIneq[0];
    auto & ineqS2  = MD.allIneq[1];
    auto & ineqS3  = MD.allIneq[2];
    auto & ineqS4  = MD.allIneq[3];

    auto & tableP     = MD.allTable[4];
    auto & tableFL    = MD.allTable[5];
    auto & tableinvFL = MD.allTable[6];

    auto & ineqP	 = MD.allIneq[4];
    auto & ineqFL    = MD.allIneq[5];
    auto & ineqinvFL = MD.allIneq[6];

    auto & matrixP = MD.allMatrix[0];
    auto & matrixFL = MD.allMatrix[1];
    auto & matrixinvFL = MD.allMatrix[2];

    //Sbox variables
    if(MI.sboxCheckType != CheckType::None){
        for(auto const & r : MI.roundOrder){
    		//Sbox layer s1,s2,s3,s4,s2,s3,s4,s1
    		for(uint i = 0; i < 8; i++){
    			vector<GRBVar> in(8);
				vector<GRBVar> out(8);
				for(uint j = 0; j < 8; j++){
				    in[j] = m.getVarByName("cx"+to_string(r)+"_"+to_string(8*i+j));
				    out[j] = m.getVarByName("y"+to_string(r)+"_"+to_string(8*i+j));
				}

				if(i == 0 || i == 7){ //S1
					if(MI.sboxCheckType == CheckType::Table)
	                    MD.allVarsTable.emplace_back(move(in), move(out), tableS1,ConstrType::Sbox);
	                else if(MI.sboxCheckType == CheckType::Ineq)
	                    MD.allVarsIneq.emplace_back(move(in), move(out), ineqS1,ConstrType::Sbox);
				}
				else if(i == 1 || i == 4){//S2
					if(MI.sboxCheckType == CheckType::Table)
	                    MD.allVarsTable.emplace_back(move(in), move(out), tableS2,ConstrType::Sbox);
	                else if(MI.sboxCheckType == CheckType::Ineq)
	                    MD.allVarsIneq.emplace_back(move(in), move(out), ineqS2,ConstrType::Sbox);
				}
				else if(i == 2 || i == 5){//S3
					if(MI.sboxCheckType == CheckType::Table)
	                    MD.allVarsTable.emplace_back(move(in), move(out), tableS3,ConstrType::Sbox);
	                else if(MI.sboxCheckType == CheckType::Ineq)
	                    MD.allVarsIneq.emplace_back(move(in), move(out), ineqS3,ConstrType::Sbox);
				}
				else{//S4
					if(MI.sboxCheckType == CheckType::Table)
	                    MD.allVarsTable.emplace_back(move(in), move(out), tableS4,ConstrType::Sbox);
	                else if(MI.sboxCheckType == CheckType::Ineq)
	                    MD.allVarsIneq.emplace_back(move(in), move(out), ineqS4,ConstrType::Sbox);
				}
    		}
        }
    }

    //MC variables

    if(MI.MCCheckType != CheckType::None){
        for(auto const & r : MI.roundOrder){
        	
    		//Lin layer in F function
            vector<GRBVar> in(64);
            vector<GRBVar> out(64);
            for(uint i = 0; i < 64; i++){
                in[i] = m.getVarByName("y"+to_string(r)+"_"+to_string(i));
                out[i] = m.getVarByName("z"+to_string(r)+"_"+to_string(i));
            }
            if(MI.MCCheckType == CheckType::Table){
                MD.allVarsTable.emplace_back(move(in), move(out), tableP,ConstrType::Lin);
            }
            else if(MI.MCCheckType == CheckType::Ineq){
                MD.allVarsIneq.emplace_back(move(in), move(out), ineqP,ConstrType::Lin);
            }
            else if(MI.MCCheckType == CheckType::Matrix){
                MD.allVarsMatrix.emplace_back(move(in), move(out), matrixP);
            }

	        if((r+MI.startingRound)%6 == 0){
	        	//FL layer
	        	vector<GRBVar> inFL(64);
	            vector<GRBVar> outFL(64);
	            for(uint i = 0; i < 64; i++){
	                inFL[i] = m.getVarByName("inFL"+to_string(r)+"_"+to_string(i));
	                outFL[i] = m.getVarByName("x"+to_string(r+1)+"_"+to_string(i));
	            }
	            if(MI.MCCheckType == CheckType::Table)
	                MD.allVarsTable.emplace_back(move(inFL), move(outFL), tableFL,ConstrType::Lin);
	            else if(MI.MCCheckType == CheckType::Ineq)
	                MD.allVarsIneq.emplace_back(move(inFL), move(outFL), ineqFL,ConstrType::Lin);
	            else if(MI.MCCheckType == CheckType::Matrix)
	                MD.allVarsMatrix.emplace_back(move(inFL), move(outFL), matrixFL);

	            //invFL
	            vector<GRBVar> ininvFL(64);
	            vector<GRBVar> outinvFL(64);
	            for(uint i = 0; i < 64; i++){
	                ininvFL[i] = m.getVarByName("inFL"+to_string(r)+"_"+to_string(i+64));
	                outinvFL[i] = m.getVarByName("x"+to_string(r+1)+"_"+to_string(i+64));
	            }
	            if(MI.MCCheckType == CheckType::Table)
	                MD.allVarsTable.emplace_back(move(ininvFL), move(outinvFL), tableinvFL,ConstrType::Lin);
	            else if(MI.MCCheckType == CheckType::Ineq)
	                MD.allVarsIneq.emplace_back(move(ininvFL), move(outinvFL), ineqinvFL,ConstrType::Lin);
	            else if(MI.MCCheckType == CheckType::Matrix)
	                MD.allVarsMatrix.emplace_back(move(ininvFL), move(outinvFL), matrixinvFL);
	        }
        }
    }

    return MD;
}

uint64_t globalCtrCheckCamellia = 0;
uint64_t globalCtrSolCamellia = 0;
uint64_t globalCtrConstrSboxCamellia = 0;
uint64_t globalCtrConstrLinCamellia = 0;
uint64_t globalCtrConstrSSBCamellia = 0;
bool existTrailCamellia(GRBModel & m,
                   modelData const & MD,
                   vector<uint8_t> const & input,
                   vector<uint8_t> const & output,
                   bool const useCallback,
                   uint const maxNumberConstr){
/*
    Return true if there is a trail from #input to #output using #model for Camellia
    - #MD contains the variables to be checked
    - if #useCallback is true then the solving using CB + LC, otherwise repeated solving
    - #maxNumberConstr is the max number of contraints added for each solution found
*/

    string outputvarPrefix = "x"+to_string(MD.rMax)+"_";
    for(uint i = 0; i < 128; i++){
        m.addConstr(m.getVarByName("x0_"+to_string(i)) == input[i]);
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
        globalCtrCheckCamellia++;
        globalCtrSolCamellia += cb.ctrSol;
        globalCtrConstrSboxCamellia += cb.ctrConstrSbox;
		globalCtrConstrLinCamellia += cb.ctrConstrLin;
		globalCtrConstrSSBCamellia += cb.ctrConstrSSB;

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

bool existTrailCamellia(modelInfo & MI,
                   vector<uint8_t> const & input,
                   vector<uint8_t> const & output){
/*
    Return true if there is a trail from #input to #output for Camellia
    #MI should contain the necessary information to generate the model and parametrize the solving
*/

    GRBModel m = getModelCamellia(MI);
    auto MD = getVariablesCamellia(m, MI);

    return existTrailCamellia(m,MD,input,output,MI.useCallback,MI.maxNumberConstr);
}

void genSboxFileCamellia(){
	//Generate the tables and ineq files for Camellia

	vector<uint32_t> S1(
		{0x70,0x82,0x2c,0xec,0xb3,0x27,0xc0,0xe5,0xe4,0x85,0x57,0x35,0xea,0x0c,0xae,0x41,
		 0x23,0xef,0x6b,0x93,0x45,0x19,0xa5,0x21,0xed,0x0e,0x4f,0x4e,0x1d,0x65,0x92,0xbd,
		 0x86,0xb8,0xaf,0x8f,0x7c,0xeb,0x1f,0xce,0x3e,0x30,0xdc,0x5f,0x5e,0xc5,0x0b,0x1a,
		 0xa6,0xe1,0x39,0xca,0xd5,0x47,0x5d,0x3d,0xd9,0x01,0x5a,0xd6,0x51,0x56,0x6c,0x4d,
		 0x8b,0x0d,0x9a,0x66,0xfb,0xcc,0xb0,0x2d,0x74,0x12,0x2b,0x20,0xf0,0xb1,0x84,0x99,
		 0xdf,0x4c,0xcb,0xc2,0x34,0x7e,0x76,0x05,0x6d,0xb7,0xa9,0x31,0xd1,0x17,0x04,0xd7,
		 0x14,0x58,0x3a,0x61,0xde,0x1b,0x11,0x1c,0x32,0x0f,0x9c,0x16,0x53,0x18,0xf2,0x22,
		 0xfe,0x44,0xcf,0xb2,0xc3,0xb5,0x7a,0x91,0x24,0x08,0xe8,0xa8,0x60,0xfc,0x69,0x50,
		 0xaa,0xd0,0xa0,0x7d,0xa1,0x89,0x62,0x97,0x54,0x5b,0x1e,0x95,0xe0,0xff,0x64,0xd2,
		 0x10,0xc4,0x00,0x48,0xa3,0xf7,0x75,0xdb,0x8a,0x03,0xe6,0xda,0x09,0x3f,0xdd,0x94,
		 0x87,0x5c,0x83,0x02,0xcd,0x4a,0x90,0x33,0x73,0x67,0xf6,0xf3,0x9d,0x7f,0xbf,0xe2,
		 0x52,0x9b,0xd8,0x26,0xc8,0x37,0xc6,0x3b,0x81,0x96,0x6f,0x4b,0x13,0xbe,0x63,0x2e,
		 0xe9,0x79,0xa7,0x8c,0x9f,0x6e,0xbc,0x8e,0x29,0xf5,0xf9,0xb6,0x2f,0xfd,0xb4,0x59,
		 0x78,0x98,0x06,0x6a,0xe7,0x46,0x71,0xba,0xd4,0x25,0xab,0x42,0x88,0xa2,0x8d,0xfa,
		 0x72,0x07,0xb9,0x55,0xf8,0xee,0xac,0x0a,0x36,0x49,0x2a,0x68,0x3c,0x38,0xf1,0xa4,
		 0x40,0x28,0xd3,0x7b,0xbb,0xc9,0x43,0xc1,0x15,0xe3,0xad,0xf4,0x77,0xc7,0x80,0x9e});

	cout << "S1" << endl;
	genAndSaveTableIneq(vector<vector<uint32_t>>({S1}),8,8,"./checkData/CamelliaS1Sbox",true);

	vector<uint32_t> S2(
		{0xe0,0x05,0x58,0xd9,0x67,0x4e,0x81,0xcb,0xc9,0x0b,0xae,0x6a,0xd5,0x18,0x5d,0x82,
		 0x46,0xdf,0xd6,0x27,0x8a,0x32,0x4b,0x42,0xdb,0x1c,0x9e,0x9c,0x3a,0xca,0x25,0x7b,
		 0x0d,0x71,0x5f,0x1f,0xf8,0xd7,0x3e,0x9d,0x7c,0x60,0xb9,0xbe,0xbc,0x8b,0x16,0x34,
		 0x4d,0xc3,0x72,0x95,0xab,0x8e,0xba,0x7a,0xb3,0x02,0xb4,0xad,0xa2,0xac,0xd8,0x9a,
		 0x17,0x1a,0x35,0xcc,0xf7,0x99,0x61,0x5a,0xe8,0x24,0x56,0x40,0xe1,0x63,0x09,0x33,
		 0xbf,0x98,0x97,0x85,0x68,0xfc,0xec,0x0a,0xda,0x6f,0x53,0x62,0xa3,0x2e,0x08,0xaf,
		 0x28,0xb0,0x74,0xc2,0xbd,0x36,0x22,0x38,0x64,0x1e,0x39,0x2c,0xa6,0x30,0xe5,0x44,
		 0xfd,0x88,0x9f,0x65,0x87,0x6b,0xf4,0x23,0x48,0x10,0xd1,0x51,0xc0,0xf9,0xd2,0xa0,
		 0x55,0xa1,0x41,0xfa,0x43,0x13,0xc4,0x2f,0xa8,0xb6,0x3c,0x2b,0xc1,0xff,0xc8,0xa5,
		 0x20,0x89,0x00,0x90,0x47,0xef,0xea,0xb7,0x15,0x06,0xcd,0xb5,0x12,0x7e,0xbb,0x29,
		 0x0f,0xb8,0x07,0x04,0x9b,0x94,0x21,0x66,0xe6,0xce,0xed,0xe7,0x3b,0xfe,0x7f,0xc5,
		 0xa4,0x37,0xb1,0x4c,0x91,0x6e,0x8d,0x76,0x03,0x2d,0xde,0x96,0x26,0x7d,0xc6,0x5c,
		 0xd3,0xf2,0x4f,0x19,0x3f,0xdc,0x79,0x1d,0x52,0xeb,0xf3,0x6d,0x5e,0xfb,0x69,0xb2,
		 0xf0,0x31,0x0c,0xd4,0xcf,0x8c,0xe2,0x75,0xa9,0x4a,0x57,0x84,0x11,0x45,0x1b,0xf5,
		 0xe4,0x0e,0x73,0xaa,0xf1,0xdd,0x59,0x14,0x6c,0x92,0x54,0xd0,0x78,0x70,0xe3,0x49,
		 0x80,0x50,0xa7,0xf6,0x77,0x93,0x86,0x83,0x2a,0xc7,0x5b,0xe9,0xee,0x8f,0x01,0x3d});

	cout << "S2" << endl;
	genAndSaveTableIneq(vector<vector<uint32_t>>({S2}),8,8,"./checkData/CamelliaS2Sbox",true);

	vector<uint32_t> S3(
		{0x38,0x41,0x16,0x76,0xd9,0x93,0x60,0xf2,0x72,0xc2,0xab,0x9a,0x75,0x06,0x57,0xa0,
		 0x91,0xf7,0xb5,0xc9,0xa2,0x8c,0xd2,0x90,0xf6,0x07,0xa7,0x27,0x8e,0xb2,0x49,0xde,
		 0x43,0x5c,0xd7,0xc7,0x3e,0xf5,0x8f,0x67,0x1f,0x18,0x6e,0xaf,0x2f,0xe2,0x85,0x0d,
		 0x53,0xf0,0x9c,0x65,0xea,0xa3,0xae,0x9e,0xec,0x80,0x2d,0x6b,0xa8,0x2b,0x36,0xa6,
		 0xc5,0x86,0x4d,0x33,0xfd,0x66,0x58,0x96,0x3a,0x09,0x95,0x10,0x78,0xd8,0x42,0xcc,
		 0xef,0x26,0xe5,0x61,0x1a,0x3f,0x3b,0x82,0xb6,0xdb,0xd4,0x98,0xe8,0x8b,0x02,0xeb,
		 0x0a,0x2c,0x1d,0xb0,0x6f,0x8d,0x88,0x0e,0x19,0x87,0x4e,0x0b,0xa9,0x0c,0x79,0x11,
		 0x7f,0x22,0xe7,0x59,0xe1,0xda,0x3d,0xc8,0x12,0x04,0x74,0x54,0x30,0x7e,0xb4,0x28,
		 0x55,0x68,0x50,0xbe,0xd0,0xc4,0x31,0xcb,0x2a,0xad,0x0f,0xca,0x70,0xff,0x32,0x69,
		 0x08,0x62,0x00,0x24,0xd1,0xfb,0xba,0xed,0x45,0x81,0x73,0x6d,0x84,0x9f,0xee,0x4a,
		 0xc3,0x2e,0xc1,0x01,0xe6,0x25,0x48,0x99,0xb9,0xb3,0x7b,0xf9,0xce,0xbf,0xdf,0x71,
		 0x29,0xcd,0x6c,0x13,0x64,0x9b,0x63,0x9d,0xc0,0x4b,0xb7,0xa5,0x89,0x5f,0xb1,0x17,
		 0xf4,0xbc,0xd3,0x46,0xcf,0x37,0x5e,0x47,0x94,0xfa,0xfc,0x5b,0x97,0xfe,0x5a,0xac,
		 0x3c,0x4c,0x03,0x35,0xf3,0x23,0xb8,0x5d,0x6a,0x92,0xd5,0x21,0x44,0x51,0xc6,0x7d,
		 0x39,0x83,0xdc,0xaa,0x7c,0x77,0x56,0x05,0x1b,0xa4,0x15,0x34,0x1e,0x1c,0xf8,0x52,
		 0x20,0x14,0xe9,0xbd,0xdd,0xe4,0xa1,0xe0,0x8a,0xf1,0xd6,0x7a,0xbb,0xe3,0x40,0x4f});

	cout << "S3" << endl;
	genAndSaveTableIneq(vector<vector<uint32_t>>({S3}),8,8,"./checkData/CamelliaS3Sbox",true);

	vector<uint32_t> S4(
		{0x70,0x2c,0xb3,0xc0,0xe4,0x57,0xea,0xae,0x23,0x6b,0x45,0xa5,0xed,0x4f,0x1d,0x92,
		 0x86,0xaf,0x7c,0x1f,0x3e,0xdc,0x5e,0x0b,0xa6,0x39,0xd5,0x5d,0xd9,0x5a,0x51,0x6c,
		 0x8b,0x9a,0xfb,0xb0,0x74,0x2b,0xf0,0x84,0xdf,0xcb,0x34,0x76,0x6d,0xa9,0xd1,0x04,
		 0x14,0x3a,0xde,0x11,0x32,0x9c,0x53,0xf2,0xfe,0xcf,0xc3,0x7a,0x24,0xe8,0x60,0x69,
		 0xaa,0xa0,0xa1,0x62,0x54,0x1e,0xe0,0x64,0x10,0x00,0xa3,0x75,0x8a,0xe6,0x09,0xdd,
		 0x87,0x83,0xcd,0x90,0x73,0xf6,0x9d,0xbf,0x52,0xd8,0xc8,0xc6,0x81,0x6f,0x13,0x63,
		 0xe9,0xa7,0x9f,0xbc,0x29,0xf9,0x2f,0xb4,0x78,0x06,0xe7,0x71,0xd4,0xab,0x88,0x8d,
		 0x72,0xb9,0xf8,0xac,0x36,0x2a,0x3c,0xf1,0x40,0xd3,0xbb,0x43,0x15,0xad,0x77,0x80,
		 0x82,0xec,0x27,0xe5,0x85,0x35,0x0c,0x41,0xef,0x93,0x19,0x21,0x0e,0x4e,0x65,0xbd,
		 0xb8,0x8f,0xeb,0xce,0x30,0x5f,0xc5,0x1a,0xe1,0xca,0x47,0x3d,0x01,0xd6,0x56,0x4d,
		 0x0d,0x66,0xcc,0x2d,0x12,0x20,0xb1,0x99,0x4c,0xc2,0x7e,0x05,0xb7,0x31,0x17,0xd7,
		 0x58,0x61,0x1b,0x1c,0x0f,0x16,0x18,0x22,0x44,0xb2,0xb5,0x91,0x08,0xa8,0xfc,0x50,
		 0xd0,0x7d,0x89,0x97,0x5b,0x95,0xff,0xd2,0xc4,0x48,0xf7,0xdb,0x03,0xda,0x3f,0x94,
		 0x5c,0x02,0x4a,0x33,0x67,0xf3,0x7f,0xe2,0x9b,0x26,0x37,0x3b,0x96,0x4b,0xbe,0x2e,
		 0x79,0x8c,0x6e,0x8e,0xf5,0xb6,0xfd,0x59,0x98,0x6a,0x46,0xba,0x25,0x42,0xa2,0xfa,
		 0x07,0x55,0xee,0x0a,0x49,0x68,0x38,0xa4,0x28,0x7b,0xc9,0xc1,0xe3,0xf4,0xc7,0x9e});

	cout << "S4" << endl;
	genAndSaveTableIneq(vector<vector<uint32_t>>({S4}),8,8,"./checkData/CamelliaS4Sbox",true);
}

void createMatrixFileCamellia(){
	//Create the matrix file ./checkData/matrixCamellia.bin

	//Slightly compressed in the code, the binary matrices are 64 x 64...
	//The values were obtained by hardcoding the binary matrix and priting it with M.printRawHex()
	//The hex values printed are actually the value of the chunks in the internal representation, so more compact code (still ugly and large)

	vector<uint64_t> val({0x0101010001010001, 
						  0x0202020002020002, 
						  0x0404040004040004, 
						  0x0808080008080008, 
						  0x1010100010100010, 
						  0x2020200020200020, 
						  0x4040400040400040, 
						  0x8080800080800080, 
						  0x0101000101000101, 
						  0x0202000202000202, 
						  0x0404000404000404, 
						  0x0808000808000808, 
						  0x1010001010001010, 
						  0x2020002020002020, 
						  0x4040004040004040, 
						  0x8080008080008080, 
						  0x0100010100010101, 
						  0x0200020200020202, 
						  0x0400040400040404, 
						  0x0800080800080808, 
						  0x1000101000101010, 
						  0x2000202000202020, 
						  0x4000404000404040, 
						  0x8000808000808080, 
						  0x0001010101010100, 
						  0x0002020202020200, 
						  0x0004040404040400, 
						  0x0008080808080800, 
						  0x0010101010101000, 
						  0x0020202020202000, 
						  0x0040404040404000, 
						  0x0080808080808000, 
						  0x0101010000000101, 
						  0x0202020000000202, 
						  0x0404040000000404, 
						  0x0808080000000808, 
						  0x1010100000001010, 
						  0x2020200000002020, 
						  0x4040400000004040, 
						  0x8080800000008080, 
						  0x0101000100010100, 
						  0x0202000200020200, 
						  0x0404000400040400, 
						  0x0808000800080800, 
						  0x1010001000101000, 
						  0x2020002000202000, 
						  0x4040004000404000, 
						  0x8080008000808000, 
						  0x0100010101010000, 
						  0x0200020202020000, 
						  0x0400040404040000, 
						  0x0800080808080000, 
						  0x1000101010100000, 
						  0x2000202020200000, 
						  0x4000404040400000, 
						  0x8000808080800000, 
						  0x0001010101000001, 
						  0x0002020202000002, 
						  0x0004040404000004, 
						  0x0008080808000008, 
						  0x0010101010000010, 
						  0x0020202020000020, 
						  0x0040404040000040, 
						  0x0080808080000080});

	Matrix M(64,64);
	for(uint i = 0; i < 64; i++)
		M.get(i)[0] = val[i];

	M.saveToFile("./checkData/matrixCamelliaP.bin");
	cout << "P : " << endl;
	M.printRawHex();

	vector<uint64_t> valFL({0x0000000180000001, 
							0x0000000200000003, 
							0x0000000400000006, 
							0x000000080000000c, 
							0x0000001000000018, 
							0x0000002000000030, 
							0x0000004000000060, 
							0x00000080000000c0, 
							0x0000010000000180, 
							0x0000020000000300, 
							0x0000040000000600, 
							0x0000080000000c00, 
							0x0000100000001800, 
							0x0000200000003000, 
							0x0000400000006000, 
							0x000080000000c000, 
							0x0001000000018000, 
							0x0002000000030000, 
							0x0004000000060000, 
							0x00080000000c0000, 
							0x0010000000180000, 
							0x0020000000300000, 
							0x0040000000600000, 
							0x0080000000c00000, 
							0x0100000001800000, 
							0x0200000003000000, 
							0x0400000006000000, 
							0x080000000c000000, 
							0x1000000018000000, 
							0x2000000030000000, 
							0x4000000060000000, 
							0x80000000c0000000, 
							0x0000000180000000, 
							0x0000000200000001, 
							0x0000000400000002, 
							0x0000000800000004, 
							0x0000001000000008, 
							0x0000002000000010, 
							0x0000004000000020, 
							0x0000008000000040, 
							0x0000010000000080, 
							0x0000020000000100, 
							0x0000040000000200, 
							0x0000080000000400, 
							0x0000100000000800, 
							0x0000200000001000, 
							0x0000400000002000, 
							0x0000800000004000, 
							0x0001000000008000, 
							0x0002000000010000, 
							0x0004000000020000, 
							0x0008000000040000, 
							0x0010000000080000, 
							0x0020000000100000, 
							0x0040000000200000, 
							0x0080000000400000, 
							0x0100000000800000, 
							0x0200000001000000, 
							0x0400000002000000, 
							0x0800000004000000, 
							0x1000000008000000, 
							0x2000000010000000, 
							0x4000000020000000, 
							0x8000000040000000});

	Matrix FL(64,64);
	for(uint i = 0; i < 64; i++)
		FL.get(i)[0] = valFL[i];

	FL.saveToFile("./checkData/matrixCamelliaFL.bin");
	cout << "FL : " << endl;
	FL.printRawHex();

	vector<uint64_t> valinvFL({0x0000000100000001, 
							   0x0000000200000002, 
							   0x0000000400000004, 
							   0x0000000800000008, 
							   0x0000001000000010, 
							   0x0000002000000020, 
							   0x0000004000000040, 
							   0x0000008000000080, 
							   0x0000010000000100, 
							   0x0000020000000200, 
							   0x0000040000000400, 
							   0x0000080000000800, 
							   0x0000100000001000, 
							   0x0000200000002000, 
							   0x0000400000004000, 
							   0x0000800000008000, 
							   0x0001000000010000, 
							   0x0002000000020000, 
							   0x0004000000040000, 
							   0x0008000000080000, 
							   0x0010000000100000, 
							   0x0020000000200000, 
							   0x0040000000400000, 
							   0x0080000000800000, 
							   0x0100000001000000, 
							   0x0200000002000000, 
							   0x0400000004000000, 
							   0x0800000008000000, 
							   0x1000000010000000, 
							   0x2000000020000000, 
							   0x4000000040000000, 
							   0x8000000080000000, 
							   0x8000000180000000, 
							   0x0000000300000001, 
							   0x0000000600000002, 
							   0x0000000c00000004, 
							   0x0000001800000008, 
							   0x0000003000000010, 
							   0x0000006000000020, 
							   0x000000c000000040, 
							   0x0000018000000080, 
							   0x0000030000000100, 
							   0x0000060000000200, 
							   0x00000c0000000400, 
							   0x0000180000000800, 
							   0x0000300000001000, 
							   0x0000600000002000, 
							   0x0000c00000004000, 
							   0x0001800000008000, 
							   0x0003000000010000, 
							   0x0006000000020000, 
							   0x000c000000040000, 
							   0x0018000000080000, 
							   0x0030000000100000, 
							   0x0060000000200000, 
							   0x00c0000000400000, 
							   0x0180000000800000, 
							   0x0300000001000000, 
							   0x0600000002000000, 
							   0x0c00000004000000, 
							   0x1800000008000000, 
							   0x3000000010000000, 
							   0x6000000020000000, 
							   0xc000000040000000});

	Matrix invFL(64,64);
	for(uint i = 0; i < 64; i++)
		invFL.get(i)[0] = valinvFL[i];

	invFL.saveToFile("./checkData/matrixCamelliainvFL.bin");
	cout << "invFL : " << endl;
	invFL.printRawHex();
}

void checkNormalDistinguisherCamellia(){
	/*
		Search for a normal (no linear combination) distinguisher over 7 rounds of Camellia
		Focus only on the already known distinguisher (last 64 bits balanced, the other are not)
		Don't use the checks with minors for MC
	*/

	vector<uint8_t> input({0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1});

	//These parameters have no effect for Camellia but are required for the modelInfo object
	string arxModel = "ARX";
	CheckType ARXCheckType = CheckType::None;
	bool smartLin = true;

	//Fixed parameters
	uint rMax = 7;
	bool firstSbox = true;
	bool lastLin = false;
	CheckType sboxCheckType = CheckType::None;
	CheckType MCCheckType = CheckType::None;
	CheckType SSBCheckType = CheckType::None;
	bool useCallback = false;
	uint maxNumberConstr = 10;
	vector<uint> roundOrder = roundOrderAscend(rMax);

	//Variable parameters
	//The simple model without checks already find the known 64 balanced bits
	vector<string> listSboxModel({"QM", "Simple"});
	vector<string> listLinModel({"CX", "Simple"});
	vector<uint> listStartinground({0,1,2,3,4,5}); //only the value of r%6 matter

	//The searchs
	for(auto const & startingRound : listStartinground){
		for(auto const & sboxModel : listSboxModel){
				for(auto const & linModel : listLinModel){

					globalCtrCheckCamellia = 0;
					globalCtrSolCamellia = 0;
					globalCtrConstrSboxCamellia = 0;
					globalCtrConstrLinCamellia = 0;
					globalCtrConstrSSBCamellia = 0;

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
						balancedBits[i] = !existTrailCamellia(MI,input,output);
						auto endBit = chrono::high_resolution_clock::now();
						if(balancedBits[i]){
							ctrBalanced++;
							durBalanced += endBit - startBit;
						}
					}
					auto end = chrono::high_resolution_clock::now();
					cout << "Starting round : " << startingRound << " - ";
					cout << " sboxModel : " << sboxModel << " - linModel : " << linModel << endl;
					cout << "sboxCheckType : " << sboxCheckType << " - MCCheckType : " << MCCheckType << endl;
					cout << ctrBalanced << " balanced bits" << endl;
					cout << chrono::duration<double>(end - start).count() << " seconds" << endl;
					cout << "Time for only balanced bits : " << chrono::duration<double>(durBalanced).count() << " second" << endl;
					cout << globalCtrCheckCamellia << " input/output pairs checked" << endl;
					cout << globalCtrSolCamellia << " solutions examined" << endl;
					cout << globalCtrConstrSboxCamellia << " constraints added for the sbox" << endl;
					cout << globalCtrConstrLinCamellia  << " constraints added for the lin layer" << endl;
					cout << globalCtrConstrSSBCamellia << " constraints added for the SSB" << endl;
			}
		}
	}
}


void timingCamellia8r(){
	vector<uint8_t> input({0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1});

	//These parameters have no effect for Camellia but are required for the modelInfo object
	string arxModel = "ARX";
	CheckType ARXCheckType = CheckType::None;
	bool smartLin = true;

	//Fixed parameters
	uint rMax = 8;
	bool firstSbox = true;
	bool lastLin = false;
	// CheckType sboxCheckType = CheckType::None;
	CheckType MCCheckType = CheckType::Matrix;
	CheckType SSBCheckType = CheckType::None;
	bool useCallback = true;
	// uint maxNumberConstr = 100;
	// vector<uint> roundOrder = roundOrderAscend(rMax);

	//Variable parameters
	//Simple model without checks find 64 balanced bits
	vector<string> listSboxModel({"QM", "Simple"});
	vector<string> listLinModel({"CX", "Simple"});
	vector<uint> listStartinground({0}); //only the value of r%6 matter
	vector<uint> listMaxNbConstr({100,10,1});
	vector<vector<uint>> listRoundOrder({roundOrderOutsideIn(rMax), roundOrderAscend(rMax), roundOrderDescend(rMax)});

	//The searchs
	for(auto const & startingRound : listStartinground){
		for(auto const & maxNumberConstr : listMaxNbConstr){
			for(auto const & roundOrder : listRoundOrder){
				for(auto const & sboxModel : listSboxModel){
					CheckType sboxCheckType = CheckType::None;
					if(sboxModel == "Simple")
						sboxCheckType = CheckType::Ineq;
					for(auto const & linModel : listLinModel){

						cout << "Starting round : " << startingRound << " - ";
						cout << "Round Order : ";
					    for(auto const & tmpr : roundOrder) cout << tmpr << " ";
					    cout << endl;
					    cout << " sboxModel : " << sboxModel << " - linModel : " << linModel << endl;
						cout << "sboxCheckType : " << sboxCheckType << " - MCCheckType : " << MCCheckType << endl;
					    cout << "maxNumberConstr : " << maxNumberConstr << endl;

					    globalCtrCheckCamellia = 0;
						globalCtrSolCamellia = 0;
			      	 	globalCtrConstrSboxCamellia = 0;
						globalCtrConstrLinCamellia = 0;
						globalCtrConstrSSBCamellia = 0;

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
							balancedBits[i] = !existTrailCamellia(MI,input,output);
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

						cout << globalCtrCheckCamellia << " input/output pairs checked" << endl;
				 		cout << globalCtrSolCamellia << " solutions examined" << endl;
				 		cout << globalCtrConstrSboxCamellia << " constraints added for the sbox" << endl;
				 		cout << globalCtrConstrLinCamellia  << " constraints added for the lin layer" << endl;
				 		cout << globalCtrConstrSSBCamellia << " constraints added for the SSB" << endl;
					}
				}
			}
		}
		cout << endl;
	}
}

void searchCamellia7r(){

	vector<uint8_t> input({0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1});

	//These parameters have no effect for Camellia but are required for the modelInfo object
	string arxModel = "ARX";
	CheckType ARXCheckType = CheckType::None;
	bool smartLin = true;

	//Fixed parameters
	uint rMax = 7;
	bool firstSbox = true;
	bool lastLin = false;
	CheckType sboxCheckType = CheckType::None;
	CheckType MCCheckType = CheckType::None;
	CheckType SSBCheckType = CheckType::None;
	bool useCallback = true;
	uint maxNumberConstr = 10;
	vector<uint> roundOrder = roundOrderAscend(rMax);

	//Variable parameters
	vector<uint> listStartinground({0,1,2,3,4,5}); //only the value of r%6 matter

	string sboxModel = "Simple";
	string linModel;

	//The searchs
	for(auto const & startingRound : listStartinground){
		if(startingRound == 0 || startingRound == 5)
			linModel = "CX";
		else
			linModel = "Simple";

		globalCtrCheckCamellia = 0;
		globalCtrSolCamellia = 0;
		globalCtrConstrSboxCamellia = 0;
		globalCtrConstrLinCamellia = 0;
		globalCtrConstrSSBCamellia = 0;

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
			balancedBits[i] = !existTrailCamellia(MI,input,output);
			auto endBit = chrono::high_resolution_clock::now();
			if(balancedBits[i]){
				ctrBalanced++;
				durBalanced += endBit - startBit;
			}
		}
		auto end = chrono::high_resolution_clock::now();
		cout << "Starting round : " << startingRound << " - ";
		cout << " sboxModel : " << sboxModel << " - linModel : " << linModel << endl;
		cout << "sboxCheckType : " << sboxCheckType << " - MCCheckType : " << MCCheckType << endl;
		cout << ctrBalanced << " balanced bits" << endl;
		cout << chrono::duration<double>(end - start).count() << " seconds" << endl;
		cout << "Time for only balanced bits : " << chrono::duration<double>(durBalanced).count() << " second" << endl;
		cout << globalCtrCheckCamellia << " input/output pairs checked" << endl;
		cout << globalCtrSolCamellia << " solutions examined" << endl;
		cout << globalCtrConstrSboxCamellia << " constraints added for the sbox" << endl;
		cout << globalCtrConstrLinCamellia  << " constraints added for the lin layer" << endl;
		cout << globalCtrConstrSSBCamellia << " constraints added for the SSB" << endl;
	}
}

void searchCamellia8r(){

	uint rMax = 8;
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
	uint maxNumberConstr = 10;
	vector<uint> roundOrder = roundOrderAscend(rMax);

	for(uint startingRound = 0; startingRound < 6; startingRound++){
		globalCtrCheckCamellia = 0;
		globalCtrSolCamellia = 0;
		globalCtrConstrSboxCamellia = 0;
		globalCtrConstrLinCamellia = 0;
		globalCtrConstrSSBCamellia = 0;
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
				balancedBits[i] = !existTrailCamellia(MI,input,output);
				auto endBit = chrono::high_resolution_clock::now();
				if(balancedBits[i]){
					ctrBalanced++;
					durBalanced += endBit - startBit;
				}
			}

			auto end = chrono::high_resolution_clock::now();
			cout << "Starting round : " << startingRound << endl;
			cout << "Constant input bit : " << indexInput << endl;
			cout << ctrBalanced << " balanced bits" << endl;
			cout << chrono::duration<double>(end - start).count() << " seconds" << endl;
			cout << "Time for only balanced bits : " << chrono::duration<double>(durBalanced).count() << " second" << endl;
			if(QUICK_LAZY_COUNT_CAMELLIA)
				break;
		}
		cout << globalCtrCheckCamellia << " input/output pairs checked" << endl;
		cout << globalCtrSolCamellia << " solutions examined" << endl;
		cout << globalCtrConstrSboxCamellia << " constraints added for the sbox" << endl;
		cout << globalCtrConstrLinCamellia  << " constraints added for the lin layer" << endl;
		cout << globalCtrConstrSSBCamellia << " constraints added for the SSB" << endl;
		if(QUICK_LAZY_COUNT_CAMELLIA)
			break;
	}

}