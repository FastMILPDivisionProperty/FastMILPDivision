#include "hight.hpp"
#include "DivTable.hpp"

using namespace std;

#define DISPLAY_GEN_OUTPUT false
#define DISPLAY_SOLVER_OUTPUT false
#define DISPLAY_NUMBER_ADDED_CONSTRAINT false

static vector<Polynome> F0(vector<Polynome> const & x) {
  unsigned const n = x.size();
  vector<Polynome> y (n);
  for (unsigned i = 0; i < n; ++i) {
    y[i] = x[(i+n-1)%n] + x[(i+n-2)%n] + x[(i+n-7)%n];
  }
  return y;
}

static vector<Polynome> F1(vector<Polynome> const & x) {
  unsigned const n = x.size();
  vector<Polynome> y (n);
  for (unsigned i = 0; i < n; ++i) {
    y[i] = x[(i+n-3)%n] + x[(i+n-4)%n] + x[(i+n-6)%n];
  }
  return y;
}

static vector<Polynome> modularAddition(vector<Polynome> const & x, vector<Polynome> const & y) {
  vector<Polynome> z (x.size());
  vector<Polynome> c (x.size());
  z[0] = x[0] + y[0];
  c[0] = x[0]*y[0];
  for (unsigned i = 1; i < x.size(); ++i) {
    z[i] = x[i] + y[i] + c[i-1];
    c[i] = (x[i] + c[i-1])*(y[i] + c[i-1]) + c[i-1];
  }
  return z;
}

vector<vector<pair<vector<uint8_t>, vector<uint8_t>>>> allInOutPairsHIGHT(unsigned i, unsigned j) { // 0 <= i < 4  - 0 <= j < 4
  vector<Polynome> x0 (8);
  for (unsigned i = 0; i < 8; ++i) x0[i] = Monome(1u << i);
  vector<Polynome> x1 (8);
  for (unsigned i = 0; i < 8; ++i) x1[i] = Monome(1u << (i+8));
  vector<Polynome> k (8);
  for (unsigned i = 0; i < 8; ++i) k[i] = Monome(1u << (i+16));

  auto f1 = F1(x0);
  for (unsigned i = 0; i < 8; ++i) f1[i] += k[i];
  auto y1 = modularAddition(f1, x1);
  vector<Polynome> anf1 (16);
  for (unsigned i = 0; i < 8; ++i) anf1[i] = x0[i];
  for (unsigned i = 0; i < 8; ++i) anf1[i+8] = y1[i];

  auto f0 = F0(x0);
  auto y0 = modularAddition(k, f0);
  for (unsigned i = 0; i < 8; ++i) y0[i] += x1[i];
  vector<Polynome> anf0 (16);
  for (unsigned i = 0; i < 8; ++i) anf0[i] = x0[i];
  for (unsigned i = 0; i < 8; ++i) anf0[i+8] = y0[i];

  vector<uint8_t> perm = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1};
  vector<uint8_t> perm_inv (16);
  for (unsigned i = 0; i < 16; ++i) perm_inv[perm[i]] = i;

  auto v_in = (i%2 == 0) ? loadTransitionsIN(anf1) : loadTransitionsIN(anf0);
  auto v_out = (j%2 == 0) ? loadTransitionsOUT(anf1) : loadTransitionsOUT(anf0);

  // cout << "v_in: " << v_in.size() << endl;
  // cout << "v_out: " << v_out.size() << endl;

  vector<vector<pair<vector<uint8_t>, vector<uint8_t>>>> res;
  for (auto pin : v_in) {
    for (auto pout : v_out) {
      vector<pair<vector<uint8_t>, vector<uint8_t>>> tmp;
      for (auto xin : pin.first) {
        for (auto xout : pout.first) {
          vector<uint8_t> in_8t (64, 1);
          vector<uint8_t> out_8t (64, 0);
          for (unsigned b = 0; b < 16; ++b) in_8t[16*i + b] = (xin >> b) & 1;
          for (unsigned b = 0; b < 16; ++b) out_8t[16*j + b] = (xout >> b) & 1;
          tmp.emplace_back(move(in_8t), move(out_8t));
        }
      }
      res.emplace_back(move(tmp));
    }
  }

  return res;

}

GRBModel getModelHIGHT(modelInfo & MI){
/*
	Get the model for HIGHT over #rMax round

	The #arxModel parameter defines how the mod addition is modelized :
	- "ARX" use the technique from Sun,Wang,Liu,Wang

	The #linModel parameter defines how the linear layer is modelized :
	- "QM" use the Quin-McCluskey algorithm as in Abdelkhalek,Sasaki,Todo,Tolba,Youssef
	- "ZR" modelize the linear layer using the technique from Zhang and Rijmen
	- "CX" modelize the linear layer using the classical copy+xor technique
	- "Simple" modelize the linear layer with the simplified constraint w(x) = w(y)
*/
    string modelName = "HIGHT_"+to_string(MI.rMax)+"r_"+MI.arxModel+"_"+MI.linModel;
    modelName += ".mps";

    if(!fileExist("./models/"+modelName)){
        string cmd = "sage ./modelGeneration/genHIGHT.sage";
        cmd += " -r " + to_string(MI.rMax);
        cmd += " -a " + MI.arxModel;
        cmd += " -m " + MI.linModel;
        auto sysstdout = exec(cmd.c_str());
        if(DISPLAY_GEN_OUTPUT) cout << sysstdout << endl;
    }

    if(!DISPLAY_SOLVER_OUTPUT)
        MI.gurobiEnv.set(GRB_IntParam_OutputFlag,0);

    GRBModel m(MI.gurobiEnv, "./models/"+modelName);
    return m;
}

modelData getVariablesHIGHT(GRBModel & m,
							modelInfo const & MI){
/*
	Read the model and return the corresponding blocks of variables
	- #MI defines various parameters for how the trails are checked (see aux_function.hpp)
*/

    modelData MD;
    MD.rMax = MI.rMax;

    MD.allTable = vector<vector<vector<uint32_t>>>(5);
    MD.allIneq = vector<vector<vector<int>>>(5);
    MD.allMatrix = vector<Matrix>(2);

    if(MI.ARXCheckType != CheckType::None){
    	cerr << "Error : ARX Check for HIGHT not implemented" << endl;
    	exit(1);
    }
    // if(ARXCheckType == CheckType::Table)
    //     MD.allTable[0] = //Table for ARX
    // else if(ARXCheckType == CheckType::Ineq)
    //     MD.allIneq[0]   = //Ineq for ARX

    if(MI.MCCheckType == CheckType::Matrix){
    	MD.allMatrix[0] = Matrix("./checkData/matrixHIGHTF0.bin");
    	MD.allMatrix[1] = Matrix("./checkData/matrixHIGHTF1.bin");
    }
    else if(MI.MCCheckType == CheckType::Table){
        MD.allTable[1] = readDivTableFromFile("./checkData/HIGHTF0_table.bin");
        MD.allTable[2] = readDivTableFromFile("./checkData/HIGHTF1_table.bin");
    }
    else if(MI.MCCheckType == CheckType::Ineq){
        MD.allIneq[1]   = readIneqFromFile("./checkData/HIGHTF0_ineq.bin");
   		MD.allIneq[2]   = readIneqFromFile("./checkData/HIGHTF1_ineq.bin");
   	}

    if(MI.SSBCheckType == CheckType::Table){
        MD.allTable[3] = readDivTableFromFile("./checkData/HIGHTSSBF0_table.bin");
        MD.allTable[4] = readDivTableFromFile("./checkData/HIGHTSSBF1_table.bin");
    }
    else if(MI.SSBCheckType == CheckType::Ineq){
        MD.allIneq[3]   = readIneqFromFile("./checkData/HIGHTSSBF0_ineq.bin");
        MD.allIneq[4]   = readIneqFromFile("./checkData/HIGHTSSBF1_ineq.bin");
    }


    auto & tableARX = MD.allTable[0];
    auto & ineqARX	= MD.allIneq[0];

    auto & tableF0  = MD.allTable[1];
    auto & ineqF0   = MD.allIneq[1];
    auto & tableF1  = MD.allTable[2];
    auto & ineqF1   = MD.allIneq[2];

    auto & tableSSBF0  = MD.allTable[3];
    auto & tableSSBF1  = MD.allTable[4];
    auto & ineqSSBF0   = MD.allIneq[3];
    auto & ineqSSBF1   = MD.allIneq[4];

    auto & matrixF0 = MD.allMatrix[0];
    auto & matrixF1 = MD.allMatrix[1];

    //ARX variables
    if(MI.ARXCheckType != CheckType::None){
        for(auto const & r : MI.roundOrder){
        	vector<GRBVar> in1(16);
        	vector<GRBVar> out1(8);
        	vector<GRBVar> in2(16);
        	vector<GRBVar> out2(8);
        	for(uint i = 0; i < 8; i++){
        		in1[i] = m.getVarByName("x"+to_string(r)+"_1_"+to_string(i));
        		in1[8+i] = m.getVarByName("y"+to_string(r)+"_0_"+to_string(i));
        		out1[i] = m.getVarByName("x"+to_string(r+1)+"_2_"+to_string(i));

        		in2[i] = m.getVarByName("x"+to_string(r)+"_5_"+to_string(i));
        		in2[8+i] = m.getVarByName("y"+to_string(r)+"_4_"+to_string(i));
        		out2[i] = m.getVarByName("x"+to_string(r+1)+"_6_"+to_string(i));
        	}

        	if(MI.ARXCheckType == CheckType::Table){
        		MD.allVarsTable.emplace_back(move(in1), move(out1), tableARX,ConstrType::ARX);
                MD.allVarsTable.emplace_back(move(in2), move(out2), tableARX,ConstrType::ARX);
            }
            else if(MI.ARXCheckType == CheckType::Ineq){
                MD.allVarsIneq.emplace_back(move(in2), move(out2), ineqARX,ConstrType::ARX);
            	MD.allVarsIneq.emplace_back(move(in1), move(out1), ineqARX,ConstrType::ARX);
            }
        }
    }

    //MC variables
    if(MI.MCCheckType != CheckType::None){
        for(auto const & r : MI.roundOrder){
        	for(uint i = 0; i < 8; i+=2){
        		vector<GRBVar> in(8);
        		vector<GRBVar> out(8);
        		for(uint j = 0; j < 8; j++){
	                in[j] = m.getVarByName("cx"+to_string(r)+"_"+to_string(i)+"_"+to_string(j));
	                out[j] = m.getVarByName("y"+to_string(r)+"_"+to_string(i)+"_"+to_string(j));
	            }
	            if(i == 0 || i == 4){
	            	if(MI.MCCheckType == CheckType::Table)
		                MD.allVarsTable.emplace_back(move(in), move(out), tableF1,ConstrType::Lin);
		            else if(MI.MCCheckType == CheckType::Ineq)
		                MD.allVarsIneq.emplace_back(move(in), move(out), ineqF1,ConstrType::Lin);
		            else if(MI.MCCheckType == CheckType::Matrix)
                		MD.allVarsMatrix.emplace_back(move(in), move(out), matrixF1);
	            }
	            else{
	            	if(MI.MCCheckType == CheckType::Table)
		                MD.allVarsTable.emplace_back(move(in), move(out), tableF0,ConstrType::Lin);
		            else if(MI.MCCheckType == CheckType::Ineq)
		                MD.allVarsIneq.emplace_back(move(in), move(out), ineqF0,ConstrType::Lin);
		            else if(MI.MCCheckType == CheckType::Matrix)
                		MD.allVarsMatrix.emplace_back(move(in), move(out), matrixF0);
	            }
        	}
        }
    }

    //SSB variables
    if(MI.SSBCheckType != CheckType::None){
    	for(auto const & r : MI.roundOrder){
    		for(uint i = 0; i < 8; i+=2){
    			vector<GRBVar> in(16);
    			vector<GRBVar> out(16);
    			for(uint j = 0; j < 8; j++){
    				in[j] = m.getVarByName("x"+to_string(r)+"_"+to_string(i)+"_"+to_string(j));
    				in[8+j] = m.getVarByName("x"+to_string(r)+"_"+to_string(i+1)+"_"+to_string(j));
    				out[j] = m.getVarByName("x"+to_string(r+1)+"_"+to_string(i+1)+"_"+to_string(j));
    				out[8+j] = m.getVarByName("x"+to_string(r+1)+"_"+to_string((i+2)%8)+"_"+to_string(j));
    			}
    			if(i == 0 || i == 4){
	            	if(MI.SSBCheckType == CheckType::Table)
		                MD.allVarsTable.emplace_back(move(in), move(out), tableSSBF1,ConstrType::SSB);
		            else if(MI.SSBCheckType == CheckType::Ineq)
		                MD.allVarsIneq.emplace_back(move(in), move(out), ineqSSBF1,ConstrType::SSB);
	            }
	            else{
	            	if(MI.SSBCheckType == CheckType::Table)
		                MD.allVarsTable.emplace_back(move(in), move(out), tableSSBF0,ConstrType::SSB);
		            else if(MI.SSBCheckType == CheckType::Ineq)
		                MD.allVarsIneq.emplace_back(move(in), move(out), ineqSSBF0,ConstrType::SSB);
	            }
    		}
    	}
    }

    return MD;
}

uint64_t globalCtrCheckHIGHT = 0;
uint64_t globalCtrSolHIGHT = 0;
uint64_t globalCtrConstrSboxHIGHT = 0;
uint64_t globalCtrConstrLinHIGHT = 0;
uint64_t globalCtrConstrSSBHIGHT = 0;
bool existTrailHIGHT(GRBModel & m,
                   modelData const & MD,
                   vector<uint8_t> const & input,
                   vector<uint8_t> const & output,
                   bool const useCallback,
                   uint const maxNumberConstr){
/*
    Return true if there is a trail from #input to #output using #model for HIGHT
    - #MD contains the variables to be checked
    - if #useCallback is true then the solving using CB + LC, otherwise repeated solving
    - #maxNumberConstr is the max number of contraints added for each solution found
*/

	uint ctr = 0;
    for(uint i = 0; i < 8; i++){
    	for(uint j = 0; j < 8; j++){
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
        globalCtrCheckHIGHT++;
        globalCtrSolHIGHT += cb.ctrSol;
        globalCtrConstrSboxHIGHT += cb.ctrConstrSbox;
		globalCtrConstrLinHIGHT += cb.ctrConstrLin;
		globalCtrConstrSSBHIGHT += cb.ctrConstrSSB;

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

bool existTrailHIGHT(modelInfo & MI,
                   vector<uint8_t> const & input,
                   vector<uint8_t> const & output){
/*
    Return true if there is a trail from #input to #output for HIGHT
    #MI should contain the necessary information to generate the model and parametrize the solving
*/

    GRBModel m = getModelHIGHT(MI);
    auto MD = getVariablesHIGHT(m, MI);

    return existTrailHIGHT(m,MD,input,output,MI.useCallback,MI.maxNumberConstr);
}

uint8_t rotl8(uint8_t const x, uint c){
	c &= 255;
	return ((x << c) | (x >> (8-c)));
}

uint8_t applyF0(uint8_t const x){
	return (rotl8(x,1) ^ rotl8(x,2) ^ rotl8(x,7));
}

uint8_t applyF1(uint8_t const x){
	return (rotl8(x,3) ^ rotl8(x,4) ^ rotl8(x,6));
}

void genMCFileHIGHT(){
	//Generate the tables and ineq files for HIGHT
	vector<uint32_t> F0(256);
	vector<uint32_t> F1(256);
	for(uint x = 0; x < 256; x++){
		F0[x] = applyF0(x);
		F1[x] = applyF1(x);
	}

	cout << "F0 : " << endl;
	genAndSaveTableIneq(vector<vector<uint32_t>>({F0}),8,8,"./checkData/HIGHTF0",true);
	cout << "F1 : " << endl;
	genAndSaveTableIneq(vector<vector<uint32_t>>({F1}),8,8,"./checkData/HIGHTF1",true);
}

void createMatrixFileHIGHT(){
	//Create the matrix files ./checkData/matrixHIGHTF0.bin and ./checkData/matrixHIGHTF1.bin

	vector<vector<uint8_t>> val0({{0,1,0,0,0,0,1,1},
								  {1,0,1,0,0,0,0,1},
								  {1,1,0,1,0,0,0,0},
								  {0,1,1,0,1,0,0,0},
								  {0,0,1,1,0,1,0,0},
								  {0,0,0,1,1,0,1,0},
								  {0,0,0,0,1,1,0,1},
								  {1,0,0,0,0,1,1,0}});
	Matrix M0(8,8);
	for(uint i = 0; i < 8; i++){
		for(uint j = 0; j < 8; j++)
			M0.set(i,j,val0[i][j]);
	}

	M0.saveToFile("./checkData/matrixHIGHTF0.bin");

	vector<vector<uint8_t>> val1({{0,0,1,0,1,1,0,0},
								  {0,0,0,1,0,1,1,0},
								  {0,0,0,0,1,0,1,1},
								  {1,0,0,0,0,1,0,1},
								  {1,1,0,0,0,0,1,0},
								  {0,1,1,0,0,0,0,1},
								  {1,0,1,1,0,0,0,0},
								  {0,1,0,1,1,0,0,0}});
	Matrix M1(8,8);
	for(uint i = 0; i < 8; i++){
		for(uint j = 0; j < 8; j++)
			M1.set(i,j,val1[i][j]);
	}

	M1.saveToFile("./checkData/matrixHIGHTF1.bin");
}

void checkNormalDistinguisherHIGHT(){
	/*
		Search for a normal (no linear combination) distinguisher over 4 rounds of HIGHT
	*/

	vector<uint8_t> input({0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1});

	//These parameters have no effect for HIGHT but are required for the modelInfo object
	string sboxModel = "QM";
	bool firstSbox = false;
	bool lastLin = true;
	bool smartLin = true;
	unsigned int startingRound = 1;
	CheckType sboxCheckType = CheckType::None;

	//Fixed parameters
	uint rMax = 18;
	string arxModel = "ARX";
	CheckType ARXCheckType = CheckType::None;
	CheckType MCCheckType = CheckType::None;
	CheckType SSBCheckType = CheckType::None;
	uint maxNumberConstr = 1;
	vector<uint> roundOrder = roundOrderOutsideIn(rMax);

	//Variable parameters
	//The simple model without checks already find the known balanced bits
	// string linModel = "Simple";
	// bool useCallback = false;
	vector<string> listLinModel({"QM", "ZR", "CX", "Simple"});
	vector<bool> listUseCB({true});
	
	for(auto const & useCallback : listUseCB){
		for(auto const & linModel : listLinModel){

			globalCtrCheckHIGHT = 0;
			globalCtrSolHIGHT = 0;
			globalCtrConstrSboxHIGHT = 0;
			globalCtrConstrLinHIGHT = 0;
			globalCtrConstrSSBHIGHT = 0;

			modelInfo MI(rMax,sboxModel,linModel,arxModel,firstSbox,lastLin,smartLin,startingRound,sboxCheckType,MCCheckType,SSBCheckType,ARXCheckType,useCallback,maxNumberConstr,roundOrder);

			auto start = chrono::high_resolution_clock::now();
			auto durBalanced = start - start;
			//Check each output bit (half balanced)
			vector<bool> balancedBits(64);
			uint ctrBalanced = 0;
			for(uint i = 0; i < 64; i++){
				vector<uint8_t> output(64,0);
				output[i] = 1;
				auto startBit = chrono::high_resolution_clock::now();
				balancedBits[i] = !existTrailHIGHT(MI,input,output);
				auto endBit = chrono::high_resolution_clock::now();
				if(balancedBits[i]){
					ctrBalanced++;
					durBalanced += endBit - startBit;
				}
			}
			auto end = chrono::high_resolution_clock::now();
			if(useCallback) cout << "Using callback - ";
			else cout << "Using iterative model - ";
			cout << "linModel : " << linModel << endl;
			cout << ctrBalanced << " balanced bits" << endl;
			cout << chrono::duration<double>(end - start).count() << " seconds" << endl;
			cout << "Time for only balanced bits : " << chrono::duration<double>(durBalanced).count() << " second" << endl;
			cout << globalCtrCheckHIGHT << " input/output pairs checked" << endl;
			cout << globalCtrSolHIGHT << " solutions examined" << endl;
			cout << globalCtrConstrSboxHIGHT << " constraints added for the sbox" << endl;
			cout << globalCtrConstrLinHIGHT  << " constraints added for the lin layer" << endl;
			cout << globalCtrConstrSSBHIGHT << " constraints added for the SSB" << endl;
		}
	}
}

void testSearchHIGHT(){

	uint rMax = 19;
	string sboxModel = "QM";
	// string linModel = "CX";
	string arxModel = "ARX";
	bool firstSbox = false;
	bool lastLin = true;
	bool smartLin = true;
	unsigned int startingRound = 1;
	CheckType sboxCheckType = CheckType::None;
	// CheckType MCCheckType = CheckType::Matrix;
	// CheckType SSBCheckType = CheckType::None;
	CheckType ARXCheckType = CheckType::None;
	// bool useCallback = false;
	// uint maxNumberConstr = 4;
	// vector<uint> roundOrder = roundOrderOutsideIn(rMax);
	auto const allTest = allInOutPairsHIGHT(0,0);

	vector<string> listLinModel({"CX", "ZR", "QM", "Simple"});
	// vector<CheckType> listMCCheckType({CheckType::Ineq, CheckType::Table, CheckType::Matrix});
	vector<CheckType> listSSBCheckType({CheckType::Ineq});
	// vector<uint> listMaxNbConstr({100,50,25,10,1});
	vector<bool> listUseCB({true});
	vector<vector<uint>> listroundOrder({roundOrderOutsideIn(rMax), roundOrderAscend(rMax), roundOrderDescend(rMax)});

	uint ctrParameterSet = 0;
	uint firstParameterSet = 0;

	for(auto const & useCallback : listUseCB){
	 for(auto const & roundOrder : listroundOrder){
	  for(auto const & linModel : listLinModel){
	   vector<CheckType> listMCCheckType({CheckType::Ineq});
	   if(linModel=="QM" || linModel=="ZR")
	   	listMCCheckType=vector<CheckType>({CheckType::None});
	   for(auto const & MCCheckType : listMCCheckType){
	    for(auto const & SSBCheckType : listSSBCheckType){
	     vector<uint> listMaxNbConstr({1});
         if(!(sboxCheckType == CheckType::None && MCCheckType == CheckType::None && SSBCheckType == CheckType::None))
          listMaxNbConstr = vector<uint>({500,100,10,1});
	     for(auto const & maxNumberConstr : listMaxNbConstr){
	     	if(ctrParameterSet >= firstParameterSet){

	     	  globalCtrCheckHIGHT = 0;
			  globalCtrSolHIGHT = 0;
	      	  globalCtrConstrSboxHIGHT = 0;
			  globalCtrConstrLinHIGHT = 0;
			  globalCtrConstrSSBHIGHT = 0;

		      modelInfo MI(rMax,sboxModel,linModel,arxModel,firstSbox,lastLin,smartLin,startingRound,sboxCheckType,MCCheckType,SSBCheckType,ARXCheckType,useCallback,maxNumberConstr,roundOrder);
		      cout << "Parameter set " << ctrParameterSet << endl;
		      if(useCallback) cout << "Using callback - ";
		      else cout << "Using iterative model - ";
		      cout << "Round Order : ";
		      for(auto const & tmpr : roundOrder) cout << tmpr << " ";
		      cout << endl;
		      cout << "linModel : " << linModel << endl;
		      cout << "MC Check : " << MCCheckType << " - SSB Check : " << SSBCheckType << " - maxNumberConstr : " << maxNumberConstr << endl;

		      auto start = chrono::high_resolution_clock::now();
		      auto haveDist = checkForDistinguisher(allTest, existTrailHIGHT, MI, false);
		      auto end = chrono::high_resolution_clock::now();
		      cout << "Time : " << chrono::duration<double>(end - start).count() << endl;
		      if(haveDist)cout << "Distinguisher !!!" << endl;
		      else cout << "No distinguisher..." << endl;
		      cout << globalCtrCheckHIGHT << " input/output pairs checked" << endl;
			  cout << globalCtrSolHIGHT << " solutions examined" << endl;
			  cout << globalCtrConstrSboxHIGHT << " constraints added for the sbox" << endl;
			  cout << globalCtrConstrLinHIGHT  << " constraints added for the lin layer" << endl;
			  cout << globalCtrConstrSSBHIGHT << " constraints added for the SSB" << endl;
		    }
		    ctrParameterSet++;
	}}}}}
}
}

//Fastest parameter set 
// Parameter set 27
// Using callback - Round Order : 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 
// linModel : CX
// MC Check : Ineq - SSB Check : Ineq - maxNumberConstr : 1

void searchDistinguisherHIGHT20r(){

	uint rMax = 18;
	string sboxModel = "QM";
	string linModel = "CX";
	string arxModel = "ARX";
	bool firstSbox = false;
	bool lastLin = true;
	bool smartLin = true;
	unsigned int startingRound = 1;
	CheckType sboxCheckType = CheckType::None;
	CheckType MCCheckType = CheckType::Ineq;
	CheckType SSBCheckType = CheckType::Ineq;
	CheckType ARXCheckType = CheckType::None;
	bool useCallback = true;
	uint maxNumberConstr = 1;
	vector<uint> roundOrder = roundOrderAscend(rMax);

	globalCtrCheckHIGHT = 0;
	globalCtrSolHIGHT = 0;
	globalCtrConstrSboxHIGHT = 0;
	globalCtrConstrLinHIGHT = 0;
	globalCtrConstrSSBHIGHT = 0;

	auto globalStart = chrono::high_resolution_clock::now();
	for(uint i = 0; i < 4; i++){
		for(uint j = 0; j < 4; j++){
			cout << "HIGHT Search i = " << i << " j = " << j << endl;
			//Get the list of input/output to test
			auto allTest = allInOutPairsHIGHT(i,j);
			//Permute to get the right order
			for(auto & setPair : allTest){
				for(auto & pair : setPair){
					vector<uint8_t> tmp(64);
					for(uint index = 0; index < 64; index++)
						tmp[(index+8)%64] = pair.first[index];
					pair.first = move(tmp);
				}
			}

			//Generate the modelInfo object
			modelInfo MI(rMax,sboxModel,linModel,arxModel,firstSbox,lastLin,smartLin,startingRound,sboxCheckType,MCCheckType,SSBCheckType,ARXCheckType,useCallback,maxNumberConstr,roundOrder);

			//Check if we can find a distinguisher
			auto start = chrono::high_resolution_clock::now();
			auto haveDist = checkForDistinguisher(allTest, existTrailHIGHT, MI, false);
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
	cout << globalCtrCheckHIGHT << " input/output pairs checked" << endl;
	cout << globalCtrSolHIGHT << " solutions examined" << endl;
	cout << globalCtrConstrSboxHIGHT << " constraints added for the sbox" << endl;
	cout << globalCtrConstrLinHIGHT  << " constraints added for the lin layer" << endl;
	cout << globalCtrConstrSSBHIGHT << " constraints added for the SSB" << endl;
}

void searchDistinguisherHIGHT21r(){

	uint rMax = 19;
	string sboxModel = "QM";
	string linModel = "CX";
	string arxModel = "ARX";
	bool firstSbox = false;
	bool lastLin = true;
	bool smartLin = true;
	unsigned int startingRound = 1;
	CheckType sboxCheckType = CheckType::None;
	CheckType MCCheckType = CheckType::Ineq;
	CheckType SSBCheckType = CheckType::Ineq;
	CheckType ARXCheckType = CheckType::None;
	bool useCallback = true;
	uint maxNumberConstr = 1;
	vector<uint> roundOrder = roundOrderOutsideIn(rMax);

	globalCtrCheckHIGHT = 0;
	globalCtrSolHIGHT = 0;
	globalCtrConstrSboxHIGHT = 0;
	globalCtrConstrLinHIGHT = 0;
	globalCtrConstrSSBHIGHT = 0;

	auto globalStart = chrono::high_resolution_clock::now();
	for(uint i = 0; i < 4; i++){
		for(uint j = 0; j < 4; j++){
			cout << "HIGHT Search i = " << i << " j = " << j << endl;
			//Get the list of input/output to test
			auto allTest = allInOutPairsHIGHT(i,j);
			//Permute to get the right order
			for(auto & setPair : allTest){
				for(auto & pair : setPair){
					vector<uint8_t> tmp(64);
					for(uint index = 0; index < 64; index++)
						tmp[(index+8)%64] = pair.first[index];
					pair.first = move(tmp);
				}
			}

			//Generate the modelInfo object
			modelInfo MI(rMax,sboxModel,linModel,arxModel,firstSbox,lastLin,smartLin,startingRound,sboxCheckType,MCCheckType,SSBCheckType,ARXCheckType,useCallback,maxNumberConstr,roundOrder);

			//Check if we can find a distinguisher
			auto start = chrono::high_resolution_clock::now();
			auto haveDist = checkForDistinguisher(allTest, existTrailHIGHT, MI, false);
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
	cout << globalCtrCheckHIGHT << " input/output pairs checked" << endl;
	cout << globalCtrSolHIGHT << " solutions examined" << endl;
	cout << globalCtrConstrSboxHIGHT << " constraints added for the sbox" << endl;
	cout << globalCtrConstrLinHIGHT  << " constraints added for the lin layer" << endl;
	cout << globalCtrConstrSSBHIGHT << " constraints added for the SSB" << endl;
}
