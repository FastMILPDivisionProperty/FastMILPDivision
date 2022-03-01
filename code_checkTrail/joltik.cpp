#include "joltik.hpp"
#include "DivTable.hpp"

using namespace std;

#define DISPLAY_GEN_OUTPUT false
#define DISPLAY_SOLVER_OUTPUT false
#define DISPLAY_NUMBER_ADDED_CONSTRAINT false

vector<vector<pair<vector<uint8_t>, vector<uint8_t>>>> allInOutPairsJoltik(unsigned pos_in, unsigned pos_out, unsigned R) { // 0 <= i < 4  - 0 <= j < 4
  vector<uint16_t> sbox = {14,4,11,2,3,8,0,9,1,10,7,15,6,12,5,13};

                          vector<vector<uint8_t>> val({{1,0,0,0,0,0,1,0,1,1,0,0,1,1,1,0},
                                          {0,1,0,0,0,0,1,1,0,0,1,0,0,0,0,1},
                                          {0,0,1,0,1,0,0,1,0,0,0,1,1,0,0,0},
                                          {0,0,0,1,0,1,0,0,1,0,0,0,1,1,0,0},
                                          {0,0,1,0,1,0,0,0,1,1,1,0,1,1,0,0},
                                          {0,0,1,1,0,1,0,0,0,0,0,1,0,0,1,0},
                                          {1,0,0,1,0,0,1,0,1,0,0,0,0,0,0,1},
                                          {0,1,0,0,0,0,0,1,1,1,0,0,1,0,0,0},
                                          {1,1,0,0,1,1,1,0,1,0,0,0,0,0,1,0},
                                          {0,0,1,0,0,0,0,1,0,1,0,0,0,0,1,1},
                                          {0,0,0,1,1,0,0,0,0,0,1,0,1,0,0,1},
                                          {1,0,0,0,1,1,0,0,0,0,0,1,0,1,0,0},
                                          {1,1,1,0,1,1,0,0,0,0,1,0,1,0,0,0},
                                          {0,0,0,1,0,0,1,0,0,0,1,1,0,1,0,0},
                                          {1,0,0,0,0,0,0,1,1,0,0,1,0,0,1,0},
                                          {1,1,0,0,1,0,0,0,0,1,0,0,0,0,0,1}});

  vector<uint16_t> linear (16, 0);
  for (unsigned l = 0; l < 16; ++l) {
    for (unsigned b = 0; b < 16; ++b) linear[l] |= uint16_t(val[l][b]) << b;
  }


  vector<uint8_t> perm_inv = {0,10,5,15,14,4,11,1,9,3,12,6,7,13,2,8}; // nibble i <-- nibble perm_inv[i]
  vector<uint8_t> perm (16); // nibble perm[i]  <-- nibble i
  for (unsigned i = 0; i < 16; ++i) perm[perm_inv[i]] = i;

  vector<uint8_t> permBit({0,1,2,3,52,53,54,55,40,41,42,43,28,29,30,31,16,17,18,19,4,5,6,7,56,57,58,59,44,45,46,47,32,33,34,35,20,21,22,23,8,9,10,11,60,61,62,63,48,49,50,51,36,37,38,39,24,25,26,27,12,13,14,15});
  vector<uint8_t> permBit_inv(64);
  for(uint i = 0; i < 64; i++)
  	permBit_inv[permBit[i]] = i;

  auto anf_s = getANF_Parallel_S(4, sbox);
  auto anf_l = convertLtoPolynome(linear);

  auto anf_los = composition(anf_l, anf_s);
  auto anf_kolos = anf_los;
  for (unsigned k = 0; k < 16; ++k) anf_kolos[k] += Monome(uint32_t(1) << (16 + k));
  auto anf_sokolos = composition(anf_s, anf_kolos);

  auto anf_sol = composition(anf_s, anf_l);

  auto v_in = loadTransitionsIN(anf_sokolos);
  cout << "v_in : " << v_in.size() << endl;
  auto v_out = (R%2 == 0) ? loadTransitionsOUT(anf_sokolos) : loadTransitionsOUT(anf_sol);
  cout << "v_out : " << v_out.size() << endl;

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

GRBModel getModelJoltik(modelInfo & MI){
/*
	Get the model for Joltik over #rMax round
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
    string modelName = "Joltik_"+to_string(MI.rMax)+"r_"+MI.sboxModel+"_"+MI.linModel;
    if(MI.lastLin) modelName += "_lastLin";
    if(!MI.firstSbox) modelName += "_noFirstSbox";
    modelName += ".mps";

    if(!fileExist("./models/"+modelName)){
        string cmd = "sage ./modelGeneration/genJoltik.sage";
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

modelData getVariablesJoltik(GRBModel & m,
	                         modelInfo const & MI){
     /*
    Read the model and return the corresponding blocks of variables
    - #rMax defines the number of rounds
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


    if(MI.MCCheckType == CheckType::Matrix)
    	MD.allMatrix[0] = Matrix("./checkData/matrixJoltik.bin");
    else if(MI.MCCheckType == CheckType::Table){
    	cerr << "Error : Table Check for Joltik MC not implemented" << endl;
    	exit(1);
        // MD.allTable[1] = //Table for MC
    }
    else if(MI.MCCheckType == CheckType::Ineq){
        MD.allIneq[1]   = readIneqFromFile("./checkData/JoltikMC_ineq.bin");
    }

    if(MI.SSBCheckType == CheckType::Table)
        MD.allTable[2] = readDivTableFromFile("./checkData/JoltikSSB_table.bin");
    else if(MI.SSBCheckType == CheckType::Ineq)
        MD.allIneq[2]   = readIneqFromFile("./checkData/JoltikSSB_ineq.bin");

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
	                vector<GRBVar> in(4);
	                vector<GRBVar> out(4);
	                for(uint j = 0; j < 4; j++){
	                    in[j] = m.getVarByName("x"+to_string(r)+"_"+to_string(4*i+j));
	                    out[j] = m.getVarByName("y"+to_string(r)+"_"+to_string(4*i+j));
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
	                    MD.allVarsTable.emplace_back(move(in), move(out), tableMC);
	                else if(MI.MCCheckType == CheckType::Ineq)
	                    MD.allVarsIneq.emplace_back(move(in), move(out), ineqMC);
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

bool existTrailJoltik(GRBModel & m,
                   modelData const & MD,
                   vector<uint8_t> const & input,
                   vector<uint8_t> const & output,
                   bool const useCallback,
                   uint const maxNumberConstr){
/*
    Return true if there is a trail from #input to #output using #model for Joltik
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

bool existTrailJoltik(modelInfo & MI,
                   vector<uint8_t> const & input,
                   vector<uint8_t> const & output){
/*
    Return true if there is a trail from #input to #output for Joltik
    #MI should contain the necessary information to generate the model and parametrize the solving
*/

    GRBModel m = getModelJoltik(MI);
    auto MD = getVariablesJoltik(m, MI);

    return existTrailJoltik(m,MD,input,output,MI.useCallback,MI.maxNumberConstr);
}

void createMatrixFileJoltik(){
	//Create the matrix file ./checkData/matrixJoltik.bin

	vector<vector<uint8_t>> val({{1,0,0,0,0,0,1,0,1,1,0,0,1,1,1,0},
								 {0,1,0,0,0,0,1,1,0,0,1,0,0,0,0,1},
								 {0,0,1,0,1,0,0,1,0,0,0,1,1,0,0,0},
								 {0,0,0,1,0,1,0,0,1,0,0,0,1,1,0,0},
								 {0,0,1,0,1,0,0,0,1,1,1,0,1,1,0,0},
								 {0,0,1,1,0,1,0,0,0,0,0,1,0,0,1,0},
								 {1,0,0,1,0,0,1,0,1,0,0,0,0,0,0,1},
								 {0,1,0,0,0,0,0,1,1,1,0,0,1,0,0,0},
								 {1,1,0,0,1,1,1,0,1,0,0,0,0,0,1,0},
								 {0,0,1,0,0,0,0,1,0,1,0,0,0,0,1,1},
								 {0,0,0,1,1,0,0,0,0,0,1,0,1,0,0,1},
								 {1,0,0,0,1,1,0,0,0,0,0,1,0,1,0,0},
								 {1,1,1,0,1,1,0,0,0,0,1,0,1,0,0,0},
								 {0,0,0,1,0,0,1,0,0,0,1,1,0,1,0,0},
								 {1,0,0,0,0,0,0,1,1,0,0,1,0,0,1,0},
								 {1,1,0,0,1,0,0,0,0,1,0,0,0,0,0,1}});
	Matrix M(16,16);
	for(uint i = 0; i < 16; i++){
		for(uint j = 0; j < 16; j++)
			M.set(i,j,val[i][j]);
	}

	M.saveToFile("./checkData/matrixJoltik.bin");
}

void computeSboxMatrixJoltik(){
	// vector<vector<uint8_t>> Mb({{1,0,0,0,0,0,1,0,1,1,0,0,1,1,1,0},
	// 							{0,1,0,0,0,0,1,1,0,0,1,0,0,0,0,1},
	// 							{0,0,1,0,1,0,0,1,0,0,0,1,1,0,0,0},
	// 							{0,0,0,1,0,1,0,0,1,0,0,0,1,1,0,0},
	// 							{0,0,1,0,1,0,0,0,1,1,1,0,1,1,0,0},
	// 							{0,0,1,1,0,1,0,0,0,0,0,1,0,0,1,0},
	// 							{1,0,0,1,0,0,1,0,1,0,0,0,0,0,0,1},
	// 							{0,1,0,0,0,0,0,1,1,1,0,0,1,0,0,0},
	// 							{1,1,0,0,1,1,1,0,1,0,0,0,0,0,1,0},
	// 							{0,0,1,0,0,0,0,1,0,1,0,0,0,0,1,1},
	// 							{0,0,0,1,1,0,0,0,0,0,1,0,1,0,0,1},
	// 							{1,0,0,0,1,1,0,0,0,0,0,1,0,1,0,0},
	// 							{1,1,1,0,1,1,0,0,0,0,1,0,1,0,0,0},
	// 							{0,0,0,1,0,0,1,0,0,0,1,1,0,1,0,0},
	// 							{1,0,0,0,0,0,0,1,1,0,0,1,0,0,1,0},
	// 							{1,1,0,0,1,0,0,0,0,1,0,0,0,0,0,1}});
	
	// vector<uint16_t> M(16);
	// for(uint i = 0; i < 16; i++){
	// 	uint16_t x = 0;
	// 	for(uint j = 0; j < 16; j++){
	// 		if(Mb[i][j] == 1)
	// 			x += (1 << j);
	// 	}
	// 	M[i] = x;
	// }

	// for(auto const & x : M)
	// 	cout << x << ",";
	// cout << endl;

	vector<uint16_t> M({29505,33986,6292,12584,14100,18476,33097,4994,16755,49796,37912,10289,5175,11336,18817,33299,});

	uint32_t bound = (1ul << 16);
	vector<uint32_t> S(bound,0);
	for(uint32_t x = 0; x < bound; x++){
		uint16_t y = 0;
		for(uint j = 0; j < 16; j++){
			if(__builtin_parity(M[j] & x))
				y |= (1 << j);
		}
		S[x] = y;
	}

	genAndSaveTableIneq(vector<vector<uint32_t>>({S}),16,16,"./JoltikMC",true);
}

void checkNormalDistinguisherJoltik(){
	/*
		Search for a normal (no linear combination) distinguisher over 5 rounds of Joltik 
	*/

	vector<uint8_t> input({1,1,1,1,0,0,0,0,0,0,0,0,1,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1});

	//These parameters have no effect for Joltik but are required for the modelInfo object
	string arxModel = "ARX";
	unsigned int startingRound = 1;
	CheckType ARXCheckType = CheckType::None;

	//Fixed parameters
	CheckType SSBCheckType = CheckType::None; //no need to check SSB for this one
	uint rMax = 5;
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
	vector<string> listLinModel({"CX", "Simple", "QM"});
	vector<bool> listUseCallback({true});
	vector<vector<uint>> listroundOrderIfNeeded({roundOrderOutsideIn(rMax), roundOrderAscend(rMax), roundOrderDescend(rMax)});
	// vector<uint> listMaxNbConstr({500,100,10,1});

	uint ctrParameterSet = 0;
	uint firstParameterSet = 5;

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
				//Check each output bit (half balanced)
				vector<bool> balancedBits(64);
				uint ctrBalanced = 0;
				for(uint i = 0; i < 64; i++){
					vector<uint8_t> output(64,0);
					output[i] = 1;
					auto startBit = chrono::high_resolution_clock::now();
					balancedBits[i] = !existTrailJoltik(MI,input,output);
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
			}
			ctrParameterSet++;
	}}}}}}}
}


void timingSearchJoltik(){
	//These parameters have no effect for Joltik but are required for the modelInfo object
	string arxModel = "ARX";
	unsigned int startingRound = 1;
	CheckType ARXCheckType = CheckType::None;

	//Fixed parameters
	CheckType SSBCheckType = CheckType::Ineq; 
	uint rMax = 4; //7 Total rounds, (SB MC SB) + (MC  3 rounds) + (SB MC SB)
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
	vector<string> listLinModel({"QM", "CX", "Simple"});
	vector<bool> listUseCallback({true});
	vector<vector<uint>> listroundOrderIfNeeded({roundOrderOutsideIn(rMax), roundOrderAscend(rMax), roundOrderDescend(rMax)});
	auto const allTest = allInOutPairsJoltik(0,0,0);

	uint ctrParameterSet = 0;
	uint firstParameterSet = 0;

	for(auto const & useCallback : listUseCallback){
	 for(auto const & linModel : listLinModel){
	  vector<CheckType> listMCCheckType({CheckType::Ineq});
	  if(linModel == "QM") //No need for checks on MC
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
				 auto haveDist = checkForDistinguisher(allTest, existTrailJoltik, MI, false);
				 auto end = chrono::high_resolution_clock::now();
	 
				 cout << "Time : " << chrono::duration<double>(end - start).count() << endl;
				 if(haveDist)cout << "Distinguisher !!!" << endl << endl;
				 else cout << "No distinguisher..." << endl << endl;
			 }
			 ctrParameterSet++;
	}}}}}}}
}