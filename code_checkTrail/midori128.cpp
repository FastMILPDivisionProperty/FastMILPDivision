#include "midori128.hpp"
#include "DivTable.hpp"

using namespace std;

#define DISPLAY_GEN_OUTPUT false
#define DISPLAY_SOLVER_OUTPUT true
#define DISPLAY_NUMBER_ADDED_CONSTRAINT true



vector<vector<pair<vector<uint8_t>, vector<uint8_t>>>> allInOutPairsMidori128(unsigned i, unsigned j) {
  vector<uint16_t> SSb0 = {
  0x11, 0x10, 0x51, 0x50, 0xb4, 0x30, 0xf4, 0x70, 0x59, 0x58, 0x19, 0x18, 0xfc, 0x78, 0xbc, 0x38,
  0x01, 0x00, 0x13, 0x12, 0xa4, 0x20, 0xb6, 0x32, 0x0b, 0x0a, 0x1b, 0x1a, 0xae, 0x2a, 0xbe, 0x3a,
  0x15, 0x31, 0x55, 0x71, 0xb5, 0x35, 0xf5, 0x75, 0x5d, 0x79, 0x1d, 0x39, 0xfd, 0x7d, 0xbd, 0x3d,
  0x05, 0x21, 0x17, 0x33, 0xa5, 0x25, 0xb7, 0x37, 0x0f, 0x2b, 0x1f, 0x3b, 0xaf, 0x2f, 0xbf, 0x3f,
  0x4b, 0x4a, 0x5b, 0x5a, 0xee, 0x6a, 0xfe, 0x7a, 0x49, 0x48, 0x41, 0x40, 0xec, 0x68, 0xe4, 0x60,
  0x03, 0x02, 0x53, 0x52, 0xa6, 0x22, 0xf6, 0x72, 0x09, 0x08, 0x43, 0x42, 0xac, 0x28, 0xe6, 0x62,
  0x4f, 0x6b, 0x5f, 0x7b, 0xef, 0x6f, 0xff, 0x7f, 0x4d, 0x69, 0x45, 0x61, 0xed, 0x6d, 0xe5, 0x65,
  0x07, 0x23, 0x57, 0x73, 0xa7, 0x27, 0xf7, 0x77, 0x0d, 0x29, 0x47, 0x63, 0xad, 0x2d, 0xe7, 0x67,
  0x95, 0xb0, 0xd5, 0xf0, 0x94, 0x90, 0xd4, 0xd0, 0xdd, 0xf8, 0x9d, 0xb8, 0xdc, 0xd8, 0x9c, 0x98,
  0x85, 0xa0, 0x97, 0xb2, 0x84, 0x80, 0x96, 0x92, 0x8f, 0xaa, 0x9f, 0xba, 0x8e, 0x8a, 0x9e, 0x9a,
  0x91, 0xb1, 0xd1, 0xf1, 0x14, 0x34, 0x54, 0x74, 0xd9, 0xf9, 0x99, 0xb9, 0x5c, 0x7c, 0x1c, 0x3c,
  0x81, 0xa1, 0x93, 0xb3, 0x04, 0x24, 0x16, 0x36, 0x8b, 0xab, 0x9b, 0xbb, 0x0e, 0x2e, 0x1e, 0x3e,
  0xcf, 0xea, 0xdf, 0xfa, 0xce, 0xca, 0xde, 0xda, 0xcd, 0xe8, 0xc5, 0xe0, 0xcc, 0xc8, 0xc4, 0xc0,
  0x87, 0xa2, 0xd7, 0xf2, 0x86, 0x82, 0xd6, 0xd2, 0x8d, 0xa8, 0xc7, 0xe2, 0x8c, 0x88, 0xc6, 0xc2,
  0xcb, 0xeb, 0xdb, 0xfb, 0x4e, 0x6e, 0x5e, 0x7e, 0xc9, 0xe9, 0xc1, 0xe1, 0x4c, 0x6c, 0x44, 0x64,
  0x83, 0xa3, 0xd3, 0xf3, 0x06, 0x26, 0x56, 0x76, 0x89, 0xa9, 0xc3, 0xe3, 0x0c, 0x2c, 0x46, 0x66};


  vector<uint16_t> SSb1 = {
  0x88, 0x8a, 0x4b, 0xcb, 0xac, 0xae, 0x6f, 0xef, 0x80, 0x82, 0x43, 0xc3, 0x94, 0x96, 0x57, 0xd7,
  0xa8, 0xaa, 0x6b, 0xeb, 0x8c, 0x8e, 0x4f, 0xcf, 0x98, 0x9a, 0x5b, 0xdb, 0x9c, 0x9e, 0x5f, 0xdf,
  0xb4, 0xb6, 0x77, 0xf7, 0xa4, 0xa6, 0x67, 0xe7, 0x90, 0x92, 0x53, 0xd3, 0x84, 0x86, 0x47, 0xc7,
  0xbc, 0xbe, 0x7f, 0xff, 0xa0, 0xa2, 0x63, 0xe3, 0xb8, 0xba, 0x7b, 0xfb, 0xb0, 0xb2, 0x73, 0xf3,
  0xca, 0xc8, 0x4a, 0x0a, 0xee, 0xec, 0x6e, 0x2e, 0xc2, 0xc0, 0x42, 0x02, 0xd6, 0xd4, 0x56, 0x16,
  0xea, 0xe8, 0x6a, 0x2a, 0xce, 0xcc, 0x4e, 0x0e, 0xda, 0xd8, 0x5a, 0x1a, 0xde, 0xdc, 0x5e, 0x1e,
  0xf6, 0xf4, 0x76, 0x36, 0xe6, 0xe4, 0x66, 0x26, 0xd2, 0xd0, 0x52, 0x12, 0xc6, 0xc4, 0x46, 0x06,
  0xfe, 0xfc, 0x7e, 0x3e, 0xe2, 0xe0, 0x62, 0x22, 0xfa, 0xf8, 0x7a, 0x3a, 0xf2, 0xf0, 0x72, 0x32,
  0x08, 0x89, 0x09, 0x8b, 0x2c, 0xad, 0x2d, 0xaf, 0x00, 0x81, 0x01, 0x83, 0x14, 0x95, 0x15, 0x97,
  0x28, 0xa9, 0x29, 0xab, 0x0c, 0x8d, 0x0d, 0x8f, 0x18, 0x99, 0x19, 0x9b, 0x1c, 0x9d, 0x1d, 0x9f,
  0x34, 0xb5, 0x35, 0xb7, 0x24, 0xa5, 0x25, 0xa7, 0x10, 0x91, 0x11, 0x93, 0x04, 0x85, 0x05, 0x87,
  0x3c, 0xbd, 0x3d, 0xbf, 0x20, 0xa1, 0x21, 0xa3, 0x38, 0xb9, 0x39, 0xbb, 0x30, 0xb1, 0x31, 0xb3,
  0x49, 0xc9, 0x48, 0x0b, 0x6d, 0xed, 0x6c, 0x2f, 0x41, 0xc1, 0x40, 0x03, 0x55, 0xd5, 0x54, 0x17,
  0x69, 0xe9, 0x68, 0x2b, 0x4d, 0xcd, 0x4c, 0x0f, 0x59, 0xd9, 0x58, 0x1b, 0x5d, 0xdd, 0x5c, 0x1f,
  0x75, 0xf5, 0x74, 0x37, 0x65, 0xe5, 0x64, 0x27, 0x51, 0xd1, 0x50, 0x13, 0x45, 0xc5, 0x44, 0x07,
  0x7d, 0xfd, 0x7c, 0x3f, 0x61, 0xe1, 0x60, 0x23, 0x79, 0xf9, 0x78, 0x3b, 0x71, 0xf1, 0x70, 0x33};



  vector<uint16_t> SSb2 = {
  0x44, 0xc3, 0x47, 0x43, 0x40, 0xc0, 0xc2, 0x42, 0x54, 0xd3, 0x57, 0x53, 0x50, 0xd0, 0xd2, 0x52,
  0x3c, 0xbb, 0x3f, 0x3b, 0x38, 0xb8, 0xba, 0x3a, 0x7c, 0xfb, 0x7f, 0x7b, 0x78, 0xf8, 0xfa, 0x7a,
  0x74, 0xf3, 0x77, 0x73, 0x70, 0xf0, 0xf2, 0x72, 0x64, 0xe3, 0x67, 0x63, 0x60, 0xe0, 0xe2, 0x62,
  0x34, 0xb3, 0x37, 0x33, 0x30, 0xb0, 0xb2, 0x32, 0x14, 0x93, 0x17, 0x13, 0x10, 0x90, 0x92, 0x12,
  0x04, 0x83, 0x07, 0x03, 0x00, 0x80, 0x82, 0x02, 0x4c, 0xcb, 0x4f, 0x4b, 0x48, 0xc8, 0xca, 0x4a,
  0x0c, 0x8b, 0x0f, 0x0b, 0x08, 0x88, 0x8a, 0x0a, 0x5c, 0xdb, 0x5f, 0x5b, 0x58, 0xd8, 0xda, 0x5a,
  0x2c, 0xab, 0x2f, 0x2b, 0x28, 0xa8, 0xaa, 0x2a, 0x6c, 0xeb, 0x6f, 0x6b, 0x68, 0xe8, 0xea, 0x6a,
  0x24, 0xa3, 0x27, 0x23, 0x20, 0xa0, 0xa2, 0x22, 0x1c, 0x9b, 0x1f, 0x1b, 0x18, 0x98, 0x9a, 0x1a,
  0x45, 0xc7, 0x46, 0x41, 0xc4, 0xc5, 0xc6, 0xc1, 0x55, 0xd7, 0x56, 0x51, 0xd4, 0xd5, 0xd6, 0xd1,
  0x3d, 0xbf, 0x3e, 0x39, 0xbc, 0xbd, 0xbe, 0xb9, 0x7d, 0xff, 0x7e, 0x79, 0xfc, 0xfd, 0xfe, 0xf9,
  0x75, 0xf7, 0x76, 0x71, 0xf4, 0xf5, 0xf6, 0xf1, 0x65, 0xe7, 0x66, 0x61, 0xe4, 0xe5, 0xe6, 0xe1,
  0x35, 0xb7, 0x36, 0x31, 0xb4, 0xb5, 0xb6, 0xb1, 0x15, 0x97, 0x16, 0x11, 0x94, 0x95, 0x96, 0x91,
  0x05, 0x87, 0x06, 0x01, 0x84, 0x85, 0x86, 0x81, 0x4d, 0xcf, 0x4e, 0x49, 0xcc, 0xcd, 0xce, 0xc9,
  0x0d, 0x8f, 0x0e, 0x09, 0x8c, 0x8d, 0x8e, 0x89, 0x5d, 0xdf, 0x5e, 0x59, 0xdc, 0xdd, 0xde, 0xd9,
  0x2d, 0xaf, 0x2e, 0x29, 0xac, 0xad, 0xae, 0xa9, 0x6d, 0xef, 0x6e, 0x69, 0xec, 0xed, 0xee, 0xe9,
  0x25, 0xa7, 0x26, 0x21, 0xa4, 0xa5, 0xa6, 0xa1, 0x1d, 0x9f, 0x1e, 0x19, 0x9c, 0x9d, 0x9e, 0x99};



  vector<uint16_t> SSb3 = {
  0x22, 0x2b, 0x20, 0x29, 0xa2, 0xab, 0x26, 0x2f, 0x4b, 0x0b, 0x49, 0x09, 0xcb, 0x8b, 0x4f, 0x0f,
  0xb2, 0xbb, 0x34, 0x3d, 0x32, 0x3b, 0x36, 0x3f, 0xdb, 0x9b, 0x5d, 0x1d, 0x5b, 0x1b, 0x5f, 0x1f,
  0x02, 0x43, 0x00, 0x41, 0x82, 0xc3, 0x06, 0x47, 0x42, 0x03, 0x40, 0x01, 0xc2, 0x83, 0x46, 0x07,
  0x92, 0xd3, 0x14, 0x55, 0x12, 0x53, 0x16, 0x57, 0xd2, 0x93, 0x54, 0x15, 0x52, 0x13, 0x56, 0x17,
  0x2a, 0x23, 0x28, 0x21, 0xaa, 0xa3, 0x2e, 0x27, 0x6b, 0x0a, 0x69, 0x08, 0xeb, 0x8a, 0x6f, 0x0e,
  0xba, 0xb3, 0x3c, 0x35, 0x3a, 0x33, 0x3e, 0x37, 0xfb, 0x9a, 0x7d, 0x1c, 0x7b, 0x1a, 0x7f, 0x1e,
  0x62, 0x63, 0x60, 0x61, 0xe2, 0xe3, 0x66, 0x67, 0x6a, 0x4a, 0x68, 0x48, 0xea, 0xca, 0x6e, 0x4e,
  0xf2, 0xf3, 0x74, 0x75, 0x72, 0x73, 0x76, 0x77, 0xfa, 0xda, 0x7c, 0x5c, 0x7a, 0x5a, 0x7e, 0x5e,
  0xb4, 0xbd, 0x24, 0x2d, 0xb6, 0xbf, 0xa6, 0xaf, 0xdd, 0x9d, 0x4d, 0x0d, 0xdf, 0x9f, 0xcf, 0x8f,
  0xb0, 0xb9, 0x30, 0x39, 0xa0, 0xa9, 0xa4, 0xad, 0xd9, 0x99, 0x59, 0x19, 0xc9, 0x89, 0xcd, 0x8d,
  0x94, 0xd5, 0x04, 0x45, 0x96, 0xd7, 0x86, 0xc7, 0xd4, 0x95, 0x44, 0x05, 0xd6, 0x97, 0xc6, 0x87,
  0x90, 0xd1, 0x10, 0x51, 0x80, 0xc1, 0x84, 0xc5, 0xd0, 0x91, 0x50, 0x11, 0xc0, 0x81, 0xc4, 0x85,
  0xbc, 0xb5, 0x2c, 0x25, 0xbe, 0xb7, 0xae, 0xa7, 0xfd, 0x9c, 0x6d, 0x0c, 0xff, 0x9e, 0xef, 0x8e,
  0xb8, 0xb1, 0x38, 0x31, 0xa8, 0xa1, 0xac, 0xa5, 0xf9, 0x98, 0x79, 0x18, 0xe9, 0x88, 0xed, 0x8c,
  0xf4, 0xf5, 0x64, 0x65, 0xf6, 0xf7, 0xe6, 0xe7, 0xfc, 0xdc, 0x6c, 0x4c, 0xfe, 0xde, 0xee, 0xce,
  0xf0, 0xf1, 0x70, 0x71, 0xe0, 0xe1, 0xe4, 0xe5, 0xf8, 0xd8, 0x78, 0x58, 0xe8, 0xc8, 0xec, 0xcc};


  vector<vector<uint16_t>> sboxes ({SSb0, SSb1, SSb2, SSb3});

  auto in = loadTransitionsIN(getANF_S(sboxes[i%4]));
  auto out = loadTransitionsOUT(getANF_S(sboxes[j%4]));

  vector<uint8_t> permBit({0,1,2,3,4,5,6,7,56,57,58,59,60,61,62,63,112,113,114,115,116,117,118,119,72,73,74,75,76,77,78,79,40,41,42,43,44,45,46,47,16,17,18,19,20,21,22,23,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,120,121,122,123,124,125,126,127,64,65,66,67,68,69,70,71,8,9,10,11,12,13,14,15,48,49,50,51,52,53,54,55,80,81,82,83,84,85,86,87,104,105,106,107,108,109,110,111,32,33,34,35,36,37,38,39,24,25,26,27,28,29,30,31});
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
          tmp.emplace_back(move(in_8t), move(out_8t));
        }
      }
      res.emplace_back(move(tmp));
    }
  }

  return res;

}

GRBModel getModelMidori128(modelInfo & MI){
/*
	Get the model for Midori128 over #rMax round
	#lastLin defines if the linear layer is applied on the last round
	#smartLin only has effect if linModel=="Simple" and defines if MC is seen as independent Lboxes or as a single matrix
	#firstSbox defines if the first Sbox layer is applied

	The #sboxModel parameter defines how the Sbox is modelized :
	- "Hull" use the Convex Hull technique
	- "Simple" use the simplified constraint with PWL

	The #linModel parameter defines how the linear layer is modelized :
	- "Lbox" modelize the linear layer as the parallel application of the same linear Sbox (thanks to the specific form of the matrix)
	- "ZR" modelize the linear layer using the technique from Zhang and Rijmen
	- "CX" modelize the linear layer using the classical copy+xor technique
	- "Simple" modelize the linear layer with the simplified constraint w(x) = w(y)
*/
    string modelName = "Midori128_"+to_string(MI.rMax)+"r_"+MI.sboxModel+"_"+MI.linModel;
    if(MI.lastLin) modelName += "_lastLin";
    if(MI.linModel=="Simple" && MI.smartLin) modelName += "_smartLin";
    if(!MI.firstSbox) modelName += "_noFirstSbox";
    modelName += ".mps";

    if(!fileExist("./models/"+modelName)){
        string cmd = "sage ./modelGeneration/genMidori128.sage";
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

modelData getVariablesMidori128(GRBModel & m,
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
        MD.allTable[0] = vector<vector<uint32_t>>({{0},{1,2,4,8},{1,2,4,8},{1,2,4},{1,2,4,8},{1,2,4,8},{1,2,4,8},{1,2},{1,2,4,8},{2,4,8},{1,2,8},{2,9,12},{1,2,4,8},{2,4,8},{1,8},{15}});
    else if(MI.sboxCheckType == CheckType::Ineq)
        MD.allIneq[0] = vector<vector<int>>({{1, 1, 4, 1, -2, -2, -2, -2, 1, 1},
											 {0, 0, 0, 3, -1, -1, -1, -1, 1, 1},
											 {-2, -5, -3, -3, 2, 1, -1, 2, 9, 1},
											 {0, 3, 0, 0, -1, -1, -1, -1, 1, 1},
											 {-2, -2, -1, -1, 5, 5, 4, 4, 0, 1},
											 {0, -1, -1, -1, -1, -2, 3, -1, 4, 1},
											 {-1, 0, 0, -1, -1, 0, 0, 1, 2, 1},
											 {1, -1, 0, -1, 0, 0, -1, 0, 2, 1},
											 {-1, -1, 0, 0, 1, 1, 1, 0, 1, 1}});

    if(MI.smartLin){
    	if(MI.MCCheckType == CheckType::Table)
		    MD.allTable[1] = vector<vector<uint32_t>>({{0},{2,4,8},{1,4,8},{3,5,6,9,10},{1,2,8},{3,5,6,9,12},{3,5,6,10,12},{11,13,14},{1,2,4},{3,5,9,10,12},{3,6,9,10,12},{7,13,14},{5,6,9,10,12},{7,11,14},{7,11,13},{15}});
		else if(MI.MCCheckType == CheckType::Ineq)
		    MD.allIneq[1] = vector<vector<int>>({{1, 1, 1, 1, -1, -1, -1, -1, 0, 0},
												 {0, 0, 0, 1, -1, -1, -1, 0, 2, 1},
												 {0, 0, 1, 0, -1, -1, 0, -1, 2, 1},
												 {0, 1, 0, 0, -1, 0, -1, -1, 2, 1},
												 {0, -1, -1, -1, 1, 0, 0, 0, 2, 1},
												 {0, 1, 0, 1, 0, -1, 0, -1, 1, 1},
												 {0, -1, -1, 0, 0, 1, 1, 0, 1, 1},
												 {0, 0, -1, -1, 0, 0, 1, 1, 1, 1},
												 {0, -1, 0, 0, 1, 0, 1, 1, 0, 1},
												 {0, 0, -1, 0, 1, 1, 0, 1, 0, 1},
												 {0, -1, 0, -1, 0, 1, 0, 1, 1, 1},
												 {0, 0, 1, 1, 0, 0, -1, -1, 1, 1},
												 {0, 1, 1, 1, -1, 0, 0, 0, 0, 1},
												 {0, 0, 0, -1, 1, 1, 1, 0, 0, 1},
												 {0, 1, 1, 0, 0, -1, -1, 0, 1, 1}});
		else if(MI.MCCheckType == CheckType::Matrix)
			MD.allMatrix[0] = Matrix("./checkData/matrixMidori128_small.bin");
    }
    else{
    	if(MI.MCCheckType == CheckType::Matrix)
    		MD.allMatrix[0] = Matrix("./checkData/matrixMidori128.bin");
	    else if(MI.MCCheckType == CheckType::Table){
	    	cerr << "Error : Table Check for Midori128 MC (non smartLin) not implemented, probably not worth it compared to using smartLin=true" << endl;
	    	exit(1);
	        // MD.allTable[1] = //Table for MC
	    }
	    else if(MI.MCCheckType == CheckType::Ineq){
	    	cerr << "Error : Ineq Check for Midori128 MC (non smartLin) not implemented, probably not worth it compared to using smartLin=true" << endl;
	    	exit(1);
	        // MD.allIneq[1]   = //Ineq for MC
	    }
    }

    if(MI.SSBCheckType != CheckType::None){
    	cerr << "Error : SSB Check for Midori128 not implemented" << endl;
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

    static const uint8_t index[4][8] = {{0, 5, 2, 7, 4, 1, 6, 3},
										{3, 4, 5, 2, 7, 0, 1, 6},
										{2, 7, 0, 1, 6, 3, 4, 5},
										{1, 2, 7, 4, 5, 6, 3, 0}};
    //Sbox variables
	uint roundFirstSbox = 0;
    if(!MI.firstSbox) roundFirstSbox = 1;

    if(MI.sboxCheckType != CheckType::None){
        for(auto const & r : MI.roundOrder){
        	if(r >= roundFirstSbox){
	            for(uint i = 0; i < 16; i++){
	            	//Sbox in Midori128 has a specific form with two 4-bit sboxes and bit permutations depending on the index
	            	//Code is a bit ugly as a result

					vector<GRBVar> in1(4);
	            	vector<GRBVar> out1(4);
	            	vector<GRBVar> in2(4);
	            	vector<GRBVar> out2(4);
	            	for(uint j = 0; j < 4; j++){
	            		in1[j] = m.getVarByName("x"+to_string(r)+"_"+to_string(4*i+index[i%4][j]));
	            		out1[j] = m.getVarByName("x"+to_string(r)+"_"+to_string(4*i+index[i%4][j]));
	            		in2[j] = m.getVarByName("x"+to_string(r)+"_"+to_string(4*i+index[i%4][j+4]));
	            		out2[j] = m.getVarByName("x"+to_string(r)+"_"+to_string(4*i+index[i%4][j+4]));
	            	}
	            	if(MI.sboxCheckType == CheckType::Table){
	                    MD.allVarsTable.emplace_back(move(in1), move(out1), tableSbox);
	                    MD.allVarsTable.emplace_back(move(in2), move(out2), tableSbox);
	                }
	                else if(MI.sboxCheckType == CheckType::Ineq){
	                    MD.allVarsIneq.emplace_back(move(in1), move(out1), ineqSbox);
	                    MD.allVarsIneq.emplace_back(move(in2), move(out2), ineqSbox);
	                }
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
    //{ 0, 10,  5, 15} -> {0,1,2,3}
    //{14,  4, 11,  1} -> {4,5,6,7}
    //{ 9,  3, 12,  6} -> {8,9,10,11}
    //{ 7, 13,  2,  8} -> {12,13,14,15}
    static const uint8_t indexIn[4][4] = {{0,10,5,15},{14,4,11,1},{9,3,12,6},{7,13,2,8}};
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

bool existTrailMidori128(GRBModel & m,
                   modelData const & MD,
                   vector<uint8_t> const & input,
                   vector<uint8_t> const & output,
                   bool const useCallback,
                   uint const maxNumberConstr){
/*
    Return true if there is a trail from #input to #output using #model for Midori128
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

bool existTrailMidori128(modelInfo & MI,
                   vector<uint8_t> const & input,
                   vector<uint8_t> const & output){
/*
    Return true if there is a trail from #input to #output for Midori128
    #MI should contain the necessary information to generate the model and parametrize the solving
*/

    GRBModel m = getModelMidori128(MI);
    auto MD = getVariablesMidori128(m, MI);

    return existTrailMidori128(m,MD,input,output,MI.useCallback,MI.maxNumberConstr);
}

void createMatrixFileMidori128(){
	//Create the matrix file ./checkData/matrixMidori128.bin

	vector<vector<uint8_t>> val({{0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
								 {0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0},
								 {0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0},
								 {0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0},
								 {0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0},
								 {0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0},
								 {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0},
								 {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1},
								 {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
								 {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0},
								 {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0},
								 {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0},
								 {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0},
								 {0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0},
								 {0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0},
								 {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1},
								 {1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
								 {0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0},
								 {0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0},
								 {0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0},
								 {0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0},
								 {0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0},
								 {0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0},
								 {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
								 {1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
								 {0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
								 {0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
								 {0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0},
								 {0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0},
								 {0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0},
								 {0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0},
								 {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0}});
	Matrix M(32,32);
	for(uint i = 0; i < 32; i++){
		for(uint j = 0; j < 32; j++)
			M.set(i,j,val[i][j]);
	}

	M.saveToFile("./checkData/matrixMidori128.bin");

	val = vector<vector<uint8_t>>({{0,1,1,1},
								   {1,0,1,1},
								   {1,1,0,1},
								   {1,1,1,0}});
	M = Matrix(4,4);
	for(uint i = 0; i < 4; i++){
		for(uint j = 0; j < 4; j++)
			M.set(i,j,val[i][j]);
	}
	M.saveToFile("./checkData/matrixMidori128_small.bin");
}

void checkNormalDistinguisherMidori128(){
	/*
		Search for a normal (no linear combination) distinguisher over 7 rounds of Midori128
	*/

	vector<uint8_t> input({1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1});

	//These parameters have no effect for Midori128 but are required for the modelInfo object
	string arxModel = "ARX";
	unsigned int startingRound = 1;
	CheckType ARXCheckType = CheckType::None;

	//Fixed parameters
	CheckType SSBCheckType = CheckType::None; //no need to check SSB for this one
	uint rMax = 7;
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
				balancedBits[i] = !existTrailMidori128(MI,input,output);
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
