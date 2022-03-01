#ifndef H_SKINNY128
#define H_SKINNY128

#include <vector>
#include <cstdint>
#include <fstream>
#include <string>

#include "aux_function.hpp"
#include "customCallback.hpp"
#include "divtrail.hpp"

std::vector<std::vector<std::pair<std::vector<uint8_t>, std::vector<uint8_t>>>> allInOutPairsSKINNY128(unsigned i, unsigned j);

GRBModel getModelSkinny128(modelInfo & MI);
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

modelData getVariablesSkinny128(GRBModel & m,
								modelInfo const & MI);
/*
	Read the model and return the corresponding blocks of variables
	- #MI defines various parameters for how the trails are checked (see aux_function.hpp)
*/

bool existTrailSkinny128(GRBModel & model,
				   modelData const & MD,
				   std::vector<uint8_t> const & input,
				   std::vector<uint8_t> const & output,
				   bool const useCallback,
				   unsigned int const maxNumberConstr);
/*
	Return true if there is a trail from #input to #output using #model for Skinny128
	- #MD contains the variables to be checked
	- if #useCallback is true then the solving using CB + LC, otherwise repeated solving
	- #maxNumberConstr is the max number of contraints added for each solution found
*/

bool existTrailSkinny128(modelInfo & MI,
                   std::vector<uint8_t> const & input,
                   std::vector<uint8_t> const & output);
/*
	Return true if there is a trail from #input to #output for Skinny128
	#MI should contain the necessary information to generate the model and parametrize the solving
*/

void createMatrixFileSkinny128();
void genSboxFileSkinny128();
//Generate the tables and ineq files for Skinny 128

void checkNormalDistinguisherSkinny128();
	/*
		Search for a normal (no linear combination) distinguisher over 10 rounds of Skinny128
	*/

void timingTestSearchNewSkinny128();
void searchDistinguisherSkinny128();
#endif
