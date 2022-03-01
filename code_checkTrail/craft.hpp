#ifndef H_CRAFT
#define H_CRAFT

#include <vector>
#include <cstdint>
#include <fstream>
#include <string>

#include "aux_function.hpp"
#include "customCallback.hpp"
#include "divtrail.hpp"

std::vector<std::vector<std::pair<std::vector<uint8_t>, std::vector<uint8_t>>>> allInOutPairsCraft(unsigned pos_in, unsigned pos_out, unsigned R); // 0 <= i < 4  - 0 <= j < 4

GRBModel getModelCRAFT(modelInfo & MI);
/*
	Get the model for CRAFT over #rMax round
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

modelData getVariablesCRAFT(GRBModel & m,
							modelInfo const & MI);
/*
	Read the model and return the corresponding blocks of variables
	- #MI defines various parameters for how the trails are checked (see aux_function.hpp)
*/

bool existTrailCRAFT(GRBModel & model,
				   modelData const & MD,
				   std::vector<uint8_t> const & input,
				   std::vector<uint8_t> const & output,
				   bool const useCallback,
				   unsigned int const maxNumberConstr);
/*
	Return true if there is a trail from #input to #output using #model for CRAFT
	- #MD contains the variables to be checked
	- if #useCallback is true then the solving using CB + LC, otherwise repeated solving
	- #maxNumberConstr is the max number of contraints added for each solution found
*/

bool existTrailCRAFT(modelInfo & MI,
                   std::vector<uint8_t> const & input,
                   std::vector<uint8_t> const & output);
/*
	Return true if there is a trail from #input to #output for CRAFT
	#MI should contain the necessary information to generate the model and parametrize the solving
*/

void createMatrixFileCRAFT();

void checkNormalDistinguisherCRAFT();
	/*
		Search for a normal (no linear combination) distinguisher over 13 rounds of CRAFT
	*/

void timingSearchCRAFT();
void searchDistinguisher14RoundsCRAFT();
void searchDistinguisher13RoundsCRAFT();
#endif
