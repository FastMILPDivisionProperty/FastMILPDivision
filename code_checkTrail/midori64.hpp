#ifndef H_MIDORI64
#define H_MIDORI64

#include <vector>
#include <cstdint>
#include <fstream>
#include <string>

#include "aux_function.hpp"
#include "customCallback.hpp"
#include "divtrail.hpp"

std::vector<std::vector<std::pair<std::vector<uint8_t>, std::vector<uint8_t>>>> allInOutPairsMidori64(unsigned pos_in, unsigned pos_out, unsigned R);

GRBModel getModelMidori64(modelInfo  & MI);
/*
	Get the model for Midori64 over #rMax round
	#lastLin defines if the linear layer is applied on the last round
	#smartLin only has effect if linModel=="Simple" and defines if MC is seen as independent Lboxes or as a single matrix

	The #sboxModel parameter defines how the Sbox is modelized :
	- "Hull" use the Convex Hull technique
	- "Simple" use the simplified constraint with PWL

	The #linModel parameter defines how the linear layer is modelized :
	- "Lbox" modelize the linear layer as the parallel application of the same linear Sbox (thanks to the specific form of the matrix)
	- "ZR" modelize the linear layer using the technique from Zhang and Rijmen
	- "CX" modelize the linear layer using the classical copy+xor technique
	- "Simple" modelize the linear layer with the simplified constraint w(x) = w(y)
*/

modelData getVariablesMidori64(GRBModel & m,
							   modelInfo const & MI);
/*
	Read the model and return the corresponding blocks of variables
	- #MI defines various parameters for how the trails are checked (see aux_function.hpp)
*/

bool existTrailMidori64(GRBModel & model,
				   modelData const & MD,
				   std::vector<uint8_t> const & input,
				   std::vector<uint8_t> const & output,
				   bool const useCallback,
				   unsigned int const maxNumberConstr);
/*
	Return true if there is a trail from #input to #output using #model for Midori64
	- #MD contains the variables to be checked
	- if #useCallback is true then the solving using CB + LC, otherwise repeated solving
	- #maxNumberConstr is the max number of contraints added for each solution found
*/

bool existTrailMidori64(modelInfo & MI,
                   std::vector<uint8_t> const & input,
                   std::vector<uint8_t> const & output);
/*
	Return true if there is a trail from #input to #output for Midori64
	#MI should contain the necessary information to generate the model and parametrize the solving
*/

void createMatrixFileMidori64();

void checkNormalDistinguisherMidori64();
	/*
		Search for a normal (no linear combination) distinguisher over 7 rounds of Midori64 
	*/
#endif
