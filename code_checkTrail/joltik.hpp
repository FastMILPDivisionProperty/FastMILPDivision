#ifndef H_JOLTIK
#define H_JOLTIK

#include <vector>
#include <cstdint>
#include <fstream>
#include <string>

#include "aux_function.hpp"
#include "customCallback.hpp"
#include "divtrail.hpp"

std::vector<std::vector<std::pair<std::vector<uint8_t>, std::vector<uint8_t>>>> allInOutPairsJoltik(unsigned pos_in, unsigned pos_out, unsigned R);

GRBModel getModelJoltik(modelInfo & MI);
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

modelData getVariablesJoltik(GRBModel & m, 
							 modelInfo const & MI);
/*
	Read the model and return the corresponding blocks of variables
	- #rMax defines the number of rounds
	- #MI defines various parameters for how the trails are checked (see aux_function.hpp)
*/

bool existTrailJoltik(GRBModel & model,
				   modelData const & MD,
				   std::vector<uint8_t> const & input,
				   std::vector<uint8_t> const & output,
				   bool const useCallback,
				   unsigned int const maxNumberConstr);
/*
	Return true if there is a trail from #input to #output using #model for Joltik
	- #MD contains the variables to be checked
	- if #useCallback is true then the solving using CB + LC, otherwise repeated solving
	- #maxNumberConstr is the max number of contraints added for each solution found
*/

bool existTrailJoltik(modelInfo & MI,
                   std::vector<uint8_t> const & input,
                   std::vector<uint8_t> const & output);
/*
	Return true if there is a trail from #input to #output for Joltik
	#MI should contain the necessary information to generate the model and parametrize the solving
*/

void createMatrixFileJoltik();
void computeSboxMatrixJoltik();

void checkNormalDistinguisherJoltik();
	/*
		Search for a normal (no linear combination) distinguisher over 5 rounds of Joltik 
	*/
void timingSearchJoltik();
#endif