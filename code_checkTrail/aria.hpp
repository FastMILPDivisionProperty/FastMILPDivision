#ifndef H_ARIA
#define H_ARIA

#include <vector>
#include <cstdint>
#include <fstream>
#include <string>

#include "aux_function.hpp"
#include "customCallback.hpp"
#include "divtrail.hpp"

std::vector<std::vector<std::pair<std::vector<uint8_t>, std::vector<uint8_t>>>> allInOutPairsARIA(unsigned pos_in, unsigned pos_out, unsigned R);

GRBModel getModelARIA(modelInfo & MI);
/*
	Get the model for ARIA over #rMax round
	#lastLin defines if the linear layer is applied on the last round
	#startingRound defines the index of the first round
	#firstSbox defines if the first Sbox layer is applied

	The #sboxModel parameter defines how the Sbox is modelized :
	- "Simple" use the simplified constraint with PWL

	The #linModel parameter defines how the linear layer is modelized :
	- "CX" modelize the linear layer using the classical copy+xor technique
	- "Simple" modelize the linear layer with the simplified constraint w(x) = w(y)
*/

modelData getVariablesARIA(GRBModel & m,
						   modelInfo const & MI);
/*
	Read the model and return the corresponding blocks of variables
	- #MI defines various parameters for how the trails are checked (see aux_function.hpp)
*/

bool existTrailARIA(GRBModel & model,
				   modelData const & MD,
				   std::vector<uint8_t> const & input,
				   std::vector<uint8_t> const & output,
				   bool const useCallback,
				   unsigned int const maxNumberConstr);
/*
	Return true if there is a trail from #input to #output using #model for ARIA
	- #MD contains the variables to be checked
	- if #useCallback is true then the solving using CB + LC, otherwise repeated solving
	- #maxNumberConstr is the max number of contraints added for each solution found
*/

bool existTrailARIA(modelInfo & MI,
                   std::vector<uint8_t> const & input,
                   std::vector<uint8_t> const & output);
/*
	Return true if there is a trail from #input to #output for ARIA
	#MI should contain the necessary information to generate the model and parametrize the solving
*/

void genSboxFileARIA();
//Generate the tables and ineq files for ARIA
void createMatrixFileARIA();
void createMatrixQMAria();

void checkNormalDistinguisherARIA();
/*
	Search for a normal (no linear combination) distinguisher over 4 rounds of ARIA (new)
*/

void testSearchARIA();
void searchDistinguisher4RoundsARIA();
void searchDistinguisher5RoundsARIA();
#endif
