#ifndef H_HIGHT
#define H_HIGHT

#include <vector>
#include <cstdint>
#include <fstream>
#include <string>

#include "aux_function.hpp"
#include "customCallback.hpp"
#include "divtrail.hpp"

std::vector<std::vector<std::pair<std::vector<uint8_t>, std::vector<uint8_t>>>> allInOutPairsHIGHT(unsigned pos_in, unsigned pos_out);

GRBModel getModelHIGHT(modelInfo & MI);
/*
	Get the model for HIGHT over #rMax round

	The #arxModel parameter defines how the mod addition is modelized :
	- "ARX" use the technique from Sun,Wang,Liu,Wang

	The #linModel parameter defines how the linear layer is modelized :
	- "Lbox" modelize the linear layer as the parallel application of the same linear Sbox
	  Less precise than usual as it does not use the reduced table as for a 4-bit sbox
	- "ZR" modelize the linear layer using the technique from Zhang and Rijmen
	- "CX" modelize the linear layer using the classical copy+xor technique
	- "Simple" modelize the linear layer with the simplified constraint w(x) = w(y)
*/

modelData getVariablesHIGHT(GRBModel & m,
							modelInfo const & MI);
/*
	Read the model and return the corresponding blocks of variables
	- #MI defines various parameters for how the trails are checked (see aux_function.hpp)
*/

bool existTrailHIGHT(GRBModel & model,
				   modelData const & MD,
				   std::vector<uint8_t> const & input,
				   std::vector<uint8_t> const & output,
				   bool const useCallback,
				   unsigned int const maxNumberConstr);
/*
	Return true if there is a trail from #input to #output using #model for HIGHT
	- #MD contains the variables to be checked
	- if #useCallback is true then the solving using CB + LC, otherwise repeated solving
	- #maxNumberConstr is the max number of contraints added for each solution found
*/

bool existTrailHIGHT(modelInfo & MI,
                   std::vector<uint8_t> const & input,
                   std::vector<uint8_t> const & output);
/*
	Return true if there is a trail from #input to #output for HIGHT
	#MI should contain the necessary information to generate the model and parametrize the solving
*/

uint8_t rotl8(uint8_t const x, uint c);

uint8_t applyF0(uint8_t const x);

uint8_t applyF1(uint8_t const x);

void genMCFileHIGHT();
//Generate the tables and ineq files for HIGHT
void createMatrixFileHIGHT();

void checkNormalDistinguisherHIGHT();
	/*
		Search for a normal (no linear combination) distinguisher over 4 rounds of HIGHT
	*/
void testSearchHIGHT();
void searchDistinguisherHIGHT20r();
void searchDistinguisherHIGHT21r();
#endif
