#ifndef H_AES
#define H_AES

#include <vector>
#include <cstdint>
#include <fstream>
#include <string>
#include <chrono> 
#include <tuple> 

#include "aux_function.hpp"
#include "customCallback.hpp"
#include "divtrail.hpp"

std::vector<std::vector<std::pair<std::vector<uint8_t>, std::vector<uint8_t>>>> allInOutPairsAES(unsigned pos_in, unsigned pos_out);

GRBModel getModelAES(modelInfo & MI);
/*
	Get the model for AES over #rMax round
	#lastLin defines if the linear layer is applied on the last round
	#firstSbox defines if the first Sbox layer is applied

	The #sboxModel parameter defines how the Sbox is modelized :
	- "QM" use the Quin-McCluskey algorithm as in Abdelkhalek,Sasaki,Todo,Tolba,Youssef
	- "Simple" use the simplified constraint with PWL

	The #linModel parameter defines how the linear layer is modelized :
	- "CX" modelize the linear layer using the classical copy+xor technique
	- "Simple" modelize the linear layer with the simplified constraint w(x) = w(y)
*/

modelData getVariablesAES(GRBModel & m,
						  modelInfo const & MI);
/*
	Read the model and return the corresponding blocks of variables
	- #MI defines various parameters for how the trails are checked (see aux_function.hpp)
*/

bool existTrailAES(GRBModel & model,
				   modelData const & MD,
				   std::vector<uint8_t> const & input,
				   std::vector<uint8_t> const & output,
				   bool const useCallback,
				   unsigned int const maxNumberConstr);
/*
	Return true if there is a trail from #input to #output using #model for AES
	- #MD contains the variables to be checked
	- if #useCallback is true then the solving using CB + LC, otherwise repeated solving
	- #maxNumberConstr is the max number of contraints added for each solution found
*/

bool existTrailAES(modelInfo & MI,
                   std::vector<uint8_t> const & input,
                   std::vector<uint8_t> const & output);
/*
	Return true if there is a trail from #input to #output for AES
	#MI should contain the necessary information to generate the model and parametrize the solving
*/

void createMatrixFileAES();
void genSboxFileAES();
//Generate the tables and ineq files for AES

void checkNormalDistinguisherAES();
/*
	Search again for a known distinguisher over AES to compare execution time
*/

void timingSearchAES();

void searchDistinguisher4RoundsAES();
void searchDistinguisher5RoundsAES();
/*
	Search for a distinguisher over 5 rounds AES
	Parameters are fixed and were determined with experiments on a subset of cases
*/

#endif
