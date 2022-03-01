#ifndef H_CAMELLIA
#define H_CAMELLIA

#include <vector>
#include <cstdint>
#include <fstream>
#include <string>

#include "aux_function.hpp"
#include "customCallback.hpp"
#include "divtrail.hpp"

GRBModel getModelCamellia(modelInfo & MI);
/*
	Get the model for Camellia over #rMax round
	#startingRound defines the index of the first round

	The #sboxModel parameter defines how the Sbox is modelized :
	- "QM" use the Quin-McCluskey algorithm as in Abdelkhalek,Sasaki,Todo,Tolba,Youssef
	- "Simple" use the simplified constraint with PWL

	The #linModel parameter defines how the linear layer is modelized :
	- "CX" modelize the linear layer using the classical copy+xor technique
	- "Simple" modelize the linear layer with the simplified constraint w(x) = w(y)
*/

modelData getVariablesCamellia(GRBModel & m,
                           modelInfo const & MI);
     /*
    Read the model and return the corresponding blocks of variables
    - #MI defines various parameters for how the trails are checked (see aux_function.hpp)
    */

bool existTrailCamellia(GRBModel & m,
                   modelData const & MD,
                   std::vector<uint8_t> const & input,
                   std::vector<uint8_t> const & output,
                   bool const useCallback,
                   uint const maxNumberConstr);
/*
    Return true if there is a trail from #input to #output using #model for Camellia
    - #MD contains the variables to be checked
    - if #useCallback is true then the solving using CB + LC, otherwise repeated solving
    - #maxNumberConstr is the max number of contraints added for each solution found
*/

bool existTrailCamellia(modelInfo & MI,
                   std::vector<uint8_t> const & input,
                   std::vector<uint8_t> const & output);
/*
    Return true if there is a trail from #input to #output for Camellia
    #MI should contain the necessary information to generate the model and parametrize the solving
*/

void createMatrixFileCamellia();
//Create the matrix file ./checkData/matrixCamellia.bin

void genSboxFileCamellia();
//Generate the tables and ineq files for Camellia

void checkNormalDistinguisherCamellia();
/*
	Search for a normal (no linear combination) distinguisher over 7 rounds of Camellia
*/

void timingCamellia8r();

void searchCamellia7r();
void searchCamellia8r();

#endif