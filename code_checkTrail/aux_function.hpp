#ifndef H_AUXFUNCTION
#define H_AUXFUNCTION

#include <fstream>
#include <vector>
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
#include <cstdint>
#include <map>

#include "Matrix.hpp"
#include "gurobi_c++.h"

class Matrix;

enum class ConstrType{
	Sbox=0,
	Lin,
	SSB,
	ARX
};

struct BlockVarsTable{
	std::vector<GRBVar> input;
	std::vector<GRBVar> output;
	std::vector<std::vector<uint32_t>> const & table;
	ConstrType ctype;

	BlockVarsTable(std::vector<GRBVar> const & xinput,
				   std::vector<GRBVar> const & xoutput,
				   std::vector<std::vector<uint32_t>> const & xtable,
				   ConstrType const xctype) :
				   input(xinput),
				   output(xoutput),
				   table(xtable),
				   ctype(xctype)
				   {};
};
/*
	Structure for a block using a table
	#input are the input variables
	#output are the output variables
	#table is a reference to a table representing the valid solutions
	 e.g. v in table[u] iif u -> v is valid, using the binary representation of u and v
*/

struct BlockVarsIneq{
	std::vector<GRBVar> input;
	std::vector<GRBVar> output;
	std::vector<std::vector<int>> const & ineq;
	ConstrType ctype;

	BlockVarsIneq(std::vector<GRBVar> const & xinput,
				  std::vector<GRBVar> const & xoutput,
				  std::vector<std::vector<int>> const & xineq,
				  ConstrType const xctype) :
				  input(xinput),
				  output(xoutput),
				  ineq(xineq),
				  ctype(xctype)
				  {};
};
/*
	Structure for a block using a list of inequalities to be verified
	#input are the input variables
	#output are the output variables
	#ineq is a list of inequalities
	 one inequality is a list of integer
	 {a0, ..., an, b0, ... bn, c, s} representing
	 a0*input[0] + ... + an*input[n] + b0*ouptut[0] + ... + bn*output[n] + c op 0
	 whre op is == if s = 0 and >= if s = 1
*/

struct BlockVarsMatrix{
	std::vector<GRBVar> input;
	std::vector<GRBVar> output;
	Matrix const & mat;

	BlockVarsMatrix(std::vector<GRBVar> const & xinput,
					std::vector<GRBVar> const & xoutput,
					Matrix const & xmat) :
					input(xinput),
					output(xoutput),
					mat(xmat)
					{};

};
/*
	Structure for a block using the minor of the matrix to be verified
	#input are the input variables
	#output are the output variables
	#mat is the matrix
*/

struct modelData{
	unsigned int rMax;
	bool firstSbox;				//true if the sbox layer is applied on the first round
	bool lastLin;
	std::vector<BlockVarsTable> allVarsTable;	 //All blocks of variables checked with a table
    std::vector<BlockVarsIneq> allVarsIneq;		 //All blocks of variables checked with ineqs
    std::vector<BlockVarsMatrix> allVarsMatrix;	 //All blocks of variables checked with matrix
    std::vector<std::vector<std::vector<uint32_t>>> allTable;	 //All table used
    std::vector<std::vector<std::vector<int>>> allIneq; //All ineq used
    std::vector<Matrix> allMatrix; 						 //All matrices used

    modelData() :
    	rMax(0),
    	firstSbox(true),
    	lastLin(true),
    	allVarsTable(),
    	allVarsIneq(),
    	allVarsMatrix(),
    	allTable(),
    	allIneq(),
    	allMatrix()
    	{};

    // modelData(unsigned int const xrMax,
    // 		  bool const xlastLin,
    // 		  std::vector<BlockVarsTable> xallVarsTable,
			 //  std::vector<BlockVarsIneq> xallVarsIneq,
			 //  std::vector<BlockVarsMatrix>
			 //  std::vector<std::vector<std::vector<uint32_t>>> xallTable,
			 //  std::vector<std::vector<std::vector<int>>> xallIneq) :
    // 		  rMax(xrMax),
    // 		  lastLin(xlastLin),
    // 		  allVarsTable(xallVarsTable),
    // 		  allVarsIneq(xallVarsIneq),
    // 		  allTable(xallTable),
    // 		  allIneq(xallIneq)
    // 		  {};
};
/*
	Structure to keep various data about the model
	Especially so that the tables/inequalities do not go out of scope
	This is tied to one specific model (unlike the modelInfo struct which is just information about the model generation
*/

enum class CheckType{
	None=0,
	Table,
	Ineq,
	Matrix
};
std::ostream& operator<<(std::ostream& out, const CheckType value);
/*
	Enum to define how a given block of variable is checked :
	- None does not check anything
	- Table uses a table check
	- Ineq uses ineqautions checks
	- Matrix uses the minor of the MC matrix for checking
*/

struct modelInfo{
	unsigned int rMax;			//Number of rounds
	std::string sboxModel;		//How the sbox is modelized
	std::string linModel;		//How the linear layer is modelized
	std::string arxModel;		//How the modular addition is modelized
	bool firstSbox;				//true if the sbox layer is applied on the first round
	bool lastLin;				//true if the linear layer is applied on the last round
	bool smartLin;				//true if MC needs to be considered as independent lboxes for simple modelization
	unsigned int startingRound; //Index of the first round
	CheckType sboxCheckType;	//How to check Sboxes variables
	CheckType MCCheckType;		//How to check MC variables
	CheckType SSBCheckType;		//How to check SSB variables
	CheckType ARXCheckType; 	//How to check ARX variables
	bool useCallback;			//true if using a callback for solving
	uint maxNumberConstr;		//Maximum number of constraints added on each solution found
	std::vector<unsigned int> roundOrder; //Order in which the rounds are checked
										  //Should contains integer form 0 to rMax-1
	GRBEnv gurobiEnv;					  //Unique Gurobi environment

	modelInfo(unsigned int const xrMax,
			  std::string const & xsboxModel,
			  std::string const & xlinModel,
			  std::string const & xarxModel,
			  bool const xfirstSbox,
			  bool const xlastLin,
			  bool const xsmartLin,
			  unsigned int xstartingRound,
			  CheckType const xsboxCheckType,
			  CheckType const xMCCheckType,
			  CheckType const xSSBCheckType,
			  CheckType const xARXCheckType,
			  bool const xuseCallback,
			  uint const xmaxNumberConstr,
			  std::vector<unsigned int> xroundOrder) :
			  rMax(xrMax),
			  sboxModel(xsboxModel),
			  linModel(xlinModel),
			  arxModel(xarxModel),
			  firstSbox(xfirstSbox),
			  lastLin(xlastLin),
			  smartLin(xsmartLin),
			  startingRound(xstartingRound),
			  sboxCheckType(xsboxCheckType),
			  MCCheckType(xMCCheckType),
			  SSBCheckType(xSSBCheckType),
			  ARXCheckType(xARXCheckType),
			  useCallback(xuseCallback),
			  maxNumberConstr(xmaxNumberConstr),
			  roundOrder(xroundOrder),
			  gurobiEnv()
			  {}

};
/*
	Some info about the model generation and checking
	Not tied to a specific model object
	Note that some of them only apply to some ciphers (e.g. smartLin) but are kept in the same struct for consistency and keeping the same prototype for the existTrail* functions
	- arxModel and ARXCheckType only have effect for HIGHT
	- smartLin only has effect for Skinny and Midori (if linModel="Simple")
	- startingRound only has effect for ARIA
	- SSBCheckType has no effect for ARIA as it does not have super sboxes
*/

bool fileExist(std::string const & filename);

std::string exec(const char* cmd);

uint16_t apply16bitMC16bitState(std::vector<uint16_t> const & M,
								uint16_t const x);

uint16_t apply4bitSbox16bitState(std::vector<uint16_t> const & S, 
							  uint16_t const x);

std::vector<uint> roundOrderAscend(uint const rMax);
	//Return a round order from 0 to rMax-1

std::vector<uint> roundOrderDescend(uint const rMax);
	//Return a round order from rMax-1 to 0

std::vector<uint> roundOrderOutsideIn(uint const rMax);
	/*
	Return a round order 0 rMax-1 1 rMax-2 etc., e.g.
	rMax = 7
	0 6 1 5 2 4 3

	rMax = 8
	0 7 1 6 2 5 3 4
	*/

std::vector<uint> roundOrderInsideOut(uint const rMax);
	/*
	Return a round order inverse of OutsideIn, start in the middle, end at begin/end
	*/

typedef bool (*t_existTrail)(modelInfo & , std::vector<uint8_t> const &, std::vector<uint8_t> const &);
bool checkForDistinguisher(std::vector<std::vector<std::pair<std::vector<uint8_t>, std::vector<uint8_t>>>> const & allSetPairs,
						   t_existTrail existTrail,
						   modelInfo & MI,
						   bool const fastSearch = true);

//Special overload to handle cases with a lot of (large) sets of pairs, e.g. LED
bool checkForDistinguisher(std::vector<std::vector<std::pair<std::vector<uint8_t>, std::vector<uint8_t>>>> const & allSetPairs,
						   t_existTrail existTrail,
						   modelInfo & MI,
						   std::map<std::pair<std::vector<uint8_t>, std::vector<uint8_t>>, bool> & checkedPairs,
						   bool const fastSearch = true);

void printVec(std::vector<uint8_t> const & v,
			  unsigned int const wordSize,
			  bool pythonArrayFormat = false);

#endif