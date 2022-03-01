#ifndef H_CUSTOMCALLBACK
#define H_CUSTOMCALLBACK

#include <vector>
#include <iostream>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <tuple>
#include <chrono>
#include <map> 

#include "gurobi_c++.h"
#include "aux_function.hpp"



class CustomCallback: public GRBCallback{

	public:
		std::vector<BlockVarsTable> const & varsTable; //variables with a table check
		std::vector<BlockVarsIneq> const & varsIneq;   //variables with an ineq check
		std::vector<BlockVarsMatrix> const & varsMatrix; //variables with a matrix check
		uint64_t maxNumberConstr;					   //Maximum number of contraints added in a single callback
		uint64_t ctrGlobalConstr;						 	   //counter for the number of constraints added
		uint64_t ctrSol;
		// uint64_t ctrNonInvMinor;
		// std::map<uint64_t, uint64_t> histoNonInvMinors;
		uint64_t ctrConstrSbox;
		uint64_t ctrConstrLin;
		uint64_t ctrConstrSSB;

		CustomCallback(std::vector<BlockVarsTable> const & xvarsTable,
					   std::vector<BlockVarsIneq> const & xvarsIneq,
					   std::vector<BlockVarsMatrix> const & xvarsMatrix,
					   uint64_t const xmaxNumberConstr);

	protected:
    	void callback();
};

std::tuple<bool,uint64_t,uint64_t> solveWithoutCallback(GRBModel & m,
						  std::vector<BlockVarsTable> const & varsTable,
						  std::vector<BlockVarsIneq> const & varsIneq,
						  std::vector<BlockVarsMatrix> const & xvarsMatrix,
						  uint64_t const maxNumberConstr);

#endif