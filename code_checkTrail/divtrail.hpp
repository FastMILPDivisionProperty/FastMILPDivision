#ifndef H_DIVTRAIL
#define H_DIVTRAIL
#include <array>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdint>
#include <unordered_set>
#include <utility>
#include <map>
#include <algorithm>

#include <boost/dynamic_bitset.hpp>

#include "aux_function.hpp"

std::vector<boost::dynamic_bitset<>> 
truthTableFromSbox(std::vector<uint32_t> const & S,
				   uint const nIn,
				   uint const nOut);
/*
	Return the truth tables for S over #nIn bits input and #nOut bits output
*/

std::vector<boost::dynamic_bitset<>> 
precomputeMask(uint const nbBits);
//Compute the masks used to help with the Moebius transform computation for a function over #nbBits

boost::dynamic_bitset<> 
computeANF(boost::dynamic_bitset<> const & TT, 
		   std::vector<boost::dynamic_bitset<>> const & masks,
		   uint const nbBits);
/*
	From the truth table #TT of a function over #nbBits, compute the corresponding ANF using Moebius transform
	#masks is an array of precomputed masks to help the computation
	This is MUCH faster than not precomputing the masks
*/

void printMonomial(unsigned int u, unsigned int const nbBits);
//Print the monomial corresponding to x^u with #nbBits variables

void printANF(boost::dynamic_bitset<> const & anf,
			  unsigned int const nbBits);
//Print the polynomial corresponding to the ANF with #nbBits variables

bool greatervec(uint64_t const a, uint64_t const b);
/*
Return true if a[i] >= b[i] for all i
	   false otherwise
*/

std::vector<uint32_t> reduceDiv(std::unordered_set<uint32_t> const & s);
std::vector<uint32_t> reduceDiv(std::vector<uint32_t> const & s);
//Reduce the division property set s to remove all redundant vectors
//Also convert it to a vector

std::vector<std::vector<uint32_t>> 
divTrailTable(std::vector<std::vector<boost::dynamic_bitset<>>> const & allTT,
  			  uint const nIn,
  			  uint const nOut,
  			  bool const reduceTable=false);
std::vector<std::vector<uint32_t>> 
divTrailTablev2(std::vector<std::vector<boost::dynamic_bitset<>>> const & allTT,
				uint const nIn,
				uint const nOut,
				bool const reduceTable=false);
/*
	Return the division property table of the functions represented by all the truth tables in #allTT
	#nIn is the input size (in bits)
	#nOut is the output size (in bits)
	#reduceTable defines if we remove redundant vectors
				 We get less inequalities if we don't

	v2 does not use a set to maintain unicity and instead directly insert in a sorted vector (without duplicates).
	Better memory efficiency but not sure about time complexity due to reallocations
*/

void uniqueInsertInSortedVector(std::vector<uint32_t> & T,
								uint32_t const x);

void printDivTable(std::vector<std::vector<uint32_t>> const & D,
				   uint const nIn,
				   uint const nOut);

std::vector<std::pair<uint32_t, uint32_t>> 
forbidDiv(std::vector<std::vector<uint32_t>> const & div, 
		  unsigned n);

void saveDivTableToFile(std::vector<std::vector<uint32_t>> const & T,
						std::string const & filename);
//Save a division table #T to the file #filename

std::vector<std::vector<uint32_t>> readDivTableFromFile(std::string const & filename);
//Read the division table #T from the file #filename

std::vector<std::vector<int>> 
convertToIneq(std::vector<std::pair<uint32_t, uint32_t>> const & T,
			  uint const nbBit);
/*
	Convert the output #T over #nbBit bits from the forbidDiv function to inequalities
		 one inequality is a list of integer
 		 {a0, ..., an, c, s} representing
 		 a0*x[0] + ... + an*x[n] + c op 0
 		 whre op is == if s = 0 and >= if s = 1
 		 In this function, we will only have s=1 (so >=)
*/

void printIneq(std::vector<int> const & ineq);
void printIneqSageFormat(std::vector<int> const & ineq);

void saveIneqToFile(std::vector<std::vector<int>> const & listIneq,
					std::string const & filename);
//Save a list of inequations to the file #filename

std::vector<std::vector<int>> readIneqFromFile(std::string const & filename);
//Read a list of inequations from the file #filename

void genAndSaveTableIneq(std::vector<std::vector<uint32_t>> const & allS,
						 uint const nIn,
						 uint const nOut,
						 std::string const & namePrefix,
						 bool const printSageIneq);
/*
	Generate the table and inequations for the family of functions given by #allS
	Each function must be #nIn bits input and #nOut bits output
	#If printSageIneq is true, also print the ineq in sage format
	File will be named #namePrefix+"_table.bin" and #namePrefix+"_ineq.bin"
*/

void readASCIITableAndConvert(std::string const & filename,
							  std::string const & namePrefix,
							  unsigned int nIn,
							  unsigned int nOut);
//Read a table saved in ASCII format from the DivTable class, convert it into a table binary file and generate the file for inequalities

#endif