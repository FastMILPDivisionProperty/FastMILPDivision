#ifndef H_MATRIX
#define H_MATRIX

#include <vector>
#include <cstdint>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <utility>

#include "aux_function.hpp"

using pairRowColum = std::vector<std::pair<std::vector<unsigned int>, std::vector<unsigned int>>>;

class Matrix{

	public:
		std::vector<std::vector<uint64_t>> rows;
		unsigned int nrows;
		unsigned int ncols;

		Matrix();
		//Default constructor, 0 x 0 matrix

		Matrix(unsigned int r,
			   unsigned int c);
		//Build a r x c zero matrix

		Matrix(std::string const & filename);
		//Read the Matrix from a file (created using saveToFile)

		std::vector<uint64_t> & get(unsigned int i);
		std::vector<uint64_t>  const & get(unsigned int i) const;
		//Get the i-th row
		//No overload on [] operator for consistency as we cannot use [][] to get a single coefficient

		unsigned int get(unsigned int const i,
				 		 unsigned int const j) const;
		//Get the [i][j] coefficient, NOT a reference to the coefficient so no [][] operator

		void set(unsigned int const i,
				 unsigned int const j,
				 unsigned int const x);
		//Set coefficient [i][j] to 0 if x=0, 1 otherwise

		bool isInvertibleMinor(std::vector<unsigned int> const & r,
							   std::vector<unsigned int> const & c) const;
		//Return true if the minor defined by the indexed in #r (for the rows) and #c (for the columns) is invertible, false otherwise

		std::pair<bool, pairRowColum>
		isInvertibleMinorWithAugment(std::vector<unsigned int> const & r,
					   				 std::vector<unsigned int> const & c) const;

		void print();
		//print the matrix in a binary format

		void printRawHex();
		//Print the content of this.rows in hex
		//Not an accurante representation as there will be bit reordering and too many zeroes (possibly)

		void saveToFile(std::string const & filename);
		//Save the matrix to a binary file



};

#endif