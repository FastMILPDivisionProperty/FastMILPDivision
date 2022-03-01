#include "Matrix.hpp"

using namespace std;
typedef unsigned int uint;

// static clock_t totalClock = 0;
// static uint64_t nbMinor = 0;

Matrix::Matrix() : 
	rows(),
	nrows(0),
	ncols(0)
	{}
	//Default constructor, 0 x 0 matrix

Matrix::Matrix(unsigned int r,
			   unsigned int c) :
			   rows(),
			   nrows(r),
			   ncols(c)
			   {
	//Build a r x c zero matrix
	if(c%64 == 0)
		rows = vector<vector<uint64_t>>(r,vector<uint64_t>(c/64,0));
	else
		rows = vector<vector<uint64_t>>(r,vector<uint64_t>(c/64 + 1,0));
}

Matrix::Matrix(string const & filename){
	//Read the Matrix from a file (created using saveToFile)
	if(!fileExist(filename)){
		cerr << "Matrix file " << filename << " does not exists, please generate the file" << endl;
		exit(1);
	}

	auto file = fstream(filename, ios::in | ios::binary);
	//First 32-bit chunk is the number of rows
	file.read((char*)&nrows, 4);
	//Second 32-bit chunk is the number of columns
	file.read((char*)&ncols, 4);

	if(ncols%64 == 0)
		rows = vector<vector<uint64_t>>(nrows,vector<uint64_t>(ncols/64,0));
	else
		rows = vector<vector<uint64_t>>(nrows,vector<uint64_t>(ncols/64 + 1,0));

	//Now the content of rows
	for(auto & r : rows)
		file.read((char*)&r[0], 8*r.size());
	
	file.close();
}

vector<uint64_t> & Matrix::get(uint i){
	return rows[i]; 
}
vector<uint64_t>  const & Matrix::get(uint i) const{
	return rows[i];
}
	//Get the i-th row
	//No overload on [] operator for consistency as we cannot use [][] to get a single coefficient

uint Matrix::get(uint const i,
		 		 uint const j) const{
	//Get the [i][j] coefficient, NOT a reference to the coefficient so no [][] operator
	if(rows[i][j/64] & (uint64_t(1) << (j%64))) return 1;
	else return 0;
}

void Matrix::set(uint const i,
		 uint const j,
		 uint const x){
	//Set coefficient [i][j] to 0 if x=0, 1 otherwise
	if(x == 0)
		rows[i][j/64] = rows[i][j/64] & (uint64_t(-1) ^ (uint64_t(1) << (j%64)));
	else
		rows[i][j/64] = rows[i][j/64] | (uint64_t(1) << (j%64));
}

bool Matrix::isInvertibleMinor(std::vector<unsigned int> const & r,
					   std::vector<unsigned int> const & c) const{
	//Return true if the minor defined by the indexed in #r (for the rows) and #c (for the columns) is invertible, false otherwise
	//Important : we need the minor to be *invertible*, not just full rank, so we need r.size() == c.size()
	//empty matrix is considered invertible

	if(r.size() != c.size()) return false;
	if(r.size() == 0) return true; //empty matrix 0 x 0

	//Build a temporary matrix with r.size() rows
	//Keep the same number of columns as #this to make it more efficient
	//We will just be smart about the pivots
	Matrix M(r.size(), this->ncols);
	for(uint i = 0; i < r.size(); i++){
		M.get(i) = this->get(r[i]);
	}

	//Noy try to echelonize on each column defined by #c
	//If not possible then the minor is not invertible
	for(uint i = 0; i < c.size(); i++){
		uint col = c[i];

		//Find the pivot row
		//If none, not invertible
		int indexPivot = -1;
		for(uint j = i; j < M.nrows; j++){
			if(M.get(j,col) == 1){
				indexPivot = j;
				break;
			}
		}
		//If no pivot found, not invertible
		if(indexPivot == -1) return false;

		//Found the pivot, swap it to row i
		swap(M.get(i), M.get(indexPivot));

		//Echelonize
		for(uint j = i+1; j < M.nrows; j++){
			if(M.get(j,col) == 1){
				for(uint k = 0; k < M.get(j).size(); k++)
					M.get(j)[k] ^= M.get(i)[k];
			}
		}
	}

	//If we get here, we found a pivot for each column required, so the minor is invertible
	return true;
}

pair<bool, pairRowColum> Matrix::isInvertibleMinorWithAugment(std::vector<unsigned int> const & r,
					   std::vector<unsigned int> const & c) const{
	//Return true if the minor defined by the indexed in #r (for the rows) and #c (for the columns) is invertible, false otherwise
	//Important : we need the minor to be *invertible*, not just full rank, so we need r.size() == c.size()
	//empty matrix is considered invertible

	if(r.size() != c.size()) return make_pair(false, pairRowColum());
	if(r.size() == 0) return make_pair(true, pairRowColum()); //empty matrix 0 x 0

	//Build a temporary matrix with r.size() rows augmented with identity
	Matrix M(r.size(), 2*c.size());
	for(uint i = 0; i < r.size(); i++){
		for(uint j = 0; j < c.size(); j++){
			M.set(i,j,this->get(r[i],c[j]));
		}
		M.set(i, c.size()+i, 1);
	}

	// cout << "initial matrix" << endl;
	// M.print();
	// cout << endl;

	//Now try to echelonize on each column of this first half
	uint nextPivot = 0;
	uint rankDeficiency = 0;
	for(uint i = 0; i < c.size(); i++){

		//Find the pivot row
		//If none, not invertible
		int indexPivot = -1;
		for(uint j = nextPivot; j < M.nrows; j++){
			if(M.get(j,i) == 1){
				indexPivot = j;
				// cout << "for i = " << i << " pivot : ";
				// for(uint k = 0; k < M.ncols; k++)
				// 	cout << M.get(j,k);
				// cout << endl;
				break;
			}
		}
		//If no pivot found, not invertible
		if(indexPivot == -1){
			// cout << "No pivot for " << i << " nextpivot = " << nextPivot << endl;
			// cout << "current :" << endl;
			// M.print();
			rankDeficiency++;
		}
		else{
			//Found the pivot, swap it to row nextPivot
			swap(M.get(nextPivot), M.get(indexPivot));

			//Echelonize
			for(uint j = nextPivot+1; j < M.nrows; j++){
				if(M.get(j,i) == 1){
					for(uint k = 0; k < M.get(j).size(); k++)
						M.get(j)[k] ^= M.get(nextPivot)[k];
				}
			}
			nextPivot++;
			// cout << "After ech on i = " << i << endl;
			// M.print();
			// for(uint j = 0; j < M.ncols; j++){
			// 	if(j == i) cout << "^";
			// 	else cout << " ";
			// }
			// cout << endl;
		}
	}

	// cout << "echelonized matrix" << endl;
	// M.print();
	// cout << endl;

	if(rankDeficiency == 0){
		//Invertible minor
		return make_pair(true, pairRowColum());
	}
	else{
		//Not invertible, the numbre of rows equal to 0 is exactly rankDeficiency, and they're the last rows
		//Extract the linear combinations that sum to 0
		pairRowColum ret;
		for(uint i = M.nrows - rankDeficiency; i < M.nrows; i++){
			vector<uint> setRow;
			vector<uint> setCol;
			Matrix tmp(1, this->ncols);
			//Sum the full rows corresponding to the 0 lin. comb. in the minor
			for(uint j = 0; j < c.size(); j++){
				if(M.get(i, c.size() + j) == 1){
					setRow.emplace_back(r[j]);
					for(uint k = 0; k < tmp.get(0).size(); k++){
						tmp.get(0)[k] ^= this->get(r[j])[k];
					}
				}
			}
			//Check which columns have a 1 coeff
			for(uint j = 0; j < tmp.ncols; j++){
				if(tmp.get(0,j) == 1)
					setCol.emplace_back(j);
			}
			ret.emplace_back(make_pair(move(setRow), move(setCol)));
		}
		return make_pair(false, ret);
	}
}

void Matrix::print(){
	//print the matrix in a binary format
	for(uint i = 0; i < nrows; i++){
		for(uint j = 0; j < ncols; j++)
			cout << get(i,j);
		cout << endl;
	}
}

void Matrix::printRawHex(){
		//Print the content of this.rows in hex
		//Not an accurante representation as there will be bit reordering and too many zeroes (possibly)

	for(auto const & r : rows){
		for(auto const & c : r)
			cout << "0x" << std::setfill('0') << std::setw(16) << hex << c << dec << " ";
		cout << endl;
	}
}

void Matrix::saveToFile(std::string const & filename){
	//Save the matrix to a binary file
	auto file = fstream(filename, ios::out | ios::binary);
	//First 32-bit chunk is the number of rows
	file.write((char*)&nrows, 4);
	//Second 32-bit chunk is the number of columns
	file.write((char*)&ncols, 4);

	//Now write the content of rows
	for(auto const & r : rows)
		file.write((char*)&r[0], 8*r.size());
	
	file.close();
}