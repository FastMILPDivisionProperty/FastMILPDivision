#include "divtrail.hpp"

using namespace std;
using namespace boost;

vector<dynamic_bitset<>> 
truthTableFromSbox(vector<uint32_t> const & S,
				   uint const nIn,
				   uint const nOut){
	/*
		Return the truth tables for S over #nIn bits input and #nOut bits output
	*/

	uint64_t bound = (1ULL << nIn);
	vector<dynamic_bitset<>> TT(nOut, dynamic_bitset<>(bound,0));
	for(uint64_t x = 0; x < bound; x++){
		uint32_t y = S[x];
		for(uint i = 0; i < nOut; i++){
			TT[i][x] = y & 1;
			y >>= 1;
		}
	}

	return TT;
}

vector<dynamic_bitset<>> precomputeMask(uint const nbBits){
	/*
	Compute the masks used to help with the Moebius transform computation for a function over #nbBits (input)
	*/
	vector<dynamic_bitset<>> masks(nbBits);
	uint nbBlock = 1;

	for(uint i = nbBits; i > 0; i--){
		string s;
		for(uint j = 0; j < nbBlock; j++)
			s = string(1ULL << (i-1), '0') + string(1ULL << (i-1), '1') + s;
		masks[i-1] = dynamic_bitset<>(s);
		nbBlock *= 2;
	}

	return masks;
}

dynamic_bitset<> computeANF(dynamic_bitset<> const & TT, 
							vector<dynamic_bitset<>> const & masks,
							uint const nbBits){
	/*
	From the truth table #TT of a function over #nbBits, compute the corresponding ANF using Moebius transform
	#masks is an array of precomputed masks to help the computation
	This is MUCH faster than not precomputing the masks
	*/

	dynamic_bitset<> anf(TT);
	for(uint i = nbBits; i > 0; i--){
		anf ^= ((anf & masks[i-1]) << (1ULL << (i-1)));
	}
	
	return anf;
}

void printMonomial(uint u, uint const nbBits){
	/*
	Print the monomial corresponding to x^u with #nbBits variables
	*/
	if(u == 0) cout << 1;
	else{
		bool firstVar = true;
		for(uint i = 0; i < nbBits; i++){
			string varname = "x" + to_string(i);
			if(u & 1){
				if(!firstVar) cout << "*" << varname;
				else{
					firstVar = false;
					cout << varname;
				}
			}
			u >>= 1;
		}
	}
}

void printANF(dynamic_bitset<> const & anf,
			  uint const nbBits){
	/*
	Print the polynomial corresponding to the ANF with #nbBits variables
	*/
	bool firstMon = true;
	for(uint64_t i = 0; i < (1ULL << nbBits); i++){
		if(anf[i]){
			if(!firstMon) cout << " + ";
			else firstMon = false;
			printMonomial(i,nbBits);
		}
	}
}

bool greatervec(uint64_t const a, uint64_t const b){
	/*
	Return true if a[i] >= b[i] for all i
		   false otherwise
	*/
	if(__builtin_popcountll(a) < __builtin_popcountll(b))
		return false;

	//At this point we know that weight(a) >= weight(b)
	//Just need to check if a[i] = 1 whenever b[i] = 1
	//This is the same as b == a & b
	return b == (a & b);
}

vector<uint32_t> reduceDiv(unordered_set<uint32_t> const & s){
	//Reduce the division property set s to remove all redundant vectors
	//Also convert it to a vector
	vector<uint32_t> t;
	for(auto const & v : s){
		bool great = false;
		for(auto const & v1 : s){
			if(v != v1 && greatervec(v,v1)){
				great = true;
				break;
			}
		}
		if(!great)
			t.emplace_back(v);
	}
	return t;
}

vector<uint32_t> reduceDiv(vector<uint32_t> const & s){
	//Reduce the division property set s to remove all redundant vectors
	//Also convert it to a vector
	vector<uint32_t> t;
	for(auto const & v : s){
		bool great = false;
		for(auto const & v1 : s){
			if(v != v1 && greatervec(v,v1)){
				great = true;
				break;
			}
		}
		if(!great)
			t.emplace_back(v);
	}
	return t;
}

vector<vector<uint32_t>> divTrailTable(vector<vector<dynamic_bitset<>>> const & allTT,
									   uint const nIn,
									   uint const nOut,
									   bool const reduceTable){
/*
	Return the division property table of the functions represented by all the truth tables in #allTT
	Each table entry will be sorted
	#nIn is the input size (in bits)
	#nOut is the output size (in bits)
	#reduceTable defines if we remove redundant vectors
				 We get less inequalities if we don't
*/

	uint64_t boundIn = (1ULL << nIn);
	uint64_t boundOut = (1ULL << nOut);

	//Use an unordered_set to guarantee unicity and convert it to a vector at the end
	vector<unordered_set<uint32_t>> T(boundIn);

	//Initializator for dynamic_bitset
	dynamic_bitset<> init1(string(boundIn, '1'));

	//Mask precomputation for the ANF computation
	auto const maskANF = precomputeMask(nIn);

	//Precompute all monomials of degree <= nIn
	vector<vector<uint32_t>> mon(nIn+1);
	for(uint64_t m = 0; m < boundIn; m++)
		mon[__builtin_popcountll(m)].emplace_back(m);
	
	//Compute the table
	for(auto const & TT : allTT){
		for(uint64_t u = 0; u < boundOut; u++){
			//Compute anf of f^u
			uint tmpu = u;
			dynamic_bitset<> fu(init1);
			for(uint i = 0; i < nOut; i++){
				if(tmpu & 1) fu &= TT[i];
				tmpu >>= 1;
			}
			auto anfyv = computeANF(fu, maskANF, nIn);

			//Populate the table with the monomials from f^u
			for(uint64_t m = 0; m < boundIn; m++){
				if(anfyv[m] == 1)
					T[m].emplace(u);
			}
		}
	}

	vector<vector<uint32_t>> finalT(boundIn);
	//Last transition (from 11...1) is not modified by the next part so convert immediately
	finalT[boundIn-1] = vector<uint32_t>(T[boundIn-1].begin(), T[boundIn-1].end());
	sort(finalT[boundIn-1].begin(), finalT[boundIn-1].end());

	for(int d = nIn-1; d >= 0; d--){
		for(auto const & m : mon[d]){
			auto & T_m = T[m];
			for(auto const & m2 : mon[d+1]){
				if(greatervec(uint64_t(m2), uint64_t(m))){
					for(auto const & v : T[m2])
						T_m.emplace(v);
				}
			}
			//Reduce if needed and convert to vector
			if(reduceTable)
				finalT[m] = reduceDiv(T[m]);
			else //Just convert to vector
				finalT[m] = vector<uint32_t>(T[m].begin(), T[m].end());
			//Sort the entry
			sort(finalT[m].begin(), finalT[m].end());
		}
	}

	return finalT;
}

void uniqueInsertInSortedVector(vector<uint32_t> & T,
								uint32_t const x){
	auto it = lower_bound(T.begin(), T.end(), x);
	//If it = T.end() then element does not exists in T but avoid dereferencing it
	if(it == T.end() || *it != x)
		T.insert(it,x);
}

vector<vector<uint32_t>> divTrailTablev2(vector<vector<dynamic_bitset<>>> const & allTT,
									   uint const nIn,
									   uint const nOut,
									   bool const reduceTable){
/*
	Return the division property table of the functions represented by all the truth tables in #allTT
	Each table entry will be sorted
	#nIn is the input size (in bits)
	#nOut is the output size (in bits)
	#reduceTable defines if we remove redundant vectors
				 We get less inequalities if we don't

	this version does not use a set to maintain unicity and instead directly insert in a sorted vector (without duplicates)
*/

	uint64_t boundIn = (1ULL << nIn);
	uint64_t boundOut = (1ULL << nOut);

	vector<vector<uint32_t>> finalT(boundIn);

	//Initializator for dynamic_bitset
	dynamic_bitset<> init1(string(boundIn, '1'));

	//Mask precomputation for the ANF computation
	auto const maskANF = precomputeMask(nIn);

	//Precompute all monomials of degree <= nIn
	vector<vector<uint32_t>> mon(nIn+1);
	for(uint64_t m = 0; m < boundIn; m++)
		mon[__builtin_popcountll(m)].emplace_back(m);
	
	//Compute the table
	for(auto const & TT : allTT){
		for(uint64_t u = 0; u < boundOut; u++){
			//Compute anf of f^u
			uint tmpu = u;
			dynamic_bitset<> fu(init1);
			for(uint i = 0; i < nOut; i++){
				if(tmpu & 1) fu &= TT[i];
				tmpu >>= 1;
			}
			auto anfyv = computeANF(fu, maskANF, nIn);

			//Populate the table with the monomials from f^u
			for(uint64_t m = 0; m < boundIn; m++){
				if(anfyv[m] == 1)
					uniqueInsertInSortedVector(finalT[m],u);
			}
		}
	}

	for(int d = nIn-1; d >= 0; d--){
		for(auto const & m : mon[d]){
			auto & T_m = finalT[m];
			for(auto const & m2 : mon[d+1]){
				if(greatervec(uint64_t(m2), uint64_t(m))){
					for(auto const & v : finalT[m2])
						uniqueInsertInSortedVector(T_m,v);
				}
			}
			//Reduce if needed
			if(reduceTable)
				finalT[m] = reduceDiv(finalT[m]);
		}
	}

	return finalT;
}

void printDivTable(vector<vector<uint32_t>> const & D,
				   uint const nIn,
				   uint const nOut){

	for(uint64_t u = 0; u < D.size(); u++){
		uint64_t tmpu = u;
		for(uint i = 0; i < nIn; i++){
			if(tmpu & 1) cout << '1';
			else cout << '0';
			tmpu >>= 1;
		}

		cout << " : ";
		for(auto const & v : D[u]){
			uint64_t tmpv = v;
			for(uint i = 0; i < nOut; i++){
				if(tmpv & 1) cout << '1';
				else cout << '0';
				tmpv >>= 1;
			}
			cout << " ";
		}
		cout << endl;
	}
}

vector<pair<uint32_t, uint32_t>> 
forbidDiv(vector<vector<uint32_t>> const & div, 
		  unsigned n){
  vector<vector<uint32_t>> v (1u << n);
  #pragma omp parallel for schedule(dynamic)
  for (uint32_t x = 0; x < 1u << n; ++x) {
    auto const & dx = div[x];
    auto & vx = v[x];
    for (uint32_t y = 1u << n; y-- != 0; ) {
      if (none_of(dx.begin(), dx.end(), [y](auto z){return (y & z) == z;})){ // x --> y is impossible
        if (none_of(vx.begin(), vx.end(), [y](auto z){return (y & z) == y;})) vx.emplace_back(y); // y is maximal
      }
    }
  }

  map<uint32_t, vector<uint32_t>> my_map;
  for (uint32_t x = 0; x < 1u << n; ++x) {
    auto const & vx = v[x];
    for (auto y : vx) {
      auto & my = my_map[y]; // x --> y is impossible and y is maximal
      if (none_of(my.begin(), my.end(), [x](auto z){return (x & z) == z;})) my.emplace_back(x); // x is minimal
    }
  }

  uint32_t const & mask = (1u << n) - 1;

  vector<pair<uint32_t, uint32_t>> res;
  for (auto const & p : my_map) {
    for (uint32_t x : p.second) { // x --> p.first is impossible
      res.emplace_back((mask & ~x) | ((mask & p.first) << n), x);
    }
  }
  return res;
}

void saveDivTableToFile(vector<vector<uint32_t>> const & T,
						string const & filename){
	//Save a division table #T to the file #filename

	auto file = fstream(filename, ios::out | ios::binary);

	//First 32-bit chunk, size of the table
	uint32_t Tsize = T.size();
	file.write((char*)&Tsize, 4);
	for(auto const & v : T){
		//Size of v
		uint32_t vsize = v.size();
		file.write((char*)&vsize, 4);
		//Content of v, exploit that content of a vector is contiguous
		file.write((char*)&v[0], 4*vsize);
	}
	file.close();
}

vector<vector<uint32_t>> readDivTableFromFile(string const & filename){
	//Read the division table #T from the file #filename
	if(!fileExist(filename)){
		cerr << "File " << filename << " does not exists, please generate the file" << endl;
		exit(1);
	}
	auto file = fstream(filename, ios::in | ios::binary);

	//First 32-bit chunk, size of the table
	uint32_t Tsize = 0;
	file.read((char*)&Tsize, 4);

	vector<vector<uint32_t>> T(Tsize);

	for(auto & v : T){
		//Size of v
		uint32_t vsize = 0;
		file.read((char*)&vsize, 4);
		v = vector<uint32_t>(vsize);
		//Content of v, exploit that content of a vector is contiguous
		file.read((char*)&v[0], 4*v.size());
	}
	file.close();

	return T;
}

vector<vector<int>> 
convertToIneq(vector<pair<uint32_t, uint32_t>> const & T,
			  uint const nbBit){
/*
	Convert the output #T over #nbBit bits from the forbidDiv function to inequalities
		 one inequality is a list of integer
 		 {a0, ..., an, c, s} representing
 		 a0*x[0] + ... + an*x[n] + c op 0
 		 whre op is == if s = 0 and >= if s = 1
 		 In this function, we will only have s=1 (so >=)
*/

	//Each ineq is nbBit+2 int, #nbBit for each term + constant + op
	vector<vector<int>> listIneq(T.size(), vector<int>(nbBit+2,0)); 
	for(uint i = 0; i < T.size(); i++){
		auto const & uv = T[i];
		uint32_t u = uv.first;
		uint32_t v = uv.second;

		//inequality will be ... >= 1, so to put in canonical form >= 0 the constant is -1
		int c = -1; 

		auto & ineq = listIneq[i];
		for(uint b = 0; b < nbBit; b++){
			if(!(u & 1)){
				ineq[b] = 1; //Term is x[b]
				if(v & 1){ //Term is (1 - x[b])
					ineq[b] = -1;
					c++;
				}
			}

			u >>= 1;
			v >>= 1;
		}

		ineq[nbBit] = c; //Constant term
		ineq[nbBit+1] = 1; //op is >=
	}

	return listIneq;
}

void printIneq(vector<int> const & ineq){
	uint nbBit = ineq.size()-2;
	bool firstTerm = true;
	for(uint i = 0; i < nbBit; i++){
		if(ineq[i] != 0){
			if(!firstTerm) cout << " + ";
			else firstTerm = false;
			cout << ineq[i] << "*x" << i;
		}
	}

	cout << " + " << ineq[nbBit];
	if(ineq[nbBit+1] == 0)
		cout << " == 0" << endl;
	else
		cout << " >= 0" << endl;
}

void printIneqSageFormat(vector<int> const & ineq){
	uint nbBit = ineq.size()-2;
	bool firstTerm = true;
	cout << "[";
	for(uint i = 0; i < nbBit; i++){
		if(!firstTerm) cout << ", ";
		else firstTerm = false;
		cout << ineq[i];
	}

	cout << ", " << ineq[nbBit];
	if(ineq[nbBit+1] == 0)
		cout << ", \"==\"]," << endl;
	else
		cout << ", \">=\"]," << endl;
}

void saveIneqToFile(vector<vector<int>> const & listIneq,
					string const & filename){
	//Save a list of inequations to the file #filename

	auto file = fstream(filename, ios::out | ios::binary);

	//First 32-bit chunk, size of the table
	uint32_t listIneq_size = listIneq.size();
	file.write((char*)&listIneq_size, 4);
	for(auto const & ineq : listIneq){
		//Size of ineq
		uint32_t ineq_size = ineq.size();
		file.write((char*)&ineq_size, 4);
		//Content of ineq, exploit that content of a vector is contiguous
		file.write((char*)&ineq[0], sizeof(int)*ineq_size);
	}
	file.close();
}

vector<vector<int>> readIneqFromFile(string const & filename){
	//Read a list of inequations from the file #filename
	if(!fileExist(filename)){
		cerr << "File " << filename << " does not exists, please generate the file" << endl;
		exit(1);
	}

	auto file = fstream(filename, ios::in | ios::binary);

	//First 32-bit chunk, size of the table
	uint32_t listIneq_size = 0;
	file.read((char*)&listIneq_size, 4);
	vector<vector<int>> listIneq(listIneq_size);

	for(auto & ineq : listIneq){
		//Size of ineq
		uint32_t ineq_size = ineq.size();
		file.read((char*)&ineq_size, 4);
		ineq = vector<int>(ineq_size);
		//Content of ineq, exploit that content of a vector is contiguous
		file.read((char*)&ineq[0], sizeof(int)*ineq_size);
	}
	file.close();

	return listIneq;
}

void genAndSaveTableIneq(vector<vector<uint32_t>> const & allS,
						 uint const nIn,
						 uint const nOut,
						 string const & namePrefix,
						 bool const printSageIneq){
	/*
		Generate the table and inequations for the family of functions given by #allS
		Each function must be #nIn bits input and #nOut bits output
		#If printSageIneq is true, also print the ineq in sage format
		File will be named #namePrefix+"_table.bin" and #namePrefix+"_ineq.bin"
	*/

	vector<vector<dynamic_bitset<>>> allTT(allS.size());
	#pragma omp parallel for schedule(static)
	for(uint i = 0; i < allS.size(); i++)
		allTT[i] = truthTableFromSbox(allS[i],nIn,nOut);
	
	auto D = divTrailTablev2(allTT,nIn,nOut);
	// printDivTable(D,nIn,nOut);
	saveDivTableToFile(D, namePrefix+"_table.bin");
	auto Tuv = forbidDiv(D,nIn);
	auto listIneq = convertToIneq(Tuv,nIn+nOut);
	saveIneqToFile(listIneq, namePrefix+"_ineq.bin");

	if(printSageIneq){
		for(auto const & ineq : listIneq)
			printIneqSageFormat(ineq);
	}
}

void readASCIITableAndConvert(string const & filename,
							  string const & namePrefix,
							  uint const nIn,
							  uint const nOut){
	ifstream infile(filename);
	std::string temp;

	std::vector<std::vector<uint32_t>> table;
	while (std::getline(infile, temp)) {
		std::istringstream buffer(temp);
		std::vector<uint32_t> line((std::istream_iterator<uint32_t>(buffer)), std::istream_iterator<uint32_t>());
		table.emplace_back(std::move(line));
	}
	saveDivTableToFile(table, namePrefix+"_table.bin");
 	auto Tuv = forbidDiv(table,nIn);
 	auto listIneq = convertToIneq(Tuv,nIn+nOut);
 	cout << listIneq.size() << " inequalities" << endl;
 	saveIneqToFile(listIneq, namePrefix+"_ineq.bin");
}