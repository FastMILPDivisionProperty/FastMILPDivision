#include "aux_function.hpp"

using namespace std;

typedef unsigned int uint;

std::ostream& operator<<(std::ostream& out, const CheckType value){
    static std::array<std::string,4> names({"None","Table","Ineq","Matrix"});
    return out << names[static_cast<int>(value)];
}

bool fileExist(string const & filename){
    ifstream f(filename.c_str());
    return f.good();
}

string exec(const char* cmd){
//Exec cmd and grab the stdout
    array<char, 128> buffer;
    string result;
    std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (!feof(pipe.get())) {
        if (fgets(buffer.data(), 128, pipe.get()) != nullptr)
            result += buffer.data();
    }
    return result;
}

uint16_t apply16bitMC16bitState(vector<uint16_t> const & M,
								uint16_t const x){
	uint16_t y = 0;
	for(uint j = 0; j < 16; j++){
		if(__builtin_parity(M[j] & x))
			y |= (1 << j);
	}
	return y;
}

uint16_t apply4bitSbox16bitState(vector<uint16_t> const & S, 
							  uint16_t const x){
	return (S[(x & 0xF)]) | 
		   (S[(x >> 4) & 0xF] << 4) |
		   (S[(x >> 8) & 0xF] << 8) |
		   (S[(x >> 12) & 0xF] << 12);
}

vector<uint> roundOrderAscend(uint const rMax){
	//Return a round order from 0 to rMax-1
	vector<uint> roundOrder(rMax);
	for(uint i = 0; i < rMax; i++)
		roundOrder[i] = i;
	return roundOrder;
}

vector<uint> roundOrderDescend(uint const rMax){
	//Return a round order from rMax-1 to 0
	vector<uint> roundOrder(rMax);
	for(uint i = 0; i < rMax; i++)
		roundOrder[i] = rMax-1-i;
	return roundOrder;
}

vector<uint> roundOrderOutsideIn(uint const rMax){
	/*
	Return a round order 0 rMax-1 1 rMax-2 etc., e.g.
	rMax = 7
	0 6 1 5 2 4 3

	rMax = 8
	0 7 1 6 2 5 3 4
	*/
	vector<uint> roundOrder(rMax);
	for(uint i = 0; i < rMax/2; i++){
		roundOrder[2*i] = i;
		roundOrder[2*i+1] = rMax-1-i;
	}
	if(rMax%2 == 1)
		roundOrder[rMax-1] = rMax/2; //special case if rMax odd
	return roundOrder;
}

vector<uint> roundOrderInsideOut(uint const rMax){
	/*
	Return a round order inverse of OutsideIn, start in the middle, end at begin/end
	*/
	vector<uint> tmp = roundOrderOutsideIn(rMax);
	vector<uint> roundOrder(rMax);
	for(uint i = 0; i < rMax; i++)
		roundOrder[i] = tmp[rMax-1-i];
	return roundOrder;
}

bool checkForDistinguisher(vector<vector<pair<vector<uint8_t>, vector<uint8_t>>>> const & allSetPairs,
						   t_existTrail existTrail,
						   modelInfo & MI,
						   bool const fastSearch){

	//Prepare a map to avoid checking duplicates
	map<pair<vector<uint8_t>, vector<uint8_t>>, bool> checkedPairs;
	auto checkedPairs_end = checkedPairs.end();
	bool ret = false;

	//Each set of pair allows to test for a distinguisher
	uint indexSet = 0;
	uint ctrDist = 0;
	for(auto const & setPair : allSetPairs){
		indexSet++;
		bool isTrail = false;
		uint ctr = 0;
		//First check if an already computed pair is in the table
		for(auto const & pair : setPair){
			if(checkedPairs.find(pair) != checkedPairs_end){
				isTrail = checkedPairs[pair];
				if(isTrail)
					break;
			}
		}
		if(!isTrail){ //Already known pairs have no trails, check everything
			for(auto const & pair : setPair){
				// cout << "Input  : ";
				// for(auto const & x : pair.first)
				// 	cout << uint(x);
				// cout << endl << "Output : ";
				// for(auto const & x : pair.second)
				// 	cout << uint(x);
				// cout << endl;
				ctr++;
				// cout << "Checking pair " << ctr << "/"  << setPair.size() << "      \r" << flush;
				// cout << "input  : "; printVec(pair.first,4,true);
				// cout << "output : "; printVec(pair.second,4,true);
				if(checkedPairs.find(pair) == checkedPairs_end){
					// cout << "MILP solving..." << endl;
					isTrail = existTrail(MI, pair.first, pair.second);
					checkedPairs[pair] = isTrail;
				}
				else{
					// cout << "Already computed" << endl;
					isTrail = checkedPairs[pair];
				}

				if(isTrail)
					break;
			}
		}
		// cout << "Set " << indexSet << "/" << allSetPairs.size() << ", " << ctr << "/" << setPair.size() << " pairs checked, " << endl;
		if(!isTrail){
			cout << "Set " << indexSet << "/" << allSetPairs.size() << ", " << ctr << "/" << setPair.size() << " pairs checked, " << endl;
			// for(auto const & pair : setPair){
			// 	cout << "i : ";
			// 	for(auto const & tmp : pair.first)
			// 		cout << int(tmp);
			// 	cout << endl << "o : ";
			// 	for(auto const & tmp : pair.second)
			// 		cout << int(tmp);
			// 	cout << endl;
			// }
			cout << "We have a distinguisher !!!" << endl;
			ctrDist++;
			ret = true;
			if(fastSearch) return true;
		}
		// else
			// cout << "no distinguisher" << endl;
	}
	cout << "Found " << ctrDist << " distinguishers" << endl;
	// return false;
	return ret;
}

//Special overload to handle cases with a lot of (large) sets of pairs, e.g. LED
bool checkForDistinguisher(vector<vector<pair<vector<uint8_t>, vector<uint8_t>>>> const & allSetPairs,
						   t_existTrail existTrail,
						   modelInfo & MI,
						   map<pair<vector<uint8_t>, vector<uint8_t>>, bool> & checkedPairs,
						   bool const fastSearch){

	//Prepare a map to avoid checking duplicates
	auto checkedPairs_end = checkedPairs.end();
	bool ret = false;

	//Each set of pair allows to test for a distinguisher
	uint indexSet = 0;
	uint ctrDist = 0;
	for(auto const & setPair : allSetPairs){
		indexSet++;
		bool isTrail = false;
		uint ctr = 0;
		//First check if an already computed pair is in the table
		for(auto const & pair : setPair){
			if(checkedPairs.find(pair) != checkedPairs_end){
				isTrail = checkedPairs[pair];
				if(isTrail)
					break;
			}
		}
		if(!isTrail){ //Already known pairs have no trails, check everything
			for(auto const & pair : setPair){
				ctr++;
				if(checkedPairs.find(pair) == checkedPairs_end){
					isTrail = existTrail(MI, pair.first, pair.second);
					checkedPairs[pair] = isTrail;
				}
				else{
					isTrail = checkedPairs[pair];
				}

				if(isTrail)
					break;
			}
		}
		if(!isTrail){
			cout << "Set " << indexSet << "/" << allSetPairs.size() << ", " << ctr << "/" << setPair.size() << " pairs checked, " << endl;
			cout << "We have a distinguisher !!!" << endl;
			ctrDist++;
			ret = true;
			if(fastSearch) return true;
		}
		// else
			// cout << "no distinguisher" << endl;
	}
	// cout << "Found " << ctrDist << " distinguishers" << endl;
	return ret;
}

void printVec(std::vector<uint8_t> const & v,
			  unsigned int const wordSize,
			  bool pythonArrayFormat){
	if(pythonArrayFormat){
		cout << "[";
		for(auto const & tmp : v)
			cout << uint(tmp) << ",";
		cout << "]" << endl;
	}
	else{
		for(uint i = 0; i < v.size(); i++){
			cout << uint(v[i]);
			if(wordSize > 0 && (i+1)%wordSize == 0)
				cout << " ";
		}
		cout << endl;
	}
}