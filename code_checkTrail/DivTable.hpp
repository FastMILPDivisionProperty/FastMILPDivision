#ifndef DIVTABLE_HPP
#define DIVTABLE_HPP

#include <cstdint>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>
#include <set>
#include <map>

#include "Polynome.hpp"

std::vector<std::vector<uint16_t>> invDivtable(std::vector<std::vector<uint16_t>> const &);
std::vector<std::vector<uint16_t>> getDivisionTable(std::vector<Polynome> const &anf);
std::vector<std::vector<uint16_t>> getDivisionTable(std::vector<std::vector<Polynome>> const &anf);
std::vector<std::vector<uint16_t>> getDivisionTableNibble(std::vector<std::vector<uint16_t>> const & div);
uint16_t convertToNibbles(uint16_t b);

class DivTable {
public:
    DivTable(std::vector<Polynome> const & anf) : table (getDivisionTable(anf)) {};
    DivTable(std::vector<std::vector<uint16_t>> t) : table (std::move(t)) {};
    DivTable(std::ifstream & input_file) {
      std::string temp;

      while (std::getline(input_file, temp)) {
        std::istringstream buffer(temp);
        std::vector<uint16_t> line((std::istream_iterator<uint16_t>(buffer)), std::istream_iterator<uint16_t>());
        table.emplace_back(move(line));
      }
      for (auto & v : table) {
        sort(v.begin(), v.end(), [](auto const & a, auto const & b){
          auto const & na = __builtin_popcount(a);
          auto const & nb = __builtin_popcount(b);
          return na < nb || (na == nb && a < b);
        });
      }
    }

    unsigned size() const {return table.size();};

    template<typename T>
    std::vector<uint16_t> const & transitionsFrom(T i) const {return table[i];};

    void save(std::ofstream & output_file) const {
      for (auto const & v : table) {
        for (auto x : v) output_file << x << " ";
        output_file << "\n";
      }
    };

private:
    std::vector<std::vector<uint16_t>> table;
};

std::vector<std::pair<std::vector<uint16_t>, unsigned>> getTransitionsOut(std::vector<Polynome> const &anf);
std::vector<std::pair<std::vector<uint16_t>, unsigned>> getTransitionsIn(std::vector<Polynome> const &anf);

void saveTransitions(std::ofstream & output_file, std::vector<std::pair<std::vector<uint16_t>, unsigned>> const & vt);
std::vector<std::pair<std::vector<uint16_t>, unsigned>> loadTransitions(std::ifstream & input_file);

std::vector<std::pair<std::vector<uint16_t>, unsigned>> loadTransitionsIN(std::vector<Polynome> const & anf);
std::vector<std::pair<std::vector<uint16_t>, unsigned>> loadTransitionsOUT(std::vector<Polynome> const & anf);

uint16_t convertToNibbles(uint16_t b);
uint64_t convertToNibbles(uint64_t b);

void printDivisionTable(DivTable const & divtable);


#endif
