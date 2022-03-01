#ifndef POLYNOME_HPP
#define POLYNOME_HPP

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <vector>
#include <set>

typedef uint32_t Monome;

class Polynome {
public:
  Polynome() = default;
  Polynome(Monome a) : v_mon(1, a){};
  Polynome(std::vector<uint32_t> v, bool sorted = false) : v_mon(std::move(v)) {
    if (!sorted) std::sort(v_mon.begin(), v_mon.end());
  };

  auto nbTerms() const { return v_mon.size(); };
  std::set<unsigned> getVariables() const;

  //Polynome set(unsigned n, uint8_t val) const;

  auto begin() const { return v_mon.begin(); };
  auto end() const { return v_mon.end(); };

  auto cbegin() const { return v_mon.cbegin(); };
  auto cend() const { return v_mon.cend(); };

  Monome first() const {return v_mon.front();};
  Monome last() const {return v_mon.back();};

  bool contains(Monome u) const {
    return std::binary_search(v_mon.begin(), v_mon.end(), u);
  };

  std::vector<Monome> getTopMonomes() const;

  Polynome &operator+=(Monome a) {
    if (v_mon.empty() || v_mon.back() < a) v_mon.emplace_back(a);
    else *this = *this + a;
    return *this;
  };
  Polynome &operator+=(Polynome const &a) {
    if (!a.isNull()) *this = *this + a;
    return *this;
  };
  Polynome &operator*=(Monome a) {
    *this = *this * a;
    return *this;
  };
  Polynome &operator*=(Polynome const &a) {
    if (a.nbTerms() != 1 || a.first() != 0) *this = *this * a;
    return *this;
  };

  bool isNull() const {return v_mon.empty();};

  friend Polynome operator+(Polynome const &a, Polynome const &b);
  friend Polynome operator+(Polynome const &a, Monome b);
  friend Polynome operator+(Monome a, Polynome const &b) { return b + a; };

  friend Polynome operator*(Polynome const &a, Polynome const &b);
  friend Polynome operator*(Polynome const &a, Monome b) { return b * a; };
  friend Polynome operator*(Monome a, Polynome const &b) {return Polynome(a) * b;};

  friend bool operator==(Polynome const &a, Polynome const &b) {
    return a.v_mon == b.v_mon;
  };
  friend bool operator<(Polynome const &a, Polynome const &b) {
    return a.v_mon < b.v_mon;
  };

  friend std::vector<Polynome> getANF_S(std::vector<uint16_t> const &sbox);

  friend std::ostream &operator<<(std::ostream &, Polynome const &);

private:
  std::vector<uint32_t> v_mon;
  //
  // auto begin() { return v_mon.begin(); };
  // auto end() { return v_mon.end(); };
};

unsigned degree(Monome m);

template<typename T, typename U>
std::vector<uint32_t> merge(T ia, T na, U ib, U nb) {
  if (ia == na) return std::vector<uint32_t> (ib, nb);
  if (ib == nb) return std::vector<uint32_t> (ia, na);
  std::vector<uint32_t> res (std::distance(ia, na) + std::distance(ib, nb));
  auto it0 = res.begin();
  auto it = it0;
  for (;;) {
    if (*ia == *ib) {
      ++ia; ++ib;
      if (ia == na) {it = std::copy(ib, nb, it); break;}
      if (ib == nb) {it = std::copy(ia, na, it); break;}
    }
    else {
      if (*ia < *ib) {
        *it = *ia; ++it;
        if (++ia == na) {it = std::copy(ib, nb, it); break;}
      }
      else {
        *it = *ib; ++it;
        if (++ib == nb) {it = std::copy(ia, na, it); break;}
      }
    }
  }
  res.resize(std::distance(it0, it));
  return res;
};

Polynome composition(Polynome const & p, std::vector<Polynome> const & v); // return p(v0, v1, ...)
std::vector<Polynome> composition(std::vector<Polynome> const & p, std::vector<Polynome> const & v); // return p0(v0, v1, ...), ...
Polynome convertLtoPolynome(uint16_t l); // "1010" -> b + d,  "1111" -> a + b + c + d
std::vector<Polynome> convertLtoPolynome(std::vector<uint16_t> const & l);

std::vector<Polynome> getANF_S(std::vector<uint16_t> const &sbox);
std::vector<Polynome> getANF_Parallel_S(unsigned n, std::vector<uint16_t> const & sbox);

Polynome product(std::vector<Polynome> const &v, uint16_t u);
Polynome product(std::vector<Polynome> const &v, uint32_t u);

void printAsMonome(Monome x);

#endif
