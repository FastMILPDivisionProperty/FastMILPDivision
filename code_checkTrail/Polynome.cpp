#include <vector>
#include <algorithm>
#include <map>
#include <execution>

#include "Polynome.hpp"

using namespace std;

unsigned degree(Monome m) {return __builtin_popcount(m);}

Polynome operator+(Polynome const & a, Polynome const & b) {
  Polynome res;
  res.v_mon = merge(a.begin(), a.end(), b.begin(), b.end());
  return res;
}

Polynome operator+(Polynome const & a, Monome b) {
    Polynome res;
    res.v_mon = merge(a.begin(), a.end(), &b, &b + 1);
    return res;
}

uint32_t getMask(uint32_t const a) {
  return uint32_t(1) << (31 -__builtin_clz(a));
}

template<typename T, typename U>
vector<uint32_t> mergeAddMask(T ia, T na, uint32_t const mask, U ib, U nb) {
  if (ia == na) return vector<uint32_t> (ib, nb);
  if (ib == nb) {
    vector<uint32_t> res (distance(ia, na));
    transform(ia, na, res.begin(), [mask](auto const & a){return a^mask;});
    return res;
  }
  vector<uint32_t> res (distance(ia, na) + distance(ib, nb));
  auto it = res.begin();
  auto ia_mask = *ia ^ mask;
  for (;;) {
    if (ia_mask == *ib) {
      ++ia; ++ib;
      if (ia == na) {it = copy(ib, nb, it); break;}
      if (ib == nb) {it = transform(ia, na, it, [mask](auto const & a){return a^mask;}); break;}
      ia_mask = *ia ^ mask;
    }
    else {
      if (ia_mask < *ib) {
        *it = ia_mask; ++ia; ++it;
        if (ia == na) {it = copy(ib, nb, it); break;}
        ia_mask = *ia ^ mask;
      }
      else {
        *it = *ib; ++ib; ++it;
        if (ib == nb) {it = transform(ia, na, it, [mask](auto const & a){return a^mask;}); break;}
      }
    }
  }
  res.resize(distance(res.begin(), it));
  return res;
}

template<typename T, typename U, typename V>
auto mergeXOR(T ia, T na, U ib, U nb, V it) {
  if (ia == na) return copy(ib, nb, it);
  if (ib == nb) return copy(ia, na, it);
  for (;;) {
    if (*ia == *ib) {
      ++ia; ++ib;
      if (ia == na) return copy(ib, nb, it);
      if (ib == nb) return copy(ia, na, it);
    }
    else {
      if (*ia < *ib) {
        *it = *ia; ++ia; ++it;
        if (ia == na) return copy(ib, nb, it);
      }
      else {
        *it = *ib; ++ib; ++it;
        if (ib == nb) return copy(ia, na, it);
      }
    }
  }
}

template<typename T, typename U, typename V>
auto mergeXOR_addMask(T ia, T na, U ib, U nb, uint32_t mask, V it) {
  if (ia == na) return transform(ib, nb, it, [mask](auto const & x){return x | mask;});
  if (ib == nb) return transform(ia, na, it, [mask](auto const & x){return x | mask;});
  for (;;) {
    if (*ia == *ib) {
      ++ia; ++ib;
      if (ia == na) return transform(ib, nb, it, [mask](auto const & x){return x | mask;});
      if (ib == nb) return transform(ia, na, it, [mask](auto const & x){return x | mask;});
    }
    else {
      if (*ia < *ib) {
        *it = *ia | mask; ++ia; ++it;
        if (ia == na) return transform(ib, nb, it, [mask](auto const & x){return x | mask;});
      }
      else {
        *it = *ib | mask; ++ib; ++it;
        if (ib == nb) return transform(ia, na, it, [mask](auto const & x){return x | mask;});
      }
    }
  }
}

template <typename T, typename U>
vector<uint32_t> multiplication(T ab, T ae, U bb, U be) {
  if (ab == ae || bb == be) return vector<uint32_t> ();
  if (*ab == 0) {auto const & p = multiplication(bb, be, next(ab), ae); return merge(bb, be, p.begin(), p.end());}
  if (*bb == 0) {auto const & p = multiplication(ab, ae, next(bb), be); return merge(ab, ae, p.begin(), p.end());}
  auto const mask = getMask(*prev(ae) | *prev(be));
  if (mask < 128) {
      unsigned const n = mask << 1;
      vector<uint8_t> tmp (n, 0);
      for (; ab != ae; ++ab) {for (auto itb = bb; itb != be; ++itb) tmp[*ab | *itb] ^= 1;}
      vector<uint32_t> res (n);
      auto beg = res.begin(); auto it = beg;
      for (unsigned i = 0; i < n; ++i) {if (tmp[i] != 0) {*it = i; ++it;}}
      res.resize(distance(beg, it));
      return res;
  }
  auto ita = lower_bound(ab, ae, mask);
  auto itb = lower_bound(bb, be, mask);
  vector<uint32_t> p;
  if (ita == ae) {
    auto const & pb_qb = mergeAddMask(itb, be, mask, bb, itb);
    p = multiplication(ab, ae, pb_qb.begin(), pb_qb.end());
  }
  else {
    auto const & pa_qa = mergeAddMask(ita, ae, mask, ab, ita);
    if (itb == be) p = multiplication(pa_qa.begin(), pa_qa.end(), bb, be);
    else {
      auto const & pb_qb = mergeAddMask(itb, be, mask, bb, itb);
      p = multiplication(pa_qa.begin(), pa_qa.end(), pb_qb.begin(), pb_qb.end());
    }
  }
  auto qaqb = multiplication(ab, ita, bb, itb);
  unsigned const qaqb_size = qaqb.size();
  if (qaqb_size == 0) {
      for (auto & x : p) x |= mask;
      return p;
  }
  else {
    qaqb.resize(2*qaqb_size + p.size());
    auto beg = qaqb.begin();
    auto it0 = beg + qaqb_size;
    it0 = mergeXOR_addMask(p.begin(), p.end(), beg, it0, mask, it0);
    qaqb.resize(distance(beg, it0));
    return qaqb;
  }
}



Polynome operator*(Polynome const & a, Polynome const & b) {
  return Polynome(multiplication(a.begin(), a.end(), b.begin(), b.end()), true );
}

Polynome product(vector<Polynome> const & v, uint16_t u) {
  Polynome res (0); // res = 1
  for (unsigned i = 0; i < 16; ++i) {
    if (((u >> i) & 1) != 0) res *= v[i];
  }
  return res;
}

Polynome product(vector<Polynome> const & v, uint32_t u) {
  Polynome res (0); // res = 1
  for (unsigned i = 0; i < 32; ++i) {
    if (((u >> i) & 1) != 0) res *= v[i];
  }
  return res;
}

vector<Monome> Polynome::getTopMonomes() const {
  vector<Monome> res;

  auto const rend = v_mon.rend();
  for (auto rit = v_mon.rbegin(); rit != rend; ++rit) {
    if (none_of(res.begin(), res.end(), [&rit](Monome x){return (*rit & x) == *rit;})) res.emplace_back(*rit);
  }
  return res;
}

void printAsMonome(Monome x) {
    if (x == 0) cout << 1;
    else {
        for (unsigned i = 0; i < 26; ++i) {
            if (((x >> i) & 1) != 0) cout << (char) ('a' + i);
        }
        for (unsigned i = 26; i < 32; ++i) {
            if (((x >> i) & 1) != 0) cout << (char) ('A' + (i-26));
        }
    }
}

ostream & operator<<(ostream & os, Polynome const & p) {
  bool first = true;
  for (auto x : p) {
    if (first) first = false;
    else os << " + ";
    if (x == 0) os << 1;
    else {
      for (unsigned i = 0; i < 26; ++i) {
        if (((x >> i) & 1) != 0) os << (char) ('a' + i);
      }
      for (unsigned i = 26; i < 32; ++i) {
        if (((x >> i) & 1) != 0) os << (char) ('A' + (i-26));
      }
    }
  }
  if (first) cout << 0;
  return os;
}

set<unsigned> Polynome::getVariables() const {
  set<unsigned> res;
  unsigned n = 32;
  vector<unsigned> vars (n);
  for (unsigned i = 0; i < n; ++i) vars[i] = i;
  for (auto m : v_mon) {
    unsigned i = 0;
    while (i < n) {
      if (((m >> vars[i]) & 1) != 0) {
        res.emplace(vars[i]);
        vars[i] = vars[--n];
      }
      else ++i;
    }
  }
  return res;
}



vector<Polynome> getANF_S(vector<uint16_t> const & sbox) {
    unsigned two_to_n = sbox.size();
    unsigned n = 0;
    while ((1u << n) != two_to_n) ++n;

    vector<Polynome> anf (n);
    vector<vector<uint16_t>> m (two_to_n, vector<uint16_t>(two_to_n+1, 0));
    // init m
    for (unsigned x = 0; x < two_to_n; ++x) {
        for (unsigned i = 0; i < two_to_n; ++i) {
            if ((x & i) == i) m[x][i] = 1;
            else m[x][i] = 0;
        }
        m[x][two_to_n] = sbox[x];
    }
    // gauss
    for (unsigned c = 0; c < two_to_n; ++c) {
        unsigned l = c;
        while (m[l][c] == 0) ++l;
        swap(m[l], m[c]);
        for (++l; l < two_to_n; ++l) {
            if (m[l][c] != 0) {
                for (unsigned k = c; k <= two_to_n; ++k) m[l][k] ^= m[c][k];
            }
        }
        for (l = 0; l < c; ++l) {
            if (m[l][c] != 0) {
                for (unsigned k = c; k <= two_to_n; ++k) m[l][k] ^= m[c][k];
            }
        }
    }
    // ANF
    for (unsigned b = 0; b < n; ++b) {
        for (unsigned i = 0; i < two_to_n; ++i) {
            if (((m[i][two_to_n] >> b) & 1) == 1) anf[b].v_mon.emplace_back(i);
        }
    }

    return anf;
  }

vector<Polynome> getANF_Parallel_S(unsigned n, vector<uint16_t> const & sbox) {
  auto const & anf_s = getANF_S(sbox);
  auto const n_s = anf_s.size();
  vector<Polynome> anf (n*n_s);
  for (unsigned i = 0; i < n; ++i) {
    for (unsigned j = 0; j < n_s; ++j) {
      auto & p = anf[i*n_s + j];
      auto const & tmp = anf_s[j];
      for (auto m : tmp) p += (m << i*n_s);
    }
  }
  return anf;
}

Polynome composition(Polynome const & p, vector<Polynome> const & v) { // return p(v0, v1, ...)
  Polynome res; //res = 0
  for (auto m : p) res += product(v, m);
  return res;
}

vector<Polynome> composition(vector<Polynome> const & p, vector<Polynome> const & v) { // return p0(v0, v1, ...), ...
  vector<Polynome> res;
  for (auto const & pp : p) res.emplace_back(composition(pp, v));
  return res;
}

Polynome convertLtoPolynome(uint16_t l) {
  Polynome res;
  for (unsigned i = 0; i < 16; ++i) if (((l >> i) & 1) != 0) res += Monome(uint16_t(1) << i);
  return res;
}

vector<Polynome> convertLtoPolynome(vector<uint16_t> const & l) {
  vector<Polynome> res;
  for (auto line : l) res.emplace_back(convertLtoPolynome(line));
  return res;
}
