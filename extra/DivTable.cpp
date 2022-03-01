#include <algorithm>
#include <execution>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <vector>
#include <future>

#include "DivTable.hpp"
#include "Polynome.hpp"

using namespace std;

vector<vector<uint16_t>> invDivtable(vector<vector<uint16_t>> const &divtable) {
  vector<vector<uint16_t>> res(1u << 16);
  for (unsigned i = 0; i < divtable.size(); ++i) {
    unsigned n_div = divtable[i].size();
    for (unsigned j = 0; j < n_div; ++j) {
      res[divtable[i][j]].emplace_back(i);
    }
  }
  return res;
}

void addToDivtable(vector<Monome> const & top, uint16_t const u, unsigned const i, vector<vector<uint16_t>> & divtable) {
  if (none_of(top.begin(), top.end(), [i](uint16_t mon) { return (mon & i) == i; })) return;
  auto & v = divtable[i];
  if (any_of(v.begin(), v.end(), [u](uint16_t x){return (x & u) == x;})) return;
  v.emplace_back(u);
}

auto fillAllProd_MapToSet(map<Monome, Polynome> const & map_c) {
  auto fn = [](Polynome const& pi, Polynome const & pj){return pi.last() > pj.last();};
  set<Polynome, bool(*)(Polynome const &, Polynome const &)> s (fn);
  for (auto const & p : map_c) {
    auto pp = p.second;
    auto it = s.begin();
    auto const & end = s.end();
    while (!pp.isNull()) {
      if (it == end || it->last() < pp.last()) {s.emplace(move(pp)); break;}
      if (it->last() == pp.last()) pp += *it;
      ++it;
    }
  }
  return s;
}

void fillAllProd(vector<vector<bool>> & allprod,  vector<Polynome> const &anf) {
  auto const & n = anf.size();
  auto const & n0 = n/2;
  auto const & n1 = n - n0;
  //uint32_t const & mask = (1u << n) - 1;

  vector<map<Monome, Polynome>> vmap0 (1u << n0);
  vector<vector<Polynome>> vs0 (1u << n0);
  #pragma omp parallel for schedule(dynamic)
  for (uint32_t c = 0; c < 1u << n0; ++c) {
    auto const & p = product(anf, c);
    auto & map_c = vmap0[c];
    for (auto m : p) map_c[m & 0xFFFF0000u] += Monome(m & 0xFFFF);
    auto const & s = fillAllProd_MapToSet(map_c);
    vs0[c] = vector<Polynome> (s.begin(), s.end());
  }

  vector<map<Monome, Polynome>> vmap1 (1u << n1);
  vector<vector<Polynome>> vs1 (1u << n1);
  #pragma omp parallel for schedule(dynamic)
  for (uint32_t c = 0; c < 1u << n1; ++c) {
    auto const & p = product(anf, c << n0);
    auto & map_c = vmap1[c];
    for (auto m : p) map_c[m & 0xFFFF0000u] += Monome(m & 0xFFFF);
    auto const & s = fillAllProd_MapToSet(map_c);
    vs1[c] = vector<Polynome> (s.begin(), s.end());
  }

  uint64_t cpt = 0;

  vector<uint32_t> vc (1u << n);
  for (uint32_t i = 0; i < 1u << n; ++i) vc[i] = i;
  sort(vc.begin(), vc.end(), [](auto const & i, auto const & j){
    auto const & bi = __builtin_popcount(i);
    auto const & bj = __builtin_popcount(j);
    return bi > bj || (bi == bj && i < j);
  });

  #pragma omp parallel for schedule(dynamic)
  for (uint32_t i = 0; i < 1u << n; ++i) {
    auto & v01 = allprod[vc[i]];
    auto const & c0 = vc[i] & ((1u << n0) - 1);
    auto const & c1 = vc[i] >> n0;
    auto const & map_c0 = vmap0[c0];
    auto const & map_c1 = vmap1[c1];
    set<Monome> keys;
    for (auto const & p0 : map_c0) {
      for (auto const & p1 : map_c1) keys.emplace(p0.first | p1.first);
    }
    if (keys.size() == map_c0.size()*map_c1.size()) { // keys are like independent
      auto const & s0 = vs0[c0];
      auto const & s1 = vs1[c1];
      for (auto const & p0 : s0) {
        for (auto const & p1 : s1) {
          auto const & pp = p0*p1;
          for (auto const & m : pp) v01[m] = true;
        }
      }
    }
    else {
      for (auto const & k : keys) {
        Polynome pp0;
        for (auto const & p0 : map_c0) if ((k & p0.first) == p0.first) pp0 += p0.second;
        Polynome pp1;
        for (auto const & p1 : map_c1) if ((k & p1.first) == p1.first) pp1 += p1.second;
        auto const & pp = pp0*pp1;
        for (auto const & m : pp) v01[m] = true;
      }
    }
    #pragma omp critical
    {
      cout << "\r" << ++cpt << "/" << (1u << n) << flush;
    }
  }
  cout << " - done" << endl;
}


vector<vector<uint16_t>> getDivisionTable(vector<Polynome> const &anf) {
  unsigned const n = anf.size();
  unsigned const n_mon = 1u << n;
  vector<vector<bool>> allprod (n_mon, vector<bool>(n_mon, false));
  fillAllProd(allprod, anf);
  vector<vector<uint16_t>> divtable(n_mon);
  divtable[0].emplace_back(0);
  divtable[n_mon - 1].emplace_back(n_mon - 1);
  map<unsigned, vector<uint16_t>> map_u;
  for (unsigned u = 1; u < n_mon; ++u) map_u[degree(u)].emplace_back(u);
  for (auto const & p : map_u) {
    auto const & v_u = p.second;
    for (auto const & u : v_u) {
      vector<Monome> top;
      auto const & v = allprod[u];
      for (unsigned j = n_mon; j-- != 0; ) {
        if (v[j] && none_of(top.begin(), top.end(), [j](unsigned x){return (x & j) == j;})) top.emplace_back(j);
      }
      for_each(execution::par, top.begin(), top.end(), [&divtable, u](auto const & m){
        auto & vv = divtable[m];
        if (none_of(vv.begin(), vv.end(), [u](uint16_t x){return (x & u) == x;})) vv.emplace_back(u);
      });
    }
  }
  for (unsigned d = n-1; d-- != 0; ) {
    auto const & v_m = map_u[d];
    for_each(execution::par, v_m.begin(), v_m.end(), [&divtable, n](auto const & m){
      auto & div_m = divtable[m];
      set<uint16_t> to_add;
      for (unsigned i = 0; i < n; ++i) {
        auto const & mm = m | (uint16_t(1) << i);
        if (mm == m) continue;
        auto const & div_mm = divtable[mm];
        to_add.insert(div_mm.begin(), div_mm.end());
      }
      for (auto const & u : to_add) {
        unsigned j = 0;
        unsigned n_m = div_m.size();
        while (j < n_m) {
          uint16_t const & inter = u & div_m[j];
          if (inter == div_m[j]) break;
          if (inter == u) div_m[j] = div_m[--n_m];
          else ++j;
        }
        if (j == n_m) {
          div_m.resize(n_m);
          div_m.emplace_back(u);
        }
      }
    });
  }
  return divtable;
}


void addToSasLinearSystem(Polynome p, set<Polynome> & s) {
  if (p.isNull()) return;
  auto it = s.begin();
  auto const end = s.end();
  while (it != end && p.first() >= it->first()) { // gaussian elimination
    if (p.first() == it->first()) {
      p += *it;
      if (p.isNull()) return;
    }
    ++it;
  }
  s.emplace(move(p));
}

void fillAllProd_In(vector<set<Polynome>> & allprod,  vector<Polynome> const &anf) {
  auto const & n = anf.size();
  auto const & n0 = n/2;
  auto const & n1 = n - n0;
  uint32_t const & mask = (1u << n) - 1;

  vector<map<Monome, Polynome>> vmap0 (1u << n0);
  vector<vector<Polynome>> vs0 (1u << n0);
  #pragma omp parallel for schedule(dynamic)
  for (uint32_t c = 0; c < 1u << n0; ++c) {
    auto const & p = product(anf, c);
    auto & map_c = vmap0[c];
    for (auto m : p) map_c[m & 0xFFFF0000u] += Monome(m & 0xFFFF);
    auto const & s = fillAllProd_MapToSet(map_c);
    vs0[c] = vector<Polynome> (s.begin(), s.end());
  }

  vector<map<Monome, Polynome>> vmap1 (1u << n1);
  vector<vector<Polynome>> vs1 (1u << n1);
  #pragma omp parallel for schedule(dynamic)
  for (uint32_t c = 0; c < 1u << n1; ++c) {
    auto const & p = product(anf, c << n0);
    auto & map_c = vmap1[c];
    for (auto m : p) map_c[m & 0xFFFF0000u] += Monome(m & 0xFFFF);
    auto const & s = fillAllProd_MapToSet(map_c);
    vs1[c] = vector<Polynome> (s.begin(), s.end());
  }

  uint64_t cpt = 0;

  vector<uint32_t> vc (1u << n);
  for (uint32_t i = 0; i < 1u << n; ++i) vc[i] = i;
  sort(vc.begin(), vc.end(), [](auto const & i, auto const & j){
    auto const & bi = __builtin_popcount(i);
    auto const & bj = __builtin_popcount(j);
    return bi > bj || (bi == bj && i < j);
  });

  #pragma omp parallel for schedule(dynamic)
  for (uint32_t i = 0; i < 1u << n; ++i) {
    auto & v01 = allprod[vc[i]];
    auto const & c0 = vc[i] & ((1u << n0) - 1);
    auto const & c1 = vc[i] >> n0;
    auto const & map_c0 = vmap0[c0];
    auto const & map_c1 = vmap1[c1];
    set<Monome> keys;
    for (auto const & p0 : map_c0) {
      for (auto const & p1 : map_c1) keys.emplace(p0.first | p1.first);
    }
    if (keys.size() == map_c0.size()*map_c1.size()) { // keys are like independent
      auto const & s0 = vs0[c0];
      auto const & s1 = vs1[c1];
      for (auto const & p0 : s0) {
        for (auto const & p1 : s1) {
          auto const & pp = p0*p1;
          Polynome tmp;
          for (unsigned b = n; b-- != 0; ) {
            Monome const & m = (~(1u << b)) & mask;
            if (pp.contains(m)) tmp += m;
          }
          addToSasLinearSystem(move(tmp), v01);
        }
      }
    }
    else {
      for (auto const & k : keys) {
        Polynome pp0;
        for (auto const & p0 : map_c0) if ((k & p0.first) == p0.first) pp0 += p0.second;
        Polynome pp1;
        for (auto const & p1 : map_c1) if ((k & p1.first) == p1.first) pp1 += p1.second;
        auto const & pp = pp0*pp1;
        Polynome tmp;
        for (unsigned b = n; b-- != 0; ) {
          Monome const & m = (~(1u << b)) & mask;
          if (pp.contains(m)) tmp += m;
        }
        addToSasLinearSystem(move(tmp), v01);
      }
    }
    #pragma omp critical
    {
      cout << "\r" << ++cpt << "/" << (1u << n) << flush;
    }
  }
  cout << " - done" << endl;
}

bool isBetterIn(pair<vector<uint16_t>, unsigned> const & p1, pair<vector<uint16_t>, unsigned> const & p2) {
  if (p1.first == p2.first) return __builtin_popcount(p1.second) <= __builtin_popcount(p2.second);
  return all_of(p1.first.begin(), p1.first.end(), [&p2](auto u1){
    return any_of(p2.first.begin(), p2.first.end(), [u1](auto u2){return (u1 & u2) == u2;});
  });
}

vector<pair<vector<uint16_t>, unsigned>> getTransitionsIn(vector<Polynome> const &anf) {
  unsigned const n = anf.size();
  unsigned const n_mon = 1u << n;

  vector<set<Polynome>> allprod (n_mon);
  fillAllProd_In(allprod, anf);

  vector<uint16_t> v_u;
  vector<uint16_t> unremovable;
  for (unsigned u = 0; u < n_mon; ++u) if (!allprod[u].empty()) {
    if (allprod[u].size() == n) unremovable.emplace_back(u);
    else v_u.emplace_back(u);
  }
  sort(v_u.begin(), v_u.end(), [](auto u1, auto u2){
    auto const n1 = degree(u1);
    auto const n2 = degree(u2);
    return n1 < n2 || (n1 == n2 && u1 < u2);
  });

  vector<pair<vector<uint16_t>, unsigned>> res;

  #pragma omp parallel for
  for (unsigned l = 1; l < n_mon; ++l) {
    map<uint16_t, uint8_t> F;
    for (unsigned i = 0; i < n; ++i) {
      uint16_t x = (~(uint16_t(1) << i)) & (n_mon - 1);
      F[x] = (l >> i) & 1;
    }
    vector<uint16_t> v;
    for (auto u : v_u) {
      if (any_of(v.begin(), v.end(), [u](auto const & uu){return (u & uu) == uu;})) continue;
      auto const & sys = allprod[u];
      if (any_of(sys.begin(), sys.end(), [&F](auto const & p){
        uint8_t sum = 0;
        for (auto const & m : p) sum ^= F[m];
        return sum != 0;
      })) v.emplace_back(u);
    }
    for (auto const & u : unremovable) if (none_of(v.begin(), v.end(), [u](auto const & uu){return (u & uu) == uu;})) v.emplace_back(u);
    pair<vector<uint16_t>, unsigned> my_pair (move(v), l);
    #pragma omp critical
    {
      if (none_of(res.begin(), res.end(), [&my_pair](auto const & x){return isBetterIn(x, my_pair);})) {
        unsigned i = 0;
        unsigned n_res = res.size();
        while (i < n_res) {
          if (isBetterIn(my_pair, res[i])) res[i] = res[--n_res];
          else ++i;
        }
        res.resize(n_res);
        res.emplace_back(move(my_pair));
      }
    }
  }
  return res;
}



bool isBetterOut(pair<vector<uint16_t>, unsigned> const & p1, pair<vector<uint16_t>, unsigned> const & p2) {
  if (p1.first == p2.first) return __builtin_popcount(p1.second) <= __builtin_popcount(p2.second);
  return all_of(p1.first.begin(), p1.first.end(), [&p2](auto u1){
    return any_of(p2.first.begin(), p2.first.end(), [u1](auto u2){return (u1 & u2) == u1;});
  });
}

vector<pair<vector<uint16_t>, unsigned>> getTransitionsOut(vector<Polynome> const &anf) {
  unsigned const n = anf.size();

  vector<pair<vector<uint16_t>, unsigned>> res;
  vector<unsigned> v_l;
  for (unsigned l = 1; l < 1u << n; ++l) v_l.emplace_back(l); // sort linear combinations by the number of non zero coefficients
  sort(v_l.begin(), v_l.end(), [](auto l1, auto l2){
    auto const n1 = __builtin_popcount(l1);
    auto const n2 = __builtin_popcount(l2);
    return n1 < n2 || (n1 == n2 && l1 < l2);
  });

  #pragma omp parallel for schedule(dynamic)
  for (unsigned i = 0; i < v_l.size(); ++i) {
    auto const & l = v_l[i];
    Polynome pl;
    for (unsigned j = 0; j < n; ++j) if (((l >> j) & 1) != 0) pl += anf[j];
    set<uint16_t> sl;
    for (auto const & m : pl) sl.emplace(m & 0xFFFF);
    vector<uint16_t> top;
    auto const rend = sl.rend();
    for (auto rit = sl.rbegin(); rit != rend; ++rit) {
      if (none_of(top.begin(), top.end(), [&rit](uint16_t const & x){return (*rit & x) == *rit;})) top.emplace_back(*rit);
    }
    pair<vector<uint16_t>, unsigned> x (move(top), l);
    #pragma omp critical
    {
      if (none_of(res.begin(), res.end(), [&x](auto const & y){return isBetterOut(y, x);})) {
        unsigned j = 0;
        unsigned n_res = res.size();
        while (j < n_res) {
          if (isBetterOut(x, res[j])) res[j] = res[--n_res];
          else ++j;
        }
        res.resize(n_res);
        res.emplace_back(move(x));
      }
    }
  }
  return res;
}

void saveTransitions(ofstream & output_file, vector<pair<vector<uint16_t>, unsigned>> const & vt) {
  for (auto const & p : vt) output_file << p.second << " ";
  output_file << "\n";
  for (auto const & p : vt) {
    for (auto const & m : p.first) output_file << m << " ";
    output_file << "\n";
  }
}

vector<pair<vector<uint16_t>, unsigned>> loadTransitions(ifstream & input_file) {
  vector<pair<vector<uint16_t>, unsigned>> res;
  string temp;
  getline(input_file, temp);
  istringstream buffer(temp);
  vector<unsigned> linear ((istream_iterator<unsigned>(buffer)), istream_iterator<unsigned>());
  for (auto const & l : linear) {
    getline(input_file, temp);
    istringstream buffer2(temp);
    vector<uint16_t> comb ((istream_iterator<uint16_t>(buffer2)), istream_iterator<uint16_t>());
    res.emplace_back(move(comb), l);
  }
  return res;
}



void printDivisionTable(DivTable const &divtable) {
  unsigned const n_mon = divtable.size();
  unsigned n = 0;
  while ((1u << n) != n_mon)
    ++n;

  vector<unsigned> v;
  v.reserve(n_mon);
  for (unsigned i = 0; i < n_mon; ++i)
    v.emplace_back(i);
  sort(v.begin(), v.end(), [](unsigned i, unsigned j) {
    return __builtin_popcount(i) < __builtin_popcount(j) ||
           (__builtin_popcount(i) == __builtin_popcount(j) && i < j);
  });
  for (auto i : v) {
    for (unsigned j = 0; j < n; ++j)
      cout << ((i >> j) & 1);
    cout << " | ";
    for (auto k : divtable.transitionsFrom(i)) {
      for (unsigned j = 0; j < n; ++j)
        cout << ((k >> j) & 1);
      cout << " - ";
    }
    cout << endl;
  }
}

uint16_t convertToNibbles(uint16_t b) {
  uint16_t res = __builtin_popcount(b & 0xF);
  for (unsigned i = 1; i < 4; ++i) {
    res |= __builtin_popcount((b >> 4*i) & 0xF) << 3*i;
  }
  return res;
}

uint64_t convertToNibbles(uint64_t b) {
  uint64_t res = __builtin_popcountll(b & 0xF);
  for (unsigned i = 1; i < 16; ++i) {
    res |= uint64_t(__builtin_popcountll((b >> 4*i) & 0xF)) << 3*i;
  }
  return res;
}


vector<vector<uint16_t>> getDivisionTableNibble(vector<vector<uint16_t>> const & div) {
  vector<set<uint16_t>> tmp (1u << 12);
  for (unsigned i = 0; i < 1u << 16; ++i) {
    auto const & div_i = div[i];
    auto & s = tmp[convertToNibbles(uint16_t(i))];
    for (auto x : div_i) s.emplace(convertToNibbles(x));
  }
  vector<vector<uint16_t>> res;
  for (auto const & s : tmp) res.emplace_back(vector<uint16_t> (s.begin(), s.end()));
  return res;
}

bool removeUnusedTrails_tmp(vector<vector<uint8_t>> const & words, unsigned const c, uint16_t const pm, DivTable const & div) {
  for (auto const & m0 : words[4*c]) {
    uint16_t const & tmp0 = m0 | pm;
    for (auto const & m1 : words[4*c + 1]) {
      uint16_t const & tmp1 = (m1 << 4) | tmp0;
      for (auto const & m2 : words[4*c + 2]) {
        uint16_t const & tmp2 = (m2 << 8) | tmp1;
        for (auto const & m3 : words[4*c + 3]) {
          uint16_t const & tmp = (m3 << 12) | tmp2;
          if (!div.transitionsFrom(tmp).empty()) return true;
        }
      }
    }
  }
  return false;
}

void removeUnusedTrailsIN(uint64_t in64, vector<vector<DivTable>> & vdiv, vector<uint8_t> const & perm) {
  uint16_t in[16];
  for (unsigned i = 0; i < 16; ++i) in[i] = 1u << ((in64 >> 4*i) & 0xF);

  for (unsigned r = 0; r < vdiv.size(); ++r) {
    uint16_t in_p [16];
    for (unsigned i = 0; i < 16; ++i) {in_p[perm[i]] = in[i]; in[i] = 0;}
    for (unsigned c = 0; c < 4; ++c) {
      for (unsigned i = 0; i < 1u << 16; ++i) {
        unsigned j = 0;
        while (j < 4 && ((uint16_t(1) << ((i >> 4*j) &0xF)) & in_p[4*c + j]) != 0) ++j;
        auto & ti = vdiv[r][c].table[i];
        if (j != 4) ti.clear();
        for (auto const & m : ti) {
          for (unsigned k = 0; k < 4; ++k) in[k + 4*c] |= uint16_t(1) << ((m >> 4*k) & 0xF);
        }
      }
    }
  }
}

void removeUnusedTrailsOUT(uint64_t out64, vector<vector<DivTable>> & vdiv, vector<uint8_t> const & perm) {
  uint64_t p_out64 = 0;
  for (unsigned i = 0; i < 16; ++i) p_out64 |= ((out64 >> 4*perm[i]) & 0xF) << 4*i;

  auto const & R = vdiv.size();
  for (unsigned c = 0; c < 4; ++c) {
    auto & t = vdiv[R-1][c].table;
    uint16_t const & check = (p_out64 >> 16*c) & 0xFFFF;
    #pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < 1u << 16; ++i) {
      bool flag = any_of(t[i].begin(), t[i].end(), [check](auto const & m){return (m & check) == m;});
      t[i].clear();
      if (flag) t[i].emplace_back(0);
    }
  }
  for (unsigned r = R-1; r-- != 0; ) {
    bool next = false;
    uint16_t out[16];
    for (unsigned c = 0; c < 4; ++c) {
      auto const & vdiv_c = vdiv[r][c];
      for (unsigned i = 0; i < 1u << 16; ++i) {
        auto const & t = vdiv_c.transitionsFrom(i);
        for (auto const & m : t) {
          for (unsigned k = 0; k < 4; ++k) out[k + 4*c] |= uint16_t(1) << ((m >> 4*k) & 0xF);
        }
      }
    }
    uint16_t out_p[16];
    for (unsigned i = 0; i < 16; ++i) {out_p[perm[i]] = out[i]; out[i] = 0;}
    for (unsigned c = 0; c < 4; ++c) {
      auto const & vdiv_c = vdiv[r+1][c];
      for (unsigned i = 0; i < 1u << 16; ++i) {
        auto const & t = vdiv_c.transitionsFrom(i);
        if (t.empty()) continue;
        unsigned j = 0;
        while (j < 4 && ((uint16_t(1) << ((i >> 4*j) &0xF)) & out_p[4*c + j]) != 0) ++j;
        if (j != 4) continue;
        for (unsigned k = 0; k < 4; ++k) out[k + 4*c] |= uint16_t(1) << ((i >> 4*k) & 0xF);
      }
    }
    for (unsigned c = 0; c < 4; ++c) {
      auto & tc = vdiv[r][c].table;
      for (unsigned i = 0; i < 1u << 16; ++i) {
        auto & ti = tc[i];
        unsigned j = 0;
        unsigned n = ti.size();
        while (j < n) {
          auto & tij = ti[j];
          unsigned k = 0;
          while (k < 4 && ((uint16_t(1) << ((tij >> 4*k) &0xF)) & out[perm[4*c + k]]) != 0) ++k;
          if (k == 4) ++j;
          else tij = ti[--n];
        }
        ti.resize(n);
        if (n == 0) next = true;
      }
    }
    if (!next) return;
  }
}

void removeUnusedTrails(uint64_t in64, uint64_t out64, vector<vector<DivTable>> & vdiv, vector<uint8_t> const & perm) {
  removeUnusedTrailsIN(in64, vdiv, perm);
  removeUnusedTrailsOUT(out64, vdiv, perm);
}

// void removeUnusedTrails(uint64_t out64, vector<vector<DivTable>> & vdiv, vector<uint8_t> const & perm) {
//   uint64_t p_out64 = 0;
//   for (unsigned i = 0; i < 16; ++i) p_out64 |= ((out64 >> 4*perm[i]) & 0xF) << 4*i;
//
//   auto const & R = vdiv.size();
//   {
//     for (unsigned c = 0; c < 4; ++c) {
//       auto & t = vdiv[R-1][c].table;
//       uint16_t const & check = (p_out64 >> 16*c) & 0xFFFF;
//       #pragma omp parallel for schedule(dynamic)
//       for (unsigned i = 0; i < 1u << 16; ++i) {
//         bool flag = any_of(t[i].begin(), t[i].end(), [check](auto const & m){return (m & check) == m;});
//         t[i].clear();
//         if (flag) t[i].emplace_back(0);
//       }
//     }
//   }
//   for (unsigned r = R-1; r-- != 0; ) {
//     bool next = false;
//     for (unsigned c = 0; c < 4; ++c) {
//       vector<uint8_t> nib (16);
//       for (unsigned i = 0; i < 16; ++i) nib.emplace_back(i);
//       vector<vector<uint8_t>> words (16, nib);
//       for (unsigned i = 0; i < 4; ++i) {
//         auto const & j = perm[i + 4*c];
//         words[j].resize(1);
//       }
//
//       vector<uint8_t> valid (1u << 16, 0);
//
//       vector<vector<uint16_t>> vpm (1u << 16);
//       for (unsigned m = 0; m < 1u << 16; ++m) {
//         vector<uint16_t> pm (4, 0);
//         for (unsigned i = 0; i < 4; ++i) {
//           auto const & j = perm[i + 4*c];
//           pm[j/4] |= ((m >> 4*i) & 0xF) << 4*(j%4);
//         }
//         vpm[m] = move(pm);
//       }
//
//       vector<set<uint16_t>> v_hassol (4);
//       for (unsigned cc = 0; cc < 4; ++cc) {
//         set<uint16_t> s;
//         for (auto const & pm : vpm) s.emplace(pm[cc]);
//         vector<uint16_t> vs (s.begin(), s.end());
//         auto const & vs_size = vs.size();
//         #pragma omp parallel for
//         for (unsigned i = 0; i < vs_size; ++i) {
//           auto const & m0 = vs[i];
//           if (removeUnusedTrails_tmp(words, cc, m0, vdiv[r+1][cc])) {
//             #pragma omp critical
//             {
//               v_hassol[cc].emplace(m0);
//             }
//           }
//         }
//       }
//
//       #pragma omp parallel for
//       for (unsigned m = 0; m < 1u << 16; ++m) {
//         auto const & pm = vpm[m];
//         unsigned cc = 0;
//         while (cc < 4 && v_hassol[cc].count(pm[cc]) != 0) ++cc;
//         if (cc == 4) valid[m] = 1;
//       }
//
//       uint64_t cpt = 0;
//       auto & t = vdiv[r][c].table;
//       #pragma omp parallel for
//       for (unsigned i = 0; i < 1u << 16; ++i) {
//         auto & ti = t[i];
//         unsigned n = ti.size();
//         unsigned j = 0;
//         while (j < n) {
//           if (valid[ti[j]] == 0) ti[j] = ti[--n];
//           else ++j;
//         }
//         #pragma omp critical
//         {
//           if (n == 0) next = true;
//           cpt += ti.size() - n;
//         }
//         ti.resize(n);
//       }
//
//       //cout << r << " - " << c << ": " << cpt << endl;
//
//     }
//     if (!next) {/*cout << " " << r << endl;*/ return;}
//   }
// }

vector<uint16_t> getAllNibblesOutput(vector<DivTable> & div) {
  vector<uint16_t> out (16);
  for (unsigned c = 0; c < 4; ++c) { // get all possible nibbles at output of step R-2
    auto const & vc = div[c];
    for (unsigned i = 0; i < 1u << 16; ++i) {
      auto const & ti = vc.transitionsFrom(i);
      for (auto const & m : ti) {
        for (unsigned j = 0; j < 4; ++j) out[4*c + j] |= uint16_t(1) << ((m >> 4*j) & 0xF);
      }
    }
  }
  return out;
}

// void removeUnnecessaryTrails(vector<vector<DivTable>> & vdiv, vector<uint8_t> const & perm) {
//   vector<map<uint64_t, bool>> better (4);
//   auto const & R = vdiv.size();
//
//   auto const & out = getAllNibblesOutput(vdiv[R-2]); // get all possible nibbles at output of step R-2
//   vector<vector<uint16_t>> p_out (4, vector<uint16_t> (4));
//   for (unsigned i = 0; i < 16; ++i) p_out[perm[i]/4][perm[i]%4] = getCorrespondingVector(out[i]);
//
//   for (unsigned c = 0; c < 4; ++c) {
//     auto & div = vdiv[R-2][c].table;
//     for (unsigned i = 0; i < 1u << 16; ++i) {
//       auto & ti = div[i];
//       vector<uint16_t> new_ti;
//       new_ti.reserve(ti.size());
//       for (auto const & m : ti) {
//         unsigned j = 0;
//         unsigned n = new_ti.size();
//         while (j < n) {
//
//         }
//       }
//     }
//   }
//
//
//
//
//
// }

void removeUnnecessaryTrails(vector<vector<DivTable>> & vdiv, vector<uint8_t> const & perm) {
  vector<vector<vector<bool>>> better (4, vector<vector<bool>>(1u << 16, vector<bool>(1u << 16, false)));

  auto const & R = vdiv.size();
  for (unsigned c = 0; c < 4; ++c) {
    #pragma omp parallel for
    for (unsigned m0 = 0; m0 < 1u << 16; ++m0) {
      auto & bcm0 = better[c][m0];
      for (unsigned m1 = 0; m1 < 1u << 16; ++m1) {
        bcm0[m1] = (m0 & m1) == m0;
       }
    }
  }

  for (unsigned r = R-1; r-- != 0; ) {
    vector<vector<vector<bool>>> better_tmp (4, vector<vector<bool>>(1u << 16, vector<bool>(1u << 16, false)));

    uint16_t out[16];
    for (unsigned c = 0; c < 4; ++c) {
      auto const & vdiv_c = vdiv[r][c];
      for (unsigned i = 0; i < 1u << 16; ++i) {
        auto const & t = vdiv_c.transitionsFrom(i);
        for (auto const & m : t) {
          for (unsigned k = 0; k < 4; ++k) out[k + 4*c] |= uint16_t(1) << ((m >> 4*k) & 0xF);
        }
      }
    }
    vector<vector<uint8_t>> words_origin (16);
    for (unsigned i = 0; i < 16; ++i) {
      auto & t = words_origin[perm[i]];
      for (unsigned j = 0; j < 16; ++j) if (((out[i] >> j) & 1) != 0) t.emplace_back(j);
    }

    for (unsigned c = 0; c < 4; ++c) {

      vector<vector<uint16_t>> all_words (4);
      {
        //vector<uint8_t> nib (16);
        //for (unsigned i = 0; i < 16; ++i) nib.emplace_back(i);
        //vector<vector<uint8_t>> words (16, nib);
        vector<vector<uint8_t>> words (words_origin);
        for (unsigned i = 0; i < 4; ++i) {
          auto const & j = perm[i + 4*c];
          words[j].resize(1);
        }
        for (unsigned cc = 0; cc < 4; ++cc) {
          for (auto const & m0 : words[4*cc]) {
            for (auto const & m1 : words[4*cc + 1]) {
              uint16_t const & tmp1 = (m1 << 4) | m0;
              for (auto const & m2 : words[4*cc + 2]) {
                uint16_t const & tmp2 = (m2 << 8) | tmp1;
                for (auto const & m3 : words[4*cc + 3]) {
                  all_words[cc].emplace_back((m3 << 12) | tmp2);
                }
              }
            }
          }
        }
      }


      vector<vector<uint16_t>> vpm (1u << 16, vector<uint16_t> (4, 0));
      for (unsigned m = 0; m < 1u << 16; ++m) {
        auto & pm = vpm[m];
        for (unsigned i = 0; i < 4; ++i) {
          auto const & j = perm[i + 4*c];
          pm[j/4] |= ((m >> 4*i) & 0xF) << 4*(j%4);
        }
      }

      vector<vector<unsigned>> v_better (4);
      for (unsigned cc = 0; cc < 4; ++cc) {
        auto const & div = vdiv[r+1][cc];
        auto const & bcc = better[cc];
        set<uint16_t> s;
        for (auto const & pm : vpm) s.emplace(pm[cc]);
        for (unsigned m0 : s) {
          for (unsigned m1 : s) {
            if ((m0 & m1) == m0 || all_of(execution::par, all_words[cc].begin(), all_words[cc].end(), [&div, m0, m1, &bcc](auto const & m){
              auto const & v0 = div.transitionsFrom(m | m0);
              auto const & v1 = div.transitionsFrom(m | m1);
              for (auto const & x1 : v1) {
                if (none_of(v0.begin(), v0.end(), [x1, &bcc](auto const & x0){return bcc[x0][x1];})) return false;
              }
              return true;
            })) v_better[cc].emplace_back((m0 << 16) | m1);
          }
        }
      }

      auto & bc = better_tmp[c];

      #pragma omp parallel for
      for (unsigned m0 = 0; m0 < 1u << 16; ++m0) {
        auto const & pm0 = vpm[m0];
        auto & bcm0 = bc[m0];
        for (unsigned m1 = 0; m1 < 1u << 16; ++m1) {
          auto const & pm1 = vpm[m1];
          unsigned cc = 0;
          while (cc < 4 && (binary_search(v_better[cc].begin(), v_better[cc].end(), (((unsigned) pm0[cc]) << 16) | pm1[cc]))) ++cc;
          if (cc == 4) bcm0[m1] = true;
        }
      }

      // vector<vector<uint16_t>> eq;
      // for (unsigned i = 0; i < 1u << 16; ++i) {
      //   bool flag = false;
      //   for (auto & v : eq) {
      //     if (bc[i][v[0]] && bc[v[0]][i]) {
      //       v.emplace_back(i);
      //       flag = true;
      //       break;
      //     }
      //   }
      //   if (!flag) eq.emplace_back(vector<uint16_t> (1, i));
      // }

      auto & table = vdiv[r][c].table;

      // if (eq.size() != (1u << 16)) {
      //   vector<uint16_t> map_eq (1u << 16);
      //   for (auto const & v : eq) {
      //     for (auto const & x : v) map_eq[x] = v[0];
      //   }
      //   for (unsigned i = 0; i < 1u << 16; ++i) {
      //     auto & ti = table[i];
      //     for (auto & x : ti) x = map_eq[x];
      //   }
      // }

      uint64_t cpt = 0;


      #pragma omp parallel for schedule(dynamic)
      for (unsigned i = 0; i < 1u << 16; ++i) {
        auto & ti = table[i];
        vector<uint16_t> tmp;
        tmp.reserve(ti.size());
        for (auto const & m : ti) {
          if (any_of(tmp.begin(), tmp.end(), [&bc, m](auto const & mm){return bc[mm][m];})) continue;
          unsigned j = 0;
          unsigned n = tmp.size();
          while (j < n) {
            if (bc[m][tmp[j]]) tmp[j] = tmp[--n];
            else ++j;
          }
          tmp.resize(n);
          tmp.emplace_back(m);
        }
        #pragma omp critical
        {
          cpt += ti.size() - tmp.size();
        }
        ti = move(tmp);
      }
      //cout << r << " - " << c << ": " << cpt << endl;
    }
    better = move(better_tmp);
  }
}
