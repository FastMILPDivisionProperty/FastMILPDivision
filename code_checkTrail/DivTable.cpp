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
      uint16_t x = ~(uint16_t(1) << i);
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


uint64_t hashANF(vector<Polynome> const & anf) {
  uint64_t seed = anf.size();
  for (auto const & p : anf) {
    seed ^= p.nbTerms() + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    for (auto const & x : p) seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }
  return seed;
}

DivTable loadDivTable(vector<Polynome> const & anf) {
  auto const & h = hashANF(anf);
  string name = "tables/table_" + to_string(h) + ".txt";
  ifstream infile(name);
  if (infile.good()) return DivTable (infile);
  cout << "table does not exist, computing it (may take a long time)" << endl;
  DivTable div (anf);
  cout << "\r" << " - computed" << endl;
  ofstream outfile(name);
  div.save(outfile);
  return div;
}

vector<pair<vector<uint16_t>, unsigned>> loadTransitionsIN(vector<Polynome> const & anf) {
  auto const & h = hashANF(anf);
  string name = "tables/tableIN_" + to_string(h) + ".txt";
  ifstream infile(name);
  if (infile.good()) return loadTransitions(infile);
  cout << "table does not exist, computing it (may take a long time)" << endl;
  auto s_in = getTransitionsIn(anf);
  cout << "\r" << " - computed" << endl;
  ofstream outfile(name);
  saveTransitions(outfile, s_in);
  return s_in;
}

vector<pair<vector<uint16_t>, unsigned>> loadTransitionsOUT(vector<Polynome> const & anf) {
  auto const & h = hashANF(anf);
  string name = "tables/tableOUT_" + to_string(h) + ".txt";
  ifstream infile(name);
  if (infile.good()) return loadTransitions(infile);
  auto s_out = getTransitionsOut(anf);
  ofstream outfile(name);
  saveTransitions(outfile, s_out);
  return s_out;
}

void printTransitions(vector<pair<vector<uint16_t>, unsigned>> const & s) {
  for (auto const & my_pair : s) {
    cout << " - " << (unsigned) my_pair.second << ": ";
    for (auto const & m : my_pair.first) {printAsMonome(m); cout << " ";}
    cout << endl;
  }
}
