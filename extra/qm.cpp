#include "DivTable.hpp"
#include "Polynome.hpp"


using namespace std;

vector<pair<uint32_t, uint32_t>> forbidDiv(DivTable const & div, unsigned n) {
  vector<vector<uint32_t>> v (1u << n);
  #pragma omp parallel for schedule(dynamic)
  for (uint32_t x = 0; x < 1u << n; ++x) {
    auto const & dx = div.transitionsFrom(x);
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
      //res.emplace_back((mask & ~x) | ((mask & p.first) << n), x);
      res.emplace_back(x, p.first);
    }
  }
  return res;
}

void simplifyImp(vector<pair<uint32_t, uint32_t>> const & imp, unsigned n){
  map<uint64_t,unsigned> histo;
  map<unsigned,set<uint64_t>> mapImp;
  map<uint64_t, set<unsigned>> mapInv;
  for (unsigned i_imp = 0; i_imp < imp.size(); ++i_imp) {
    auto const & p = imp[i_imp];
    for (unsigned i = 0; i < n; ++i) {
      if (((1ull << i) & p.first) != 0) {
        uint64_t x = ((1ull << i) ^ p.first) | (p.second << n);
        mapImp[i_imp].emplace(x);
        histo[x] += 1;
      }
    }
    for (unsigned i = 0; i < n; ++i) {
      if (((1ull << i) & ~p.second) != 0) {
        uint64_t x = p.first | ((p.second | (1ull << i))  << n);
        mapImp[i_imp].emplace(x);
        histo[x] += 1;
      }
    }
  }

  for (auto const & p : mapImp) {
    for (auto const & x : p.second) mapInv[x].emplace(p.first);
  }

  unsigned cpt = 0;

  map<unsigned, set<uint64_t>> histoInv;
  for (auto const & p : histo) histoInv[p.second].emplace(p.first);

  map<uint64_t, uint64_t> finalA;
  vector<pair<uint64_t, bool>> myorder;

  while (!histo.empty()) {
    cout << "\r" << histo.size() << "    " << flush;
    set<uint64_t> toprocess;

    for (auto c : histoInv[1]) {
      histo.erase(c);
      auto it_inv = mapInv.find(c);
      auto x = *it_inv->second.begin();
      auto & s = mapImp[x];
      s.erase(c);
      if (s.size() == 1) {
        toprocess.emplace(*s.begin());
      }
      else if (s.size() == 0) {
        ++cpt;
        mapImp.erase(x);
      }
      mapInv.erase(it_inv);
    }
    histoInv.erase(1);
    for (auto c : histoInv[0]) {
      histo.erase(c);
    }
    histoInv.erase(0);

    bool flag_forced = !toprocess.empty();

    if (toprocess.empty()) {
      if (histo.empty()) break;
      auto x = *(histoInv.rbegin()->second.begin());
      histoInv.rbegin()->second.erase(x);
      if (histoInv.rbegin()->second.empty()) histoInv.erase(histoInv.rbegin()->first);
      toprocess.emplace(x);
    }

    for (auto x : toprocess) {
      finalA[x];
      myorder.emplace_back(x, flag_forced);
      auto it_inv = mapInv.find(x);
      if (it_inv == mapInv.end()) continue;
      ++cpt;
      for (auto i : it_inv->second) {
        auto it_imp = mapImp.find(i);
        for (auto y : it_imp->second) {
          auto & hy = histo[y];
          histoInv[hy].erase(y);
          if (histoInv[hy].empty()) histoInv.erase(hy);
          hy -= 1;
          histoInv[hy].emplace(y);
          if (x != y) mapInv[y].erase(i);
        }
        mapImp.erase(it_imp);
      }
      mapInv.erase(it_inv);
    }
    //cout << " >= " << mymax << endl;
  }
  cout << endl;

  cout << "cpt: " << cpt << endl;

  for (unsigned i_imp = 0; i_imp < imp.size(); ++i_imp) {
    auto const & p = imp[i_imp];
    for (unsigned i = 0; i < n; ++i) {
      if (((1ull << i) & p.first) != 0) {
        uint64_t x = ((1ull << i) ^ p.first) | (p.second << n);
        auto it = finalA.find(x);
        if (it != finalA.end()) {
          it->second |= (1ull << i);
        }
      }
    }
    for (unsigned i = 0; i < n; ++i) {
      if (((1ull << i) & ~p.second) != 0) {
        uint64_t x = p.first | ((p.second | (1ull << i))  << n);
        auto it = finalA.find(x);
        if (it != finalA.end()) {
          it->second |= (1ull << (i+n));
        }
      }
    }
  }

  for (auto q : myorder) {
    auto const & p = *finalA.find(q.first);
    cout << q.second << ": ";
    auto coef = __builtin_popcountll(p.second);
    bool plus = false;
    cout << coef << "*(";
    for (unsigned i = 0; i < n; ++i) {
      if (((p.first >> i) & 1) != 0) {
        if (plus) cout << " + ";
        else plus = true;
        cout << "1 - u_" << i;
      }
    }
    for (unsigned i = 0; i < n; ++i) {
      if (((p.first >> (i+n)) & 1) == 0) {
        if (plus) cout << " + ";
        else plus = true;
        cout << "v_" << i;
      }
    }
    cout << ")";
    for (unsigned i = 0; i < n; ++i) {
      if (((p.second >> i) & 1) != 0) {
        if (plus) cout << " + ";
        else plus = true;
        cout << "1 - u_" << i;
      }
    }
    for (unsigned i = 0; i < n; ++i) {
      if (((p.second >> (i+n)) & 1) != 0) {
        if (plus) cout << " + ";
        else plus = true;
        cout << "v_" << i;
      }
    }
    cout << " >= " << coef << endl;
  }

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
