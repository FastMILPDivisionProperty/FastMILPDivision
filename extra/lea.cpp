#include <vector>
#include <cstdint>


#include "/opt/gurobi902/linux64/include/gurobi_c++.h"

void modularAddition(GRBModel & model, vector<GRBVar> & x, vector<GRBVar> & y, vector<GRBVar> & z) { // add inequalities z = x + y mod 2^n
  unsigned n = x.size();
  vector<GRBVar> a (n), b (n), c (n), d (n), e (n);
  for (unsigned i = 0; i < n; ++i) {
    a[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
    b[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
    c[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
    d[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
    e[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
  }
  for (unsigned i = 0; i < n; ++i) model.addConstr(a[i] + 2*b[i] == x[i] + y[i]);
  for (unsigned i = 0; i < n; ++i) model.addConstr(z[i] == c[i] + d[i] + e[i]);
  for (unsigned i = 0; i < n; ++i) {
    GRBLinExpr le = 0;
    for (unsigned j = i + 1; j < n; ++j) le += e[j] - b[j];
    model.addConstr(le >= b[i]);
  }
  for (unsigned i = 0; i < n; ++i) {
    GRBLinExpr le = 0;
    for (unsigned j = i + 1; j < n; ++j) le += e[j] - b[j] + d[j];
    model.addConstr(le >= a[i] + b[i] - c[i]);
  }
}

void testLEA(unsigned R) {
  GRBEnv env = GRBEnv(true);
  env.set("LogFile", "mip1.log");
  env.start();

  GRBModel model = GRBModel(env);
  model.set(GRB_IntParam_OutputFlag, 0);

  vector<vector<GRBVar>> x (R+1, vector<GRBVar> (128));
  vector<vector<GRBVar>> y0 (R, vector<GRBVar> (128));
  vector<vector<GRBVar>> y1 (R, vector<GRBVar> (128));
  vector<vector<GRBVar>> z (R, vector<GRBVar> (96));
  vector<vector<GRBVar>> h (R, vector<GRBVar> (96));

  for (unsigned r = 0; r < R; ++r) {
    for (unsigned b = 0; b < 128; ++b) x[r][b] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
    for (unsigned b = 0; b < 128; ++b) y0[r][b] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
    for (unsigned b = 0; b < 128; ++b) y1[r][b] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
    for (unsigned b = 0; b < 96; ++b) h[r][b] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);

    for (unsigned i = 0; i < 3; ++i) {
      for (unsigned b = 0; b < 32; ++b) model.addConstr(y0[r][32*i + b] + y1[r][32*i + b] >= x[r][32*i+b]);
    }
    for (unsigned b = 0; b < 32; ++b) {
      model.addConstr(y0[r][32*3 + b] == x[r][32*3+b]);
      model.addConstr(y1[r][32*3 + b] == 0);
    }

    for (unsigned i = 0; i < 3; ++i) {
      vector<GRBVar> in1 (32), in2 (32), out (32);
      for (unsigned b = 0; b < 32; ++b) in1[b] = y1[r][32*i + b];
      for (unsigned b = 0; b < 32; ++b) in2[b] = y0[r][32*(i+1) + b];
      for (unsigned b = 0; b < 32; ++b) out[b] = h[r][32*i + b];
      modularAddition(model, in1, in2, out);
    }
  }
  for (unsigned b = 0; b < 128; ++b) x[R][b] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
  for (unsigned r = 0; r < R; ++r) {
    for (unsigned b = 0; b < 32; ++b) model.addConstr(x[r+1][96+b] == y0[r][b]);
    for (unsigned b = 0; b < 32; ++b) model.addConstr(x[r+1][64+b] == h[r][64+ ((b+3)%32)]);
    for (unsigned b = 0; b < 32; ++b) model.addConstr(x[r+1][32+b] == h[r][32+ ((b+5)%32)]);
    for (unsigned b = 0; b < 32; ++b) model.addConstr(x[r+1][((b+9)%32)] == h[r][b]);
  }

  GRBLinExpr e_out = 0;
  for (unsigned b = 0; b < 128; ++b) e_out += x[R][b];
  model.addConstr(e_out == 1);

  map<unsigned, vector<vector<int>>> imp;

  vector<int> tab_in = {0,0,0,0};
  int data = 127;
  while (data > 0) {
    bool flag = true;
    GRBLinExpr e_in = 0;
    for (unsigned b = 0; b < 128; ++b) e_in += x[0][b];
    auto c0 = model.addConstr(e_in == data);
    for (tab_in[0] = 32; tab_in[0] >= 0 && flag; --tab_in[0]) {
      for (tab_in[1] = 32; tab_in[1] >= 0 && flag; --tab_in[1]) {
        for (tab_in[2] = 32; tab_in[2] >= 0 && flag; --tab_in[2]) {
          for (tab_in[3] = 32; tab_in[3] >= 0 && flag; --tab_in[3]) {
            if (tab_in[0] + tab_in[1] + tab_in[2] + tab_in[3] == data) {
              GRBLinExpr e = 0;
              for (unsigned i = 0; i < 4; ++i) {
                for (int b = 0; b < 32 - tab_in[i]; ++b) e += x[0][b + 32*i];
              }
              auto c1 = model.addConstr(e == 0);
              for (unsigned b = 0; b < 128; ++b) {
                auto & all = imp[b];
                if (any_of(all.begin(), all.end(), [&tab_in](auto & v){return tab_in[0] <= v[0] && tab_in[1] <= v[1] && tab_in[2] <= v[2] && tab_in[3] <= v[3];})) continue;
                auto c2 = model.addConstr(x[R][b] == 1);
                model.optimize();
                if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
                  if (flag) cout << "(" << tab_in[0] << "," << tab_in[1] << "," << tab_in[2] << ", " << tab_in[3] << "): ";
                  flag = false;
                  cout << " " << b;
                }
                else all.emplace_back(tab_in);
                model.remove(c2);
              }
              if (!flag) cout << endl;
              model.remove(c1);
            }
          }
        }
      }
    }
    if (flag) break;
    model.remove(c0);
    --data;
  }



  // for (unsigned in = 0; in < 128; ++in) {
  //   if ((in <= 4 && in >= 0) || (in <= 36 && in >= 32)) model.addConstr(x[0][in] == 0);
  //   else model.addConstr(x[0][in] == 1);
  // }
  //
  // for (unsigned b = 0; b < 128; ++b) {
  //   if (b == 59) model.addConstr(x[R][b] == 1);
  //   else model.addConstr(x[R][b] == 0);
  // }
  //
  // model.optimize();
  //
  // if (model.get(GRB_IntAttr_Status) != GRB_INFEASIBLE) {
  //   for (unsigned r = 0; r <= R; ++r) {
  //     for (unsigned b = 0; b < 128; ++b) {
  //       if (b != 0 && b%32 == 0) cout << " | ";
  //       cout << (int) (x[r][b].get(GRB_DoubleAttr_X) + 0.5) << " ";
  //     }
  //     cout << endl;
  //   }
  //   cout << " --------------- " << endl;
  //   for (unsigned r = 0; r < R; ++r) {
  //     for (unsigned i = 0; i < 3; ++i) {
  //       cout << "x: ";
  //       for (unsigned b = 0; b < 32; ++b) {
  //         if (b != 0 && b%32 == 0) cout << " | ";
  //         cout << (int) (y1[r][32*i + b].get(GRB_DoubleAttr_X) + 0.5) << " ";
  //       }
  //       cout << endl;
  //       cout << "y: ";
  //       for (unsigned b = 0; b < 32; ++b) {
  //         if (b != 0 && b%32 == 0) cout << " | ";
  //         cout << (int) (y0[r][32*(i+1) + b].get(GRB_DoubleAttr_X) + 0.5) << " ";
  //       }
  //       cout << endl;
  //       cout << "z: ";
  //       for (unsigned b = 0; b < 32; ++b) {
  //         if (b != 0 && b%32 == 0) cout << " | ";
  //         cout << (int) (h[r][32*i + b].get(GRB_DoubleAttr_X) + 0.5) << " ";
  //       }
  //       cout << endl;
  //       cout << " -------------- " << endl;
  //     }
  //   }
  // }


}

int main(int argc, char const *argv[]) {

  testLEA(8);

  return 0;
}
