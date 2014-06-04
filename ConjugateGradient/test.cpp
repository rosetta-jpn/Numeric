#include <bits/stdc++.h>
#include "Matrix.h"
#include "LU.h"
#include "GaussianElimination.h"
#include "CG.h"
#include "Jacobi.h"
#include "SOR.h"
#include "GMRES.h"
using namespace std;

int main(){
  cin.tie(0);
  ios_base::sync_with_stdio(0);
  Matrix A({{3,1,1},{1,3,1},{1,1,3}});
  Matrix b({{0},{4},{6}});
  auto x = GMRES_method(A,b);
  cout << x;
  cout << "\n";
  return 0;
}

