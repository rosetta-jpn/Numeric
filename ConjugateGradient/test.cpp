#include <bits/stdc++.h>
#include "Matrix.h"
#include "LU.h"
#include "GaussianElimination.h"
#include "CG.h"
#include "Jacobi.h"
#include "SOR.h"
using namespace std;

int main(){
  cin.tie(0);
  ios_base::sync_with_stdio(0);
  Matrix A({{3,1,1},{1,3,1},{1,1,3}});
  vector <double> b = {0,4,6};
  auto x = SOR_method(A,b);
  for(auto a :x){
    cout << a << " ";
  }
  cout << "\n";
  return 0;
}

