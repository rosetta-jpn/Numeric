#include <bits/stdc++.h>
#include "Matrix.h"
#include "QR.h"

using namespace std;

typedef long long ll;


int main(){
  cin.tie(0);
  ios_base::sync_with_stdio(0);

  Matrix A({{1,2,3},{4,5,6},{7,8,9}});
  auto QR = QR_decomp(A);//(Q,R)
  cout << QR.first;
  cout << endl;
  cout << QR.second;
  cout << endl;
  cout << QR.first * QR.second << endl;
  Matrix B({{1,2,3},{4,5,6},{7,8,9}});
  auto lam = QR_method(B);
  for(auto x:lam){
    cout << x << " ";
  }
  cout << endl;
  
  return 0;
}

