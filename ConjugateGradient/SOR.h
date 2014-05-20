#ifndef SOR_H
#define SOR_H
#include "Matrix.h"
#include <vector>
#define SOR_EPS 1e-9
#define SOR_OMEGA 1.5 // default omega 1.5
bool SORisfinish(std::vector <double>&x ,std::vector <double> &x_){
  for(int i = 0;i < x.size();i++){
    if(std::abs(x[i] -x_[i]) > SOR_EPS){
      return false;
    }
  }
  return true;
}

std::vector <double> SOR_method(Matrix &A,std::vector <double> &b,bool force){
  if(!A.isSuperdiagonalAngle()){
    if(force){
    std::cerr << "Warning Matrix is not superdiagonal angle\n";
    }else{
      std::cerr << "Matrix is not superdiagonal angle\n";
      exit(-1);
    }
  }
  std::vector <double> x(A.column,0);
  std::vector <double> x_(A.column);
  
  do{
    std::swap(x_,x);
    for(int i = 0;i < x.size();i++){
      double next_xi = b[i];
      for(int j = 0;j < i;j++){
        next_xi -= A.matrix[i][j] * x[j];
      }
      for(int j = i+1;j < x.size();j++){
        next_xi -= A.matrix[i][j] * x_[j];
      }
      x[i] = next_xi / A.matrix[i][i];
      x[i] = x_[i] + SOR_OMEGA * (x[i] - x_[i]);
    }
  }while(!SORisfinish(x_,x));
  return x;
}
std::vector <double> SOR_method(Matrix &A,std::vector <double> &b){
  return SOR_method(A,b,false);
}
#endif
