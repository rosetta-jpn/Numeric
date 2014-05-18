#ifndef JACOBI_H
#define JACOBI_H
#include "Matrix.h"
#include <vector>
#include <algorithm>
#include <iostream>
#define JACOBI_EPS 1e-9
bool jacobi_isfinish(std::vector <double> x1,std::vector <double> x2){
  for(int i = 0; i < x1.size();i++){
    if(std::abs(x1[i] -x2[i]) > JACOBI_EPS){
      return false;
    }
  }
  return true;
}

std::vector <double> Jacobi_method(Matrix &A,std::vector <double> b,bool force){
if(!A.isSuperdiagonalAngle()){
    if(force){
    std::cerr << "Warning Matrix is not superdiagonal angle\n";
    }else{
      std::cerr << "Matrix is not superdiagonal angle\n";
      exit(-1);
    }
  }
  
  std::vector <double> x(A.row,0);
  std::vector <double> x_(A.row,0);
  for(int i = 0;i < A.row;i++){
    double ax = 0;
    for(int j = 0; j < A.column;j++){
      if(i != j){
        ax += A.matrix[i][j] * x_[j];
      }
    }
    x[i] = (b[i] - ax) / A.matrix[i][i];
  }
  
  while(!jacobi_isfinish(x,x_)){
    std::swap(x,x_);
    for(int i = 0;i < A.row;i++){
      double ax = 0;
      for(int j = 0; j < A.column;j++){
        if(i != j){
          ax += A.matrix[i][j] * x_[j];
        }
      }
      x[i] = (b[i] - ax) / A.matrix[i][i];
    }
  }
  return x;
}
std::vector <double> Jacobi_method(Matrix &A,std::vector <double> b){
  return Jacobi_method(A,b,false);
}
#endif
