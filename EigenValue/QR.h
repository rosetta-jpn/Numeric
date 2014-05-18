#ifndef QR_H
#define QR_H
#include <utility>
#include <vector>
#include <algorithm>
#include "Matrix.h"
#define QR_EPS 1e-9
bool isfinish(std::vector <double> r,std::vector <double> r_){
  for(int i = 0; i < r.size();i++){
    if(std::abs(r[i]-r_[i]) > QR_EPS){
      return false;
    }
  }
  return true;
}


std::pair <Matrix,Matrix> QR_decomp(Matrix A){//Householder transformation
  if(!(A.row > A.column)){
    std::cerr << "row size <= column size\n";
    exit(-1);
  }
  //still
}

std::vector <double> QR_method(Matrix A){
  std::vector <double> r(A.row);
  for(int i = 0;i < A.row;i++){
    r[i] = A.matrix[i][i];
  }
  std::vector <double> r_(A.row);
  std::pair <Matrix,Matrix> QR = QR_decomp(A);
  A = QR.second * QR.first;
  for(int i = 0;i < A.row;i++){
    r_[i] = A.matrix[i][i];
  }
  
  while(!isfinish(r,r_)){
    std::swap(r_,r);
    QR = QR_decomp(A);
    A = QR.second * QR.first;
    for(int i = 0;i < A.row;i++){
      r[i] = A.matrix[i][i];
    }
  }
  return r;
}
#endif
