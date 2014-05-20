#ifndef QR_H
#define QR_H
#include <utility>
#include <vector>
#include <algorithm>
#include <cmath>
#include "Matrix.h"
#define QR_EPS 1e-9
#define DEBUG std::cerr << "DEBUG" << __LINE__ << "\n";

bool isfinish(std::vector <double> r,std::vector <double> r_){
  for(int i = 0; i < r.size();i++){
    if(std::abs(r[i]-r_[i]) > QR_EPS){
      return false;
    }
  }
  return true;
}

std::pair <Matrix,Matrix> QR_decomp(Matrix& A){//Householder transformation
  if(!(A.row > A.column)){
    std::cerr << "row size <= column size\n";
    exit(-1);
  }
  
  Matrix R(A.row,A.column),Q(A.row,A.row);
  R = A;
  Q.Identity();
  for(int k = 0; k < A.column;k++){
    std::vector <double> x(R.row-k);
    std::vector <double> u(R.row-k);
    std::vector <double> ut_R(R.column-k,0);
    std::vector <double> Q_u(Q.row,0);
    double x_len = 0,u_len = 0;
    for(int i = k; i < R.row;i++){
      x[i-k] = R.matrix[i][k];
      x_len += x[i-k] * x[i-k];
    }
    x_len = std::sqrt(x_len);
    u[0] = x[0] + x_len;
    u_len = u[0] * u[0];
    for(int i = 1;i < R.row-k;i++){
      u[i] = x[i];//x - y y = 0(i != 0)
      u_len += u[i] * u[i];
    }

    for(int i = k;i < R.row;i++){
      for(int j = k; j < R.column;j++){
        ut_R[j-k] += R.matrix[i][j] * u[i-k];
      }
    }
    
    for(int i = k;i < R.row;i++){
      for(int j = k;j < R.column;j++){
        R.matrix[i][j] = R.matrix[i][j] - 2 * u[i-k] * ut_R[j-k] / u_len;
      }
    }
    for(int i = 0;i < Q.row;i++){
      for(int j = k; j < Q.column;j++){
        Q_u[i] += Q.matrix[i][j] * u[j-k];
      }
    }
    for(int i = 0;i < Q.row;i++){
      for(int j = k; j < Q.column;j++){
        Q.matrix[i][j] = Q.matrix[i][j] - 2 * Q_u[i] * u[j-k] / u_len;
      }
    }
  }
 
  return std::make_pair(Q,R);
}

std::vector <double> QR_method(Matrix A){
  std::vector <double> r(A.row);
  for(int i = 0;i < A.row;i++){
    r[i] = A.matrix[i][i];
  }
  std::vector <double> r_(A.row);
  
  do{
    std::swap(r_,r);
    std::pair <Matrix,Matrix> QR = QR_decomp(A);
    A = QR.second * QR.first;
    for(int i = 0;i < A.row;i++){
      r[i] = A.matrix[i][i];
    }
  }while(!isfinish(r,r_));
  return r;
}
#endif
