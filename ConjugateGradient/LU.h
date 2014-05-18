#ifndef LU_H
#define LU_H
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>
#include "Matrix.h"

std::pair <Matrix,Matrix> LU_decomp(Matrix A){//no pivot selection
  if(!A.isSquare()){
    std::cerr << "Matrix row size != Matrix column size\n";
    std::cerr << A;
    exit(-1);
  }
  Matrix L(A.row,A.column),U(A.row,A.column);
  for(int k = 0;k < A.row;k++){
    L.matrix[k][k] = 1;
    double u = U.matrix[k][k] = A.matrix[k][k];
    for(int i = k+1;i < A.row;i++){
      L.matrix[i][k] = A.matrix[i][k] / u;
      U.matrix[i][k] = 0;
    }
    for(int j = k+1; j < A.column;j++){
      L.matrix[k][j] = 0;
      U.matrix[k][j] = A.matrix[k][j];
    }
    
    for(int i = k+1;i < A.row;i++){
      for(int j = k+1;j < A.column;j++){
        A.matrix[i][j] -= L.matrix[i][k] * U.matrix[k][j];
      }
    }
  }
  return std::make_pair(L,U);
}

std::pair <std::pair <Matrix,Matrix>,std::vector <int> > LU_pivotdecomp(Matrix A){//pivot selection
  if(!A.isSquare()){
    std::cerr << "Matrix row size != Matrix column size\n";
    std::cerr << A;
    exit(-1);
  }
  
  Matrix L(A.row,A.column),U(A.row,A.column);
  std::vector <int> b_index(A.column);
  for(int k = 0;k < A.column;k++){
    b_index[k] = k;
  }
  
  for(int k = 0; k < A.row;k++){
    
    double u_max = A.matrix[k][k];
    int u_maxindex = k;
    
    for(int i = k+1;i < A.row;i++){
      if(u_max < A[i][k]){
        u_max = A.matrix[i][k];
        u_maxindex = i;
      }
    }
    std::swap(b_index[k],b_index[u_maxindex]);
    std::swap(A.matrix[k],A.matrix[u_maxindex]);
    std::swap(L.matrix[k],L.matrix[u_maxindex]);
    L.matrix[k][k] = 1;
    double u = U.matrix[k][k] = A.matrix[k][k];
    for(int i = k+1;i < A.row;i++){
      L.matrix[i][k] = A.matrix[i][k] / u;
    }
    for(int j = k+1;j < A.column;j++){
      U.matrix[k][j] = A.matrix[k][j];
    }
    for(int i = k+1; i < A.row;i++){
      for(int j = k+1;j < A.column;j++){
        A.matrix[i][j] -= L.matrix[i][k] * U.matrix[k][j];
      }
    }
  }
  return std::make_pair(std::make_pair(L,U),b_index);
}


std::vector <double> LU_solver(Matrix A,std::vector <double> b,bool pivot_flag){
  if(A.row != b.size()){
    std::cerr << "Matrix row size != vector size\n";
    std::cerr << A;
    exit(-1);
  }else if(!A.isSquare()){
    std::cerr << "Matrix row size != Matrix column size\n";
    std::cerr << A;
    exit(-1);
  }
  
  std::vector <double> x(A.column);
  std::vector <double> y(A.column);
  Matrix L(A.row,A.column),U(A.row,A.column);
  if(pivot_flag){
    std::pair <std::pair <Matrix,Matrix>,std::vector <int> > LUb = LU_pivotdecomp(A);
    
    L = LUb.first.first;
    U = LUb.first.second;
    std::vector <double> newb(b.size());
    for(int i = 0;i < b.size();i++){
      newb[i] = b[LUb.second[i]];
    }
    std::swap(newb,b);

  }else{
    std::pair <Matrix,Matrix> LU = LU_decomp(A);
    L = LU.first;
    U = LU.second;
  }

  for(int i = 0; i < A.row;i++){
    double other = 0;
    for(int j = 0;j < i;j++){
      other += L.matrix[i][j] * y[j];
    }
    b[i] -= other;
    if(L.matrix[i][i] == 0){
      std::cerr << "LU solver error, zero division\n";
      std::cerr << L << "\n";
      std::cerr << U << "\n";
      exit(-1);
    }
    y[i] = b[i] / L.matrix[i][i];
  }

  for(int i = A.row-1;i >= 0;i--){
    double other = 0;
    for(int j = i+1;j < A.column;j++){
      other += U.matrix[i][j] * x[j];
    }
    y[i] -= other;
    
    if(U.matrix[i][i] == 0){
      std::cerr << "LU solver error, zero division\n";
      std::cerr << L << "\n" << U << "\n";
      exit(-1);
    }
    x[i] = y[i] / U.matrix[i][i];
  }
  return x;
}








#endif
