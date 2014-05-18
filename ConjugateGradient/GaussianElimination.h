#ifndef GAUSSIANELIMINATION_H
#define GAUSSIANELIMINATION_H
#include <iostream>
#include <vector>
#include <algorithm>
#include "Matrix.h"

std::vector <double> GaussianElimination(Matrix A,std::vector <double> b){//pivot selection
  if(!A.isSquare()){
    std::cerr << "row size != column size\n";
    std::cerr << A;
    exit(-1);
  }

  for(int k = 0; k < A.row;k++){
    double u_max = A.matrix[k][k];
    double u_maxindex = k;
    for(int i = k+1; i < A.row;i++){
      if(u_max < A.matrix[i][k]){
        u_max = A.matrix[i][k];
        u_maxindex = i;
      }
    }
    std::swap(A.matrix[k],A.matrix[u_maxindex]);
    std::swap(b[k],b[u_maxindex]);
    if(u_max == 0){
      //Ax=b no solution
      std::cerr << "Ax=b no solution";
      exit(-1);
    }
    double u = A.matrix[k][k];
    A.matrix[k][k] = 1;
    for(int j = k+1;j < A.column;j++){
      A.matrix[k][j] = A.matrix[k][j] / u;
    }
    b[k] /= u;

    for(int i = k+1;i < A.row;i++){
      double u = A.matrix[i][k];
      A.matrix[i][k] = 0;
      for(int j = k+1;j < A.column;j++){
        A.matrix[i][j] -= u * A.matrix[k][j];
      }
      b[i] -= u * b[k];
    }
  }
  std::vector <double> x(A.column);
  for(int i = A.row-1;i>=0;i--){
    double other = 0;
    for(int j = i+1;j < A.column;j++){
      other += A.matrix[i][j] * x[j];
    }
    b[i] -= other;
    x[i] = b[i] / A.matrix[i][i];
  }
  return x;
}


#endif
