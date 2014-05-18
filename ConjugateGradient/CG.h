#ifndef CG_H
#define CG_H
#include "Matrix.h"
#include <iostream>
#include <vector>
#define CG_EPS 1e-8

bool cg_isfinish(Matrix &r){
  for(int i = 0 ; i < r.row;i++){
    if(r.matrix[i][0] > CG_EPS){
      return false;
    }
  }
  return true;
}

std::vector <double> ConjugateGradient(Matrix &A,std::vector <double> &_b,bool force){
  if(!A.isSquare()){
    std::cerr << "A is not Sqaure\n";
    std::cerr << A;
    exit(-1);
  }else if(!A.isSymmetric()){
    std::cerr << "A is not Symmetric\n";
    std::cerr << A;
    exit(-1);
  }
  //Definiteness judged by superdiagonal angle
  if(!A.isSuperdiagonalAngle()){
    if(force){
    std::cerr << "Warning Matrix can't be judged Definiteness by superdiagonal angle\n";
    }else{
      std::cerr << "Matrix can't be judged Definiteness by superdiagonal angle\n";
      exit(-1);
    }
  }


  Matrix b(A.row,1);//b : vector -> Matrix
  for(int i = 0 ;i < A.row;i++){
    b.matrix[i][0] = _b[i];
  }
  Matrix x(A.row,1);//init x is 0
  Matrix r(A.row,1); Matrix r_(A.row,1); // r,r_ swaps
  Matrix p(A.row,1);
  
  p = r = b - A*x;

  double alpha;
  double beta;

  while(!cg_isfinish(r)){
    alpha = (r.t() * r).matrix[0][0] / (p.t() * A * p).matrix[0][0];
    r_ = r - (alpha * A * p);
    x = x + alpha * p;
    beta = (r_.t() * r_).matrix[0][0] / (r.t() * r).matrix[0][0];
    p = r_ + beta * p;
    std::swap(r_,r);
  }    
  // x : Matrix -> vector
  std::vector <double> x_(x.row);
  for(int i = 0; i < x.row;i++){
    x_[i] = x.matrix[i][0];
  }
  return x_;
}

std::vector <double> ConjugateGradient(Matrix &A,std::vector <double> &b){
  return ConjugateGradient(A,b,false);
}





#endif
