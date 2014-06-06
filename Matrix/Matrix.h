#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <algorithm>
#include <vector>
#include <string>

class Matrix{
  friend std::ostream& operator<<(std::ostream& os,const Matrix& A);
  friend std::istream& operator<<(std::istream& os,const Matrix& A);
 public:
  std::vector <std::vector <double> > matrix;
  int column;
  int row;
  Matrix();
  Matrix(int n,int m){
    matrix = std::vector <std::vector <double> > (n,std::vector <double>(m,0));
    row = n;
    column = m;
  }
  Matrix(std::vector <std::vector <double> > A){
    matrix = A;
    row = A.size();
    column = A[0].size();
  }
 
  void Identity(){
    for(int i = 0; i < std::max(this->row,this->column);i++){
      this->matrix[i][i] = 1;
    }
  }

  Matrix operator+(const Matrix &A){
    int n = this->row;
    int m = this->column;
    if(!isMatrix_nxm(A,n,m)){
      Matrix_calc_error(*this,A,"sum");
    }
    Matrix ret(n,m);
    for(int i = 0; i < n ; i++){
      for(int j = 0;j < m;j++){
        ret.matrix[i][j] = this->matrix[i][j] + A.matrix[i][j];
      }
    }
    return ret;
  }
  
  Matrix& operator+=(const Matrix &A){
    int n = this->row;
    int m = this->column;
    if(!isMatrix_nxm(A,n,m)){
      Matrix_calc_error(*this,A,"sum");
    }
    for(int i = 0; i < n ; i++){
      for(int j = 0;j < m;j++){
        this->matrix[i][j] += A.matrix[i][j];
      }
    }
    return *this;
  }

  Matrix operator-(const Matrix &A){
    int n = this->row;
    int m = this->column;
    Matrix ret(n,m);
    if(!isMatrix_nxm(A,n,m)){
      Matrix_calc_error(*this,A,"sub");
    }
    for(int i = 0; i < n ;i++){
      for(int j = 0; j < m ;j++){
        ret.matrix[i][j] = this->matrix[i][j] - A.matrix[i][j];
      }
    }
    return ret;
  }

  Matrix operator-=(const Matrix &A){
    int n = this->row;
    int m = this->column;
    Matrix ret(n,m);
    if(!isMatrix_nxm(A,n,m)){
      Matrix_calc_error(*this,A,"sub");
    }
    for(int i = 0;i < n;i++){
      for(int j = 0; j < m;j++){
        this->matrix[i][j] -= A.matrix[i][j];
      }
    }
    return *this;
  }
  
  Matrix operator*(Matrix &A){
    int n = this->row;
    int m = this->column;
    int l = A.column;
    if(this->column != A.row){
      Matrix_calc_error(*this,A,"mul");
    }
    Matrix A_t(A.column,A.row);
    A_t = A.t();
    Matrix ret(n,l);
    for(int i = 0; i < n;i++){
      for(int j = 0; j < l;j++){
        for(int k = 0;k < m;k++){
          ret.matrix[i][j] += this->matrix[i][k] * A_t.matrix[j][k];
        }
      }
    }
    return ret;
  }

  Matrix operator*(std::vector <double> &b){
    if(b.size() != this->column){
      Matrix_calc_error(*this,"vector mul");
    }
    Matrix ret(this->row,1);
    for(int i = 0; i < this->row;i++){
      for(int j = 0; j < this->column;j++){
        ret.matrix[i][0] += this->matrix[i][j] * b[j];
      }
    }
    return ret;
  }


  Matrix operator*=(Matrix &A){
    int n = this->row;
    int m = this->column;
    int l = A.column;
    if(this->column != A.row){
      Matrix_calc_error(*this,A,"mul");
    }
    Matrix A_t(A.column,A.row);
    Matrix ret(n,l);
    A_t = A.t();
    for(int i = 0; i < n;i++){
      for(int j = 0; j < l;j++){
        for(int k = 0; k < m;k++){
          ret.matrix[i][j] += this->matrix[i][k] * A_t.matrix[j][k];
        }
      }
    }
    *this = ret;
    return *this;
  }

  Matrix operator/(const double x){
    int n = this->row;
    int m = this->column;
    Matrix ret(n,m);
    for(int i = 0; i < n;i++){
      for(int j = 0;j < m;j++){
        ret.matrix[i][j] = this->matrix[i][j] / x;
      }
    }
    return ret;
  }
  
  Matrix operator/=(const double x){
    int n = this->row;
    int m = this->column;
    Matrix ret(n,m);
    for(int i = 0; i < n;i++){
      for(int j = 0; j < m;j++){
        ret.matrix[i][j] = this->matrix[i][j] / x;
      }
    }
    return ret;
  }

  Matrix operator=(const Matrix &A){
    int n = this->row;
    int m = this->column;
    if(!isMatrix_nxm(A,n,m)){
      Matrix_calc_error(*this,A,"substitue");
    }
    for(int i = 0;i < n;i++){
      for(int j = 0;j < m;j++){
        this->matrix[i][j] = A.matrix[i][j];
      }
    }
    return *this;
  }
  std::vector <double> operator[](const int n){
    return this->matrix[n];
  }
  
  Matrix t(){//trans
    int n = this->row;
    int m = this->column;
    Matrix ret(m,n);
    for(int i = 0; i < m;i++){
      for(int j = 0; j < n;j++){
        ret.matrix[i][j] = this->matrix[j][i];
      }
    }
    return ret;
  }
  
  double norm2(){//for(n x 1) matrix
    double ret = 0;
    if(!isMatrix_nxm(*this,this->row,1)){
      Matrix_calc_error(*this,"norm2");
    }
    for(int i = 0; i < this->row;i++){
      ret += this->matrix[i][0] * this->matrix[i][0];
    }
    return std::sqrt(ret);
  }

  bool isMatrix_nxm(const Matrix &A,int n,int m){//is n x m matrix
    return n == A.row && m == A.column;
  }
  
  bool isSquare(){
    return this->row == this->column;
  }
  
  bool isSymmetric(){//is Symmetirc?
    for(int i = 0; i < this->row;i++){
      for(int j = i+1;j < this->column;j++){
        if(this->matrix[i][j] != this->matrix[j][i]){
          return false;
        }
      }
    }
    return true;
  }
  
  bool isSuperdiagonalAngle(){//is Super diagonal angle?
    for(int i = 0;i < this->row;i++){
      double a = 0;
      for(int j = 0;j < this->column;j++){
        if(i != j){
          a += abs(this->matrix[i][j]);
        }
      }
      if(!(abs(this->matrix[i][i]) >= a)){
        return false;
      }
    }
    return true;
  }
  void Matrix_calc_error(const Matrix &A,const Matrix &B,std::string operate){//operate error 
    std::cerr << "Matrix " << operate << " Error\n";
    std::cerr << A;
    std::cerr << std::endl;
    std::cerr << B;
    exit(-1);
  }
   void Matrix_calc_error(const Matrix &A,std::string operate){//operate error 
    std::cerr << "Matrix " << operate << " Error\n";
    std::cerr << A;
    std::cerr << std::endl;
    exit(-1);
  }
};

std::istream& operator>>(std::istream& is,Matrix &A){
  int n = A.row;
  int m = A.column;
  for(int i = 0 ; i < n;i++){
    for(int j = 0; j < m;j++){
      is >> A.matrix[i][j];
    }
  }
  return is;
}
 

std::ostream& operator<<(std::ostream& os,const Matrix &A){
  int n = A.row;
  int m = A.column;
  for(int i = 0;i < n;i++){
    for(int j = 0; j < m;j++){
      if(j) os << " ";
      os << std::fixed << std::setprecision(10) << A.matrix[i][j];
    }
    os << std::endl;
  }
  return os;
}

Matrix operator*(double a,Matrix &A){
  int n = A.row;
  int m = A.column;
  Matrix ret(n,m);
  for(int i = 0 ; i < n;i++){
    for(int j = 0; j < m;j++){
      ret.matrix[i][j] = a * A.matrix[i][j];
    }
  }
  return ret;
}

#endif
