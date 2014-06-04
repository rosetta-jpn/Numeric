#ifndef GMRES_H
#define GMRES_H
#define GMRES_EPS 1e-12

Matrix Arnoldi_method(Matrix &A, Matrix &b, int n){// A(m x m) , b (m x 1) => (m x n)
  Matrix ret(A.row,n);
  Matrix omega(A.row,1);
  Matrix v(b.row,b.column);
  double h = 0;
  v = b;
  for(int j = 0; j < n;j++){
    omega = A * v;
    for(int i = 0; i < j;i++){
      h = (v.t() * omega).matrix[0][0];
      omega -= h * v;
    }
    h = omega.norm2();
    if(std::abs(h) < GMRES_EPS){
      break;
    }
    for(int i = 0;i < b.row;i++){
      ret.matrix[i][j] = v.matrix[i][0];
    }
    v = omega / h;
  }
  return ret;
}


Matrix GMRES_method(Matrix &A, Matrix &b){
  if(std::abs(b.norm2() - 1) > GMRES_EPS){ //to regular
    double h = 0;
    for(int i = 0; i < b.row;i++){
      h += b.matrix[i][0] * b.matrix[i][0];
    }
    b /= h;
    A /= h;
  }
  Matrix ret(A.row,1);
  return ret;
}

#endif
