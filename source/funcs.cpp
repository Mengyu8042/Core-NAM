// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/SVD>
using namespace Eigen;


// [[Rcpp::export]]
SEXP eigenMultSolveSp(MatrixXd X, SparseMatrix<double> XXt, VectorXd Y){
  MatrixXd XtX = XXt * X;
  VectorXd XtY = XXt * Y;
  VectorXd C = XtX.completeOrthogonalDecomposition().solve(XtY);
  
  return Rcpp::wrap(C);
}


// [[Rcpp::export]]
SEXP eigenTraceHatSp(MatrixXd X, SparseMatrix<double> XXt, int n){
  MatrixXd XtX = XXt * X;
  MatrixXd W = XtX.completeOrthogonalDecomposition().pseudoInverse() * XXt;
  double trA = 0;
  for (int i = 0; i < n; ++i) {
    VectorXd Xi = X.row(i);
    double Ai = Xi.dot(W.col(i));
    trA += Ai;
  }
  
  return Rcpp::wrap(trA);
}
