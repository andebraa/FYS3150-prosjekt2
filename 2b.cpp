//#include "Functions.h"
#include <iostream>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;

//Definere funk
double max_offdiag(mat, int&, int&);

void rotation(mat, mat, int, int, int);

/*
void unit_test(mat a, int n){
       mat g = a.t() * a.i();
       cout<<g;
}
*/



int main(int argc, char* argv[]){
  int n = atoi(argv[1]);

  mat A = mat(n,n);
  double h = 1.0/((double) n);
  double hh = h*h;
  mat R(n,n,fill::eye);

  A.fill(0.0);
  for (int i=0; i< n-1; i++){
    A(i,i) = 2/hh;
    A(i,i+1) = -1/hh;
    A(i+1,i) = -1/hh;
  }
  A(n-1,n-1) = 2/hh;

  int k, l;
  double max_value;

  max_value = max_offdiag(A, k, l);

  double eps = 1e-8;
  int max_iterat = 1e6;
  int iterat = 0;

  while (max_value > eps && iterat < max_iterat){
    max_value = max_offdiag(A, k, l);
    rotation(A, R, k, l, n);

    iterat = 1 + iterat;
  }

  A.print("A = ");
  return 0;

}


//Fylle funksjonene
double max_offdiag(mat A, int &k, int &l){
  mat B = A;

  B.diag().zeros();

  double max_value;
  uvec max_index;


  if (abs(B.max()) < abs(B.min()) ){
    max_value = abs(B.min());
    max_index = ind2sub(size(B),B.index_min());

  }
  else{
    max_value = abs(B.max());
    max_index = ind2sub(size(B),B.index_max());
  }

  k = max_index[0];
  l = max_index[1];

  return max_value;
}


void rotation(mat A,mat R, int k, int l, int n){
  //Her er indekseringen basert på at jeg lager en matrise med for-loop

  double t, tau, s, c;

  //Bestemmer sin og cos
  if (A(k,l) != 0.0){
  tau = (A(l,l) - A(k,k))/(2*A(k,l));
  if (tau > 0){
    t = 1.0/(tau + sqrt(1 + tau*tau));
  }
  else{
    t = -1.0/(tau - sqrt(1 + tau*tau));
  }
  c = 1.0/sqrt(1+t*t);
  s = t*c;
}
  else{
    c = 1.0;
    s = 0;
  }

    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    //Endrer verdien til den største ikke-diag element
    a_kk = A(k,k);
    a_ll = A(l,l);
    A(k,k) = a_kk*c*c - 2*A(k,l)*c*s + a_ll*s*s;
    A(l,l) = a_ll*c*c + 2*A(k,l)*c*s + a_kk*s*s;

    A(k,l) = 0;
    A(l,k) = 0;
    //Endrer de resterende elementene
    for (int i = 0; i < n; i++){
      if (i != l && i != k){
        a_ik = A(i,k);
        a_il = A(i,l);
        A(i,k) = a_ik*c - a_il*s;
        A(k,i) = A(i,k);
        A(i,l) = a_il*c + a_ik*s;
        A(l,i) = A(i,l);
      }
      r_ik = R(i,k);
      r_il = R(i,k);
      R(i,k) = c*r_ik - s*r_il;
      R(i,l) = c*r_il + s*r_ik;
    }

}
