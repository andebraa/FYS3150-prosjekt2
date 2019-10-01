#include <iostream>
#include <cmath>
#include <armadillo>
#include <fstream>
#include <string>
#include <time.h>
#include <math.h>


using namespace std;
using namespace arma;

void d_func(unsigned int n, double rhomax);
void b_func();
void e_func(double omega);
bool maxval_test();
bool eigval_test(mat b, mat r, vec eigenval, mat eigenvec, uword n, int rhomax, double omega);
void rotate(mat &b,mat &r, uword &k, uword &l, unsigned int n);
void anal(uword n, vec &ana);
void anal_e(uword n, vec &ana);
void anal_b(uword n, vec &ana);
void findmax(mat b, uword &l, uword &k, double &maxval);



int main(){

    /*
    for(int n = 5; n<505; n+=100){
        for(int rhomax = 3; rhomax < 50; rhomax += 5){
            d_func(n,rhomax);
        }

    }
    */

//    d_func(10, 15);



    double omega_values [4] = {0.01, 0.5, 1.0, 5.0};
    for(int i =0; i<4; i++){
        e_func(omega_values[i]);
    }


    //b_func();



    //calling maxval test
    bool returnvalue_maxval = maxval_test();
    if (returnvalue_maxval == false){
        cout<<"fail in maxval_test"<<endl;
    }

}



void findmax(mat b, uword &l, uword &k, double &maxval){
    mat a = b;
    a.diag().zeros();

    uword max = a.index_max();
    uword min = a.index_min();


    if(abs(a(min)) > abs(a(max))){
        maxval= a(min);
        uvec max_ind = ind2sub(size(a), min);
        k = max_ind[0];
        l = max_ind[1];
    }
    else{
        maxval = abs(a(max));
        uvec max_ind = ind2sub(size(a), max);
        k = max_ind[0];
        l = max_ind[1];
    }

    if(k == 0 && l == 0){
        k = 1;
        l = 0;

    }

}




void b_func(){

    cout<<"--------------------"<<endl;
    uword l = 0;
    uword k = 0;

    double maxval=100;
    double eps = 1E-13;
    unsigned int n = 100;
    double rhomin = 0;
    double rhomax = 1;
    double h = (rhomax - rhomin)/n;

    double d = 2/(h*h);
    double a = -1/(h*h);

    //filling matrix
    Mat<double> b(n,n, fill::zeros);
    for(unsigned int i=0;i< n-1; i++){
        b(i,i) = d;
        b(i+1,i) =a;
        b(i,i+1) = a;

    }
    b(n-1,n-1) = d;




    int count = 0;

    //these values are used to test if future eigvals are similar
    vec eigval(n);
    mat eigvec(n,n);
    mat r = eye(n,n);
    eig_sym(eigval, eigvec, b);


    sort(eigval);


    time_t start_bfunc, finish_bfunc;

    start_bfunc = clock();

    while (count < 1000000000 && abs(maxval) > eps){

        findmax(b, l, k, maxval);


        rotate(b, r, k, l, n);
        //cout<<b<<endl<<endl;

        count++;


    }

    finish_bfunc = clock();
    double time_used = (double) (finish_bfunc-start_bfunc)/(CLOCKS_PER_SEC);

    cout<<"time used bfunc \n";
    cout<<time_used<<endl;

    cout<<"---------------------------------------------------------------- \n";
    cout<<count<<endl;
    cout<<"---------------------------------------------------------------- \n";






    //calling eigenvalue test
    double omega = 0; //Omega is not used in this part of the function, but i need it to be an input variable in eigval_test
    bool returnvalue_eigenvec = eigval_test(b, r, eigval, eigvec, n, rhomax, omega);
    if (returnvalue_eigenvec == false){
        cout<<"fail in eigval test"<<endl;
    }






    }

void d_func(unsigned int n, double rhomax){

    uword l = 0;
    uword k = 0;

    double maxval=100;
    double eps = 1E-13;

    double rho_min = 0;


    double h = (rhomax - rho_min)/n;

    double d = (2/(h*h));
    double e = -1/(h*h);

    Mat<double> G(n,n, fill::zeros);
    for(unsigned int i=0;i< n-1; i++){
        G(i,i) = d + (i+1)*(i+1)*h*h;
        G(i+1,i) =e;
        G(i,i+1) = e;

    }
    G(n-1,n-1) = d + (n)*(n)*h*h;

    int count = 0;
    mat r = eye(n,n);

    vec eigval(n);
    mat eigvec(n,n);

    time_t start_arma, stop_arma;
    start_arma = clock();

    eig_sym(eigval, eigvec, G);

    stop_arma = clock();

    double time_used_arma = (double) (stop_arma - start_arma)/(CLOCKS_PER_SEC);
    cout<<"time used arma \n";
    cout<<time_used_arma<<endl;

    while (count < 100000000 && abs(maxval) > eps){

        findmax(G, l, k, maxval);
        //cout<<G<<endl;

        rotate(G, r, k, l, n);
        //cout<<b<<endl<<endl;


        count++;
    }





    sort(eigval);


    //calling eigenvalue test
    double omega = 0; //Omega is not used in this part of the function, but i need it to be an input variable in eigval_test
    bool returnvalue_eigenvec = eigval_test(G, r, eigval, eigvec, n, rhomax, omega);
    if (returnvalue_eigenvec == false){
        cout<<"fail in eigval test"<<endl;
        }



}

void rotate(mat &b, mat &r, uword &k, uword &l, unsigned int n){
    double s, c;

    if (b(k,l) != 0.0){


        double tau = ((b(l,l) - b(k,k))/(2*b(k,l)));
        double t;
        if(tau>=0){
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        }
        else{
            t = -1.0/(-tau + sqrt( 1.0 + tau*tau));
        }
        c = 1.0/(sqrt(1 + t*t));
        s  = c*t;


    }
    else{
       return;
    }

    double b_kk = b(k,k);
    double b_ll = b(l,l);
    double b_ik, b_il, r_il, r_ik;

    b(k,k) = c*c*b_kk - 2.0*c*s*b(k,l) + s*s*b_ll;
    b(l,l) = s*s*b_kk + 2.0*c*s*b(k,l) + c*c*b_ll;
    b(k,l) = 0;
    b(l,k) = 0;

    for (uword i=0; i<n; i++){
        if( i != k && i != l){
            b_ik = b(i,k);
            b_il = b(i,l);

            b(i,k) = c*b_ik - s*b_il;
            b(k,i) = b(i,k);
            b(i,l) = c*b_il + s*b_ik;
            b(l,i) = b(i,l);
        }
        r_ik = r(i,k);
        r_il = r(i,l);
        r(i,k) = c*r_ik - s*r_il;
        r(i,l) = c*r_il + s*r_ik;


    }





}







bool maxval_test(){
    uword n = 3;
    Mat<double> g(n,n, fill::zeros);
    g(2,2) = 3;
    g(0,2) = -5;

    double maxval_t = 0;
    uword l_t =0;
    uword k_t =0;
    findmax(g, l_t, k_t, maxval_t);
    double tol = 1E-8;


    if (maxval_t - -5.0 < tol && l_t - 2.0 < tol && k_t - 0.0 < tol){
        return true;
        }
    else{
        return false;
    }
}

bool eigval_test(mat b, mat r, vec eigval, mat eigvec, uword n, int rhomax, double omega){
    //now only checks the final values, I.E, the diagonal elements in finished jacobie matrix

    vec new_eigval(n);
    mat new_eigvec(n,n);

    time_t start_arma, stop_arma;
    start_arma = clock();

    eig_sym(new_eigval, new_eigvec, b);

    stop_arma = clock();

    double time_used_arma = (double) (stop_arma - start_arma)/(CLOCKS_PER_SEC);
    cout<<"time used arma \n";
    cout<<time_used_arma<<endl;


    //opening file for storing results
    string s11 = "eigvec";
    string s1 = "res";
    string s2 = to_string(n);
    string s3 = ".txt";
    string s4 = to_string(rhomax);
    string s5 = to_string(omega);

    ofstream myfile;

//    myfile.open(s1+s2+s4+s5+s3); //typ res10200,5.txt  //For plotting the error of different rhomax and n

//    myfile<<"n= "+s2+ " "+"rhomax= "+s4+ + " " +"omega =" + s5+"\n";

    //Finding the smallest eigenvalue and extracting the indices so i can write the corresponding eigenvector to file and plot it

    uword min = b.index_min();
    uvec max_ind = ind2sub(size(b), min);
    int k = max_ind[0];
    int l = max_ind[1];

    //For writing the eigenvectors to a seperate file the eigenvector of the smallest eigenvalue

    myfile.open(s11+s2+s4+s5+s3); //for writing eigenvectors to file and plotting
    myfile<<"n= "+s2+ " "+"rhomax= "+s4+ + " " +"omega= " + s5+"\n";

    //The line below should write the eigenvector of the smallest eigenvalue to file
    myfile<<r.col(l);
    myfile.close();


    double tol = 1E-10;
    bool res = true;
    vec f = sort(b.diag());
    sort(eigval);

//    cout<<f<<endl;
//    cout<<eigval<<endl;

    //this checks if the eigenvalues are close enough to what they were before rotation.
    for(int i = 0; i<n-1; i++){

        if (abs(f(i) - eigval(i)) < tol){
            res= true;
        }
        else{
            cout<<b.diag()<<endl;
//            cout<<b(i,i)<<endl;
            cout<<eigval<<endl;

            res = false;
            return res;
        }
    }

    //comapring the analytical and numerical results for D
    //and writing it to file


    vec ana(n);
    anal(n, ana);
    vec diff(n);

    cout<<"--------f---------\n";
    cout<<f<<endl;
    cout<<"-----anal----------\n";
    cout<<ana<<endl;
    cout<<"---------armadillo----------\n";
    cout<<new_eigval<<endl;
    cout<<"----------diff-------------\n";
    diff = abs(f-ana);
    cout<<diff<<endl;
//    myfile<<diff;






    return res;
}


void anal(uword n, vec & b ){
    //returns a vector of analytical eigenvectors
    for(unsigned int i=0; i<n; i++){
        b(i) = 3 + 4*i;
    }
}
void anal_b(uword n, vec &b){
    double rhomin = 0;
    double rhomax = 1;
    double h = (rhomax - rhomin)/n;

    double d = 2/(h*h);
    double a = -(1/(h*h));
    for (unsigned int i = 1; i<n+1;i++){
        b(i-1) = d + 2*a * cos((i*M_PI)/(n+1));
    }
}

void anal_e(uword n, vec &b){
    //return a vector of analytical eigenvalues, but for two particle system
}

void e_func(double omega){
    uword l = 0;
    uword k = 0;

    double maxval=100;
    double eps = 1E-10;
    unsigned int n = 10;
    double rhomin = 0;
    double rhomax = 15;
    double h = (rhomax - rhomin)/n;

    double rho;


    double a = -1/(h*h);

    //filling matrix
    Mat<double> P(n,n, fill::zeros);
    for(unsigned int i=0;i< n-1; i++){
        rho = (i+1)*h;
        P(i,i) = omega*rho*rho + (1/rho);
        P(i+1,i) =a;
        P(i,i+1) = a;

    }
    double end = h*n;
    P(n-1,n-1) = omega*end*end + (1/end);

    cout << P <<endl;


    int count = 0;

    mat r = eye(n,n);
    while (count < 100000 && abs(maxval) > eps){

        findmax(P, l, k, maxval);
        //cout<<G<<endl;

        rotate(P, r, k, l, n);
        //cout<<b<<endl<<endl;


        count++;
    }

    vec eigval(n);
    mat eigvec(n,n);
    eig_sym(eigval, eigvec, P);



    sort(eigval);


    //calling eigenvalue test
    bool returnvalue_eigenvec = eigval_test(P, r, eigval, eigvec, n, rhomax, omega);
    if (returnvalue_eigenvec == false){
        cout<<"fail in eigval test"<<endl;
        }


}



// sjekke symetrien ved fabs((a_kl - a_lk)/a_kl) < epsilon
// sjekke ortognalitet av egenvektorer?


