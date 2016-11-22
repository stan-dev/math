
#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <nngp_log.hpp>

#include <iostream>
#include <fstream>

using Eigen::Dynamic;
using Eigen::Matrix;
using std::vector;

int main() {
    
    const int N = 200;
    const int M = 8;
    const int P = 2;
    
    
    /*--- inputting X Y neardist neardistM and nearind ---*/
    
    std::ifstream y_file ("./data/y.txt");
    Matrix<double,Dynamic,1> y(N, 1);
    if(y_file.is_open()){
        std::cout << "able to open file y.txt"
        << std::endl;
        for (int i = 0; i < N; i++ ){
                y_file >> y(i, 0);
        }
        y_file.close();
    }
    else std::cout << "Unable to open file y.txt"
        << std::endl;
    
    
    std::ifstream X_file ("./data/X.txt");
    Matrix<double,Dynamic, P> X(N, P);
    if(X_file.is_open()){
        std::cout << "able to open file X.txt"
        << std::endl;
        for(int j = 0; j < P; j++){
            for (int i = 0; i < N; i++ ){
                X_file >> X(i,j);
            }
        }
        X_file.close();
    }
    else std::cout << "Unable to open file X.txt"
        << std::endl;
    
    std::ifstream ndist_file ("./data/neardist.txt");
    Matrix<double, Dynamic, M> neardist((N-1), M);
    if(ndist_file.is_open()){
        std::cout << "able to open file neardist.txt"
        << std::endl;
        for(int i = 0; i < (N-1); i++){
            for (int j = 0; j < M; j++ ){
                ndist_file >> neardist(i,j);
            }
        }
        ndist_file.close();
    }
    else std::cout << "Unable to open file neardist.txt"
        << std::endl;
    
    std::ifstream ndistM_file ("./data/neardistM.txt");
    Matrix<double, Dynamic, M*M> neardistM( (N-1) , M*M);
    if(ndistM_file.is_open()){
        std::cout << "able to open file neardistM.txt"
        << std::endl;
        for(int i = 0; i < (N-1); i++){
            for (int j = 0; j < M*M; j++ ){
                ndistM_file >> neardistM(i,j);
            }
        }
        ndistM_file.close();
    }
    else std::cout << "Unable to open file neardistM.txt"
        << std::endl;
    
    std::ifstream distind_file ("./data/nearind.txt");
    Matrix<int, Dynamic, M>  nearind( (N-1), M);
    if(distind_file.is_open()){
        std::cout << "able to open file nearind.txt"
        << std::endl;
        for(int i = 0; i < (N-1); i++){
            for (int j = 0; j < M; j++ ){
                distind_file >> nearind(i,j);
            }
        }
        distind_file.close();
    }
    else std::cout << "Unable to open file nearind.txt"
        << std::endl;
    
    
    
    /*---  set parameters  ---*/
    
    double sigma = 1.0;
    double tau = sqrt(0.1);
    double phi = 12.0;
    Matrix<double,Dynamic,1> beta(2,1);
    beta << 1.0, 5.0;
    
    /*---  calculate the nngp  ---*/
    
    std::cout << "nngp_log(Y, X, beta, sigma, tau, phi,"
    "neardist, neardistM, nearind)="
    << stan::math::nngp_log(y, X, beta, sigma, tau, phi,
                            neardist, neardistM, nearind)
    << std::endl;
    
    
    return 0;
}









































