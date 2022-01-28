#ifndef EIGEN_BINGHAM
#define EIGEN_BINGHAM

#include <Eigen/Dense>
#include <functional>

/*
    Methods: Generate a Bingham Distribution,
    Fitting, Mixture Models, Hypersphere Methods,
    No Meshing Functions or whatever...
    Lookup Table Construction as Raw Binary!!
    Probably slow as hell...
    Potentially Faster Computation of Partition Function
    Gradient Computation of F
    A lot of stuff in libbingham is extra crap that isn't needed.
*/

namespace Bingham {

using namespace Eigen;
using namespace std;

using T = double;                      //Floating Point Precision
const size_t DIM = 3;                 //Dimensionality of a Point
const T SPI = 2.0 / sqrt(M_PI);

const T EPSILON = 1E-8;
const size_t MAXFACT = 10000;
const size_t ITERATION_MULT = 10;
const size_t MIN_ITERATIONS = 10;

typedef Matrix<T, DIM, DIM>     matD; //DIM x DIM Matrix
typedef Matrix<T, Dynamic, DIM> matX; //DIM x X Matrix
typedef Matrix<T, DIM, 1>       vecD; //DIM x 1 Vector

/*
================================================================================
          Computation of the Hypergeometric Function of Matrix Argument
================================================================================
*/

// Log Factorials and Factorials

T lfact(size_t k){

  static T logf[MAXFACT];   //lfact lookuptable
  static bool first = true;

  if(first){
    first = false;
    logf[0] = 0;
    for (size_t i = 1; i < MAXFACT; i++)
      logf[i] = log(i) + logf[i-1];
  }

  return logf[k];

}

T fact(size_t k){ //Factorial of k
  return exp(lfact(k));
}

// Surface Area of N-Sphere

T NSphereArea(size_t D){
  if(D == 0) return 2;
  if(D == 1) return 2*M_PI;
  if(D == 2) return 4*M_PI;
  if(D == 3) return 2*M_PI*M_PI;
  return (2*M_PI/((T)D-1))*NSphereArea(D-2);
}

// Hypergeometric Functions of Matrix Argument

template<size_t D>
T C1F1(T a, T b, vecD Z, size_t iter);

template<>
T C1F1<1>(T a, T b, vecD Z, size_t iter){

  T g, F = 0;
  vecD logZ = Z.array().log();

  for(size_t i = 0; i < iter; i++){
    g = 0;
    g += lgamma( i + a );
    g -= lgamma( i + b );
    g += i * logZ(DIM-1);
    g -= lfact(i);
    //if ((i > Z(0)) && exp(g) < EPSILON * F)
  //    break;
    F += exp(g);
  }

  F /= 2;
  if(DIM % 2 == 0) return F * SPI;
  return F;

}

template<>
T C1F1<2>(T a, T b, vecD Z, size_t iter){

  T g, F = 0;
  vecD logZ = Z.array().log();

  for(size_t i = 0; i < iter; i++){
    for(size_t j = 0; j < iter; j++){
      g = 0;
      g += lgamma( i + a ) + lgamma( j + a );
      g -= lgamma( i + j + b );
      g += i * logZ(DIM-2) + j * logZ(DIM-1);
      g -= lfact(i) + lfact(j);
  //    if ((i > Z(0) || j > Z(1)) /*&& exp(g) < EPSILON * F*/)
  //      break;
      F += exp(g);
    }
  }

  F /= 2 * sqrt(M_PI);
  if(DIM % 2 == 0) return F * SPI;
  return F;

}

template<>
T C1F1<3>(T a, T b, vecD Z, size_t iter){

  T g, F = 0;
  vecD logZ = Z.array().log();

  for(size_t i = 0; i < iter; i++){
    for(size_t j = 0; j < iter; j++){
      for(size_t k = 0; k < iter; k++){
        g = 0;
        g += lgamma( i + a ) + lgamma( j + a ) + lgamma( k + a );
        g -= lgamma( i + j + k + b);
        g += i * logZ(DIM-3) + j * logZ(DIM-2) + k * logZ(DIM-1);
        g -= lfact(i) + lfact(j) + lfact(k);
  //      if ((i > Z(0) || j > Z(1) || k > Z(2)) && exp(g) < EPSILON * F)
  //        break;
        F += exp(g);
      }
    }
  }

  F /= 2 * M_PI;
  if(DIM % 2 == 0) return F * SPI;
  return F;

}

template<>
T C1F1<4>(T a, T b, vecD Z, size_t iter){

  T g, F = 0;
  vecD logZ = Z.array().log();

  for(size_t i = 0; i < iter; i++){
    for(size_t j = 0; j < iter; j++){
      for(size_t k = 0; k < iter; k++){
        for(size_t l = 0; l < iter; l++){
          g = 0;
          g += lgamma( i + a ) + lgamma( j + a ) + lgamma( k + a ) + lgamma( l + a );
          g -= lgamma( i + j + k + l + b);
          g += i * logZ(DIM-4) + j * logZ(DIM-3) + k * logZ(DIM-2) + l * logZ(DIM-1);
          g -= lfact(i) + lfact(j) + lfact(k) + lfact(l);
    //      if ((i > Z(0) || j > Z(1) || k > Z(2)) && exp(g) < EPSILON * F)
    //        break;
          F += exp(g);
        }
      }
    }
  }

  F /= 2 * M_PI * sqrt(M_PI);
  if(DIM % 2 == 0) return F * SPI;
  return F;

}

T MAX(T a, T b){
  return (a > b)?a:b;
}

/*
================================================================================
                            Distribution Wrapper
================================================================================
*/

struct Bingham {

  size_t size = 0;                    //Number of Samples
  matD C = matD::Identity();          //Underlying Gaussian Covariance
  matD M = matD::Identity();          //Orientation Matrix (Orthogonal)
  vecD Z = vecD::Zero();              //Concentration Matrix (Diagonal)
  T F = 1.0f;                         //Normalization Constant

  T d = 0.0f;
  //Constructor from Data

  Bingham(){}

  Bingham(matX& X){
    compute(X);
  }

  const void compute(matX& X){

    size = X.rows();
    if(size >= DIM){

      C = X.adjoint()*X;							   //Covariance Matrix
      C = C/(size-1);                    //Normalize

      SelfAdjointEigenSolver<matD> EIG(C);        //Eigensolver
      M = EIG.eigenvectors();                     //Eigenvectors
      Z = EIG.eigenvalues();                      //Eigenvalues

    }

    F = normalize();

  }

  //Normalization Constant

  T normalize(){

    /*
        Note: We subtract the smallest eigenvalue from all eigenvalues and shift everything!
        Then: The sphere dimensionality is actually only 2.
    */

    /*
        Ok: The relation appears to hold true, therefore I assume the formula is valid.
    */

    /*
        Basically: For a quaternion, they are placed on the unit hypersphere in 3 space.
        This means I set DIM to 4 and I have to shift everything s.t. the smallest eigenvalue becomes zero.
        Since they are ordered increasing, everything else will be positive automatically.
        I just need to make sure that the scaling factors are correct in that case.
    */

    //Re-Shift the Set!

    const T shift = 0.0f;//-Z(0);
    vecD NZ = Z.array() + shift;

    //Compute the Hypergeometric Function of Matrix Argument

    size_t iter = 20;

    T FZ;

    //Truncated Computation

    if(DIM >= 1 && fabs(NZ(DIM-1)) < EPSILON)
      FZ = 1;

    else if(DIM >= 2 && fabs(NZ(DIM-2)) < EPSILON)
      FZ = C1F1<1>(0.5, DIM*0.5, NZ, iter);

    else if(DIM >= 3 && fabs(NZ(DIM-3)) < EPSILON)
      FZ = C1F1<2>(0.5, DIM*0.5, NZ, iter);

    else if(DIM >= 4 && fabs(NZ(DIM-4)) < EPSILON)
      FZ = C1F1<3>(0.5, DIM*0.5, NZ, iter);

    //Otherwise: Direct Computation

    else FZ = C1F1<DIM>(0.5, DIM*0.5, NZ, iter);

    return NSphereArea(DIM-1)*exp(-shift)*FZ;

  }

  //Probability Density Function

  T P(vecD x){
    return exp(x.transpose()*(M*Z.asDiagonal()*M.transpose())*x)/F;
  }

};

};

#endif
