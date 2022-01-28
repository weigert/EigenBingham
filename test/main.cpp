#include <iostream>
#include "../source/bingham.h"

using namespace std;

int main( int argc, char* args[] ) {

//  srand(time(NULL));

  Bingham::Bingham B;

  /*

  Bingham::matX X = Bingham::matX(6, 3);
  X(0, 0) = 1;
  X(0, 1) = 0;
  X(0, 2) = 0;
  X(1, 0) = 0;
  X(1, 1) = 1;
  X(1, 2) = 0;
  X(2, 0) = 0;
  X(2, 1) = 0;
  X(2, 2) = 1;

  X(3, 0) = -1;
  X(3, 1) = 0;
  X(3, 2) = 0;
  X(4, 0) = 0;
  X(4, 1) = -1;
  X(4, 2) = 0;
  X(5, 0) = 0;
  X(5, 1) = 0;
  X(5, 2) = -1;

  */

/*
  const size_t N = 500;
  Bingham::matX X = Bingham::matX(N, 3);
  for(size_t i = 0; i < N; i++){

    Bingham::vecD r = Bingham::vecD::Random();
    if(r(0) < 0) r(0) = -r(0);
    if(r(1) < 0) r(1) = -r(1);
    if(r(2) < 0) r(2) = -r(2);
    r = r.normalized();

    X(i, 0) = r(0);
    X(i, 1) = r(1);
    X(i, 2) = r(2);

  }
*/

  //B.compute(X);
  B.Z << 0.00001, 0.00001, 0.00001;//, 1;
  B.F = B.normalize();


  cout<<"Concentration: "<<B.Z<<endl;
  cout<<"Normalization: "<<B.F<<endl;

  Bingham::vecD x;
  x << 0, 0, 1;//, 0, 0;
  cout<<"Probability: "<<B.P(x)<<endl;
  cout<<B.F<<endl;
  cout<<B.Z<<endl;
  /*
  cout<<"Bingham Distribution: "<<endl;

//  cout<<"Samples: "<<X<<endl;
  cout<<"Covariance: "<<B.C<<endl;
  cout<<"Orientation: "<<B.M<<endl;
  cout<<"Normalization: "<<B.F<<endl;
  */

  return 0;
}
