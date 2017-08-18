#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(x);
  DATA_VECTOR(y);
  DATA_VECTOR(n);
  PARAMETER(a);
  PARAMETER(b);
  PARAMETER(log_gamma); Type gamma = 1 + 1e-10 + exp(log_gamma); ADREPORT(gamma);
  
  Type nll = 0;
  for (int i=0; i<y.size(); i ++) {
    Type eta = a + b*x(i);
    Type p     = 1/(1+exp(-eta));
    Type alpha = p      *(n(i) - 1)/(gamma - 1) + 1e-100;
    Type beta  = (1 - p)*(n(i) - 1)/(gamma - 1) + 1e-100;
    // log of Beta binomial point mass function
    nll -=
      lfactorial(n(i)) - lfactorial(y(i)) - lfactorial(n(i) - y(i))                     \
      + lgamma(alpha + y(i)) + lgamma(beta + n(i) - y(i)) - lgamma(alpha + beta + n(i)) \
      - lgamma(alpha)        - lgamma(beta)               + lgamma(alpha + beta);
  }
  return nll;
}
