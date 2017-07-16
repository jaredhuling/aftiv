
#include "estimating_equations.h"

using namespace Rcpp;
using namespace RcppEigen;

bool compare(int a, int b, double* data)
{
    return data[a]<data[b];
}


Eigen::VectorXd rcppRev(Eigen::VectorXd x) {
  using Eigen::Map;
  using Eigen::VectorXd;
  typedef Eigen::VectorXd   Vd;
  Vd revX(x);
  
  revX.reverse();

  return revX;
} 

Eigen::VectorXd cumsum_rev(Eigen::VectorXd x){
  using Eigen::VectorXd;
  typedef Eigen::VectorXd   Vd;
  Vd xx(x);
  const int n(xx.size());
  xx.reverseInPlace();
  std::partial_sum(xx.data(), xx.data() + n, xx.data());
  xx.reverseInPlace();
  return xx;    // compute the result vector and return it
}
/* fix this
const double step_func(double newx, Eigen::VectorXd x, Eigen::VectorXd y) {
  using namespace Rcpp;
  using namespace std;
  
  try {
    using Rcpp::NumericVector;  
  
    //const NumericVector xx(x);
    const NumericVector yy(y);
    const double new_x(as<double>(newx));
    std::vector<double> xx = Rcpp::as<std::vector<double> >(x);
    std::vector<double>::iterator low;
    low = lower_bound(xx.begin(), xx.end(), new_x);
    
    //return yy.lower_bound(new_x);
    
    return wrap(yy[low - xx.begin()]);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}

*/
/*
RcppExport SEXP approx1000(SEXP vv, SEXP xx, SEXP yy)
{
  using namespace Rcpp;
  using namespace std;
  
  try {
    // Approximate  y(v),  given (x,y)[i], i = 0,..,n-1 
    int i, j, ij;
    
    const double v(Rcpp::as<double>(vv));
    
    std::vector<double> x = Rcpp::as<std::vector<double> >(xx);
    std::vector<double> y = Rcpp::as<std::vector<double> >(yy);

    int n(x.size());

    i = 0;
    j = n - 1;

    // handle out-of-domain points
    if(v < x[i]) return wrap(y[i]);
    if(v > x[j]) return wrap(y[j]);

    // find the correct interval by bisection 
    while(i < j - 1) { // x[i] <= v <= x[j] 
	ij = (i + j)/2; // i+1 <= ij <= j-1 
	if(v < x[ij]) j = ij; else i = ij;
	// still i < j 
    }
    // provably have i == j-1 

    // interpolation

    if(v == x[j]) return wrap(y[j]);
    if(v == x[i]) return wrap(y[i]);
    // impossible: if(x[j] == x[i]) return y[i];

	  return wrap(y[i] + (y[j] - y[i]) * ((v - x[i])/(x[j] - x[i])));
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}

RcppExport SEXP cumsum_rev_cp(SEXP x) {
  using namespace Rcpp;
  using namespace RcppEigen;
  using namespace std;
  
  try {
    using Eigen::VectorXd;
    
    typedef Eigen::Map<VectorXd> MapVecd;
    const MapVecd xx(as<MapVecd>(x));
    VectorXd xx2(xx);
    const int n(xx.size());
    xx2.reverseInPlace();
    std::partial_sum(xx2.data(), xx2.data() + n, xx2.data());
    xx2.reverseInPlace();
    return wrap(xx2);    // compute the result vector and return it
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}
*/

// port function from R
RcppExport SEXP AFTScorePre(SEXP XX,        // design matrix
		       SEXP delta_vec) 
{
  using namespace Rcpp;
  using namespace RcppEigen;
  using namespace std;
  
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using Eigen::Lower;
    using Rcpp::NumericVector;
    using Rcpp::IntegerVector;
    typedef Eigen::Map<VectorXd> MapVecd;
    
    const Eigen::Map<MatrixXd> X(as<Map<MatrixXd> >(XX));
    const MapVecd delta(as<MapVecd>(delta_vec));
    const int nobs(X.rows());
    const int nvars(X.cols());
    //int index[nobs];
    
    //VectorXd err(VectorXd::Zero(nobs));
    VectorXd at_risk_X_terms(VectorXd::Zero(nobs));
    VectorXd at_risk_terms(VectorXd::Zero(nobs));
    VectorXd X_col(VectorXd::Zero(nobs));
    VectorXd to_sum(VectorXd::Zero(nobs));
    VectorXd ret_vec(VectorXd::Zero(nvars));
    
    //err = log_t_vec - X * beta_vec;
    
    //std::iota(index.begin(), index.end(), 0); // fill index with {0,1,2,...} This only needs to happen once
    //std::sort(index.begin(), index.end(), std::bind(compare,  _1, _2, err ));
    
    //std::iota(at_risk_terms.data() + nobs, at_risk_terms.data(), 1);
    for (int i = nobs; i > 0; i--) 
    {
      at_risk_terms[nobs - i] = i;
    }
    at_risk_terms /= nobs;

    //err = err[index];
    //log_t_vec = log_t_vec[index];
    //X = X(index,_ );
    //delta_vec = delta_vec[index];
    
    for (int i = 0; i < nvars; i++) 
    {
      X_col = X.col(i);
      at_risk_X_terms = cumsum_rev(X_col) / nobs;

      to_sum = delta.array() * (at_risk_terms.array() * X_col.array() - at_risk_X_terms.array()) / sqrt(nobs);

      ret_vec[i] = to_sum.sum();
    }
    
    return wrap(ret_vec);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
  
}