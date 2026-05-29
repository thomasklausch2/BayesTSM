#include <Rcpp.h>
using namespace Rcpp;

// ---------- helpers ----------

inline bool is_na_or_nan(const double x) {
  return NumericVector::is_na(x) || ISNAN(x);
}

inline double cdf_loglog(const double x, const double lambda, const double gamma) {
  if (is_na_or_nan(x) || is_na_or_nan(lambda) || is_na_or_nan(gamma)) return NA_REAL;
  if (x == R_PosInf) return 1.0;
  if (x == R_NegInf || x <= 0.0) return 0.0;

  const double z = std::pow(x / lambda, gamma);
  return z / (1.0 + z);
}

inline double quantile_loglog(const double p, const double lambda, const double gamma) {
  if (is_na_or_nan(p) || is_na_or_nan(lambda) || is_na_or_nan(gamma)) return NA_REAL;
  if (p <= 0.0) return 0.0;
  if (p >= 1.0) return R_PosInf;

  return lambda * std::pow(p / (1.0 - p), 1.0 / gamma);
}

inline int dist_to_code(const std::string& dist) {
  if (dist == "exp")       return 1;
  if (dist == "gamma")     return 2;
  if (dist == "weibull")   return 3;
  if (dist == "loglog")    return 4;
  if (dist == "lognormal") return 5;

  stop("Unsupported dist. Implemented: exp, gamma, weibull, loglog, lognormal.");
  return -1;
}

inline double cdf_one(
    const double x,
    const double par1,
    const double par2,
    const int dist_code
) {
  switch (dist_code) {
  case 1: // exp: par1 = rate
    return R::pexp(x, par1, /*lower_tail=*/1, /*log_p=*/0);

  case 2: // gamma: par1 = shape, par2 = rate
    return R::pgamma(x, par1, 1.0 / par2, /*lower_tail=*/1, /*log_p=*/0);

  case 3: // weibull: par1 = scale, par2 = shape
    return R::pweibull(x, par2, par1, /*lower_tail=*/1, /*log_p=*/0);

  case 4: // loglog: par1 = lambda, par2 = gamma
    return cdf_loglog(x, par1, par2);

  case 5: // lognormal: par1 = meanlog, par2 = sdlog
    return R::plnorm(x, par1, par2, /*lower_tail=*/1, /*log_p=*/0);

  default:
    stop("Unsupported dist_code in cdf_one().");
  }

  return NA_REAL;
}

inline double quantile_one(
    const double p,
    const double par1,
    const double par2,
    const int dist_code
) {
  switch (dist_code) {
  case 1: // exp
    return R::qexp(p, par1, /*lower_tail=*/1, /*log_p=*/0);

  case 2: // gamma
    return R::qgamma(p, par1, 1.0 / par2, /*lower_tail=*/1, /*log_p=*/0);

  case 3: // weibull
    return R::qweibull(p, par2, par1, /*lower_tail=*/1, /*log_p=*/0);

  case 4: // loglog
    return quantile_loglog(p, par1, par2);

  case 5: // lognormal
    return R::qlnorm(p, par1, par2, /*lower_tail=*/1, /*log_p=*/0);

  default:
    stop("Unsupported dist_code in quantile_one().");
  }

  return NA_REAL;
}

//' Draw latent S values
//'
//' @param par Numeric matrix of parameters.
//' @param dist Distribution name.
//' @param C Numeric vector.
//' @param d Integer vector.
//' @param L Numeric vector.
//' @param R Numeric vector.
//' @param tol Tolerance.
//' @return Numeric vector.
//' @noRd
// [[Rcpp::export]]
NumericVector pst_S_cpp(
  const NumericMatrix& par,
  const std::string& dist,
  const NumericVector& X,
  const IntegerVector& d,
  const NumericVector& L,
  const NumericVector& R,
  const double tol = 1e-8
) {
  RNGScope scope;

  const int n = d.size();

  if (par.nrow() != n) stop("par must have n rows.");
  if (X.size()   != n) stop("X must have length n.");
  if (L.size()   != n) stop("L must have length n.");
  if (R.size()   != n) stop("R must have length n.");
  if (par.ncol() < 1)  stop("par must have at least one column.");

  const int dist_code = dist_to_code(dist);
  if (dist_code != 1 && par.ncol() < 2) {
    stop("par must have at least two columns for this distribution.");
  }

  NumericVector out(n, NA_REAL);

  for (int i = 0; i < n; ++i) {
    const double par1 = par(i, 0);
    const double par2 = (par.ncol() >= 2) ? par(i, 1) : NA_REAL;

    if (is_na_or_nan(par1) || (dist_code != 1 && is_na_or_nan(par2))) {
      out[i] = NA_REAL;
      continue;
    }

    if (d[i] == 1) {
      // Equivalent in distribution to rdist(...):
        // sample from the full distribution via inverse CDF.
      const double u = R::runif(0.0, 1.0);
      out[i] = quantile_one(u, par1, par2, dist_code);
      continue;
    }

    if (d[i] == 2) {
      // a = R - X, b = Inf
      if (is_na_or_nan(R[i]) || is_na_or_nan(X[i])) {
        out[i] = NA_REAL;
        continue;
      }

      const double a = R[i] - X[i];
      const double b = R_PosInf;

      const double cdf_a = cdf_one(a, par1, par2, dist_code);
      const double cdf_b = cdf_one(b, par1, par2, dist_code);

      if (is_na_or_nan(cdf_a) || is_na_or_nan(cdf_b)) {
        out[i] = NA_REAL;
        continue;
      }

      const double width = cdf_b - cdf_a;

      if (width < tol) {
        out[i] = a;
      } else {
        const double u = cdf_a + R::runif(0.0, 1.0) * width;
        out[i] = quantile_one(u, par1, par2, dist_code);
      }
      continue;
    }

    if (d[i] == 3) {
      // a = 0, b = min(R - X, R - L)
      if (is_na_or_nan(R[i]) || is_na_or_nan(X[i]) || is_na_or_nan(L[i])) {
        out[i] = NA_REAL;
        continue;
      }

      const double a = 0.0;
      const double b = std::min(R[i] - X[i], R[i] - L[i]);  // min is essential

      const double cdf_a = cdf_one(a, par1, par2, dist_code);
      const double cdf_b = cdf_one(b, par1, par2, dist_code);

      if (is_na_or_nan(cdf_a) || is_na_or_nan(cdf_b)) {
        out[i] = NA_REAL;
        continue;
      }

      const double width = cdf_b - cdf_a;

      if (width < tol) {
        out[i] = a;
      } else {
        const double u = cdf_a + R::runif(0.0, 1.0) * width;
        out[i] = quantile_one(u, par1, par2, dist_code);
      }
      continue;
    }

    stop("d must contain only 1, 2, or 3.");
  }

  return out;
}
