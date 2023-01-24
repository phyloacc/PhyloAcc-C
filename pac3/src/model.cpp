#include <algorithm>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

void stuff_Phi_matrix(arma::mat& Phi, double phi_a, double phi_b, double phi_ac) {
  Phi(0, 0) = 1 - phi_a; Phi(0, 1) = phi_a;
  Phi(1, 1) = 1 - phi_b; Phi(1, 2) = phi_b;
  // was Phi(2, 2) = 1 under Dollo approach but now:
  Phi(2, 1) = phi_ac;
  Phi(2, 2) = 1 - phi_ac;
}

double logaddexp(const arma::rowvec& xs) {
  double max_x = arma::max(xs);
  double y = log(arma::sum(arma::exp(xs - max_x))) + max_x;
  NumericVector yy(1); yy(0) = y;
  if (Rcpp::any(Rcpp::is_nan(yy))) { xs.print("*LOGADDEXP XS*"); }
  return y;
}

int sample_ulp(const arma::rowvec& log_ps) {
  double max_lps = arma::max(log_ps);
  auto log_ps_ = log_ps - max_lps;
  arma::rowvec ps = arma::exp(log_ps_);
  ps = ps / arma::sum(ps); // note / is elem wise division
  double r = R::runif(0.0, 1.0);
  arma::rowvec cs = cumsum(ps);
  for (int i = 0; i < cs.n_elem; i++) { if (r < cs(i)) return i; }
  Rcpp::stop("invalid sample");
}

// Set each hidden node to midpoint of its children.
void initialize_Y(arma::vec& Y, const arma::Mat<int>& edge) {
  for (int e = 0; e < edge.n_rows; e += 2) {
    int i, j, k;
    i = edge(e,     0);
    j = edge(e,     1); // edges always come in pairs
    k = edge(e + 1, 1);
    Y(i) = 0.5*(Y(j) + Y(k)); // TODO: distance weight this?
  }
}

// Probabilities at leaves should not change during MCMC: set them up once.
void init_dp_Xup(
    arma::mat& lp_Xup,
    int s, // site we consider
    const int (&acgt)[4],
    int L,
    const arma::Mat<int>& X)
  {
  // leaf probabilities based on data
  for (int l = 0; l < L; l++) {
    if (4 == X(l,s)) {
      for (int x : acgt) { lp_Xup(l,x) = log(0.25); } // treat gaps/unknown as uninformative
    } else {
      // this code will be common case
      for (int x : acgt) {
        if (x == X(l,s)) { lp_Xup(l,x) = log(1); } else { lp_Xup(l,x) = log(0.0); }
      }
    }
  }
}

// Joe F's pruning algorithm.
void change_aware_dp_Xup(
    arma::mat& lp_Xup,
    int s, // site we consider
    const int (&acgt)[4],
    int L,
    int E,
    const arma::Mat<int>& X, // todo: remove X parameter as not used
    const arma::Mat<int>& edge,
    const arma::Col<int>& Z,
    arma::cube * P_nuctrans[3],
    const std::vector<bool>& update_node)
  {
  // node probabilities based on traversal leaf to root
  for (int e = 0; e < E; e += 2) {
    int i = edge(e,     0);
    if (!update_node[i]) { continue; }
    int j = edge(e,     1); // edges always come in pairs
    int k = edge(e + 1, 1);
    arma::mat p_ij = P_nuctrans[Z[j]]->slice(e    ); // arma::expmat(Q * r[Z[j]] * distance(e    ));
    arma::mat p_ik = P_nuctrans[Z[k]]->slice(e + 1); // arma::expmat(Q * r[Z[k]] * distance(e + 1));
    for (int x : acgt) { // for subtree root possibilities
      double p_subj = arma::dot(p_ij.row(x), arma::exp(lp_Xup.row(j)));
      double p_subk = arma::dot(p_ik.row(x), arma::exp(lp_Xup.row(k)));
      lp_Xup(i,x) = log(p_subj) + log(p_subk);
    }
  }
}

double mini_perturb() { return R::rnorm(0, 0.01); }

double Big_prob = 0.40, // setting lower from 0.5 (e.g. .4 seems to give better Rhat)
  Perturb_Y_sd  = 0.10, // put this .05 from .1
  Perturb_lv_sd = 0.75,
  Perturb_lb_sd = 0.75;

void perturb_lvby(double& lv, double& lb2, double& lb3, arma::vec& Y, int R,
                  const int& E,
                  const arma::Mat<int>& edge,
                  const int N) {
  if (R::rbinom(1, Big_prob)) {
    int n = floor(R::runif(R, N));
    Y[n] += R::rnorm(0, Perturb_Y_sd);
  } else {
    double q = R::runif(0, 1);
    if (q < 0.33) {
      double x = R::rnorm(0, Perturb_lv_sd);
      lv  += x + mini_perturb();
      lb2 -= x + mini_perturb(); // 0.5 * x ?
      lb3 -= x + mini_perturb();
    } else if (q < 0.66) {
      double x = R::rnorm(0, Perturb_lb_sd);
      lv  -= x + mini_perturb();  
      lb3 -= x + mini_perturb();
      lb2 += x + mini_perturb();
    } else {
      double x = R::rnorm(0, Perturb_lb_sd);
      lv  -= x + mini_perturb();
      lb2 -= x + mini_perturb();
      lb3 += x + mini_perturb();
    }
  }
}

// [[Rcpp::export]]
Rcpp::List Sample_n(
    const int n_iter, /* N.B. everything passed in is assumed 0 index */
    const int L, // n tips/leaves
    const int R, // root node
    arma::Mat<int> edge, // !! can I just pass ape tree?
    arma::vec distance,
    arma::Mat<int> X, // alignment + buffer
    arma::vec Y, // trait + buffer
    arma::mat transition_matrix, // nucleotide rate matrix
    arma::Row<double> pi_distro, // stationary nucleotide distro. at root
    double prior_r2_a, //=  5.0; // gamma hyperprior on con rate (from C++)
    double prior_r2_b, //=  0.04;
    double prior_r3_a, //= 10.0; // gamma hyperprior on acc rate (from C++), in paper p. 1098 this is 15
    double prior_r3_b, //=  0.2; // in paper p. 1098 this is 0.1
    double prior_zr_1, // P(z_R=1) i.e. root is in neu/bg state
    double prior_gc_a, // prior beta a to gain conservation
    double prior_gc_b, //       beta b
    double prior_lc_a, // prior beta a to lose conservation
    double prior_lc_b, //       beta b
    double prior_ac_a, // no-Dollo assumption: prior beta a to back-switch to con from acc
    double prior_ac_b, // no-Dollo assumption: prior beta b to back-switch to con from acc
    double prior_root_u,    // prior ancestral trait at root: mean
    double prior_root_sd,   //                                sd
/* 17-12-2022 if set uniform prior on variance (below), next 2 lv args used only to set initial value */
    double prior_lv_u,      // prior on log of mean/sd of bg/con/acc variance multipliers
    double prior_lv_sd,
    double prior_lb2_u,
    double prior_lb2_sd,
    double prior_lb3_u,
    double prior_lb3_sd,
    bool do_Y, // set do_Y to false to ignore trait
    bool uniform_prior_lv, // if true, use uniform prior instead of above lognormal one
    const int inner_iter, // n steps within Metropolis loop for trait
    bool verbose) // print out status of iterations etc.
{
  const int acgt[4] { 0, 1, 2, 3 };
  
  const int N = 2*L - 1;  // total node count
  const int S = X.n_cols; // n sites
  const int E = edge.n_rows; if(E != N - 1) Rcpp::stop("E != N - 1"); // n edges
  
  arma::Col<int> Z(N); // conservation state 0, 1, 2 = neu, con, acc
  Z.ones();
  Z[R] = 1;
  
  arma::rowvec prior_ZR = {{prior_zr_1, 1.0 - prior_zr_1, 0}}; // prior for Z leading to root
  
  arma::vec r = {{1.0, 0.1, 0.7}};
  arma::vec r_(3); // convention is trailing _ means Metropolis proposals go here
  
  arma::mat Q = transition_matrix;
  arma::rowvec PI_D = pi_distro;
   
  double phi_a = R::rbeta(prior_gc_a, prior_gc_b); // alpha
  double phi_b = R::rbeta(prior_lc_a, prior_lc_b); // beta
  double phi_ac= R::rbeta(prior_ac_a, prior_ac_b); // non-Dollo prob acc to con switch

  arma::mat Phi(3, 3, arma::fill::zeros);
  stuff_Phi_matrix(Phi, phi_a, phi_b, phi_ac);
  
  arma::Mat<int> z_trace(n_iter, N); z_trace.fill(arma::datum::nan);
  arma::mat      r_trace(n_iter, 3); r_trace.fill(arma::datum::nan);
  arma::mat    phi_trace(n_iter, 3); phi_trace.fill(arma::datum::nan);

  double lv  = R::rnorm(prior_lv_u, prior_lv_sd);   // log base var
  double lb2 = R::rnorm(prior_lb2_u, prior_lb2_sd); // log con multiplier
  double lb3 = R::rnorm(prior_lb3_u, prior_lb3_sd); // log acc multiplier
  double lv_;
  double lb2_;
  double lb3_;
  
  initialize_Y(Y, edge); // Y supplied as argument
  arma::vec Y_(Y.n_rows); Y_ = Y; // after this only fill nodes not leaves
  
  arma::mat Y_trace(n_iter, N); Y_trace.fill(arma::datum::nan);
  arma::mat v_trace(n_iter, 3); v_trace.fill(arma::datum::nan);

  // nucleotide transition matrix caches, indexed by edge.
  arma::cube cache_P_neux(4, 4, E); // never changes
  arma::cube cache_P_con1(4, 4, E);
  arma::cube cache_P_con2(4, 4, E);
  arma::cube cache_P_acc1(4, 4, E);
  arma::cube cache_P_acc2(4, 4, E);

  // swap these pointers whenever we pick new rates
  arma::cube * P_nuctrans [3] = { &cache_P_neux, &cache_P_con1, &cache_P_acc1 };
  arma::cube * P_nuctrans_[3] = { &cache_P_neux, &cache_P_con2, &cache_P_acc2 };
  
  // setup current cache
  for (int e = 0; e < E; e++) {
    P_nuctrans[0]->slice(e) = arma::expmat(Q * r[0] * distance(e));
    P_nuctrans[1]->slice(e) = arma::expmat(Q * r[1] * distance(e));
    P_nuctrans[2]->slice(e) = arma::expmat(Q * r[2] * distance(e));
  }
  
  // cache for pruning algorithm
  arma::cube cache_lp_Xups1(N, 4, S, arma::fill::zeros); // prob subtree at i if nuc(i) = j, indexed by site
  arma::cube cache_lp_Xups2(N, 4, S, arma::fill::zeros);
  
  // swap these pointers whenever we pick new rates
  arma::cube * lp_Xups  = &cache_lp_Xups1;
  arma::cube * lp_Xups_ = &cache_lp_Xups2;
  
  // initialize leaves in both caches
  for (int s = 0; s < S; s++) {
    arma::mat& lp_Xup  = lp_Xups->slice(s);
    init_dp_Xup(lp_Xup , s, acgt, L, X);
    arma::mat& lp_Xup_ = lp_Xups_->slice(s); // todo: replace with a copy
    init_dp_Xup(lp_Xup_, s, acgt, L, X);
  }
  
  // arrays to avoid unnecessary updating between sampling Z and deciding between r, r_
  std::vector<bool> update_node(N, false);
  const std::vector<bool> update_everything(N, true);
  
  // initialize internal values in current cache only
  for (int s = 0; s < S; s++) {
    arma::mat& lp_Xup = lp_Xups->slice(s);
    // use DP to calculate prob of alignment above node i
    change_aware_dp_Xup(lp_Xup, s, acgt, L, E, X, edge, Z, P_nuctrans, update_everything);
  }
  
  for (int iter = 0; iter < n_iter; iter++) {
    if ((0 == iter % 100) && verbose) { Rcout << "iter: " << iter << "\n"; }
    if (do_Y) {
      // calculate trait likelihood under existing situation
      arma::vec varmults { exp(lv ), exp(lv  + lb2 ), exp(lv  + lb3 ) };
      double lly_M1  = log(1);
      for (int e = E-1; e >= 0; e--) {
          int i = edge(e, 0);
          int j = edge(e, 1);
          double vij  = varmults (Z(j)) * distance(e);
          double lp_move  = R::dnorm(Y (j), Y (i), sqrt(vij ), true);
          lly_M1  += lp_move;
      }
      double lp_lv = 0.0;
      if (!uniform_prior_lv) {
        lp_lv   = R::dnorm(lv,  prior_lv_u, prior_lv_sd, true); // PG 17-12-2022 now uniform prior on v is possible
      }

      double lp_lb2  = R::dnorm(lb2 , prior_lb2_u, prior_lb2_sd, true);
      double lp_lb3  = R::dnorm(lb3 , prior_lb3_u, prior_lb3_sd, true);
      double lp_root = R::dnorm(Y (R), prior_root_u, prior_root_sd, true);
      // final llhood existing
      double ll_denom_1 = lly_M1 + lp_lv + lp_lb2 + lp_lb3 + lp_root;
      // propose, calculate, perhaps accept and update      
      for (int ii = 0; ii < inner_iter; ii++) {
        // -Metropolis Y
        lv_ = lv; lb2_ = lb2; lb3_ = lb3; Y_.rows(R,N-1) = Y.rows(R,N-1);
        perturb_lvby(lv_, lb2_, lb3_, Y_, R, E, edge, N);
        arma::vec varmults_ { exp(lv_), exp(lv_ + lb2_), exp(lv_ + lb3_) };
        double lly_M1_ = log(1);
        for (int e = E-1; e >= 0; e--) {
          int i = edge(e, 0);
          int j = edge(e, 1);
          double vij_ = varmults_(Z(j)) * distance(e);
          double lp_move_ = R::dnorm(Y_(j), Y_(i), sqrt(vij_), true);
          lly_M1_ += lp_move_;
        }
        // y prior part
        double lp_lv_ = 0.0;
        if (!uniform_prior_lv) {
          lp_lv_ = R::dnorm(lv_, prior_lv_u, prior_lv_sd, true); // PG 17-12-2022 now uniform prior on v is possible
        }
        double lp_lb2_ = R::dnorm(lb2_, prior_lb2_u, prior_lb2_sd, true);
        double lp_lb3_ = R::dnorm(lb3_, prior_lb3_u, prior_lb3_sd, true);
        double lp_root_ = R::dnorm(Y_(R), prior_root_u, prior_root_sd, true);
        // final llhood proposed
        lly_M1_ = lly_M1_ + lp_lv_ + lp_lb2_ + lp_lb3_ + lp_root_;
        double ll_num_1 = lly_M1_;
        double log_mh_ratio1 = ll_num_1 - ll_denom_1;
        double ap1 = std::min(1.0, exp(log_mh_ratio1));
        if (R::runif(0, 1) < ap1) {
            Y.rows(R,N-1) = Y_.rows(R,N-1);
            lv = lv_; lb2 = lb2_; lb3 = lb3_;
            ll_denom_1 = ll_num_1;
        }
      }
    }
    // -Gibbs X  
    for (int s = 0; s < S; s++) {
      arma::mat& lp_Xup = lp_Xups->slice(s);
/*    pruning for this step is conducted before main loop entry, and at rate selection step
      dp_Xup(lp_Xup, s, acgt, L, E, X, edge, Z, P_nuctrans); // use DP to calculate prob of alignment above node i */
      // sample hidden X, from root to tip
      arma::rowvec tmp = arma::log(PI_D) + lp_Xup.row(R);
      X(R,s) = sample_ulp(tmp); // kick of by using PI
      for (int e = E-1; e >= 0; e--) {
        int j = edge(e, 1); // visit the 'to end' of edges
        if (j < R) continue;
        int i = edge(e, 0);
        int xi = X(i,s); // parent nucl
        arma::mat p_ij = P_nuctrans[Z[j]]->slice(e); // arma::mat p_ij = arma::expmat(Q * r[Z[j]] * distance(e));
        arma::rowvec lps = arma::log(p_ij.row(xi)) + lp_Xup.row(j);
        X(j,s) = sample_ulp(lps);
      }
    }
    // tip to root, prob sequence under Z=1,2,3 scenarios
    arma::mat log_emission(N, 3); log_emission.fill(arma::datum::nan);
    for (int e = 0; e < E; e++) {
      int i = edge(e, 0); // parent
      int j = edge(e, 1); // child
      arma::mat lp_ij_neu = arma::log(P_nuctrans[0]->slice(e)); // arma::log(arma::expmat(Q * r[0] * distance(e))); // transition scenario
      arma::mat lp_ij_con = arma::log(P_nuctrans[1]->slice(e)); // arma::log(arma::expmat(Q * r[1] * distance(e))); // "
      arma::mat lp_ij_acc = arma::log(P_nuctrans[2]->slice(e)); // arma::log(arma::expmat(Q * r[2] * distance(e))); // "

      double lp_neu = 0; // accumulators
      double lp_con = 0;
      double lp_acc = 0;

      for (int s = 0; s < S; s++) { // Z acts across whole alignment
        if (4 == X(j,s)) {
          if (j >= R) { stop("-/N should only appear at tips."); }
          ; // gaps/unknown do not provide information
        } else {
          int xi = X(i,s);
          int xj = X(j,s);
          lp_neu += lp_ij_neu(xi,xj);
          lp_con += lp_ij_con(xi,xj);
          lp_acc += lp_ij_acc(xi,xj);
        }
      }
      
      log_emission(j,0) = lp_neu;
      log_emission(j,1) = lp_con;
      log_emission(j,2) = lp_acc;   
      
      // begin trait part
      if (do_Y) {
        double var_neu = distance(e) * exp(lv      );
        double var_con = distance(e) * exp(lv + lb2);
        double var_acc = distance(e) * exp(lv + lb3);
  
        double lp_move_neu = R::dnorm(Y(j), Y(i), sqrt(var_neu), true);
        double lp_move_con = R::dnorm(Y(j), Y(i), sqrt(var_con), true);
        double lp_move_acc = R::dnorm(Y(j), Y(i), sqrt(var_acc), true);

        log_emission(j,0) += lp_move_neu;
        log_emission(j,1) += lp_move_con;
        log_emission(j,2) += lp_move_acc;
      }
      // end trait part
    }
    for (int z = 0; z < 3; z++) {
      log_emission(R,z) = log(1);
    }

    arma::mat lp_Zup(N, 3); lp_Zup.fill(arma::datum::nan);
    lp_Zup.rows(0,L) = log_emission.rows(0,L);
    
    for (int e = 0; e < E; e+=2) {
      int i, j, k;
      i = edge(e,     0);
      j = edge(e,     1); // edges always come in pairs
      k = edge(e + 1, 1);
      for (int z = 0; z < 3; z++) { // subtree probs when z(i) = z
        arma::rowvec lp_subj = arma::log(Phi.row(z)) + lp_Zup.row(j);
        arma::rowvec lp_subk = arma::log(Phi.row(z)) + lp_Zup.row(k);      
        lp_Zup(i,z) = log_emission(i,z) + logaddexp(lp_subj) + logaddexp(lp_subk);
      }
    }
    // about to sample Z, record i where Z for i -> j changed
    std::fill(update_node.begin(), update_node.end(), false);
    // sample Z root to tip
    arma::rowvec lp_root_z = lp_Zup.row(R) + arma::log(prior_ZR);
    Z[R] = sample_ulp(lp_root_z);
    // sample all hidden nodes from root to tip
    for (int e=E-1; e>=0; e--) {
      int i = edge(e, 0); // parent
      int j = edge(e, 1); // child, we visit 'to end' of edges
      arma::rowvec lp_jz = lp_Zup.row(j) + arma::log(Phi.row(Z[i]));
      int old_Z = Z[j]; // to see if updated
      Z[j] = sample_ulp(lp_jz); // update
      if (old_Z != Z[j]) { update_node[i] = true; } // Did update occur? If so, parent in DP will be updated.
    }
    // Which nodes to update in DP based on changed Z? This includes parents of parents etc.
    for (int e = 0; e < E; e++) {
      int i = edge(e, 0);
      int j = edge(e, 1);
      if (update_node[j]) { update_node[i] = true; }
    }
    // sample r
    double ll  = 0; // under current rates
    for (int s = 0; s < S; s++) {
      arma::mat& lp_Xup = lp_Xups->slice(s);
      // use DP to calculate prob of alignment above node i
      change_aware_dp_Xup(lp_Xup, s, acgt, L, E, X, edge, Z, P_nuctrans, update_node);
      ll += logaddexp(arma::log(PI_D) + lp_Xup.row(R));
    }
    // propose new rates under MH scheme
    const double prop_sigma = 0.1;
    r_(0) = r(0);
    r_(1) = exp(log(r[1]) + R::rnorm(0, prop_sigma));
    r_(2) = exp(log(r[2]) + R::rnorm(0, prop_sigma));
    // fill transition probabilities cache using proposed rates
    for (int e = 0; e < E; e++) {
      P_nuctrans_[1]->slice(e) = arma::expmat(Q * r_(1) * distance(e));
      P_nuctrans_[2]->slice(e) = arma::expmat(Q * r_(2) * distance(e));
    }
    // calculate prob tree under proposed rates
    double ll_ = 0; // under proposed rates
    for (int s = 0; s < S; s++) {
      arma::mat& lp_Xup = lp_Xups_->slice(s);
      // use DP to calculate prob of alignment above node i
      change_aware_dp_Xup(lp_Xup, s, acgt, L, E, X, edge, Z, P_nuctrans_, update_everything);
      ll_ += logaddexp(arma::log(PI_D) + lp_Xup.row(R));
    }
    // Metropolis acm of nucleotide substitution rates
    // N.B. these priors would be written e.g. dgamma(x, prior_r3_a, scale=prior_r3_b)
    // where scale = 1/rate and rate is default R parameterization, but from Cpp
    // call sale is the 2nd argument, so no reciprocal is needed.
    const double lp_rates  = R::dgamma(r [1], prior_r2_a, prior_r2_b, true) + R::dgamma(r [2], prior_r3_a, prior_r3_b, true);
    const double lp_rates_ = R::dgamma(r_[1], prior_r2_a, prior_r2_b, true) + R::dgamma(r_[2], prior_r3_a, prior_r3_b, true);
    double log_mh_ratio = ll_ + lp_rates_ - ll - lp_rates;

    const double ap = std::min(1.0, exp(log_mh_ratio));
    if (R::runif(0, 1) < ap) {
      r[0] = r_[0]; r[1] = r_[1]; r[2] = r_[2]; // update rates
      std::swap(P_nuctrans, P_nuctrans_); // use new cache on next iter
      std::swap(lp_Xups, lp_Xups_); // use new cache on next iter
    }
    
    // Sample alpha/beta for transition matrix
    int n11 = 0;
    int n12 = 0;
    int n22 = 0;
    int n23 = 0;
    int n32 = 0; // no-Dollo assumption
    int n33 = 0; // no-Dollo assumption
    for (int e=0; e<E; e++) {
      int Z_i = Z(edge(e,0));
      int Z_j = Z(edge(e,1));
      n11 += Z_i==0 && Z_j==0 ? 1 : 0;
      n12 += Z_i==0 && Z_j==1 ? 1 : 0;
      n22 += Z_i==1 && Z_j==1 ? 1 : 0;
      n23 += Z_i==1 && Z_j==2 ? 1 : 0;
      n32 += Z_i==2 && Z_j==1 ? 1 : 0; // no-Dollo assumption
      n33 += Z_i==2 && Z_j==2 ? 1 : 0; // no-Dollo assumption
    }
    phi_a = R::rbeta(n12 + prior_gc_a, n11 + prior_gc_b);
    phi_b = R::rbeta(n23 + prior_lc_a, n22 + prior_lc_b);
    phi_ac= R::rbeta(n32 + prior_ac_a, n33 + prior_ac_b);
    stuff_Phi_matrix(Phi, phi_a, phi_b, phi_ac); // no-Dollo assumption
    // record state
    Y_trace.row(iter) = Y.t();
    v_trace(iter, 0) = lv;
    v_trace(iter, 1) = lb2; 
    v_trace(iter, 2) = lb3; 
    z_trace.row(iter) = Z.t();
    r_trace.row(iter) = r.t();
    phi_trace(iter,0) = phi_a;
    phi_trace(iter,1) = phi_b;
    phi_trace(iter,2) = phi_ac; // no-Dollo assumption
  }
  return Rcpp::List::create(Rcpp::Named("z_trace") = z_trace,
                            Rcpp::Named("r_trace") = r_trace,
                            Rcpp::Named("Phi_trace") = phi_trace,
                            Rcpp::Named("y_trace") = Y_trace,
                            Rcpp::Named("v_trace") = v_trace);
}

