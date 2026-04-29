#include <RcppEigen.h>
#include <vector>
#include <algorithm>
#include <numeric>

using namespace Rcpp;
using namespace Eigen;


struct LogrankResult {
    Eigen::VectorXd observed;
    Eigen::VectorXd expected;
    Eigen::MatrixXd variance;
    Eigen::VectorXd n_by_group;
};


LogrankResult fast_logrank_internal(const double* time, const int* event, const int* group, int n, int K) {
    // Build sort order: ascending time, events before censorings at ties
    std::vector<int> ord(n);
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](int a, int b) {
        if (time[a] != time[b]) return time[a] < time[b];
        return event[a] > event[b]; // events first
    });

    LogrankResult res;
    res.observed = Eigen::VectorXd::Zero(K);
    res.expected = Eigen::VectorXd::Zero(K);
    res.variance = Eigen::MatrixXd::Zero(K, K);
    res.n_by_group = Eigen::VectorXd::Zero(K);

    // Count per-group totals (= initial at-risk counts)
    Eigen::VectorXd atrisk = Eigen::VectorXd::Zero(K);

    for (int i = 0; i < n; i++) {
        atrisk[group[i] - 1] += 1.0;
        res.n_by_group[group[i] - 1] += 1.0;
    }

    // Walk through sorted observations, processing at each distinct event time
    std::vector<int> cens_groups;
    cens_groups.reserve(n);
    Eigen::VectorXd d_k = Eigen::VectorXd::Zero(K);
    int i = 0;
    while (i < n) {
        d_k.setZero();
        int n_cens = 0;
        cens_groups.clear();        

        // Count events and censored at this time point per group
        double t_cur = time[ord[i]];
        int j = i;
        // First pass: events at t_cur (they come first due to sort order)
        while (j < n && time[ord[j]] == t_cur && event[ord[j]] == 1) {
            d_k[group[ord[j]] - 1] += 1.0;
            j++;
        }
        // Second pass: censored at t_cur
        while (j < n && time[ord[j]] == t_cur) {
            cens_groups.push_back(group[ord[j]] - 1);
            n_cens++;
            j++;
        }

        double d_total = d_k.sum();
        double n_total = atrisk.sum();

        // Accumulate O, E, V only if there are events at this time
        if (d_total > 0 && n_total > 0) {
            res.observed += d_k;
            res.expected += atrisk * (d_total / n_total);

            // Hypergeometric variance contribution
            // V_kl = d * (n - d) / (n^2 * (n - 1)) * (n_k * I(k==l) * n - n_k * n_l)
            if (n_total > 1) {
                double factor = d_total * (n_total - d_total) /
                                (n_total * n_total * (n_total - 1.0));
                for (int k = 0; k < K; k++) {
                    for (int l = k; l < K; l++) {
                        double v;
                        if (k == l) {
                            v = factor * atrisk[k] * (n_total - atrisk[k]);
                        } else {
                            v = -factor * atrisk[k] * atrisk[l];
                        }
                        res.variance(k, l) += v;
                        if (k != l) res.variance(l, k) += v;
                    }
                }
            }
        }

        // Remove events and censored from at-risk set
        for (int k = 0; k < K; k++) {
            atrisk[k] -= d_k[k];
        }
        for (int c = 0; c < n_cens; c++) {
            atrisk[cens_groups[c]] -= 1.0;
        }

        i = j;
    }
    return res;
}





// [[Rcpp::export]]
List hly_sample (NumericMatrix S_ij, NumericVector left_times, IntegerVector group_var_s, int H) {
    int n = S_ij.nrow();
    int m = S_ij.ncol();
    int K = *std::max_element(group_var_s.begin(), group_var_s.end());

    Eigen::VectorXd U_sum = Eigen::VectorXd::Zero(K);
    Eigen::MatrixXd UUt_sum = Eigen::MatrixXd::Zero(K, K);
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(K, K);

    std::vector<double> imputed_times(n);
    std::vector<int> events(n, 1);
    std::vector<int> groups(group_var_s.begin(), group_var_s.end());


    for (int h = 0; h < H; h++) {
        if (h % 1000 == 0) Rcpp::checkUserInterrupt();
        for (int i = 0; i < n; i++) {
            double q = R::runif(0.0, 1.0);
            int j = 0;
            // increment j until S_ij(i, j) > q, then impute time as left_times[j]
            while (j < m && S_ij(i, j) <= q) j++;
            if (j == m) j = m - 1;
            imputed_times[i] = left_times[j];
        }
        LogrankResult lr_res = fast_logrank_internal(imputed_times.data(), events.data(), groups.data(), n, K);
        Eigen::VectorXd u_h = lr_res.observed - lr_res.expected;
        U_sum += u_h;
        UUt_sum += u_h * u_h.transpose();
        V += lr_res.variance / H;
    }

    Eigen::VectorXd U_mean = U_sum / H;
    Eigen::MatrixXd V_s = V - (UUt_sum - H * U_mean * U_mean.transpose()) / (H - 1);

    return List::create(
        Named("U") = U_mean,
        Named("V") = V_s
    );
}



// [[Rcpp::export]]
List fast_logrank(NumericVector time, IntegerVector event, IntegerVector group) {
    int n = time.size();
    int K = *std::max_element(group.begin(), group.end());
    LogrankResult res = fast_logrank_internal(time.begin(), event.begin(), group.begin(), n, K);
    return List::create(
        Named("observed") = res.observed,
        Named("expected") = res.expected,
        Named("variance") = res.variance,
        Named("n_by_group") = res.n_by_group
    );
}
