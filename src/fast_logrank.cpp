#include <RcppEigen.h>
#include <vector>
#include <algorithm>
#include <numeric>

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
List fast_logrank(NumericVector time, IntegerVector event, IntegerVector group) {
    int n = time.size();

    // Determine number of groups (expect 1..K)
    int K = *std::max_element(group.begin(), group.end());

    // Build sort order: ascending time, events before censorings at ties
    std::vector<int> ord(n);
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](int a, int b) {
        if (time[a] != time[b]) return time[a] < time[b];
        return event[a] > event[b]; // events first
    });

    // Accumulators
    Eigen::VectorXd observed = Eigen::VectorXd::Zero(K);
    Eigen::VectorXd expected = Eigen::VectorXd::Zero(K);
    Eigen::MatrixXd variance = Eigen::MatrixXd::Zero(K, K);
    Eigen::VectorXd n_by_group = Eigen::VectorXd::Zero(K);

    // Count per-group totals (= initial at-risk counts)
    Eigen::VectorXd atrisk(K);
    atrisk.setZero();
    for (int i = 0; i < n; i++) {
        atrisk[group[i] - 1] += 1.0;
        n_by_group[group[i] - 1] += 1.0;
    }

    // Walk through sorted observations, processing at each distinct event time
    int i = 0;
    while (i < n) {
        int idx = ord[i];
        double t_cur = time[idx];

        // Count events and censorings at this time point per group
        Eigen::VectorXd d_k = Eigen::VectorXd::Zero(K);
        int n_cens = 0;
        std::vector<int> cens_groups;

        int j = i;
        // First pass: events at t_cur (they come first due to sort order)
        while (j < n && time[ord[j]] == t_cur && event[ord[j]] == 1) {
            d_k[group[ord[j]] - 1] += 1.0;
            j++;
        }
        // Second pass: censorings at t_cur
        while (j < n && time[ord[j]] == t_cur) {
            cens_groups.push_back(group[ord[j]] - 1);
            n_cens++;
            j++;
        }

        double d_total = d_k.sum();
        double n_total = atrisk.sum();

        // Accumulate O, E, V only if there are events at this time
        if (d_total > 0 && n_total > 0) {
            observed += d_k;

            for (int k = 0; k < K; k++) {
                expected[k] += atrisk[k] * d_total / n_total;
            }

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
                        variance(k, l) += v;
                        if (k != l) variance(l, k) += v;
                    }
                }
            }
        }

        // Remove events and censorings from at-risk set
        for (int k = 0; k < K; k++) {
            atrisk[k] -= d_k[k];
        }
        for (int c = 0; c < n_cens; c++) {
            atrisk[cens_groups[c]] -= 1.0;
        }

        i = j;
    }

    // Convert Eigen objects to R
    NumericVector r_obs(K), r_exp(K), r_n(K);
    for (int k = 0; k < K; k++) {
        r_obs[k] = observed[k];
        r_exp[k] = expected[k];
        r_n[k] = n_by_group[k];
    }

    NumericMatrix r_var(K, K);
    for (int k = 0; k < K; k++) {
        for (int l = 0; l < K; l++) {
            r_var(k, l) = variance(k, l);
        }
    }

    return List::create(
        Named("observed") = r_obs,
        Named("expected") = r_exp,
        Named("variance") = r_var,
        Named("n_by_group") = r_n
    );
}
