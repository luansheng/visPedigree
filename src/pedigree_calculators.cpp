#include <RcppArmadillo.h>
#include <vector>
#include <map>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// ============================================================================
// OpenMP Thread Control Functions
// ============================================================================
// These functions provide runtime control over the number of OpenMP threads
// used in parallel calculations. They are thread-safe and handle the case
// where OpenMP is not available.
// ============================================================================

// OpenMP thread controls
// [[Rcpp::export]]
int cpp_set_num_threads(int threads) {
#ifdef _OPENMP
    int prev = omp_get_max_threads();
    if (threads > 0) {
        omp_set_num_threads(threads);
    }
    return prev;
#else
    return 1;
#endif
}

// [[Rcpp::export]]
bool cpp_openmp_available() {
#ifdef _OPENMP
    return true;
#else
    return false;
#endif
}

// ============================================================================
// Calculate Inbreeding Coefficients and Dii (Diagonal of D⁻¹)
// ============================================================================
// Algorithm: Meuwissen & Luo (1992) - Path tracing with caching
// Complexity: O(n²) in worst case, O(n) for simple pedigrees
// Parallelization: Not implemented due to recursive dependencies
// ============================================================================
// Calculate Inbreeding and Dii using Meuwissen & Luo (1992)
// [[Rcpp::export]]
List cpp_calculate_inbreeding(IntegerVector sire, IntegerVector dam) {
    int n = sire.size();
    if (dam.size() != n) {
        stop("sire and dam vectors must have the same length");
    }
    
    NumericVector f(n);
    NumericVector dii(n);
    std::vector<double> a(n, 1.0); // a_ii = 1 + f_i
    
    // Work vector L and visited tracking
    std::vector<double> L(n, 0.0);
    std::vector<int> visited;
    visited.reserve(2048);
    
    // Cache for parent pair results: key is (min(s,d), max(s,d))
    std::map<std::pair<int, int>, double> pair_cache;

    for (int i = 0; i < n; ++i) {
        int s = sire[i] - 1;
        int d = dam[i] - 1;

        if (s >= n || d >= n) {
            stop("Parent index out of bounds");
        }

        // 1. Full-sib optimization
        if (i > 0 && sire[i] == sire[i-1] && dam[i] == dam[i-1]) {
            f[i] = f[i-1];
            dii[i] = dii[i-1];
            a[i] = a[i-1];
            continue; 
        }
        
        // 2. Compute dii
        double fs = (s >= 0) ? f[s] : 0.0;
        double fd = (d >= 0) ? f[d] : 0.0;
        if (s >= 0 && d >= 0) {
            dii[i] = 0.5 - 0.25 * (fs + fd);
        } else if (s >= 0 || d >= 0) {
            dii[i] = 0.75 - 0.25 * (s >= 0 ? fs : fd);
        } else {
            dii[i] = 1.0;
        }
        
        // 3. Compute f[i] using path trace
        if (s >= 0 && d >= 0) {
            std::pair<int, int> parents = (s < d) ? std::make_pair(s, d) : std::make_pair(d, s);
            auto it = pair_cache.find(parents);
            if (it != pair_cache.end()) {
                a[i] = it->second;
                f[i] = a[i] - 1.0;
            } else {
                visited.clear();
                L[i] = 1.0;
                visited.push_back(i);
                
                double sum_L2d = 0.0;
                // Meuwissen & Luo algorithm: trace back
                for (int k = i; k >= 0; --k) {
                    if (L[k] == 0) continue;
                    
                    int sk = sire[k] - 1;
                    int dk = dam[k] - 1;
                    
                    if (sk >= 0) {
                        if (L[sk] == 0) visited.push_back(sk);
                        L[sk] += 0.5 * L[k];
                    }
                    if (dk >= 0) {
                        if (L[dk] == 0) visited.push_back(dk);
                        L[dk] += 0.5 * L[k];
                    }
                    sum_L2d += L[k] * L[k] * dii[k];
                }
                a[i] = sum_L2d;
                f[i] = a[i] - 1.0;
                pair_cache[parents] = a[i];

                for (int idx : visited) L[idx] = 0.0;
            }
        } else {
            f[i] = 0.0;
            a[i] = 1.0;
        }
    }
    
    return List::create(
        Named("f") = f, 
        Named("dii") = dii
    );
}

// ============================================================================
// Build A-Inverse Sparse Matrix (Henderson's Rules)
// ============================================================================
// Algorithm: Henderson (1976) direct construction of A⁻¹
// Complexity: O(n) - highly efficient, no matrix inversion needed
// Parallelization: Enabled for n >= 5000 with thread-local storage
// Performance: 
//   - Serial preferred for n < 5000 (avoids parallel overhead)
//   - Parallel achieves ~1.3-1.5x speedup for large pedigrees (n > 10000)
// ============================================================================
// Build A-Inverse Sparse Matrix Components (Henderson's Rules)
// [[Rcpp::export]]
List cpp_build_ainv_triplets(IntegerVector sire, IntegerVector dam, NumericVector dii) {
    int n = sire.size();
    if (dam.size() != n || dii.size() != n) {
        stop("sire, dam, and dii vectors must have the same length");
    }
    
    // Check bounds before parallel region
    for (int i = 0; i < n; ++i) {
        if (sire[i] > n || dam[i] > n) {
            stop("Parent index out of bounds");
        }
    }
    
    // Adaptive threshold: use serial for small pedigrees to avoid parallel overhead
    const int PARALLEL_THRESHOLD = 5000;
    
    if (n < PARALLEL_THRESHOLD) {
        // Serial version for small pedigrees (avoids thread creation/merge overhead)
        std::vector<int> row, col;
        std::vector<double> val;
        
        // Pre-allocate: each individual generates ~6 triplets on average
        int est_size = n * 6;
        row.reserve(est_size);
        col.reserve(est_size);
        val.reserve(est_size);
        
        for (int i = 0; i < n; ++i) {
            double alpha = 1.0 / dii[i];
            int id = i + 1; 
            int s = sire[i];
            int d = dam[i];
            
            // Diagonal element
            row.push_back(id); col.push_back(id); val.push_back(alpha);
            
            // Sire contributions
            if (s > 0) {
                int r = (id > s) ? id : s;
                int c = (id > s) ? s : id;
                row.push_back(r); col.push_back(c); val.push_back(-0.5 * alpha);
                row.push_back(s); col.push_back(s); val.push_back(0.25 * alpha);
            }
            
            // Dam contributions
            if (d > 0) {
                int r = (id > d) ? id : d;
                int c = (id > d) ? d : id;
                row.push_back(r); col.push_back(c); val.push_back(-0.5 * alpha);
                row.push_back(d); col.push_back(d); val.push_back(0.25 * alpha);
            }
            
            // Sire-Dam interaction
            if (s > 0 && d > 0) {
                int r = (s > d) ? s : d;
                int c = (s > d) ? d : s;
                row.push_back(r); col.push_back(c); val.push_back(0.25 * alpha);
            }
        }
        
        return List::create(
            Named("i") = wrap(row), 
            Named("j") = wrap(col), 
            Named("v") = wrap(val)
        );
    }
    
    // Parallel version for large pedigrees (n >= 5000)
    #ifdef _OPENMP
    int max_threads = omp_get_max_threads();
    #else
    int max_threads = 1;
    #endif

    std::vector<std::vector<int>> thread_row(max_threads);
    std::vector<std::vector<int>> thread_col(max_threads);
    std::vector<std::vector<double>> thread_val(max_threads);
    
    for(int t=0; t<max_threads; ++t) {
        thread_row[t].reserve((n / max_threads) * 10);
        thread_col[t].reserve((n / max_threads) * 10);
        thread_val[t].reserve((n / max_threads) * 10);
    }

    #pragma omp parallel
    {
        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif
        
        #pragma omp for
        for (int i = 0; i < n; ++i) {
            double alpha = 1.0 / dii[i];
            int id = i + 1; 
            int s = sire[i];
            int d = dam[i];
            
            thread_row[tid].push_back(id); thread_col[tid].push_back(id); thread_val[tid].push_back(alpha);
            
            if (s > 0) {
                int r = (id > s) ? id : s;
                int c = (id > s) ? s : id;
                thread_row[tid].push_back(r); thread_col[tid].push_back(c); thread_val[tid].push_back(-0.5 * alpha);
                thread_row[tid].push_back(s); thread_col[tid].push_back(s); thread_val[tid].push_back(0.25 * alpha);
            }
            if (d > 0) {
                int r = (id > d) ? id : d;
                int c = (id > d) ? d : id;
                thread_row[tid].push_back(r); thread_col[tid].push_back(c); thread_val[tid].push_back(-0.5 * alpha);
                thread_row[tid].push_back(d); thread_col[tid].push_back(d); thread_val[tid].push_back(0.25 * alpha);
            }
            if (s > 0 && d > 0) {
                int r = (s > d) ? s : d;
                int c = (s > d) ? d : s;
                thread_row[tid].push_back(r); thread_col[tid].push_back(c); thread_val[tid].push_back(0.25 * alpha);
            }
        }
    }
    
    int total_size = 0;
    for(int t=0; t<max_threads; ++t) total_size += thread_row[t].size();
    
    IntegerVector final_row(total_size);
    IntegerVector final_col(total_size);
    NumericVector final_val(total_size);
    
    int offset = 0;
    for(int t=0; t<max_threads; ++t) {
        std::copy(thread_row[t].begin(), thread_row[t].end(), final_row.begin() + offset);
        std::copy(thread_col[t].begin(), thread_col[t].end(), final_col.begin() + offset);
        std::copy(thread_val[t].begin(), thread_val[t].end(), final_val.begin() + offset);
        offset += thread_row[t].size();
    }
    
    return List::create(Named("i") = final_row, Named("j") = final_col, Named("v") = final_val);
}

// ============================================================================
// Calculate Additive Relationship Matrix (A)
// ============================================================================
// Algorithm: Tabular method with full-sibling optimization
// Complexity: O(n²) for dense matrix construction
// Parallelization: Not implemented
// Reason: Strong row dependencies (A[i,j] depends on A[j, parents])
//         Potential inner-loop parallelization would only yield ~1.2-1.5x speedup
//         while losing full-sib optimization benefits
// ============================================================================
// Calculate Additive Relationship Matrix (A) using Armadillo
// [[Rcpp::export]]
arma::mat cpp_calculate_A(IntegerVector sire, IntegerVector dam) {
    int n = sire.size();
    if (dam.size() != n) {
        stop("sire and dam vectors must have the same length");
    }
    
    // Protect against massive memory allocation
    if (n > 20000) {
        stop("Pedigree too large for dense A matrix calculation (n > 20000). Use sparse methods instead.");
    }
    
    arma::mat A(n, n, arma::fill::zeros);
    for (int i = 0; i < n; ++i) {
        int si = sire[i] - 1;
        int di = dam[i] - 1;
        
        if (si >= n || di >= n) {
            stop("Parent index out of bounds");
        }
        
        double fi = (si >= 0 && di >= 0) ? 0.5 * A(si, di) : 0.0;
        A(i, i) = 1.0 + fi;
        if (i > 0 && sire[i] == sire[i-1] && dam[i] == dam[i-1]) {
            if (i > 1) {
                A.row(i).subvec(0, i-2) = A.row(i-1).subvec(0, i-2);
                A.col(i).subvec(0, i-2) = A.col(i-1).subvec(0, i-2);
            }
            double val = 0.5 * ((si >= 0 ? A(i-1, si) : 0.0) + (di >= 0 ? A(i-1, di) : 0.0));
            A(i, i-1) = A(i-1, i) = val;
            continue;
        }
        for (int j = 0; j < i; ++j) {
            double val = 0.5 * ((si >= 0 ? A(j, si) : 0.0) + (di >= 0 ? A(j, di) : 0.0));
            A(i, j) = A(j, i) = val;
        }
    }
    return A;
}

// Calculate mean off-diagonal relationship for a target subset
// [[Rcpp::export]]
double cpp_mean_relationship(IntegerVector sire, IntegerVector dam, IntegerVector target_idx) {
    int n = sire.size();
    if (dam.size() != n) {
        stop("sire and dam vectors must have the same length");
    }
    if (target_idx.size() < 2) {
        return NA_REAL;
    }

    std::vector<unsigned char> is_target(n, 0);
    int n_target = 0;
    for (int k = 0; k < target_idx.size(); ++k) {
        int idx = target_idx[k] - 1;
        if (idx < 0 || idx >= n) {
            stop("target_idx contains out-of-bounds index");
        }
        if (!is_target[idx]) {
            is_target[idx] = 1;
            ++n_target;
        }
    }
    if (n_target < 2) {
        return NA_REAL;
    }

    const size_t tri_size = static_cast<size_t>(n) * static_cast<size_t>(n + 1) / 2;
    std::vector<double> A_tri(tri_size, 0.0);
    auto tri_idx = [](int i, int j) -> size_t {
        if (i < j) std::swap(i, j);
        return static_cast<size_t>(i) * static_cast<size_t>(i + 1) / 2 + static_cast<size_t>(j);
    };

    std::vector<int> seen_targets;
    seen_targets.reserve(n_target);
    double lower_sum = 0.0;

    for (int i = 0; i < n; ++i) {
        int si = sire[i] - 1;
        int di = dam[i] - 1;

        if (si >= n || di >= n) {
            stop("Parent index out of bounds");
        }

        double fi = (si >= 0 && di >= 0) ? 0.5 * A_tri[tri_idx(si, di)] : 0.0;
        A_tri[tri_idx(i, i)] = 1.0 + fi;

        for (int j = 0; j < i; ++j) {
            double val = 0.5 * ((si >= 0 ? A_tri[tri_idx(j, si)] : 0.0) +
                                (di >= 0 ? A_tri[tri_idx(j, di)] : 0.0));
            A_tri[tri_idx(i, j)] = val;
        }

        if (is_target[i]) {
            for (int tj : seen_targets) {
                lower_sum += A_tri[tri_idx(i, tj)];
            }
            seen_targets.push_back(i);
        }
    }

    return (2.0 * lower_sum) / (static_cast<double>(n_target) * static_cast<double>(n_target - 1));
}

// ============================================================================
// Calculate Dominance Relationship Matrix (D)
// ============================================================================
// Algorithm: D_ij = 0.25 * (A_si,sj * A_di,dj + A_si,dj * A_di,sj)
// Complexity: O(n²) for dense matrix construction
// Parallelization: ✅ Fully parallelized with OpenMP
// Performance: 
//   - Threshold: n > 500 (avoids overhead for small matrices)
//   - Speedup: ~1.46-1.71x on 4 cores for n > 4000
//   - Strategy: Row-level parallelism (each row independent)
//   - Schedule: dynamic with chunk size 32 for load balancing
// Note: Full-sibling optimization disabled for thread safety
// ============================================================================
// Calculate Dominance Matrix (D) using Armadillo
// Parallelized version - each row computed independently
// [[Rcpp::export]]
arma::mat cpp_calculate_D(IntegerVector sire, IntegerVector dam, const arma::mat& A) {
    int n = sire.size();
    if (dam.size() != n) {
        stop("sire and dam vectors must have the same length");
    }
    if (A.n_rows != n || A.n_cols != n) {
        stop("A matrix dimensions must match pedigree size");
    }
    
    // Protect against massive memory allocation
    if (n > 20000) {
        stop("Pedigree too large for dense D matrix calculation (n > 20000).");
    }
    
    arma::mat D(n, n, arma::fill::zeros); 
    
    // Parallel computation of D matrix
    // Each D(i,j) depends only on A values (read-only), so rows are independent
    // Threshold: n > 500 to avoid parallel overhead for small matrices
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 32) if(n > 500)
    #endif
    for (int i = 0; i < n; ++i) {
        // Diagonal: D_ii = 2 - A_ii = 1 - F_i
        D(i, i) = 2.0 - A(i, i);
        
        int si = sire[i] - 1;
        int di = dam[i] - 1;
        
        // Skip if either parent is missing (founders have D_ij = 0 for i≠j)
        if (si < 0 || di < 0) continue;
        
        // We don't throw inside OpenMP, but we can skip invalid indices
        if (si >= n || di >= n) continue;
        
        // Compute off-diagonal elements: D_ij = 0.25 * (A_si,sj * A_di,dj + A_si,dj * A_di,sj)
        // Note: Full-sib optimization disabled for thread safety
        for (int j = 0; j < i; ++j) {
            int sj = sire[j] - 1;
            int dj = dam[j] - 1;
            
            if (sj >= 0 && dj >= 0 && sj < n && dj < n) {
                double val = 0.25 * (A(si, sj) * A(di, dj) + A(si, dj) * A(di, sj));
                D(i, j) = val;
                D(j, i) = val;
            }
        }
    }
    return D;
}

// ============================================================================
// Calculate Epistatic Relationship Matrix (AA = A ⊙ A)
// ============================================================================
// Algorithm: Hadamard (element-wise) product of A with itself
// Complexity: O(n²) for dense matrices
// Parallelization: Depends on Armadillo/BLAS implementation
// Note: Modern BLAS libraries (OpenBLAS, MKL) may auto-parallelize
//       Performance is primarily limited by memory bandwidth
// ============================================================================
// Hadamard Product
// [[Rcpp::export]]
arma::mat cpp_calculate_AA(const arma::mat& A) {
    return A % A;
}

// Invert Dense Matrix (general purpose)
// [[Rcpp::export]]
arma::mat cpp_invert_dense(const arma::mat& M) {
    return arma::inv(M);
}

// Invert Symmetric Positive-Definite Matrix (optimized)
// Uses Cholesky decomposition, approximately 2x faster than general inversion
// for symmetric positive-definite matrices
// [[Rcpp::export]]
arma::mat cpp_invert_sympd(const arma::mat& M) {
    return arma::inv_sympd(M);
}

// Auto-detect matrix type and use optimal inversion method
// Checks if matrix is symmetric, then tries Cholesky (for positive-definite)
// Falls back to general LU decomposition if Cholesky fails
// [[Rcpp::export]]
arma::mat cpp_invert_auto(const arma::mat& M) {
    int n = M.n_rows;
    
    // Check if matrix is square
    if (M.n_rows != M.n_cols) {
        Rcpp::stop("Matrix must be square for inversion");
    }
    
    // Check symmetry with tolerance
    double tol = 1e-10 * n;
    bool is_symmetric = arma::approx_equal(M, M.t(), "absdiff", tol);
    
    if (is_symmetric) {
        // Try Cholesky decomposition for symmetric positive-definite matrices
        // inv_sympd() will throw error if matrix is not positive-definite
        try {
            return arma::inv_sympd(M);
        } catch(...) {
            // If Cholesky fails, fall back to general inversion
            // (matrix might be symmetric but not positive-definite)
            Rcpp::warning("Cholesky decomposition failed, using general LU inversion");
            return arma::inv(M);
        }
    } else {
        // Non-symmetric matrix: use general LU decomposition
        return arma::inv(M);
    }
}

// Solve A*x = b using Path Logic
// [[Rcpp::export]]
arma::vec cpp_solve_A(IntegerVector sire, IntegerVector dam, NumericVector dii, arma::vec b) {
    int n = sire.size();
    if (dam.size() != n || dii.size() != n || b.size() != n) {
        stop("sire, dam, dii, and b vectors must have the same length");
    }
    
    arma::vec x = b;
    for (int i = n - 1; i >= 0; --i) {
        int s = sire[i] - 1; int d = dam[i] - 1;
        if (s >= n || d >= n) stop("Parent index out of bounds");
        if (s >= 0) x[s] -= 0.5 * x[i];
        if (d >= 0) x[d] -= 0.5 * x[i];
    }
    for (int i = 0; i < n; ++i) x[i] *= dii[i];
    for (int i = 0; i < n; ++i) {
        int s = sire[i] - 1; int d = dam[i] - 1;
        if (s >= 0) x[i] += 0.5 * x[s];
        if (d >= 0) x[i] += 0.5 * x[d];
    }
    return x;
}

// Generations Top-down
// [[Rcpp::export]]
IntegerVector cpp_assign_generations_top(IntegerVector sire, IntegerVector dam, IntegerVector topo_order) {
    int n = sire.size();
    if (dam.size() != n || topo_order.size() != n) {
        stop("sire, dam, and topo_order vectors must have the same length");
    }
    
    IntegerVector gen(n, 1);
    for (int i = 0; i < n; ++i) {
        int idx = topo_order[i] - 1;
        if (idx < 0 || idx >= n) stop("topo_order index out of bounds");
        
        int s = sire[idx] - 1; int d = dam[idx] - 1;
        if (s >= n || d >= n) stop("Parent index out of bounds");
        
        int sg = (s >= 0) ? gen[s] : 0;
        int dg = (d >= 0) ? gen[d] : 0;
        if (s >= 0 || d >= 0) gen[idx] = (sg > dg ? sg : dg) + 1;
    }
    return gen;
}

// ============================================================================
// Calculate Partial Inbreeding (pF)
// ============================================================================
// Decomposes the inbreeding coefficient into contributions from specific ancestors.
// Based on the recursive path tracing property of the Meuwissen & Luo algorithm.
// ============================================================================
// [[Rcpp::export]]
NumericMatrix cpp_calculate_partial_inbreeding(IntegerVector sire, IntegerVector dam, NumericVector dii, IntegerVector ancestors) {
    int n = sire.size();
    int n_anc = ancestors.size();
    NumericMatrix pF(n, n_anc);
    
    // We use a reusable flow vector to avoid reallocation
    std::vector<double> flow(n);

    // DEBUG
    // Rcpp::Rcout << "Calculating partial inbreeding for " << n_anc << " ancestors. N=" << n << "\n";

    // For each target ancestor k
    for (int j = 0; j < n_anc; ++j) {
        int anc_idx = ancestors[j] - 1; // 0-based index
        // Rcpp::Rcout << "Ancestor " << j << " Index: " << anc_idx << "\n";
        
        // 1. Calculate gene flow from ancestor k to all descendents (Forward Pass)
        std::fill(flow.begin(), flow.end(), 0.0);
        flow[anc_idx] = 1.0;
        
        // Iterate from ancestor to end of pedigree
        // Assumes pedigree is sorted or at least parents appear before children
        // (Our tidyped guarantees partial order, but rigorous topological sort is best. 
        //  However, IndNum order usually respects birth date/generation).
        for (int i = anc_idx; i < n; ++i) {
            
            // This is "pull" logic (looking at parents), but we need "push" (parent to offspring).
            // But we don't have offspring list. 
            // So we stick to standard "pull" logic? No, pull calculates ancestor's contribution.
            // Wait, standard relationship calculation A_ik is exactly what we need.
            // A_ik = 0.5 * A_sk + 0.5 * A_dk.
            // This can be computed forward if we process i=1..N.
            
            int s = sire[i] - 1;
            int d = dam[i] - 1;
            
            if (i > anc_idx) { // Skip the ancestor itself
                double from_s = (s >= 0) ? flow[s] : 0.0;
                double from_d = (d >= 0) ? flow[d] : 0.0;
                flow[i] = 0.5 * (from_s + from_d);
            }
        }
        
        // 2. Compute partial inbreeding contribution
        // F_i(k) = 0.5 * flow[s] * flow[d] * dii[k]
        
        double dk = dii[anc_idx];
        
        for (int i = 0; i < n; ++i) {
            int s = sire[i] - 1;
            int d = dam[i] - 1;
            
            if (s >= 0 && d >= 0) {
                pF(i, j) = 0.5 * flow[s] * flow[d] * dk;
                // if (pF(i, j) > 0) Rcpp::Rcout << "  pF[" << i << "," << j << "] = " << pF(i, j) << "\n";
            }
        }
    }
    
    return pF;
}

// Generations Bottom-up
// [[Rcpp::export]]
IntegerVector cpp_assign_generations_bottom(IntegerVector sire, IntegerVector dam, IntegerVector topo_order) {
    int n = sire.size();
    if (dam.size() != n || topo_order.size() != n) {
        stop("sire, dam, and topo_order vectors must have the same length");
    }
    
    IntegerVector height(n, 0);
    for (int i = n - 1; i >= 0; --i) {
        int idx = topo_order[i] - 1;
        if (idx < 0 || idx >= n) stop("topo_order index out of bounds");
        
        int s = sire[idx] - 1; int d = dam[idx] - 1;
        if (s >= n || d >= n) stop("Parent index out of bounds");
        
        if (s >= 0) if (height[idx] + 1 > height[s]) height[s] = height[idx] + 1;
        if (d >= 0) if (height[idx] + 1 > height[d]) height[d] = height[idx] + 1;
    }
    int max_h = 0;
    for (int i = 0; i < n; ++i) if (height[i] > max_h) max_h = height[i];
    IntegerVector gen(n);
    for (int i = 0; i < n; ++i) gen[i] = max_h - height[i] + 1;
    return gen;
}

// ============================================================================
// Calculate Coancestry-based Ne (Delta c)
// ============================================================================
// [[Rcpp::export]]
double cpp_calculate_sampled_coancestry_delta(IntegerVector sire, IntegerVector dam, IntegerVector target_idx, NumericVector ecg) {
    int n = sire.size();
    if (dam.size() != n || ecg.size() != n) {
        stop("sire, dam, and ecg vectors must have the same length");
    }
    if (target_idx.size() < 2) return NA_REAL;
    std::vector<unsigned char> is_target(n, 0);
    int n_target = 0;
    for (int k = 0; k < target_idx.size(); ++k) {
        int idx = target_idx[k] - 1;
        if (idx < 0 || idx >= n) stop("target_idx contains out-of-bounds index");
        if (!is_target[idx]) { is_target[idx] = 1; ++n_target; }
    }
    if (n_target < 2) return NA_REAL;
    const size_t tri_size = static_cast<size_t>(n) * static_cast<size_t>(n + 1) / 2;
    std::vector<double> A_tri(tri_size, 0.0);
    auto tri_idx = [](int i, int j) -> size_t {
        if (i < j) std::swap(i, j);
        return static_cast<size_t>(i) * static_cast<size_t>(i + 1) / 2 + static_cast<size_t>(j);
    };
    std::vector<int> seen_targets;
    seen_targets.reserve(n_target);
    double sum_delta_c = 0.0;
    double num_pairs = 0.0;
    for (int i = 0; i < n; ++i) {
        int si = sire[i] - 1;
        int di = dam[i] - 1;
        if (si >= n || di >= n) stop("Parent index out of bounds");
        double fi = (si >= 0 && di >= 0) ? 0.5 * A_tri[tri_idx(si, di)] : 0.0;
        A_tri[tri_idx(i, i)] = 1.0 + fi;
        for (int j = 0; j < i; ++j) {
            double val = 0.5 * ((si >= 0 ? A_tri[tri_idx(j, si)] : 0.0) + (di >= 0 ? A_tri[tri_idx(j, di)] : 0.0));
            A_tri[tri_idx(i, j)] = val;
        }
        if (is_target[i]) {
            for (int tj : seen_targets) {
                double C_ij = A_tri[tri_idx(i, tj)] / 2.0;
                double g_ij = (ecg[i] + ecg[tj]) / 2.0;
                if (g_ij > 0.0) { // Note: changed from > 1.0 as Cervantes uses g_ij natively
                    double delta_c = 1.0 - std::pow(1.0 - C_ij, 1.0 / g_ij); 
                    sum_delta_c += delta_c;
                    num_pairs += 1.0;
                }
            }
            seen_targets.push_back(i);
        }
    }
    if (num_pairs == 0.0) return NA_REAL;
    return sum_delta_c / num_pairs;
}

// ============================================================================
// Founder and Ancestor Contributions (Boichard 1997)
// ============================================================================
// Algorithm:
//   Founders: single backward (gene-drop) pass, O(N).
//   Ancestors: Boichard's iterative peeling – select max-contribution ancestor,
//     sever its parental links, mask its descendants, repeat. Each iteration is
//     O(N) (backward pass + descendant scan).  Total O(K·N) for K ancestors.
// Inputs (all 1-based; 0 = unknown parent):
//   sire, dam        – integer parent vectors of length N
//   cohort_pos       – 1-based positions of reference (cohort) individuals
//   mode             – 1 = founder only, 2 = ancestor only, 3 = both
// Returns a named List whose elements depend on mode:
//   founder_idx, founder_contrib   (mode 1 or 3)
//   ancestor_idx, ancestor_contrib (mode 2 or 3)
// ============================================================================
// [[Rcpp::export]]
List cpp_pedcontrib(IntegerVector sire, IntegerVector dam,
                    IntegerVector cohort_pos, int mode) {

    const int n = sire.size();
    const int n_cohort = cohort_pos.size();
    if (dam.size() != n)
        stop("sire and dam vectors must have the same length");
    if (n_cohort == 0)
        stop("cohort_pos must not be empty");

    const double w = 1.0 / static_cast<double>(n_cohort);

    List result;

    // ------------------------------------------------------------------
    // Founder contributions (single backward pass)
    // ------------------------------------------------------------------
    if (mode == 1 || mode == 3) {
        std::vector<double> contrib(n, 0.0);
        for (int i = 0; i < n_cohort; ++i)
            contrib[cohort_pos[i] - 1] += w;

        for (int i = n - 1; i >= 0; --i) {
            if (contrib[i] == 0.0) continue;
            int s = sire[i] - 1;
            int d = dam[i] - 1;
            if (s >= 0) contrib[s] += 0.5 * contrib[i];
            if (d >= 0) contrib[d] += 0.5 * contrib[i];
        }

        // Collect founders (both parents unknown)
        std::vector<int>    f_idx;
        std::vector<double> f_val;
        for (int i = 0; i < n; ++i) {
            if (sire[i] == 0 && dam[i] == 0) {
                f_idx.push_back(i + 1);       // 1-based
                f_val.push_back(contrib[i]);
            }
        }
        result["founder_idx"]    = wrap(f_idx);
        result["founder_contrib"] = wrap(f_val);
    }

    // ------------------------------------------------------------------
    // Ancestor contributions – Boichard's iterative peeling
    // ------------------------------------------------------------------
    if (mode == 2 || mode == 3) {

        // --- BFS: discover all ancestors of the cohort ---
        std::vector<bool> visited(n, false);
        std::vector<bool> is_anc(n, false);
        std::vector<int>  bfs_queue;
        bfs_queue.reserve(n);

        for (int i = 0; i < n_cohort; ++i) {
            int pos = cohort_pos[i] - 1;
            if (!visited[pos]) {
                visited[pos] = true;
                is_anc[pos] = true; // COHORT MEMBERS Themselves are candidates!
                bfs_queue.push_back(pos);
            }
        }

        size_t head = 0;
        while (head < bfs_queue.size()) {
            int pos = bfs_queue[head++];
            int s = sire[pos] - 1;
            int d = dam[pos] - 1;
            if (s >= 0 && !visited[s]) {
                visited[s] = true;
                is_anc[s] = true;
                bfs_queue.push_back(s);
            }
            if (d >= 0 && !visited[d]) {
                visited[d] = true;
                is_anc[d] = true;
                bfs_queue.push_back(d);
            }
        }

        // Candidate ancestor list (0-based positions)
        std::vector<int> cand_idx;
        for (int i = 0; i < n; ++i)
            if (is_anc[i]) cand_idx.push_back(i);

        const int n_cand = static_cast<int>(cand_idx.size());

        if (n_cand == 0) {
            result["ancestor_idx"]    = IntegerVector(0);
            result["ancestor_contrib"] = NumericVector(0);
        } else {
            // Mutable parent arrays (will be modified by severing links)
            std::vector<int> w_sire(n), w_dam(n);
            for (int i = 0; i < n; ++i) {
                w_sire[i] = sire[i];
                w_dam[i]  = dam[i];
            }

            std::vector<bool> active(n_cand, true);
            std::vector<bool> is_sel(n, false);
            std::vector<int>    sel_idx;
            std::vector<double> sel_q;
            sel_idx.reserve(n_cand);
            sel_q.reserve(n_cand);

            // Reusable work arrays
            std::vector<double> contrib(n);
            std::vector<double> a_array(n, 0.0);

            for (;;) {
                // --- Backward pass with current (modified) parent links (Calculates Q) ---
                std::fill(contrib.begin(), contrib.end(), 0.0);
                for (int i = 0; i < n_cohort; ++i)
                    contrib[cohort_pos[i] - 1] += w;

                for (int i = n - 1; i >= 0; --i) {
                    if (contrib[i] == 0.0) continue;
                    int s = w_sire[i] - 1;
                    int d = w_dam[i]  - 1;
                    if (s >= 0) contrib[s] += 0.5 * contrib[i];
                    if (d >= 0) contrib[d] += 0.5 * contrib[i];
                }

                // --- Forward pass with updated parent links (Calculates A) ---
                std::fill(a_array.begin(), a_array.end(), 0.0);
                for (int i = 0; i < n; ++i) {
                    if (is_sel[i]) {
                        a_array[i] = 1.0;
                    } else {
                        int s = w_sire[i] - 1;
                        int d = w_dam[i]  - 1;
                        double a_s = (s >= 0) ? a_array[s] : 0.0;
                        double a_d = (d >= 0) ? a_array[d] : 0.0;
                        a_array[i] = 0.5 * a_s + 0.5 * a_d;
                    }
                }

                // --- Find best among active candidates ---
                int    best_j   = -1;
                double best_val = -1.0;
                for (int j = 0; j < n_cand; ++j) {
                    if (active[j]) {
                        int pos = cand_idx[j];
                        // Marginal contribution = q_k * (1 - a_k)
                        double p_k = contrib[pos] * std::max(0.0, 1.0 - a_array[pos]);
                        if (p_k > best_val) {
                            best_val = p_k;
                            best_j   = j;
                        }
                    }
                }
                
                if (best_val <= 1e-10) break;

                int best_pos = cand_idx[best_j];
                sel_idx.push_back(best_pos + 1);   // 1-based index
                sel_q.push_back(best_val);

                // Sever parent links of selected ancestor (to form pseudo-founder)
                w_sire[best_pos] = 0;
                w_dam[best_pos]  = 0;
                // Mark as selected so it stops receiving and begins transmitting its own a_array=1.0 downwards
                active[best_j]   = false;
                is_sel[best_pos] = true;
            }

            result["ancestor_idx"]    = wrap(sel_idx);
            result["ancestor_contrib"] = wrap(sel_q);
        }
    }

    return result;
}

// Calculate Ancestry (Forward P)
// [[Rcpp::export]]
NumericMatrix cpp_calculate_ancestry(IntegerVector sire, IntegerVector dam, NumericMatrix res_mat, IntegerVector ind_to_row) {
    int n = sire.size();
    int n_cols = res_mat.ncol();
    
    for(int i=0; i<n; ++i) {
       int s = sire[i] - 1;
       int d = dam[i] - 1;
       if (s >= 0 || d >= 0) {
           for(int k=0; k<n_cols; ++k) {
               double ps = (s >= 0) ? res_mat(s,k) : 0.0;
               double pd = (d >= 0) ? res_mat(d,k) : 0.0;
               res_mat(i, k) = 0.5 * (ps + pd);
           }
       }
    }
    return res_mat;
}
