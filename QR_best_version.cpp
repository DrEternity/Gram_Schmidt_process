#include <iostream>
#include <cmath>
#include <cstring>
#include <chrono>

using namespace std;


void print_matrix(double *A, int n) {
    cout << "[";
    for(int i = 0 ; i < n; i++) {
        cout << "[";
        for (int j = 0 ; j < n; j++) {
            cout << A[j * n + i] << ", "; 
        }
        cout << "],";
        cout << endl;
    }
    cout << "]";
    cout << endl;
}


void random_init(double *A, int n) {
    for(int i = 0; i < n * n; i++) {
        A[i] = int(rand()) % 100;
    }
}  


double Frob(double *A, double *B, int n, int m) {
    double error = 0;
    for (int i = 0 ; i < n; i++) {
        for (int j = 0 ; j < m; j++) {
            error += (A[j * n + i] - B[j * n + i]) * (A[j * n + i] - B[j * n + i]);     
        }
    }
    return pow(error, 0.5);
}


double norm(double *a, int n) {
    double ans = 0;
    for (int i = 0; i < n; i++) {
        ans += a[i] * a[i];
    }
    return pow(ans, 0.5);
}


void matmul_checker(double *A, double *B, double *C, int n) { // C = AB
    memset(C, 0, n * n * sizeof(double));
    for (int j = 0; j < n; j++) {
        for (int k = 0; k < n; k++) {
            for(int i = 0 ; i < n; i++) {
                C[j * n + i] += A[k * n + i] * B[j * n + k];
            }
        }
    }
}


double check_unitary(double *A, int n) {
    double C[n * n];
    memset(C, 0, n * n * sizeof(double));
    for (int j = 0; j < n; j++) {
        for (int k = 0; k < n; k++) {
            for(int i = 0 ; i < n; i++) {
                C[j * n + i] += A[k * n + i] * A[k * n + j];
            }
        }
        C[j * n + j] -= 1;
    }
    double error = norm(C, n * n) / pow(n, 0.5);
    return error;
}


void mult_block_A_b(double *A, double *b, double *res, int n_A, int m_A) {
    for(int i = 0 ; i < m_A; i++) {
        for(int j = 0 ; j < n_A; j++) {
            res[j] += A[i * n_A + j] * b[i];
        }
    }
}


void copy_to_ram(double *A, double *destination, int i, int j, int n_A, int m_A, int row_block, int col_block, bool transpose = false) {
    if (transpose == false) {
        for(int w = 0 ; w < col_block; w++) {
            for (int v = 0; v < row_block; v++) {
                destination[w * row_block + v] = A[(j + w) * n_A + (i + v)];
            }
        }
    } else {
        for (int v = 0; v < row_block; v++) {
            for(int w = 0; w < col_block; w++) {
                destination[w * row_block + v] = A[(i + v) * n_A + (j + w)];
            }
        }
    }
}


void mult_A_b(double *A, double *b, double *res, int n_A, int m_A, int n_b, bool transpose = false, int size_block = 32) {
    static double block_A[32 * 32]; // put in ram
    static double block_b[32];
    static double block_res[32];

    int row_A = n_A, col_A = m_A;

    if (transpose) {
        swap(n_A, m_A);
    }

    for (int i = 0 ; i < n_A - (n_A % size_block); i += size_block) {
        memset(block_res, 0, size_block * sizeof(double));

        for (int j = 0; j < m_A - (m_A % size_block); j += size_block) {
            copy_to_ram(A, block_A, i, j, row_A, col_A, size_block, size_block, transpose);
            copy_to_ram(b, block_b, j, 0, m_A, 1, size_block, 1, false);
            mult_block_A_b(block_A, block_b, block_res, size_block, size_block);    
        }

        if (m_A % size_block > 0) { // остаточные блоки с уменшенным количеством стоблцов
            copy_to_ram(A, block_A, i, m_A - (m_A % size_block), row_A, col_A, size_block, m_A % size_block, transpose);
            copy_to_ram(b, block_b, m_A - (m_A % size_block), 0, m_A, 1, m_A % size_block, 1, false);
            mult_block_A_b(block_A, block_b, block_res, size_block, m_A % size_block); 
        }

        for (int j = 0; j < size_block; j++) {
            res[i + j] = block_res[j];
        }
    }

    if (n_A % size_block > 0) { // остаточные блоки с уменьшенным количеством строк
        memset(block_res, 0, size_block * sizeof(double));
        for (int j = 0; j < m_A - (m_A % size_block); j += size_block) {
            copy_to_ram(A, block_A, n_A - (n_A % size_block), j, row_A, col_A, n_A % size_block, size_block, transpose);
            copy_to_ram(b, block_b, j, 0, m_A, 1, size_block, 1, false);
            mult_block_A_b(block_A, block_b, block_res, n_A % size_block, size_block);      
        }   
        
        if (m_A % size_block > 0) { // оставшийся маленький блок в углу
            copy_to_ram(A, block_A, n_A - (n_A % size_block), m_A - (m_A % size_block), row_A, col_A, n_A % size_block, m_A % size_block, transpose);
            copy_to_ram(b, block_b, m_A - (m_A % size_block), 0, m_A, 1, m_A % size_block, 1, false);
            mult_block_A_b(block_A, block_b, block_res, n_A % size_block, m_A % size_block);  
        }
        for (int j = 0; j < n_A % size_block; j++) {
            res[n_A - (n_A % size_block) + j] = block_res[j];
        }
    }
}


void iteration(double *a, double *Q, double *R, int n, int num_iter, int size_block = 32) {
    mult_A_b(Q, a, R + (num_iter - 1) * n, n, num_iter - 1, n, true, size_block);
    
    double tmp[n];
    mult_A_b(Q, R + (num_iter - 1) * n, tmp, n, num_iter - 1, num_iter - 1, false, size_block);
    for(int i = 0; i < n; i++) {
        tmp[i] *= -1;
        tmp[i] += a[i];
    }
    double norm_q = norm(tmp, n);
    for(int i = 0 ; i < n; i++) {
        Q[(num_iter - 1) * n + i] = tmp[i] / norm_q;
    }
    R[(num_iter - 1) * n + (num_iter - 1)] = norm_q;
}


void QR(double *A, double *Q, double *R, int n, int size_block) {
    double norm_q_1 = norm(A, n);
    R[0] = norm_q_1;
    for(int i = 0; i < n; i++) {
        Q[i] = A[i] / norm_q_1;
    }
    
    for(int i = 2 ; i <= n; i++) {
        iteration(A + (i - 1) * n , Q, R, n, i, size_block);
        // Logs:
        if (i % 50 == 0) {
            cout << "Step " << i << endl;
        }
    }
}


int main() {
    int n, block_size = 32;
    cout << "Enter size of matrix: " << endl;
    cin >> n;
    double *A = (double *) malloc(n * n * sizeof(double));
    double *C = (double *) malloc(n * n * sizeof(double));
    random_init(A, n);
    double *Q = (double *) malloc(n * n * sizeof(double));
    double *R = (double *) malloc(n * n * sizeof(double));
    

    auto begin = std::chrono::steady_clock::now();
    QR(A, Q, R, n, block_size);
    auto end = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "The time: " << elapsed_ms.count() << " ms\n";


    /*
    matmul_checker(Q, R, C, n);
    cout << "||A - QR|| = " << Frob(C, A, n, n) << endl; 
    cout << "Unitary: " << check_unitary(Q, n) << endl; 
    */

    return 0;
}
