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


void mult_block_A_B(double *A, double *B, double *res, int n_A, int m_A, int n_B, int m_B) {
    for(int k = 0; k < m_B; k++) {
        for(int i = 0; i < m_A; i++) {
            for(int j = 0; j < n_A; j++) {
                res[k * n_A + j] += A[i * n_A + j] * B[k * n_B + i];
            }
        }
    }
}


void copy_to_ram(double *A, double *destination, int i, int j, int n_A, int m_A, int row_block, int col_block, bool transpose = false) {
    // A - указатель на подматрицу возможно в большой матрице
    // i, j - положение в подматрице, которую нам дали
    // n_A, m_A - Размеры глобальной структуры матрицы A
    // row_block, col_block - размеры выделяемого блока
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


void mult_A_B(double *A, int row_A, int col_A, int n_A, int m_A, // rowA, colA - количество строк в матрице (структуре) A 
              double *B, int row_B, int col_B, int n_B, int m_B, // A - указатель на подматрицу в выделенной структуре
              double *C, int row_C, int col_C, int n_C, int m_C, // n_A, m_A - количество строк и столбцов в матрице A
              bool transpose = false, int size_block = 32) {
   
    static double block_A[32 * 32]; // put in ram
    static double block_res[32 * 32];


    if (transpose) {
        swap(n_A, m_A);
    }

    for (int i = 0 ; i < n_A - (n_A % size_block); i += size_block) {
        for (int k = 0; k < m_A - (m_A % size_block); k += size_block) {
            memset(block_A, 0, size_block * size_block * sizeof(double));
            copy_to_ram(A, block_A, i, k, row_A, col_A, size_block, size_block, transpose);
            for (int j = 0; j < m_B - (m_B % size_block); j += size_block) {
                memset(block_res, 0, size_block * size_block * sizeof(double));
                double block_B[size_block * size_block];
                copy_to_ram(B, block_B, k, j, row_B, col_B, size_block, size_block);
                mult_block_A_B(block_A, block_B, block_res, size_block, size_block, size_block, size_block);    
                 
                for(int w1 = 0; w1 < size_block; w1++) {
                    for (int w2 = 0; w2 < size_block; w2++) {
                        C[(j + w1) * row_C + (i + w2)] += block_res[w1 * size_block + w2];
                    }
                }
            }
            memset(block_res, 0, size_block * size_block * sizeof(double));
            if (m_B % size_block > 0) { // остаточные блоки с уменьшенным количеством стоблцов 
                double block_B[size_block * size_block];
                copy_to_ram(B, block_B, k, m_B - (m_B % size_block), row_B, col_B, size_block, m_B % size_block, false);
                mult_block_A_B(block_A, block_B, block_res, size_block, size_block, size_block, m_B % size_block); 
            }

            for(int w1 = 0; w1 < m_B % size_block; w1++) {
                for (int w2 = 0; w2 < size_block; w2++) {
                    C[((m_B - (m_B % size_block)) + w1) * row_C + (i + w2)] += block_res[w1 * size_block + w2];
                }
            }
        }
        if (m_A % size_block > 0) {
            memset(block_A, 0, size_block * size_block * sizeof(double));
            copy_to_ram(A, block_A, i, m_A - (m_A % size_block), row_A, col_A, size_block, m_A % size_block, transpose);
             
            for (int j = 0; j < m_B - (m_B % size_block); j += size_block) {
                memset(block_res, 0, size_block * size_block * sizeof(double));
                double block_B[size_block * size_block];
                copy_to_ram(B, block_B, m_A - (m_A % size_block), j, row_B, col_B, m_A % size_block, size_block);
                mult_block_A_B(block_A, block_B, block_res, size_block, m_A % size_block, m_A % size_block, size_block);    
                
                for(int w1 = 0; w1 < size_block; w1++) {
                    for (int w2 = 0; w2 < size_block; w2++) {
                        C[(j + w1) * row_C + (i + w2)] += block_res[w1 * size_block + w2];
                    }
                }
            }

            memset(block_res, 0, size_block * size_block * sizeof(double));
            if (m_B % size_block > 0) { // остаточные блоки с уменьшенным количеством стоблцов 
                double block_B[size_block * size_block];
                copy_to_ram(B, block_B, m_A - (m_A % size_block), m_B - (m_B % size_block), row_B, col_B, m_A % size_block, m_B % size_block, false);
                mult_block_A_B(block_A, block_B, block_res, size_block, m_A % size_block, m_A % size_block, m_B % size_block); 
            }

            for(int w1 = 0; w1 < m_B % size_block; w1++) {
                for (int w2 = 0; w2 < size_block; w2++) {
                    C[((m_B - (m_B % size_block)) + w1) * row_C + (i + w2)] += block_res[w1 * size_block + w2];
                }
            }
        }
    }

    if (n_A % size_block > 0) {
        for (int k = 0; k < m_A - (m_A % size_block); k += size_block) {
            memset(block_A, 0, size_block * size_block * sizeof(double));
            copy_to_ram(A, block_A, n_A - (n_A % size_block), k, row_A, col_A, n_A % size_block, size_block, transpose);

            for (int j = 0; j < m_B - (m_B % size_block); j += size_block) {
                memset(block_res, 0, size_block * size_block * sizeof(double));
                double block_B[size_block * size_block];
                copy_to_ram(B, block_B, k, j, row_B, col_B, size_block, size_block);
                mult_block_A_B(block_A, block_B, block_res, n_A % size_block, size_block, size_block, size_block);    
                
                for(int w1 = 0; w1 < size_block; w1++) {
                    for (int w2 = 0; w2 < n_A % size_block; w2++) {
                        C[(j + w1) * row_C + (n_A - (n_A % size_block) + w2)] += block_res[w1 * (n_A % size_block) + w2];
                    }
                }
            }
            memset(block_res, 0, size_block * size_block * sizeof(double));
            if (m_B % size_block > 0) { // остаточные блоки с уменшенным количеством стоблцов 
                double block_B[size_block * size_block];
                copy_to_ram(B, block_B, k, m_B - (m_B % size_block), row_B, col_B, size_block, m_B % size_block, false);
                mult_block_A_B(block_A, block_B, block_res, n_A % size_block, size_block, size_block, m_B % size_block); 
            }

            for(int w1 = 0; w1 < m_B % size_block; w1++) {
                for (int w2 = 0; w2 < n_A % size_block; w2++) {
                    C[((m_B - (m_B % size_block)) + w1) * row_C + (n_A - (n_A % size_block) + w2)] += block_res[w1 * (n_A % size_block) + w2];
                }
            }
        }
        if (m_A % size_block > 0) {
            memset(block_A, 0, size_block * size_block * sizeof(double));
            copy_to_ram(A, block_A, n_A - (n_A % size_block), m_A - (m_A % size_block), row_A, col_A, n_A % size_block, m_A % size_block, transpose);
            for (int j = 0; j < m_B - (m_B % size_block); j += size_block) {
                memset(block_res, 0, size_block * size_block * sizeof(double));
                double block_B[size_block * size_block];
                copy_to_ram(B, block_B, m_A - (m_A % size_block), j, row_B, col_B, m_A % size_block, size_block);
                mult_block_A_B(block_A, block_B, block_res, n_A % size_block, m_A % size_block, m_A % size_block, size_block);    
                
                for(int w1 = 0; w1 < size_block; w1++) {
                    for (int w2 = 0; w2 < n_A % size_block; w2++) {
                        C[(j + w1) * row_C + (n_A - (n_A % size_block) + w2)] += block_res[w1 * (n_A % size_block) + w2];
                    }
                }
            }
            memset(block_res, 0, size_block * size_block * sizeof(double));
            if (m_B % size_block > 0) { // остаточные блоки с уменьшенным количеством стоблцов 
                double block_B[size_block * size_block];
                copy_to_ram(B, block_B, m_A - (m_A % size_block), m_B - (m_B % size_block), row_B, col_B, m_A % size_block, m_B % size_block, false);
                mult_block_A_B(block_A, block_B, block_res, n_A % size_block, m_A % size_block, m_A % size_block, m_B % size_block); 
            }

            for(int w1 = 0; w1 < m_B % size_block; w1++) {
                for (int w2 = 0; w2 < n_A % size_block; w2++) {
                    C[((m_B - (m_B % size_block)) + w1) * row_C + (n_A - (n_A % size_block) + w2)] += block_res[w1 * (n_A % size_block) + w2];
                }
            }
        }
    }
}

void mult_A_b(double *A, double *b, double *res, int n_A, int m_A, int n_b, bool transpose = false, int size_block = 4096) {
    static double tmp_ram[4096]; /// уточнить этот момент
    memset(tmp_ram, 0, 4096 * sizeof(double));

    if (transpose) {
        for (int i = 0; i < n_A; i++) {
            tmp_ram[i] = b[i];
        }
        for (int i = 0; i < m_A; i++) {
            double sum = 0;
            for (int j = 0 ; j < n_A; j++) {
                sum += A[i * n_A + j] * tmp_ram[j];
            }
            res[i] = sum;
        }
    } else {
        for (int i = 0 ; i < m_A; i++) {
            double tmp = b[i];
            for (int j = 0; j < n_A; j++) {
                tmp_ram[j] += A[i * n_A + j] * tmp;
            }
        }
        for (int i = 0 ; i < n_A; i++) {
            res[i] = tmp_ram[i];
        }
    }
}


void iteration(double *a, double *Q, double *R, int n, int offset, int size_block = 32) {
    mult_A_b(Q, a, R + (offset - 1) * n, n, offset - 1, n, true, size_block);
    
    double tmp[n];
    mult_A_b(Q, R + (offset - 1) * n, tmp, n, offset - 1, offset - 1, false, size_block);
    for(int i = 0; i < n; i++) {
        tmp[i] *= -1;
        tmp[i] += a[i];
    }
    double norm_q = norm(tmp, n);
    for(int i = 0 ; i < n; i++) {
        Q[(offset - 1) * n + i] = tmp[i] / norm_q;
    }
    R[(offset - 1) * n + (offset - 1)] = norm_q;
}


void intraortho(double *Q, double *R, int n, int s, int size_block = 32) {
    double norm_q_1 = norm(Q, n);
    R[0] = norm_q_1;
    for(int i = 0; i < n; i++) {
        Q[i] /= norm_q_1;
    }
    for (int i = 2; i <= s; i++) {
        iteration(Q + (i - 1) * n, Q, R, n, i, size_block); 
    }
}


void const_mul(double *X, int count, double alpha) {
    for (int i = 0 ; i < count; i++) {
        X[i] *= alpha;
    }
}


void QR(double *Q, double *R, int n, int s, int size_block = 32) { 
    intraortho(Q, R, n, s, size_block); 
    for (int k = 1; k <= n / s - 1; k++) {
        mult_A_B(Q, n, n, n, k * s, 
                Q + n * s * k, n, n, n, s,
                R + n * s * k, n, n, k * s, s,
                true, size_block);
        const_mul(Q + n * s * k, n * s, -1);  
        mult_A_B(Q, n, n, n, k * s,
                R + n * s * k, n, n, k * s, s,
                Q + n * s * k, n, n, n, s,
                false, size_block);
        const_mul(Q + n * s * k, n * s, -1); 
        intraortho(Q + n * s * k, R + n * s * k + s * k, n, s, size_block);
    }
}


int main() {

       
    int n, block_size = 32;
    int s = 128;
    cout << "Enter size of matrix (<= 4096): " << endl;
    cin >> n;
    
    double *A = (double *) malloc(n * n * sizeof(double));
    random_init(A, n);
    double *Q = (double *) malloc(n * n * sizeof(double));
    memcpy(Q, A, n * n * sizeof(double));
    double *R = (double *) malloc(n * n * sizeof(double)); 
    
    
    auto begin = std::chrono::steady_clock::now();
    QR(Q, R, n, s, block_size);
    auto end = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "The time: " << elapsed_ms.count() << " ms\n";

    /*
    double *C = (double *) malloc(n * n * sizeof(double));
    matmul_checker(Q, R, C, n);
    cout << "||A - QR|| = " << Frob(C, A, n, n) << endl; 
    cout << "Unitary: " << check_unitary(Q, n) << endl; 
    */
    return 0;
}
