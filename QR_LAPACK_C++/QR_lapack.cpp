#include <iostream>
#include <cmath>
#include <cstring>
#include <chrono>
#include <vector> 


using namespace std;

void random_init(double *A, int n) {
    for(int i = 0; i < n * n; i++) {
        A[i] = int(rand()) % 100;
    }
}  


extern "C" void dgeqrf_(int* M, int* N, double* A, int* LDA, double* TAU, double* WORK, int* LDWORK, int* INFO);


int main(){

    vector <int> sizes = {128, 256, 512, 768, 1024, 1280, 1536, 1792, 2048, 2560,
                          2816, 3072, 3414, 3756, 4096};

    for (int i = 0; i < sizes.size(); i++) {

        int n = sizes[i];
        double *A = (double *) malloc(n * n * sizeof(double));
        random_init(A, n);

        int lda = n, ldwork = -1, info = 0;
        double* tau = (double *) malloc(n * sizeof(double));
        double* sample = (double *) malloc(n * sizeof(double));

        dgeqrf_(&n, &n, A, &lda, tau, sample, &ldwork, &info);
        ldwork = sample[0];

        double* work = (double *) malloc(ldwork * sizeof(double));

    	auto begin = std::chrono::steady_clock::now();
    	dgeqrf_(&n, &n, A, &lda, tau, work, &ldwork, &info);
    	auto end = std::chrono::steady_clock::now();
    	auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    	std::cout << "N=" << n <<  "  The time: " << elapsed_ms.count() << " ms\n";
    }
    return 0;
}
