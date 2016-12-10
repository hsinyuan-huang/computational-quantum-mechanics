#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <utility>
#include <map>
#include <queue>
#include <set>
#include <vector>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

using namespace std;
using namespace Eigen;

const double pi = M_PI;

double potential2D(double x, double y){
	double r = sqrt(x * x + y * y);

	//return 0.5 * r * r; // SHO
	return -600 * 0.3 / (abs(x) + 0.3) - 600 * 0.3 / (abs(y) + 0.3)
		   -600 * 0.3 / (abs(r - 0.5) + 0.3); //Weird
}

int main(int argc, char* argv[]){
	int topK = 500; //number of excited states output

	double hbar = 1; // unit reduced Plank constant
	double m = 1; // unit mass

	int n = 40; // base grid
	int N = 2 * n + 1; // Number of grid points

	double dx = 0.05;
	double dk = 2 * pi / N / dx;

	double *T_pre = (double*)malloc(sizeof(double) * (N + 2));

	for(int di = 0; di <= N; di++){
		T_pre[di] = 0;
		for(int q = 0; q <= n; q++){
			double Tv = (hbar * q * dk) * (hbar * q * dk) / (2.0 * m);
			T_pre[di] += 2.0 * cos(2.0 * pi * q * di / N) / N * Tv;
		}
	}

	MatrixXd H(N*N, N*N); // For 2D

	printf("Start constructing Hamiltonian\n");
	for(int i = -n; i <= n; i++) for(int k = -n; k <= n; k++){ // xi = dx * i, yk = dx * k
		for(int j = -n; j <= n; j++) for(int l = -n; l <= n; l++){ // xj = dx * j, yl = dx * l
			int di = i < j? j - i: i - j;
			int dl = k < l? l - k: k - l;

			double &matrix_elem = H((i + n) * N + (k + n), (j + n) * N + (l + n));
			matrix_elem = 0;

			if(dl == 0)
				matrix_elem += T_pre[di]; // Kinetic Energy (X-dir)

			if(di == 0)
				matrix_elem += T_pre[dl]; // Kinetic Energy (Y-dir)

			if(i == j && k == l) // Kronecker Delta
				H((i + n) * N + (k + n), (j + n) * N + (l + n))
					+= potential2D(dx * i, dx * k); // Potential Energy
		}
	}

	printf("Eigenvalue decomposition ...\n");
	SelfAdjointEigenSolver<MatrixXd> es(H);
	VectorXd e_val = es.eigenvalues();

	printf("Top energy values:\n");
	for(int i = 0; i < topK; i++){
		printf("%f\n", e_val(i));
	}

	FILE* out = fopen("top_eigenstates.txt", "w");
	for(int i = 0; i < topK; i++){
		VectorXd e_vec = es.eigenvectors().col(i);
		for(int a = -n; a <= n; a++){
			for(int b = -n; b <= n; b++){
				fprintf(out, "%f ", e_vec((a + n) * N + (b + n)));
			}
			fprintf(out, "\n");
		}
	}

	return 0;
}
