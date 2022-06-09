#include<mpi.h>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include <iostream>
#include<immintrin.h>
#include <windows.h>
#include<algorithm>
#include<omp.h>
#pragma comment(lib,"mpi.lib")
using namespace std;
const int N = 2000, numProcess = 9;
float a[N][N];
void init()
{
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			a[i][j] = float(rand()) / 10;
}
void SequentialAlgorithm()
{
	for (int k = 0; k < N; k++)
	{
		for (int j = k + 1; j < N; j++)
			a[k][j] /= a[k][k];
		a[k][k] = 1.0;
		for (int i = k + 1; i < N; i++) {

			for (int j = k + 1; j < N; j++)
				a[i][j] -= a[i][k] * a[k][j];
			a[i][k] = 0;
		}
	}
}

void show()
{

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++)
			cout << a[i][j] << ' ';
		cout << endl;
	}
}
int main(int argc, char *argv[])
{
	int myid;
	MPI_Status status;
	MPI_Init(NULL, NULL);
	double start, end;
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	if (myid == 0)
	{
		init();
		for (int i = 1; i < numProcess; ++i)
			for (int j = i; j < N; j += numProcess)
				MPI_Send(&a[j][0], N, MPI_FLOAT, i, j, MPI_COMM_WORLD);
	}
	else
	{
		for (int j = myid; j < N; j += numProcess)
			MPI_Recv(&a[j][0], N, MPI_FLOAT, 0, j, MPI_COMM_WORLD, &status);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();
	int i, j, k, r2;
#pragma omp parallel num_threads(threadCount),private(i,j,k,r2)
	for (k = 0; k < N; ++k)
	{
#pragma omp single
		{
			if (myid == 0)
			{
				for (j = k + 1; j < N; ++j)
					a[k][j] /= a[k][k];
				a[k][k] = 1.0;
				for (j = 1; j < numProcess; ++j)
				{
					MPI_Send(&a[k][0], N, MPI_FLOAT, j, j, MPI_COMM_WORLD);
				}
			}
			else
				MPI_Recv(&a[k][0], N, MPI_FLOAT, 0, myid, MPI_COMM_WORLD, &status);
			r2 = myid;
			while (r2 < k + 1)r2 += numProcess;
		}
#pragma omp for
		for (i = r2; i < N; i += numProcess)
		{
			for (j = k + 1; j < N; ++j)
				a[i][j] -= a[k][j] * a[i][k];
			a[i][k] = 0.0;
		}
#pragma omp single
		{
			if ((k + 1) % numProcess == myid && myid != 0)
				MPI_Send(&a[k + 1][0], N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
			if (myid == 0 && (k + 1) % numProcess != 0 && k + 1 < N)
				MPI_Recv(&a[k + 1][0], N, MPI_FLOAT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		}
		
	}
	end = MPI_Wtime();
	if (myid == 0) {
		show();
		cout << "MPIAlgorithm time " << (end - start) * 1000 << " ms" << endl;
	}
	MPI_Finalize();
	return 0;
}