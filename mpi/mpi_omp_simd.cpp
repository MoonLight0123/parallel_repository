#include<mpi.h>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<immintrin.h>
#include<algorithm>
#include<omp.h>
#pragma comment(lib,"mpi.lib")
using namespace std;
const int N = 2000, numProcess = 6, threadCount = 6;
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
			printf("%f ", a[i][j]);
		printf("\n");
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
	__m128 t0, t1, t2, t3;
	int i, j, k, r2;
	float temp1[4], temp2[4];
#pragma omp parallel num_threads(threadCount),private(i,j,k,r2,t0,t1,t2,t3,temp1,temp2)
	for (k = 0; k < N; ++k)
	{
#pragma omp single
		{
			if (myid == 0)
			{
				temp1[0] = temp1[1] = temp1[2] = temp1[3] = a[k][k];
				t0 = _mm_loadu_ps(temp1);
				for (j = k + 1; j+3 < N; j+=4)
				{
					t1 = _mm_loadu_ps(a[k] + j);
					t2 = _mm_div_ps(t1, t0);
					_mm_storeu_ps(a[k] + j, t2);
				}
				for (; j < N; j++)
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
			temp2[0] = temp2[1] = temp2[2] = temp2[3] = a[k][k];
			t0 = _mm_loadu_ps(temp2);
			for (j = k + 1; j + 3 < N; j += 4)
			{
				t1 = _mm_loadu_ps(a[k] + j);
				t2 = _mm_loadu_ps(a[i] + j);
				t3 = _mm_mul_ps(t0, t1);
				t2 = _mm_sub_ps(t2, t3);
				_mm_storeu_ps(a[i] + j, t2);
			}
			for (; j < N; j++)
				a[i][j] -= a[i][k] * a[k][j];
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
		//show();
		printf("MPIAlgorithm time  %f  ms", (end - start) * 1000);
	}
	MPI_Finalize();
	return 0;
}
