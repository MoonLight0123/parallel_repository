#include<mpi.h>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include <iostream>
#include<immintrin.h>
#include <windows.h>
#include<algorithm>
#include<omp.h>
#include<fstream>
#include<cstring>
#pragma comment(lib,"mpi.lib")
using namespace std;
const int maxN = 30000;//采用位图存储,消元子和消元行的个数之和
const int maxM = 23500;//矩阵列数
const int numProcess = 6;
int R[maxN][maxM / 32 + 1];//消元子
int E[maxN][maxM / 32 + 1];//被消元行
int Rpos[maxM];//Rpos[i]表示首项为i的消元行在R中第Rpos行
int Epos[maxN];//Epos[i]表示E第i行首项的位置
int sumR, sumE, sumEtemp;
inline void setBit(int& a, int pos, int flag)//置位为flag
{
	if (flag)
		a |= (1 << pos);
	else
		a &= ~(1 << pos);
}
inline int getBit(int& a, int pos)
{
	return (a >> pos) & 1;
}
void RInput()
{//读入消元子
	ifstream infile("消元子.txt");
	int temp;
	bool newLine = true;
	for (int i = 0; i < maxM; i++)
		Rpos[i] = -1;
	while (infile >> temp)
	{
		if (newLine)
		{
			Rpos[temp] = sumR;
			newLine = false;
		}
		int pos1, pos2;
		pos1 = temp / 32;
		pos2 = temp % 32;
		setBit(R[sumR][pos1], pos2, 1);
		infile.get();
		if (infile.peek() == '\n')
		{
			infile.get();
			sumR++;
			newLine = true;
		}
	}
	infile.close();
}
void EInput()
{//读入被消元行
	int temp;
	bool newLine = true;
	ifstream infile("被消元行.txt");
	while (infile >> temp)
	{
		if (newLine)
		{
			Epos[sumE] = temp;
			newLine = false;
		}
		int pos1, pos2;
		pos1 = temp / 32;
		pos2 = temp % 32;
		setBit(E[sumE][pos1], pos2, 1);
		infile.get();
		if (infile.peek() == '\n')
		{
			infile.get();
			sumE++;
			newLine = true;
		}
	}
	sumEtemp = sumE;
	infile.close();
}
void Output()
{//将结果输出到result.txt中
	ofstream outfile("result.txt");
	outfile << "************被消元行*************" << endl;
	for (int i = 0; i < sumEtemp; i++)
	{
		if (!E[i][0] && Epos[i] <= 31)
			continue;//空行跳过

		for (int d = Epos[i] / 32; d >= 0; d--)
		{
			for (int j = 31; j >= 0; j--)
			{
				if (getBit(E[i][d], j))
					outfile << j + d * 32 << ' ';
			}
		}
		outfile << endl;
	}
}
int main(int argc,char* argv[])
{

	int myid;
	MPI_Status status;
	MPI_Init(NULL, NULL);
	double start, end;
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	int colSize = maxM / 32 + 1;
	EInput();
	RInput();
	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();
	if(myid==0)//主从式
		for (int i = 0; i < sumEtemp; i++)
		{
			while ((E[i][0] || Epos[i] > 31) && sumE > 0)
			{
				if (Rpos[Epos[i]] != -1)
				{
					int newEpos = -1, d = Epos[i] / 32;
					for (int id = 1; id < numProcess; i++)
					{
						MPI_Send(&i, 1, MPI_INT, id, 0, MPI_COMM_WORLD);//发送行号
						MPI_Send(&d, 1, MPI_INT, id, 1, MPI_COMM_WORLD);//发送任务量
						MPI_Send(&E[i][0], colSize, MPI_INT, id, 2, MPI_COMM_WORLD);
					}
					for (int id = 1; id < numProcess; id++)
					{
						int processMax;
						MPI_Recv(&processMax, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
						newEpos = max(processMax, newEpos);
						int s = d / (numProcess - 1)*(id - 1);//start end
						int e = (id == numProcess - 1) ? d : d / (numProcess - 1)*id;
						MPI_Recv(&E[i][s], e - s + 1, MPI_INT, id, 1, MPI_COMM_WORLD, &status);
					}
					Epos[i] = newEpos;
				}
				else
				{
					Rpos[Epos[i]] = sumR;
					memcpy(R[sumR], E[i], sizeof(R[sumR]));//E[i]升级为消元子
					sumR++; sumE--;
					break;
				}
			}
		}
	else
		while (1)
		{
			int turn, d;
			MPI_Recv(&turn, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&d, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
			MPI_Recv(&E[turn], colSize, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
			int s = d / (numProcess - 1)*(myid - 1);//start end
			int e = (myid == numProcess - 1) ? d : d / (numProcess - 1)*myid;
			for (int j = e - 1; j >= s; j++)
			{
				E[turn][j] ^= R[Rpos[Epos[turn]]][j];
			}
			int newpos = -1;
			for (int j = e - 1; j >= s; j--)
			{
				if (E[turn][j] == 0)continue;
				for(int k=31;k>=0;k--)
					if (getBit(E[turn][j], k))
					{
						newpos = 32 * j + k;
						break;
					}
				if (newpos != -1)break;
			}
			MPI_Send(&newpos, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			MPI_Send(&E[turn][s], e - s + 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
		}
	MPI_Barrier(MPI_COMM_WORLD);
	end = MPI_Wtime();
	if (myid == 0)
	{
		Output();
		cout << "MPIAlgorithm time " << (end - start) * 1000 << " ms" << endl;
	}
	MPI_Finalize();
	return 0;
}