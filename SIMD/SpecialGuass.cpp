#include <iostream>
#include<fstream>
#include<cstring>
#include<sys/time.h>
#include<immintrin.h>
using namespace std;
const int maxN = 30000;//采用位图存储,消元子和消元行的个数之和
const int maxM = 23500;//矩阵列数
int R[maxN][maxM / 32 + 1];//消元子
int E[maxN][maxM / 32 + 1];//被消元行
int Rpos[maxM];//Rpos[i]表示首项为i的消元行在R中第Rpos行
int Epos[maxN];//Epos[i]表示E第i行首项的位置
int sumR, sumE, sumEtemp;
void setBit(int& a, int pos, int flag)//置位为flag
{
	if (flag)
		a |= (1 << pos);
	else
		a &= ~(1 << pos);
}
int getBit(int& a, int pos)
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
void SequentialEliminate()
{
	for (int i = 0; i < sumEtemp; i++)
	{
		while ((E[i][0] || Epos[i] > 31) && sumE > 0)
		{
			if (Rpos[Epos[i]] != -1)//对应的消元行存在
			{
				int d = Epos[i] / 32, newEpos = -1;
				for (int j = d; j >= 0; j--)
				{
					E[i][j] ^= R[Rpos[Epos[i]]][j];
					if (newEpos==-1&&E[i][j] != 0)//更新Epos[i]
					{//未更新的条件为被消元行消为0
						for (int k = 31; k >= 0; k--)
							if (getBit(E[i][j], k))
							{
								newEpos = 32 * j + k;
								break;
							}
					}
				}
				Epos[i] = newEpos;
			}
			else
			{
				Rpos[Epos[i]] = sumR;
				memcpy(R[sumR], E[i], sizeof(R[sumR]));//E[i]升级为消元子
				//memset(E[i], 0, sizeof(E[i]));//清除E[i]
				//Epos[i] = -1;
				sumR++; sumE--;
				break;
			}
		}
	}
}
void ParallelEliminate()
{
	__m128i t0, t1, t3;
	for (int i = 0; i < sumEtemp; i++)
	{
		while ((E[i][0] || Epos[i] > 31) && sumE > 0)
		{
			if (Rpos[Epos[i]] != -1)//对应的消元行存在
			{
				int d = Epos[i] / 32, newEpos = -1;
				d += (4 - d % 4);//d设置为4倍数，防止越界
				for (int j = d - 4; j >= 0; j -= 4)
				{
					t0 = _mm_loadu_si128((__m128i*)(E[i] + j));
					t1 = _mm_loadu_si128((__m128i*)(R[Rpos[Epos[i]]] + j));
					t1 = _mm_xor_si128(t0, t1);
					_mm_storeu_si128((__m128i*)(E[i] + j), t1);
					if (newEpos == -1)
					{
						t0 = _mm_xor_si128(t0, t0);
						t3 = _mm_cmpeq_epi16(t1, t0);//判断是否全部消为零
						int eq[4];
						_mm_storeu_si128((__m128i*)(eq), t3);
						for (int m = 3; m >= 0; m--)
						{
							if (eq[m] != -1 && newEpos == -1)
							{
								for (int k = 31; k >= 0; k--)
									if (getBit(E[i][j + m], k))
									{
										newEpos = 32 * (j + m) + k;
										break;
									}
							}
						}
					}
				}
				Epos[i] = newEpos;
			}
			else
			{
				Rpos[Epos[i]] = sumR;
				memcpy(R[sumR], E[i], sizeof(R[sumR]));//E[i]升级为消元子
				//memset(E[i], 0, sizeof(E[i]));//清除E[i]
				//Epos[i] = -1;
				sumR++; sumE--;
				break;
			}
		}
	}
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

//	outfile << "************消元子*************" << endl;
//	for (int i = maxM - 1; i >= 0; i--)
//	{
//		if (Rpos[i] == -1)continue;
//		for (int d = i / 32; d >= 0; d--)
//		{
//			for (int j = 31; j >= 0; j--)
//			{
//				if (getBit(R[Rpos[i]][d], j))
//					outfile << j + d * 32 << ' ';
//			}
//		}
//		outfile << endl;
//	}
}
int main()
{

    RInput();
	EInput();

    timeval *start=new timeval();
    timeval *stop=new timeval();
    double durationTime=0.0;

//    gettimeofday(start,NULL);
//    SequentialEliminate();
//    gettimeofday(stop,NULL);
//    durationTime =stop->tv_sec*1000+double(stop->tv_usec)/1000-start->tv_sec*1000-double(start->tv_usec)/1000;
//    cout << " SequentialEliminate time: " << double(durationTime) << " ms" << endl;

//	RInput();
//	EInput();


    gettimeofday(start,NULL);
    ParallelEliminate();
    gettimeofday(stop,NULL);
    durationTime =stop->tv_sec*1000+double(stop->tv_usec)/1000-start->tv_sec*1000-double(start->tv_usec)/1000;
    cout << " ParallelEliminate time: " << double(durationTime) << " ms" << endl;
    Output();
	return 0;
}
