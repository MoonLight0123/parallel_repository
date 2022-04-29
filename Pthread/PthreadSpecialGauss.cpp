#include <iostream>
#include<fstream>
#include<cstring>
#include<immintrin.h>
#include<pthread.h>
#include<semaphore.h>
#include<algorithm>
#include<windows.h>
#include<sys/time.h>
using namespace std;
struct threadParm
{
	int ID, pos;
};
const int maxN = 30000;//采用位图存储,消元子和消元行的个数之和
const int maxM = 23500;//矩阵列数
const int threadCount = 3;
int R[maxN][maxM / 32 + 1];//消元子
int E[maxN][maxM / 32 + 1];//被消元行
int Rpos[maxM];//Rpos[i]表示首项为i的消元行在R中第Rpos行
int Epos[maxN];//Epos[i]表示E第i行首项的位置
int sumR, sumE, Enum;
int workCount;//每一次消去操作中的任务总量，即需要__m128i的数量
int turn;//当前消去至被消元行的哪一行
pthread_t handle[threadCount];
threadParm tp[threadCount];
sem_t semMain, semWorkstart[threadCount],test;
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
	Enum = sumE;
	infile.close();
}
void SequentialEliminate()
{
	for (int i = 0; i < Enum; i++)
	{
		while ((E[i][0] || Epos[i] > 31) && sumE > 0)
		{
			if (Rpos[Epos[i]] != -1)//对应的消元行存在
			{
				int d = Epos[i] / 32, newEpos = -1;
				for (int j = d; j >= 0; j--)
				{
					E[i][j] ^= R[Rpos[Epos[i]]][j];
					if (newEpos == -1 && E[i][j] != 0)//更新Epos[i]
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
void* PthreadFunction(void* parm)//外层循环查找首项
{
	__m128i t0, t1;
	threadParm* p = (threadParm*)parm;
	int id = p->ID;
	while (1)
	{
		sem_wait(&semWorkstart[id]);
		if (turn >= Enum)
			return NULL;
		int workPerThread = workCount / threadCount;
		int s = workPerThread * id * 4, e;//本次消去的起始位置与结束位置
		if (id == threadCount - 1) e = workCount * 4;//最后一个线程负责边界条件
		else e = workPerThread * (id + 1) * 4;
		for (int j = e - 4; j >= s; j -= 4)
		{
			t0 = _mm_loadu_si128((__m128i*)(E[turn] + j));
			t1 = _mm_loadu_si128((__m128i*)(R[Rpos[Epos[turn]]] + j));
			t1 = _mm_xor_si128(t0, t1);
			_mm_storeu_si128((__m128i*)(E[turn] + j), t1);
		}
		int newpos = -1;
		for (int j = e - 1; j >= s; j--)
		{
			if (E[turn][j] == 0)continue;
			for (int k = 31; k >= 0; k--)
				if (getBit(E[turn][j], k))
				{
					newpos = 32 * j + k;
					break;
				}
			if (newpos != -1)break;
		}
		p->pos = newpos;
		sem_post(&semMain);

	}
}
void* threadFunction(void* parm)//内层循环查找首项
{
	__m128i t0, t1, t3;
	threadParm* p = (threadParm*)parm;
	int id = p->ID;
	while (1)
	{
		sem_wait(&semWorkstart[id]);
		if (turn >= Enum)
			return NULL;
		int workPerThread = workCount / threadCount;
		int s = workPerThread * id * 4, e;//本次消去的起始位置与结束位置
		if (id == threadCount - 1) e = workCount * 4;//最后一个线程负责边界条件
		else e = workPerThread * (id + 1) * 4;
		int firstPos = -1;
		for (int j = e - 4; j >= s; j -= 4)
		{
			t0 = _mm_loadu_si128((__m128i*)(E[turn] + j));
			t1 = _mm_loadu_si128((__m128i*)(R[Rpos[Epos[turn]]] + j));
			t1 = _mm_xor_si128(t0, t1);
			_mm_storeu_si128((__m128i*)(E[turn] + j), t1);
			if (firstPos == -1)//判断是否寻找新的E起点
			{
				t0 = _mm_xor_si128(t0, t0);
				t3 = _mm_cmpeq_epi16(t1, t0);//判断是否全部消为零
				int eq[4];
				_mm_storeu_si128((__m128i*)(eq), t3);
				for (int m = 3; m >= 0; m--)
					if (eq[m] != -1 && firstPos == -1)
						for (int k = 31; k >= 0; k--)
							if (getBit(E[turn][j + m], k))
							{
								firstPos = 32 * (j + m) + k;
								break;
							}
			}
		}
		p->pos = firstPos;
		sem_post(&semMain);
	}
}
void PthreadParallelEliminate()
{

	sem_init(&semMain, 0, 0);
	for (int id = 0; id < threadCount; id++)
		sem_init(&semWorkstart[id], 0, 0);
    sem_init(&test, 0, 0);
	for (int id = 0; id < threadCount; id++)
	{
		tp[id].ID = id;
		pthread_create(&handle[id], NULL, PthreadFunction, &tp[id]);
	}
	for (turn = 0; turn < Enum; turn++)
	{
		while ((E[turn][0] || Epos[turn] > 31) && sumE > 0)
		{
			if (Rpos[Epos[turn]] != -1)//对应的消元行存在
			{
				int d = Epos[turn] / 32;
				int newEpos = -1;
				d += (4 - d % 4);//d设置为4倍数，防止越界
				workCount = d / 4;
				for (int id = 0; id < threadCount; id++)
					sem_post(&semWorkstart[id]);
				for (int id = 0; id < threadCount; id++)
					sem_wait(&semMain);
				for (int id = 0; id < threadCount; id++)
					newEpos = max(newEpos, tp[id].pos);
				Epos[turn] = newEpos;//更新
			}
			else
			{
				Rpos[Epos[turn]] = sumR;
				memcpy(R[sumR], E[turn], sizeof(R[sumR]));//E[i]升级为消元子
				sumR++; sumE--;
				break;
			}
		}
	}
	for (int id = 0; id < threadCount; id++)
		sem_post(&semWorkstart[id]);
	for (int id = 0; id < threadCount; id++)
		pthread_join(handle[id], NULL);

	sem_destroy(&semMain);
	for (int id = 0; id < threadCount; id++)
		sem_destroy(&semWorkstart[id]);
}
void ParallelEliminate()
{
	__m128i t0, t1, t3;
	for (int i = 0; i < Enum; i++)
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
				sumR++; sumE--;
				break;
			}
		}
	}
}
void Output()
{//将结果输出到result.txt中
	string dirName = "result.txt";
	remove(dirName.c_str());
	ofstream outfile;
	outfile.open(dirName, ios::app);
	outfile << "************被消元行*************" << endl;
	for (int i = 0; i < Enum; i++)
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
int main()
{
    timeval *start=new timeval();
    timeval *stop=new timeval();
    double durationTime=0.0;


//    RInput();
//	EInput();
//    gettimeofday(start,NULL);
//    ParallelEliminate();
//    gettimeofday(stop,NULL);
//    durationTime =stop->tv_sec*1000+double(stop->tv_usec)/1000-start->tv_sec*1000-double(start->tv_usec)/1000;
//    cout << " ParallelEliminate time: " << double(durationTime) << " ms" << endl;


    RInput();
	EInput();
    gettimeofday(start,NULL);
    PthreadParallelEliminate();
    gettimeofday(stop,NULL);
    durationTime =stop->tv_sec*1000+double(stop->tv_usec)/1000-start->tv_sec*1000-double(start->tv_usec)/1000;
    cout << " PthreadParallelEliminate time: " << double(durationTime) << " ms" << endl;
	Output();
	return 0;
}
