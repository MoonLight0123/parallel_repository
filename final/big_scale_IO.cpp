#include<iostream>
#include<fstream>
#include<ostream>
#include<pthread.h>
#include<semaphore.h>
#include<algorithm>
#include<math.h>
#include <time.h>
#include <windows.h>
#include<immintrin.h>
#pragma comment(lib, "pthreadVC2.lib")
using namespace std;
const int maxN = 25000;//采用位图存储,消元子和消元行的个数之和
const int maxM = 1000000;//矩阵列数
const int inputERow = 1717;
int eNum = 10000;
const int rNum = 15000;
int R[maxN][maxM / 32 + 1];//消元子
int E[maxN][maxM / 32 + 1];//被消元行
int rPos[maxM];//rPos[i]表示首项为i的消元行在R中第rPos行
int ePos[maxN];//ePos[i]表示E第i行首项的位置
int rNumNow, eNumNow, eNumOrigin;
int gettimeofday(struct timeval *tp, void *tzp)
{
	time_t clock;
	struct tm tm;
	SYSTEMTIME wtm;
	GetLocalTime(&wtm);
	tm.tm_year = wtm.wYear - 1900;
	tm.tm_mon = wtm.wMonth - 1;
	tm.tm_mday = wtm.wDay;
	tm.tm_hour = wtm.wHour;
	tm.tm_min = wtm.wMinute;
	tm.tm_sec = wtm.wSecond;
	tm.tm_isdst = -1;
	clock = mktime(&tm);
	tp->tv_sec = clock;
	tp->tv_usec = wtm.wMilliseconds * 1000;
	return 0;
}
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
void inputR()
{
	ifstream infile("消元子.txt");
	int temp;
	bool newLine = true;
	for (int i = 0; i < maxM; i++)
		rPos[i] = -1;
	while (infile >> temp)
	{
		if (newLine)
		{
			rPos[temp] = rNumNow;
			newLine = false;
		}
		int pos1, pos2;
		pos1 = temp / 32;
		pos2 = temp % 32;
		setBit(R[rNumNow][pos1], pos2, 1);
		infile.get();
		if (infile.peek() == '\n')
		{
			infile.get();
			rNumNow++;
			newLine = true;
		}
	}
	infile.close();

	newLine = true;
	ifstream infil("被消元行.txt");
	while (infil >> temp)
	{
		if (newLine)
		{
			ePos[eNumNow] = temp;
			newLine = false;
		}
		int pos1, pos2;
		pos1 = temp / 32;
		pos2 = temp % 32;
		setBit(E[eNumNow][pos1], pos2, 1);
		infil.get();
		if (infil.peek() == '\n')
		{
			infil.get();
			eNumNow++;
			newLine = true;
		}
	}
	eNumOrigin = eNumNow;
	infil.close();
	//rPos[-1] = -1;
}
void outputE()
{
	string dirName = "result.txt";
	remove(dirName.c_str());
	ofstream outfile;
	outfile.open(dirName, ios::app);
	for (int i = 0; i < eNumOrigin; i++)
	{
		if (!E[i][0] && ePos[i] <= 31)
			continue;//空行跳过

		for (int d = ePos[i] / 32; d >= 0; d--)
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
void eliminate()
{
	for (int i = 0; i < eNumOrigin; i++)
	{
		while (rPos[ePos[i]] != -1)
		{
			int d = ePos[i] / 32, newEpos = -1;
			for (int j = d; j >= 0; j--)
			{
				E[i][j] ^= R[rPos[ePos[i]]][j];
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
			ePos[i] = newEpos;
			if (newEpos == -1)break;//该行消为空行
		}
		if (ePos[i] != -1)
		{
			rPos[ePos[i]] = rNumNow;//当空行时，该语句为rPos[-1]=rNumNow，需要跳过
			memcpy(R[rNumNow], E[i], sizeof(R[rNumNow]));//E[i]升级为消元子
			rNumNow++; eNumNow--;
		}
	}
}
void inputE()
{
	int temp;
	bool newLine = true;
	ifstream infil("被消元行.txt");
	for(int i=0;i<=inputERow&&eNum;++i,--eNum)
	{
		if (newLine)
		{
			ePos[eNumNow] = temp;
			newLine = false;
		}
		int pos1, pos2;
		pos1 = temp / 32;
		pos2 = temp % 32;
		setBit(E[eNumNow][pos1], pos2, 1);
		infil.get();
		if (infil.peek() == '\n')
		{
			infil.get();
			eNumNow++;
			newLine = true;
		}
	}
	eNumOrigin = eNumNow;
	infil.close();
}
void bigScaleIO()
{
	inputR();
	for (int k = 0; k < eNum; k += inputERow)
	{
		inputE();
		for (int i = k; i < inputERow&&i < eNum; i++)
		{
			while (rPos[ePos[i]] != -1)
			{
				int d = ePos[i] / 32, newEpos = -1;
				for (int j = d; j >= 0; j--)
				{
					E[i][j] ^= R[rPos[ePos[i]]][j];
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
				ePos[i] = newEpos;
				if (newEpos == -1)break;//该行消为空行
			}
			if (ePos[i] != -1)
			{
				rPos[ePos[i]] = rNumNow;//当空行时，该语句为rPos[-1]=rNumNow，需要跳过
				memcpy(R[rNumNow], E[i], sizeof(R[rNumNow]));//E[i]升级为消元子
				rNumNow++; eNumNow--;
			}
		}
		outputE();
	}
}
int main()
{
	timeval *tart = new timeval();
	timeval *top = new timeval();
	double urationTime = 0.0;

	timeval *start = new timeval();
	timeval *stop = new timeval();
	double durationTime = 0.0;

	gettimeofday(start, NULL);
	gettimeofday(tart, NULL);
	//eliminate();
	bigScaleIO();
	gettimeofday(top, NULL);
	gettimeofday(stop, NULL);
	//gettimeofday(stop, NULL);

	urationTime = top->tv_sec * 1000 + double(top->tv_usec) / 1000 - tart->tv_sec * 1000 - double(tart->tv_usec) / 1000;

	durationTime = stop->tv_sec * 1000 + double(stop->tv_usec) / 1000 - start->tv_sec * 1000 - double(start->tv_usec) / 1000;
	cout << " all time: " << double(durationTime) << " ms" << endl;
	cout << " SSEalgorithm time: " << double(urationTime) << " ms" << endl;
	cout << " io time: " << double(durationTime - urationTime) << " ms" << endl;

	return 0;
}