#include<iostream>
#include<fstream>
#include<ostream>
#include<pthread.h>
#include<semaphore.h>
#include<algorithm>
#include<math.h>
#pragma comment(lib, "pthreadVC2.lib")
using namespace std;
struct threadParam
{
	int threadID, threadPos;
};
const int maxN = 10000;//采用位图存储,消元子和消元行的个数之和
const int maxM = 10000;//矩阵列数
int R[maxN][maxM / 32 + 1];//消元子
int E[maxN][maxM / 32 + 1];//被消元行
int rPos[maxM];//rPos[i]表示首项为i的消元行在R中第rPos行
int ePos[maxN];//ePos[i]表示E第i行首项的位置
int rNumNow, eNumNow, eNumOrigin;
const int threadCount = 1;
pthread_t handle[threadCount];
sem_t semMain, semWorkstart[threadCount];
threadParam tp[threadCount];
int turn;//所有线程当前消去的轮次
int s, e;//消去过程的起点终点
int finish;//记录完成情况
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
void input()
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
}
void output()
{
	string dirName = "result.txt";
	remove(dirName.c_str());
	ofstream outfile;
	outfile.open(dirName, ios::app);
	outfile << "1************被消元行*************2" << endl;
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
void* pthreadFunction(void* param)
{
	threadParam* p = (threadParam*)param;
	int id = p->threadID;
	while (1)
	{
		sem_wait(&semWorkstart[id]);
		int threadTurn = id + turn;
		if (threadTurn >= eNumOrigin)
		{
			setBit(finish, id, 0);
			return NULL;
		}
		while (ePos[threadTurn] >= e && ePos[threadTurn] <= s && rPos[ePos[threadTurn]] != -1)
		{
			int d = ePos[threadTurn] / 32, newEpos = -1;
			for (int j = d; j >= 0; j--)
			{
				E[threadTurn][j] ^= R[rPos[ePos[threadTurn]]][j];
				if (newEpos == -1 && E[threadTurn][j] != 0)//更新Epos[i]
				{//未更新的条件为被消元行消为0
					for (int k = 31; k >= 0; k--)
						if (getBit(E[threadTurn][j], k))
						{
							newEpos = 32 * j + k;
							break;
						}
				}
			}
			ePos[threadTurn] = newEpos;
			if (newEpos = -1)break;
		}
		sem_post(&semMain);
		//p->threadPos = newEpos;
	}
}
void pthreadEliminateImproved()
{
	sem_init(&semMain, 0, 0);
	for (int id = 0; id < threadCount; id++)
		sem_init(&semWorkstart[id], 0, 0);
	for (int id = 0; id < threadCount; id++)
	{
		tp[id].threadID = id;
		pthread_create(&handle[id], NULL, pthreadFunction, &tp[id]);
	}
	for (turn = 0; turn < eNumOrigin; turn += threadCount)
	{
		finish = pow(2, threadCount) - 1;//采用位来记录是否完成 finish=32'b0000...0000是完成本轮消去,循环结束
		while (finish)
		{
			e = -1, s = -1;//本次消元操作的起点与终点
			for (int id = 0; id < threadCount; id++)
				if (getBit(finish, id))
					s = max(ePos[turn + id], s);
			if (s == -1)break;
			e = s;
			while (rPos[--e] != -1); e += 1;
			//s到e为本次消去的区间
			for (int id = 0; id < threadCount; id++)
				if (getBit(finish, id))
					sem_post(&semWorkstart[id]);
			for (int id = 0; id < threadCount; id++)
				if (getBit(finish, id))
					sem_wait(&semMain);
			for (int id = 0; id < threadCount; id++)//处理有数据依赖的部分
			{
				int threadTurn = id + turn;
				if (getBit(finish, id) && rPos[ePos[threadTurn]] == -1)
				{
					setBit(finish, id, 0);
					rPos[ePos[threadTurn]] = rNumNow;
					memcpy(R[rNumNow], E[threadTurn], sizeof(R[rNumNow]));//E[i]升级为消元子
					rNumNow++; eNumNow--;
				}
			}
			//if (newEpos == -1)break;//该行消为空行
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
int main()
{
	input();
	//pthreadEliminateImproved();
	eliminate();
	output();
	return 0;
}