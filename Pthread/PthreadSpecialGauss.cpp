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
const int maxN = 30000;//����λͼ�洢,��Ԫ�Ӻ���Ԫ�еĸ���֮��
const int maxM = 23500;//��������
const int threadCount = 3;
int R[maxN][maxM / 32 + 1];//��Ԫ��
int E[maxN][maxM / 32 + 1];//����Ԫ��
int Rpos[maxM];//Rpos[i]��ʾ����Ϊi����Ԫ����R�е�Rpos��
int Epos[maxN];//Epos[i]��ʾE��i�������λ��
int sumR, sumE, Enum;
int workCount;//ÿһ����ȥ�����е���������������Ҫ__m128i������
int turn;//��ǰ��ȥ������Ԫ�е���һ��
pthread_t handle[threadCount];
threadParm tp[threadCount];
sem_t semMain, semWorkstart[threadCount],test;
void setBit(int& a, int pos, int flag)//��λΪflag
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
{//������Ԫ��
	ifstream infile("��Ԫ��.txt");
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
{//���뱻��Ԫ��
	int temp;
	bool newLine = true;
	ifstream infile("����Ԫ��.txt");
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
			if (Rpos[Epos[i]] != -1)//��Ӧ����Ԫ�д���
			{
				int d = Epos[i] / 32, newEpos = -1;
				for (int j = d; j >= 0; j--)
				{
					E[i][j] ^= R[Rpos[Epos[i]]][j];
					if (newEpos == -1 && E[i][j] != 0)//����Epos[i]
					{//δ���µ�����Ϊ����Ԫ����Ϊ0
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
				memcpy(R[sumR], E[i], sizeof(R[sumR]));//E[i]����Ϊ��Ԫ��
				//memset(E[i], 0, sizeof(E[i]));//���E[i]
				//Epos[i] = -1;
				sumR++; sumE--;
				break;
			}
		}
	}
}
void* PthreadFunction(void* parm)//���ѭ����������
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
		int s = workPerThread * id * 4, e;//������ȥ����ʼλ�������λ��
		if (id == threadCount - 1) e = workCount * 4;//���һ���̸߳���߽�����
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
void* threadFunction(void* parm)//�ڲ�ѭ����������
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
		int s = workPerThread * id * 4, e;//������ȥ����ʼλ�������λ��
		if (id == threadCount - 1) e = workCount * 4;//���һ���̸߳���߽�����
		else e = workPerThread * (id + 1) * 4;
		int firstPos = -1;
		for (int j = e - 4; j >= s; j -= 4)
		{
			t0 = _mm_loadu_si128((__m128i*)(E[turn] + j));
			t1 = _mm_loadu_si128((__m128i*)(R[Rpos[Epos[turn]]] + j));
			t1 = _mm_xor_si128(t0, t1);
			_mm_storeu_si128((__m128i*)(E[turn] + j), t1);
			if (firstPos == -1)//�ж��Ƿ�Ѱ���µ�E���
			{
				t0 = _mm_xor_si128(t0, t0);
				t3 = _mm_cmpeq_epi16(t1, t0);//�ж��Ƿ�ȫ����Ϊ��
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
			if (Rpos[Epos[turn]] != -1)//��Ӧ����Ԫ�д���
			{
				int d = Epos[turn] / 32;
				int newEpos = -1;
				d += (4 - d % 4);//d����Ϊ4��������ֹԽ��
				workCount = d / 4;
				for (int id = 0; id < threadCount; id++)
					sem_post(&semWorkstart[id]);
				for (int id = 0; id < threadCount; id++)
					sem_wait(&semMain);
				for (int id = 0; id < threadCount; id++)
					newEpos = max(newEpos, tp[id].pos);
				Epos[turn] = newEpos;//����
			}
			else
			{
				Rpos[Epos[turn]] = sumR;
				memcpy(R[sumR], E[turn], sizeof(R[sumR]));//E[i]����Ϊ��Ԫ��
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
			if (Rpos[Epos[i]] != -1)//��Ӧ����Ԫ�д���
			{
				int d = Epos[i] / 32, newEpos = -1;
				d += (4 - d % 4);//d����Ϊ4��������ֹԽ��
				for (int j = d - 4; j >= 0; j -= 4)
				{
					t0 = _mm_loadu_si128((__m128i*)(E[i] + j));
					t1 = _mm_loadu_si128((__m128i*)(R[Rpos[Epos[i]]] + j));
					t1 = _mm_xor_si128(t0, t1);
					_mm_storeu_si128((__m128i*)(E[i] + j), t1);
					if (newEpos == -1)
					{
						t0 = _mm_xor_si128(t0, t0);
						t3 = _mm_cmpeq_epi16(t1, t0);//�ж��Ƿ�ȫ����Ϊ��
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
				memcpy(R[sumR], E[i], sizeof(R[sumR]));//E[i]����Ϊ��Ԫ��
				sumR++; sumE--;
				break;
			}
		}
	}
}
void Output()
{//����������result.txt��
	string dirName = "result.txt";
	remove(dirName.c_str());
	ofstream outfile;
	outfile.open(dirName, ios::app);
	outfile << "************����Ԫ��*************" << endl;
	for (int i = 0; i < Enum; i++)
	{
		if (!E[i][0] && Epos[i] <= 31)
			continue;//��������

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
