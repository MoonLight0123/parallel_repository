#include <iostream>
#include<fstream>
#include<cstring>
#include<sys/time.h>
#include<immintrin.h>
using namespace std;
const int maxN = 30000;//����λͼ�洢,��Ԫ�Ӻ���Ԫ�еĸ���֮��
const int maxM = 23500;//��������
int R[maxN][maxM / 32 + 1];//��Ԫ��
int E[maxN][maxM / 32 + 1];//����Ԫ��
int Rpos[maxM];//Rpos[i]��ʾ����Ϊi����Ԫ����R�е�Rpos��
int Epos[maxN];//Epos[i]��ʾE��i�������λ��
int sumR, sumE, sumEtemp;
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
	sumEtemp = sumE;
	infile.close();
}
void SequentialEliminate()
{
	for (int i = 0; i < sumEtemp; i++)
	{
		while ((E[i][0] || Epos[i] > 31) && sumE > 0)
		{
			if (Rpos[Epos[i]] != -1)//��Ӧ����Ԫ�д���
			{
				int d = Epos[i] / 32, newEpos = -1;
				for (int j = d; j >= 0; j--)
				{
					E[i][j] ^= R[Rpos[Epos[i]]][j];
					if (newEpos==-1&&E[i][j] != 0)//����Epos[i]
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
void ParallelEliminate()
{
	__m128i t0, t1, t3;
	for (int i = 0; i < sumEtemp; i++)
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
				//memset(E[i], 0, sizeof(E[i]));//���E[i]
				//Epos[i] = -1;
				sumR++; sumE--;
				break;
			}
		}
	}
}
void Output()
{//����������result.txt��
	ofstream outfile("result.txt");
	outfile << "************����Ԫ��*************" << endl;
	for (int i = 0; i < sumEtemp; i++)
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

//	outfile << "************��Ԫ��*************" << endl;
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
