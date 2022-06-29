#include<iostream>
#include<fstream>
#include<ostream>
#include<cstring>
using namespace std;
const int maxN = 1000;//最大消元子被消元子个数
const int maxM = 1000;//向量的最大长度
int E[maxN][maxM], R[maxN][maxM];
int eLen[maxN], rLen[maxN];//E[i]代表E第i行的元素数量,例：E[i]这一行仅E[i][0]有元素，eLen[i]=1;
int rPos[maxN];//rPos[i]表示首项为i的消元行在R中第rPos行
int rNumNow, eNumNow, eNumOrigin;
void input()
{
	ifstream infile("消元子.txt");
	int temp;
	bool newLine = true;
	for (int i = 0; i < maxM; i++)
		rPos[i] = -1;
	int count;
	while (infile >> temp)
	{
		if (newLine)
		{
			count = 0;
			rPos[temp] = rNumNow;
			newLine = false;
		}
		R[rNumNow][count++] = temp;
		infile.get();
		if (infile.peek() == '\n')
		{
			rLen[rNumNow++] = count;
			infile.get();
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
			count = 0;
			newLine = false;
		}
		E[eNumNow][count++] = temp;
		infil.get();
		if (infil.peek() == '\n')
		{
			eLen[eNumNow++] = count;
			infil.get();
			newLine = true;
		}
	}
	eNumOrigin = eNumNow;
	infil.close();
}
void output()
{//当每行遇到零时停止输出  忽视那些含有零的消元结果
	string dirName = "result.txt";
	remove(dirName.c_str());
	ofstream outfile;
	outfile.open(dirName, ios::app);
	outfile << "1************被消元行*************2" << endl;
	for (int i = 0; i < eNumOrigin; i++)
	{
		for (int d = 0; d < eLen[i]; d++)
			outfile << E[i][d] << ' ';
		outfile << endl;
	}
}
void eliminate()
{
	for (int i = 0; i < eNumOrigin; i++)
	{
		int eTemp[maxM] = { 0 }, eNew[maxM] = { 0 };//本轮消去中用来取代旧的被消元行
		memcpy(eTemp, E[i], eLen[i] * sizeof(int));
		int len = eLen[i], k;
		while (len && rPos[eTemp[0]] != -1)
		{
			k = 1, len = 0;
			int p = rPos[eTemp[0]], q = rLen[rPos[eTemp[0]]];
			for (int j = 1; j < eLen[i] || k < q; j++)
			{
				while (R[p][k] > eTemp[j])
					eNew[len++] = R[p][k++];
				
				if (R[p][k] == eTemp[j])
					k++;
				else// if (R[p][k] < eTemp[j])
					eNew[len++] = eTemp[j];
			}
			memcpy(eTemp, eNew, len * sizeof(int));
			eLen[i] = len;
		}
		eLen[i] = rLen[rNumNow] = len;
		memcpy(E[i], eNew, len * sizeof(int));
		rPos[E[i][0]] = rNumNow;
		memcpy(R[rNumNow], E[i], len * sizeof(int));
		rNumNow++; eNumNow--;
	}
}
int main()
{
	input();
	eliminate();
	output();
	return 0;
}