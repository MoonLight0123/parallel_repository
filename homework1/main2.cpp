#include<iostream>
#include<sys/time.h>
using namespace std;
const int size=4000;
int matrix[size][size];
int b[size];
int sum[size];
void init()
{
    for(int i=0;i<size;i++)
    {
        b[i]=i;
        for(int j=0;j<size;j++)
        {
            matrix[i][j]=i+j;
        }
    }
}
void col_major()
{
    for(int i=0;i<size;i++)
    {

        sum[i]=0;
        for(int j=0;j<size;j++)
        {

            sum[i]=matrix[j][i]*b[j];
        }
    }
}
void row_major()
{
    for(int i=0;i<size;i++)
        sum[i]=0;
    for(int j=0;j<size;j++)
    {

        for(int i=0;i<size;i++)
            sum[i]+=matrix[j][i]*b[j];
    }
}
int main()
{
    timeval *start=new timeval();
    timeval *stop=new timeval();
    double durationTime=0.0;
    init();
    gettimeofday(start,NULL);
    col_major();
    gettimeofday(stop,NULL);
    durationTime =stop->tv_sec*1000+double(stop->tv_usec)/1000-start->tv_sec*1000-double(start->tv_usec)/1000;
    cout << "col_major time: " << double(durationTime) << " ms" << endl;

    init();
    gettimeofday(start,NULL);
    row_major();
    gettimeofday(stop,NULL);
    durationTime =stop->tv_sec*1000+double(stop->tv_usec)/1000-start->tv_sec*1000-double(start->tv_usec)/1000;
    cout << "row_major time: " << double(durationTime) << " ms" << endl;
    return 0;
}
