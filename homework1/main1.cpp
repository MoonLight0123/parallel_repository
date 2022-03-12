#include<iostream>
#include<ctime>
#include<sys/time.h>

using namespace std;
int size =1000000;
double *a;
double sum=0;
timeval *start;
timeval *stop;
void init()
{
    a=new double[size];
    for(int i=0;i<size;i++)a[i]=i;
    sum=0;
    stop=new timeval();
    start=new timeval();
}
void destroy1()
{
    delete []a;
    delete stop,start;
    sum=0;
}
void trivialAlgorithm()
{
    for(int i=0;i<size;i++)sum+a[i];
}
void optimizationAlgorithm1()
{
    for (int m = size; m > 1; m /= 2)
        for (int i = 0; i < m / 2; i++)
            a[i ] = a[i * 2] + a[i * 2 + 1];
    sum=a[0];
}
void optimizationAlgorithm2(int n)
{
    if (n == 1)
        return;
    else
    {
        for (int i = 0; i < n / 2; i++)
            a[i]+=a[n-i-1];
        n = n / 2;
        optimizationAlgorithm2(n);
    }
}
void optimizationAlgorithm3()
{
    int sum1 = 0, sum2 = 0;
    for (int i = 0;i < size; i += 2) {
        sum1 += a[i];
        sum2 += a[i + 1];
    }
    sum = sum1 + sum2;
}
int main()
{
    timeval *start=new timeval();
    timeval *stop=new timeval();
    double durationTime=0.0;

    init();
    gettimeofday(start,NULL);
    trivialAlgorithm();
    gettimeofday(stop,NULL);
    destroy1();
    durationTime =stop->tv_sec*1000+double(stop->tv_usec)/1000-start->tv_sec*1000-double(start->tv_usec)/1000;
    cout << " trivialAlgorithm time: " << double(durationTime) << " ms" << endl;

    init();
    gettimeofday(start,NULL);
    optimizationAlgorithm1();
    gettimeofday(stop,NULL);
    destroy1();
    durationTime =stop->tv_sec*1000+double(stop->tv_usec)/1000-start->tv_sec*1000-double(start->tv_usec)/1000;
    cout << " optimizationAlgorithm1 time: " << double(durationTime) << " ms" << endl;

    init();
    gettimeofday(start,NULL);
    optimizationAlgorithm2(size-1);
    gettimeofday(stop,NULL);
    destroy1();
    durationTime =stop->tv_sec*1000+double(stop->tv_usec)/1000-start->tv_sec*1000-double(start->tv_usec)/1000;
    cout << " optimizationAlgorithm1 time: " << double(durationTime) << " ms" << endl;

    init();
    gettimeofday(start,NULL);
    optimizationAlgorithm3();
    gettimeofday(stop,NULL);
    destroy1();
    durationTime =stop->tv_sec*1000+double(stop->tv_usec)/1000-start->tv_sec*1000-double(start->tv_usec)/1000;
    cout << " optimizationAlgorithm1 time: " << double(durationTime) << " ms" << endl;



    return 0;
}
