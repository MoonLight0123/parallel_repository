#include <iostream>
#include<sys/time.h>
#include<immintrin.h>
using namespace std;
const int N=1400;
float a[N][N];
//void* MallocAlign16(size_t size)
//{
//	//错误版本
//	/*void *mem = malloc(1024 + 16);
//	void *ptr = ((char *)mem + 16) & ~(char *)0x0F;
//	free(mem);*/
//	int ptrSize = sizeof(void*);
//	char* ptr = (char*)malloc(size + 16 + ptrSize);
//	char* alignedPtr = (char*)(((size_t)ptr) + 15 & ~15);
//	if ((alignedPtr - ptr) < ptrSize)
//	{
//		alignedPtr += 16;
//	}
//	*((size_t*)(alignedPtr - ptrSize)) = (size_t)ptr;
//	return (void*)alignedPtr;
//}
//void FreeAlign16(void* ptr)
//{
//	int ptrSize = sizeof(void*);
//	free((void *) *((size_t *)(((char *)ptr) - ptrSize)));
//}
void init()
{
    for(int i=0;i<N;i++)
        for(int j=0;j<N;j++)
            a[i][j]=float(rand())/10;
}
void SequentialAlgorithm()
{
    for(int k=0;k<N;k++)
    {
        for(int j=k+1;j<N;j++)
            a[k][j]/=a[k][k];
        a[k][k]=1.0;
        for(int i=k+1;i<N;i++){

            for(int j=k+1;j<N;j++)
                a[i][j]-=a[i][k]*a[k][j];
            a[i][k]=0;
        }
    }
}
void ParallelAlgorithm()
{
    __m128 t0,t1,t2,t3;
    for(int k=0;k<N;k++)
    {
        float temp1[4]={a[k][k],a[k][k],a[k][k],a[k][k]};
        t0=_mm_loadu_ps(temp1);
        int j;
        for(j=k+1;j+3<N;j+=4)//j+4的值为下一次load操作的起始位置
        {
            t1=_mm_loadu_ps(a[k]+j);
            t2=_mm_div_ps(t1,t0);
            _mm_storeu_ps(a[k]+j,t2);
        }
        for(;j<N;j++)
            a[k][j]/=a[k][k];

        a[k][k]=1.0;
        for(int i=k+1;i<N;i++)
        {
            float temp2[4]={a[i][k],a[i][k],a[i][k],a[i][k]};
            t0=_mm_loadu_ps(temp2);
            int j;
            for(j=k+1;j+3<N;j+=4)
            {
                t1=_mm_loadu_ps(a[k]+j);
                t2=_mm_loadu_ps(a[i]+j);
                t3=_mm_mul_ps(t0,t1);
                t2=_mm_sub_ps(t2,t3);
                _mm_storeu_ps(a[i]+j,t2);
            }
            for(;j<N;j++)
                 a[i][j]-=a[i][k]*a[k][j];
            a[i][k]=0.0;
        }
    }

}
void AlignedParallelAlgorithm()
{
    __m128 t0,t1,t2,t3;
    for(int k=0;k<N;k++)
    {
        float temp1[4]={a[k][k],a[k][k],a[k][k],a[k][k]};
        t0=_mm_loadu_ps(temp1);
        int j;
        for(j=k+1;j<N;j++){
            if(((size_t)(a[k]+j))%16==0)
                break;
            a[k][j]/=a[k][k];
        }
        for(;j+3<N;j+=4)//j+4的值为下一次load操作的起始位置
        {
            t1=_mm_load_ps(a[k]+j);
            t2=_mm_div_ps(t1,t0);
            _mm_store_ps(a[k]+j,t2);
        }
        for(;j<N;j++)
            a[k][j]/=a[k][k];

        a[k][k]=1.0;
        for(int i=k+1;i<N;i++)
        {
            float temp2[4]={a[i][k],a[i][k],a[i][k],a[i][k]};
            t0=_mm_loadu_ps(temp2);
            int j;
            for(j=k+1;j<N;j++){
                   if(((size_t)(a[i]+j))%16==0)
                    break;
                 a[i][j]-=a[i][k]*a[k][j];
            }
            for(;j+3<N;j+=4)
            {
                t1=_mm_loadu_ps(a[k]+j);
                t2=_mm_load_ps(a[i]+j);
                t3=_mm_mul_ps(t0,t1);
                t2=_mm_sub_ps(t2,t3);
                _mm_store_ps(a[i]+j,t2);
            }
            for(;j<N;j++)
                 a[i][j]-=a[i][k]*a[k][j];
            a[i][k]=0.0;
        }
    }

}
void show()
{

    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++)
            cout<<a[i][j]<<' ';
        cout<<endl;
    }
}


int main()
{


    timeval *start=new timeval();
    timeval *stop=new timeval();
    double durationTime=0.0;

//    init();
//    gettimeofday(start,NULL);
//    SequentialAlgorithm();
//    gettimeofday(stop,NULL);
//    durationTime =stop->tv_sec*1000+double(stop->tv_usec)/1000-start->tv_sec*1000-double(start->tv_usec)/1000;
//    cout << " SequentialAlgorithm time: " << double(durationTime) << " ms" << endl;

    init();
    gettimeofday(start,NULL);
    ParallelAlgorithm();
    gettimeofday(stop,NULL);
    durationTime =stop->tv_sec*1000+double(stop->tv_usec)/1000-start->tv_sec*1000-double(start->tv_usec)/1000;
    cout << " ParallelAlgorithm time: " << double(durationTime) << " ms" << endl;

    init();
    gettimeofday(start,NULL);
    AlignedParallelAlgorithm();
    gettimeofday(stop,NULL);
    durationTime =stop->tv_sec*1000+double(stop->tv_usec)/1000-start->tv_sec*1000-double(start->tv_usec)/1000;
    cout << " AilgnedParallelAlgorithm time: " << double(durationTime) << " ms" << endl;
   // init();
    //gettimeofday(start,NULL);
    //SSE_gauss(a);
    //gettimeofday(stop,NULL);
    //durationTime =stop->tv_sec*1000+double(stop->tv_usec)/1000-start->tv_sec*1000-double(start->tv_usec)/1000;
    //cout << " SSE_gauss time: " << double(durationTime) << " ms" << endl;
    return 0;
}
