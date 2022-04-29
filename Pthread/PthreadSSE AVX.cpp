#include <iostream>
#include<pthread.h>
#include<semaphore.h>
#include<immintrin.h>
#include<sys/time.h>
using namespace std;
const int N = 2500, threadCount = 4;
sem_t semMain, semWorkstart[threadCount], semWorkend[threadCount];
float a[N][N];
struct threadParm
{
    int ID;
};
void* threadFunction(void* parm)
{
	__m128 t0, t1, t2, t3;
	threadParm* p=(threadParm*)parm;
	int id=p->ID;
	for (int k = 0; k < N; k++)
	{
		sem_wait(&semWorkstart[id]);

		for (int i = k + 1 + id; i < N; i += threadCount)//可考虑在此处cache优化
		{
			float temp2[4] = { a[i][k],a[i][k],a[i][k],a[i][k] };
			t0 = _mm_loadu_ps(temp2);
			int j;
			for (j = k + 1; j + 3 < N; j += 4)
			{
				t1 = _mm_loadu_ps(a[k] + j);
				t2 = _mm_loadu_ps(a[i] + j);
				t3 = _mm_mul_ps(t0, t1);
				t2 = _mm_sub_ps(t2, t3);
				_mm_storeu_ps(a[i] + j, t2);
			}
			for (; j < N; j++)
				a[i][j] -= a[i][k] * a[k][j];
			a[i][k] = 0.0;
		}
		sem_post(&semMain);
		sem_wait(&semWorkend[id]);
	}
	return NULL;
}
void init()
{
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			a[i][j] = float(rand()) / 10;
}
void SequentialAlgorithm()
{
	for (int k = 0; k < N; k++)
	{
		for (int j = k + 1; j < N; j++)
			a[k][j] /= a[k][k];
		a[k][k] = 1.0;
		for (int i = k + 1; i < N; i++) {

			for (int j = k + 1; j < N; j++)
				a[i][j] -= a[i][k] * a[k][j];
			a[i][k] = 0;
		}
	}
}
void PthreadParallelAlgorithm()
{
	sem_init(&semMain, 0, 0);
	for(int id=0;id<threadCount;id++){
	sem_init(&semWorkstart[id], 0, 0);
	sem_init(&semWorkend[id], 0, 0);
	}
	pthread_t handle[threadCount];
	threadParm tp[threadCount];
	__m128 t0, t1, t2;
	for (int id = 0; id < threadCount; id++)
    {
        tp[id].ID=id;
		pthread_create(&handle[id], NULL, threadFunction, &tp[id]);
	}
	for (int k = 0; k < N; k++)
	{
		float temp1[4] = { a[k][k],a[k][k],a[k][k],a[k][k] };
		t0 = _mm_loadu_ps(temp1);
		int j;
		for (j = k + 1; j + 3 < N; j += 4)//j+4的值为下一次load操作的起始位置
		{
			t1 = _mm_loadu_ps(a[k] + j);
			t2 = _mm_div_ps(t1, t0);
			_mm_storeu_ps(a[k] + j, t2);
		}
		for (; j < N; j++)
			a[k][j] /= a[k][k];
		a[k][k] = 1.0;
		for (int id = 0; id < threadCount; id++)
			sem_post(&semWorkstart[id]);
		for (int id = 0; id < threadCount; id++)
			sem_wait(&semMain);
		for (int id = 0; id < threadCount; id++)
			sem_post(&semWorkend[id]);
	}
	for (int id = 0; id < threadCount; id++)
		pthread_join(handle[id], NULL);
	sem_destroy(&semMain);
	for(int id=0;id<threadCount;id++){
	sem_destroy(&semWorkstart[id]);
	sem_destroy(&semWorkend[id]);
	}
}
void ParallelAlgorithm()
{
	__m128 t0, t1, t2, t3;
	for (int k = 0; k < N; k++)
	{
		float temp1[4] = { a[k][k],a[k][k],a[k][k],a[k][k] };
		t0 = _mm_loadu_ps(temp1);
		int j;
		for (j = k + 1; j + 3 < N; j += 4)//j+4的值为下一次load操作的起始位置
		{
			t1 = _mm_loadu_ps(a[k] + j);
			t2 = _mm_div_ps(t1, t0);
			_mm_storeu_ps(a[k] + j, t2);
		}
		for (; j < N; j++)
			a[k][j] /= a[k][k];

		a[k][k] = 1.0;
		for (int i = k + 1; i < N; i++)
		{
			float temp2[4] = { a[i][k],a[i][k],a[i][k],a[i][k] };
			t0 = _mm_loadu_ps(temp2);
			int j;
			for (j = k + 1; j + 3 < N; j += 4)
			{
				t1 = _mm_loadu_ps(a[k] + j);
				t2 = _mm_loadu_ps(a[i] + j);
				t3 = _mm_mul_ps(t0, t1);
				t2 = _mm_sub_ps(t2, t3);
				_mm_storeu_ps(a[i] + j, t2);
			}
			for (; j < N; j++)
				a[i][j] -= a[i][k] * a[k][j];
			a[i][k] = 0.0;
		}
	}
}
void AlignedParallelAlgorithm()
{
	__m128 t0, t1, t2, t3;
	for (int k = 0; k < N; k++)
	{
		float temp1[4] = { a[k][k],a[k][k],a[k][k],a[k][k] };
		t0 = _mm_loadu_ps(temp1);
		int j;
		for (j = k + 1; j < N; j++) {
			if (((size_t)(a[k] + j)) % 16 == 0)
				break;
			a[k][j] /= a[k][k];
		}
		for (; j + 3 < N; j += 4)//j+4的值为下一次load操作的起始位置
		{
			t1 = _mm_load_ps(a[k] + j);
			t2 = _mm_div_ps(t1, t0);
			_mm_store_ps(a[k] + j, t2);
		}
		for (; j < N; j++)
			a[k][j] /= a[k][k];

		a[k][k] = 1.0;
		for (int i = k + 1; i < N; i++)
		{
			float temp2[4] = { a[i][k],a[i][k],a[i][k],a[i][k] };
			t0 = _mm_loadu_ps(temp2);
			int j;
			for (j = k + 1; j < N; j++) {
				if (((size_t)(a[i] + j)) % 16 == 0)
					break;
				a[i][j] -= a[i][k] * a[k][j];
			}
			for (; j + 3 < N; j += 4)
			{
				t1 = _mm_loadu_ps(a[k] + j);
				t2 = _mm_load_ps(a[i] + j);
				t3 = _mm_mul_ps(t0, t1);
				t2 = _mm_sub_ps(t2, t3);
				_mm_store_ps(a[i] + j, t2);
			}
			for (; j < N; j++)
				a[i][j] -= a[i][k] * a[k][j];
			a[i][k] = 0.0;
		}
	}

}
void show()
{

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++)
			cout << a[i][j] << ' ';
		cout << endl;
	}
}


int main()
{
    timeval *start=new timeval();
    timeval *stop=new timeval();
    double durationTime=0.0;

    init();
    gettimeofday(start,NULL);
    SequentialAlgorithm();
    gettimeofday(stop,NULL);
    durationTime =stop->tv_sec*1000+double(stop->tv_usec)/1000-start->tv_sec*1000-double(start->tv_usec)/1000;
    cout << " SequentialAlgorithm time: " << double(durationTime) << " ms" << endl;

    init();
    gettimeofday(start,NULL);
    ParallelAlgorithm();
    gettimeofday(stop,NULL);
    durationTime =stop->tv_sec*1000+double(stop->tv_usec)/1000-start->tv_sec*1000-double(start->tv_usec)/1000;
    cout << " ParallelAlgorithm time: " << double(durationTime) << " ms" << endl;

    init();
    gettimeofday(start,NULL);
    PthreadParallelAlgorithm();
    gettimeofday(stop,NULL);
    durationTime =stop->tv_sec*1000+double(stop->tv_usec)/1000-start->tv_sec*1000-double(start->tv_usec)/1000;
    cout << " PthreadParallelAlgorithm time: " << double(durationTime) << " ms" << endl;


    //show();


    return 0;
}
//#include<iostream>
//#include<pthread.h>
//#include<semaphore.h>
//#pragma comment(lib, "pthreadVC2.lib")
//using namespace std;
//int thread_count = 4;
//void* Hello(void* rank)
//{
//	long my_rank = (long)rank);  /* Use long in case of 64-bit system */
//
//	printf("Hello from thread %ld of %d\n", my_rank, thread_count);
//
//	return NULL;
//
//}
//int main()
//{
//	long thread; /* Use long in case of a 64-bit system */
//	pthread_t* thread_handles;
//
//	/* Get number of threads from command line */
//
//	thread_handles = (pthread_t *)malloc(thread_count * sizeof(pthread_t));
//	for (thread = 0; thread < thread_count; thread++)
//		pthread_create(&thread_handles[thread], NULL, Hello, (void*)thread);
//
//	printf("Hello from the main thread\n");
//
//	for (thread = 0; thread < thread_count; thread++)
//		pthread_join(thread_handles[thread], NULL);
//
//	free(thread_handles);
//	return 0;
//}
