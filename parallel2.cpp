// можно параллельно вычитать из всех строк нужную

#include "stdio.h"
#include "math.h"
#include "time.h"
#include "stdlib.h"
#include "pthread.h"

typedef struct _ARGS
{
	double *mat;
	double *ans;
	int *ind;
	int k;
	int n;
	int thread_num;
	int total_threads;	
} ARGS;


void Print(int n, double *a, int *ind, int m);
double norm(int n, double *mat, double *ans, int *ind);


void enter(int choice, int n, double *a, double *b);
double f(int i, int j);
int proga(int n, double *mat, double *ans, int *ind, double eps, pthread_t *threads, int nthreads, ARGS *args);
void *func(void *pa);

void synchronize(int total_threads);
void step(double *mat, double *ans, int *ind, int k, int n, int thread_num, int total_threads);

int main()
{
	int nthreads = 4; //количество процессоров
	pthread_t *threads = (pthread_t *) malloc(nthreads * sizeof (pthread_t));
	
	printf("Enter the size of the matrix\n");
	int n;
	scanf("%d", &n);
	printf("How do you want to enter a matrix? 1 - from file, 2 - according to the formula\n");
	int choice;
	scanf("%d", &choice);
	
	double *mat = new double[n * n];
	double *mat_copy = new double[n * n];
	enter(choice, n, mat, mat_copy);
	int *ind = new int[n];
	for (int i = 0; i < n; i++)
	{
		ind[i] = i;
	}
	
	//Print(n, mat, ind, n);
	//printf("\n");
	double *ans = new double[n * n];
	for (int i = 0; i < n * n; i++)
	{
		if (i / n == i % n)
		{
			ans[i] = 1;
		}
		else
		{
			ans[i] = 0;
		}
	}
	double eps;
	printf("Enter accuracy\n");
	scanf("%lf", &eps);
	
	ARGS *args = (ARGS *)malloc(nthreads * sizeof(ARGS));
	for (int i = 0; i < nthreads; i++)
	{
		args[i].mat = mat;
		args[i].ind = ind;
		args[i].ans = ans;
		args[i].n = n;
		args[i].total_threads = nthreads;
		args[i].thread_num = i;
	}
	
	int start = clock();
	if (proga(n, mat, ans, ind, eps, threads, nthreads, args) == -1)
	{
		printf("det = 0\n");
		return 0;
	}
	int final =  clock();
	
	printf("How many elements of answer do you want to view?\n");
	int m;
	scanf("%d", &m);
	
	Print(n, ans, ind, m);
	printf("runtime = %lf seconds\n", (double)(final - start) / CLK_TCK);
	printf("norm of discrepancy = %.20lf\n", norm(n, mat_copy, ans, ind));
	
	free(threads);
	delete [] mat;
	delete [] ans;
	delete [] ind;
	return 0;
}

void synchronize(int total_threads)
{
	static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
	static int threads_in = 0;
	static int threads_out = 0;
	pthread_mutex_lock(&mutex);
	threads_in++;
	if (threads_in >= total_threads)
	{
		threads_out = 0;
		pthread_cond_broadcast(&condvar_in);
	}
	else
	{
		while (threads_in < total_threads)
		{
			pthread_cond_wait(&condvar_in, &mutex);
		}
	}
	threads_out++;
	if (threads_out >= total_threads)
	{
		threads_in = 0;
		pthread_cond_broadcast(&condvar_out);
	}
	else
	{
		while (threads_out < total_threads)
		{
			pthread_cond_wait(&condvar_out, &mutex);
		}
	}
	pthread_mutex_unlock(&mutex);
}

void Print(int n, double *a, int *ind, int m)
{
	if (m > n)
	{
		m = n;
	}
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			printf("%.10lf ", a[ind[i] * n + j]);
		}
		printf("\n");
	}
	printf("\n");
}

double norm(int n, double *mat, double *ans, int *ind)
{
	double w = 0;
	for (int i = 0; i < n; i++)
	{
		double sum = 0;
		for (int j = 0; j < n; j++)
		{
			double s = 0;
			for (int k = 0; k < n; k++)
			{
				s += mat[i * n + k] * ans[ind[k] * n + j];
			}
			if (i == j)
			{
				s -= 1;
			}
			sum += s;
			//printf("%lf ", s);
		}
	//printf("%lf\n", sum);
		if (sum > w)
		{
			w = sum;
		}
	}
	return w;
}

void enter(int choice, int n, double *a, double *b)
{
	if (choice == 1)
	{
		FILE *f = fopen("input.txt", "r");
		int n1;
		fscanf(f, "%d", &n1);
		for (int i = 0; i < n * n; i++)
		{
			fscanf(f, "%lf", &a[i]);
			b[i] = a[i];
		}
	}
	else
	{
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				a[i * n + j] = f(i, j);
				b[i * n + j] = a[i * n + j];
			}
		}
	}
}

double f(int i, int j)
{
	return fabs((double)(i - j));
}

int proga(int n, double *mat, double *ans, int *ind, double eps, pthread_t *threads, int nthreads, ARGS *args)
{
	//printf("hhhh\n");
	
	for (int k = 0; k < n; k++)
	{
		double max = fabs(mat[ind[k] * n + k]);
		int index = k;
		for (int j = k + 1; j < n; j++)
		{
			if (fabs(mat[ind[j] * n + k]) > max)
			{
				max = fabs(mat[ind[j] * n + k]);
				index = j;
			}
		}
		if (max < 0.000000000000001)
		{
			return -1;
		}
		int s = ind[k];
		ind[k] = ind[index];
		ind[index] = s;
		//Print(n, mat, ind, n);
		//printf("\n");
		double y = mat[ind[k] * n + k];
		for (int j = 0; j < n; j++)
		{
			mat[ind[k] * n + j] /= y;
			ans[ind[k] * n + j] /= y;
		}
		//printf("gggg\n");
		for (int i = 0; i < nthreads; i++)
		{
			//printf("fffff\n");
			args[i].k = k;
			if (pthread_create(threads + i, NULL, func, args + i))
			{
				printf("cannot create thread %d\n", i);
				return -1;
			}
		}
		for (int i = 0; i < nthreads; i++)
		{
			if (pthread_join(threads[i], NULL))
			{
				printf("cannot wait thread %d\n", i);
			}
		}
		//Print(n, mat, ind, n);
		//printf("\n");
	}
	return 0;
}

void *func(void *pa)
{
	ARGS *pargs = (ARGS *)pa;
	step(pargs->mat, pargs->ans, pargs->ind, pargs->k, pargs->n, pargs->thread_num, pargs->total_threads);
	return 0;
}

void step(double *mat, double *ans, int *ind, int k, int n, int thread_num, int total_threads)
{
	//Print(n, pargs->mat, pargs->ind, n);
	//printf("\n");
	for (int i = thread_num * n / total_threads; i < (thread_num + 1) * n / total_threads; i++)
	{
		if (i != ind[k])
		{
			double y = mat[i * n + k];
			for (int j = 0; j < n; j++)
			{
				mat[i * n + j] -= mat[ind[k] * n + j] * y;
				ans[i * n + j] -= ans[ind[k] * n + j] * y;
			}
		}
	}
	synchronize(total_threads);
}
