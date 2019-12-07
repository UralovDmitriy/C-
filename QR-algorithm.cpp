#include "stdio.h"
#include "math.h"
#include "time.h"

void Print_matrix(int n, double *a, int m);

void Print(int n, double *own, int m);
double norm1(int n, double *a, double *own);
double norm2(int n, double *a, double *b);


int enter(int choice, int n, double *a, double *b);
double f(int i, int j, int n);
int proga(int n, double *mat, double *own, double eps, double *sin, double * cos, double *xk);


int main()
{
	printf("Enter the size of the matrix\n");
	int n;
	scanf("%d", &n);
	printf("How do you want to enter a matrix? 1 - from file, 2 - according to the formula\n");
	int choice;
	scanf("%d", &choice);
	
	double *mat = new double[n * n];
	double *mat_copy = new double[n * n];
	double *own = new double[n];
	double *sin = new double[n - 1];
	double *cos = new double[n - 1];
	double *xk = new double[n - 1];
	if (enter(choice, n, mat, mat_copy) == -1)
	{
		printf("Enter error\n");
		return 0;
	}
	//Print_matrix(n, mat, n);
	double eps;
	printf("Enter accuracy\n");
	scanf("%lf", &eps);
	
	int start =  clock();
	if (proga(n, mat, own, eps, sin, cos, xk) == -1)
	{
		printf("Proga error\n");
		return 0;
	}
	int final =  clock();
	
	printf("How many elements of answer do you want to view?\n");
	int m;
	scanf("%d", &m);
	
	Print(n, own, m);
	printf("\nruntime = %lf seconds\n", (double)(final - start) / CLK_TCK);
	printf("norm of discrepancy by matrix trace = %.16lf\n", norm1(n, mat_copy, own));
	printf("norm of discrepancy by length as vector = %.16lf\n", norm2(n, mat_copy, mat));
	
	delete [] mat_copy;
	delete [] mat;
	delete [] own;
	delete [] sin;
	delete [] cos;
	delete [] xk;
	return 0;
}

void Print_matrix(int n, double *a, int m)
{
	if (m > n)
	{
		m = n;
	}
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			printf("%.5lf ", a[i * n + j]);
		}
		printf("\n");
	}
	printf("\n");
}

void Print(int n, double *own, int m)
{
	if (m > n)
	{
		m = n;
	}
	for (int i = 0; i < m; i++)
	{
		printf("%.16lf ", own[i]);
	}
}

double norm1(int n, double *a, double *own)
{
	double tr = 0;
	for (int i = 0; i < n; i++)
	{
		tr += a[i * n + i];
		tr -= own[i];
	}	
	return tr;
}

double norm2(int n, double *a, double *b)
{
	double l1 = 0, l2 = 0;
	for (int i = 0; i < n * n; i++)
	{
		l1 += a[i] * a[i];
		l2 += b[i] * b[i];
	}
	return sqrt(l1) - sqrt(l2);
}

int enter(int choice, int n, double *a, double *b)
{
	if (choice == 1)
	{
		FILE *f = fopen("input.txt", "r");
		if (f == NULL)
		{
			return -1;
		}
		int n1;
		if (fscanf(f, "%d", &n1) != 1)
		{
			return -1;
		}
		for (int i = 0; i < n * n; i++)
		{
			if (fscanf(f, "%lf", &a[i]) != 1)
			{
				return -1;
			}
			b[i] = a[i];
		}
	}
	else
	{
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				a[i * n + j] = f(i, j, n);
				b[i * n + j] = a[i * n + j];
			}
		}
	}
	return 0;
}

double f(int i, int j, int n)
{
	if (i < n-1 && j < n-1)
	 return 1;
	 if (i == n-1)
	 return j + 1;
	 return i+1;
}

int proga(int n, double *mat, double * own, double eps, double *sin, double *cos, double *xk)
{
	//приведение к почти треугольному виду
	for (int k = 0; k < n - 2; k++)
	{
		double Sk = 0;
		for (int j = k+2; j < n; j++)
		{
			Sk += mat[j * n + k] * mat[j * n + k];
		}
		double a1_norm = sqrt(mat[(k+1) * n + k] * mat[(k+1) * n + k] + Sk);
		xk[0] = mat[(k+1) * n + k] - a1_norm;
		double xk_norm = sqrt(xk[0] * xk[0] + Sk);
		if (xk_norm < eps)
		{
			continue;
		}
		xk_norm = 1 / xk_norm;
		xk[0] *= xk_norm;
		for (int j = 1; j < n - k - 1; j++)
		{
			xk[j] = mat[(k + j + 1) * n + k] * xk_norm;
		}

		for (int j = k; j < n; j++)
		{
			double scalar_product = 0;
			for (int i = k+1; i < n; i++)
			{
				scalar_product += xk[i - k - 1] * mat[i * n + j];
			}
			for (int i = k + 1; i < n; i++)
			{
				mat[i * n + j] = mat[i * n + j] - 2 * xk[i - k - 1] * scalar_product;
			}
		}

		for (int i = 0; i < n; i++)
		{
			double scalar_product = 0;
			for (int j = k+1; j < n; j++)
			{
				scalar_product += xk[j - k - 1] * mat[i * n + j];
			}
			for (int j = k+1; j < n; j++)
			{
				mat[i * n + j] = mat[i * n + j] - 2 * scalar_product * xk[j - k - 1];
			}
		}
	}
	printf("Reduction to almost triangular form done\n");
	//Print_matrix(n, mat, n);

	int iterations = 0;
	int q = 0;

	//QR - алгоритм
	while (1)
	{
		iterations++;
		//сдвиг для ускорения сходимости
		double sk = mat[(n-1-q) * n + n-1-q];
		for (int i = 0; i < n - q; i++)
		{
			mat[i*n + i] -= sk;
		}
		//Print_matrix(n, mat, n);
		
		//QR-разложение
		for (int k = 0; k < n - q - 1; k++)
		{
			if (fabs(mat[(k+1)*n + k]) < eps)
			{
				mat[(k+1)*n + k] = 0;
				sin[k] = 0;
				if (fabs(mat[k*n + k]) < eps)
				{
					cos[k] = 0;
					mat[k*n + k] = 0;
				}
				else
				{
					cos[k] = 1;
				}
				continue;
			}
			double norm = sqrt(mat[k*n + k] * mat[k*n + k] + mat[(k+1)*n + k] * mat[(k+1)*n + k]);
			double norm1 = 1 / norm;
			sin[k] = -mat[(k+1)*n + k] * norm1;
			cos[k] = mat[k*n + k] * norm1;
			mat[k * n + k] = norm;
			mat[(k+1) * n + k] = 0;
			for (int j = k + 1; j < n - q; j++)
			{
				double a = mat[k * n + j];
				double b = mat[(k+1) * n + j];
				mat[k * n + j] = a * cos[k] - b * sin[k];
				mat[(k+1) * n + j] = a * sin[k] + b * cos[k];
			}
		}
		//Print_matrix(n, mat, n);
				
		//нахождение RQ
		for (int k = 0; k < n - 1 - q; k++)
		{
			for (int i = 0; i <= k + 1; i++)
			{
				if (fabs(sin[k]) < eps && fabs(cos[k]) < eps)
				{
					continue;
				}
				double a = mat[i * n + k];
				double b = mat[i * n + k + 1];
				mat[i * n + k] = a * cos[k] + b * (-sin[k]);
				mat[i * n + k + 1] = a * (sin[k]) + b * cos[k];
			}
		}
		//Print_matrix(n, mat, n);
		
		//сдвиг обратно
		for (int i = 0; i < n - q; i++)
		{
			mat[i*n + i] += sk;
		}
		//Print_matrix(n, mat, n);
		
		//если обнулился mat[n][n - 1], получили одно собственное значение
		if (mat[(n-q-1)*n + n-q-2] < eps)
		{
			own[q] = mat[(n-q-1)*n + n-q-1];
			q++;
			if (q == n - 2)
			{
				double d = sqrt((mat[0] - mat[n + 1]) * (mat[0] - mat[n + 1]) + 4 * mat[1] * mat[n]);
				double c = mat[0] + mat[n + 1];
				own[n - 2] = (c + d) / 2;
				own[n - 1] = (c - d) / 2;
				printf("Number of steps in QR-algorithm = %d\n", iterations);
				return 1;
			}
		}		
	}
}
