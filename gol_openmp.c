#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>

int continue_running(int **board, int **next, long int M, long int N)
{
	long int i, j;
	int same = 1;
	long int sum = 0;
	#pragma omp parallel for collapse(2) reduction(+: sum)
	for(i = 0; i < M; i++)
	{
		for(j = 0; j < N ; j++)
		{
			sum = sum + next[i][j];
			if(board[i][j] != next[i][j])
				same = 0;
		}
	}
	if((sum > 0) && (!same)) return 1;
	else return 0;
}

void calc(int **board, int **next, long int M, long int N)
{
	long int i, j, sum, m, n;
	m = M - 1;
	n = N - 1;
	sum = 0;
	#pragma omp parallel for collapse(2) reduction(+:sum)
	for(i = 1; i < m; i++) //ta eswterika prwta
	{
		for(j = 1; j < n; j++)
		{
			sum = board[i+1][j] + board[i-1][j] + board[i][j+1] + board[i][j-1] + 
			board[i+1][j+1] + board[i+1][j-1] + board[i-1][j+1] + board[i-1][j-1];
			if((board[i][j] == 1) && (sum<2 || sum>3)) next[i][j] = 0;
			else if((board[i][j] == 0) && (sum == 3)) next[i][j] = 1;
			else next[i][j] = board[i][j];
		}
	}
	#pragma omp parallel for reduction(+:sum)
	for(j = 1; j < n; j++) //prwth k teleutaia grammh (xwris gwnia)
	{
		//prwth
		sum = board[1][j] + board[0][j+1] + board[0][j-1] + 
			board[1][j+1] + board[1][j-1] +
			board[m][j] + board[m][j-1] + board[m][j+1]; //periodicity
		if((board[0][j] == 1) && (sum<2 || sum>3)) next[0][j] = 0;
		else if((board[0][j] == 0) && (sum == 3)) next[0][j] = 1;
		else next[0][j] = board[0][j];

		//teleutaia
		sum = board[m-1][j] + board[m][j+1] + board[m][j-1] + 
			board[m-1][j+1] + board[m-1][j-1] +
			board[0][j] + board[0][j-1] + board[0][j+1]; //periodicity
		if((board[m][j] == 1) && (sum<2 || sum>3)) next[m][j] = 0;
		else if((board[m][j] == 0) && (sum == 3)) next[m][j] = 1;
		else next[m][j] = board[m][j];
	}
	#pragma omp parallel for reduction(+:sum)
	for(i = 1; i < m; i++) //prwth k teleutaia sthlh (xwris gwnia)
	{
		//prwth
		sum = board[i][1] + board[i+1][0] + board[i-1][0] + board[i+1][1] + board[i-1][1] +
			board[i][n] + board[i-1][n] + board[i+1][n]; //periodicity
		if((board[i][0] == 1) && (sum<2 || sum>3)) next[i][0] = 0;
		else if((board[i][0] == 0) && (sum == 3)) next[i][0] = 1;
		else next[i][0] = board[i][0];

		//teleutaia
		sum = board[i][n-1] + board[i+1][n] + board[i-1][n] + board[i+1][n-1] + board[i-1][n-1] +
		board[i][0] + board[i-1][0] + board[i+1][0]; //periodicity
		if((board[i][n] == 1) && (sum<2 || sum>3)) next[i][n] = 0;
		else if((board[i][n] == 0) && (sum == 3)) next[i][n] = 1;
		else next[i][n] = board[i][n];
	}
	//gwnies
	//x = 0 y = 0
	sum = board[0][1] + board[1][0] + board[1][1]  + 
		board[m][n] + board[m][0] + board[m][1] + board[0][n] + board[1][n]; //periodicity
	if((board[0][0] == 1) && (sum<2 || sum>3)) next[0][0] = 0;
	else if((board[0][0] == 0) && (sum == 3)) next[0][0] = 1;
	else next[0][0] = board[0][0];
	//x = m y = n
	sum = board[m-1][n] + board[m][n-1] + board[m-1][n-1] +
		board[m-1][0] + board[m][0] + board[0][0] + board[0][n] + board[0][n-1]; //periodicity
	if((board[m][n] == 1) && (sum<2 || sum>3)) next[m][n] = 0;
	else if((board[m][n] == 0) && (sum == 3)) next[m][n] = 1;
	else next[m][n] = board[m][n];
	//x = m y = 0
	sum = board[m-1][0] + board[m][1] + board[m-1][1] +
		board[m-1][n] + board[m][n] + board[0][n] + board[0][0] + board[0][1]; //periodicity
	if((board[m][0] == 1) && (sum<2 || sum>3)) next[m][0] = 0;
	else if((board[m][0] == 0) && (sum == 3)) next[m][0] = 1;
	else next[m][0] = board[m][0];
	//x = 0 y = n
	sum = board[1][n] + board[0][n-1] + board[1][n-1] +
		board[m][n-1] + board[m][n] + board[m][0] + board[0][0] + board[1][0]; //periodicity
	if((board[0][n] == 1) && (sum<2 || sum>3)) next[0][n] = 0;
	else if((board[0][n] == 0) && (sum == 3)) next[0][n] = 1;
	else next[0][n] = board[0][n];
}

void run(int **board, int **next, long int M, long int N, long int gens)
{
	int** temp;
	int flag = 1;
	int runs = 0;
	while(flag)
	{
		runs++;
		calc(board, next, M, N);
		if(runs >= gens)
		{
			flag = 0;
			break;
		}
		if(continue_running(board, next, M, N))
		{
			temp = board;
			board = next;
			next = temp;
		}
		else flag = 0;
	}
}

int main(int argc, char *argv[])
{
	long int M, N, gens;
	long int i, j;
	int **board, **next;
	double start_t, end_t;
	
	M = strtol(argv[1], NULL, 10);
	N = strtol(argv[2], NULL, 10);
	if(argc==4)
	{
		gens = strtol(argv[3], NULL, 10);
	}
	else
	{
		gens=999;									// THRESHOLD (yparxoun periptwseis opoy to game of life paiftei se loop opote dinetai to 999 tyxaia)
	}
	
	board = (int **)malloc(M * sizeof(int *));
	#pragma omp parallel for
    for (i = 0; i < N; i++)
    {
         board[i] = (int *)malloc(M * sizeof(int));
	}

	next = (int **)malloc(M * sizeof(int *));
	#pragma omp parallel for
    for (i = 0; i < N; i++)
    {
         next[i] = (int *)malloc(M * sizeof(int));
	}
	//randomize
	srand(time(NULL)); 
	for(i = 0; i < M; i++)
	{
		for(j = 0; j < N; j++)
		{
			board[i][j] = rand()%2;
		}
	}
	start_t = omp_get_wtime();
	run(board, next, M, N, gens);
	printf("Destroying boards.. \n");
	#pragma omp parallel for
	for(i = 0; i < M; i++)
    {
		free(board[i]);
	}
	for(i = 0; i < M; i++)
    {
		free(next[i]);
	}
	free(board);
	free(next);
	end_t = omp_get_wtime();
	printf("Ending with %f time \n", end_t - start_t);
	return 0;
}
