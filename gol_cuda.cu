#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#define threads_per_block 256

__global__ void run(int rows, int cols, int* board, int* next)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	int neighbors = 0;
	if(i < rows*cols)
	{
		if((i == 0) || (i%cols == 0))
		{
			if(i == 0)
			{
				neighbors = board[1] + board[cols-1] + board[cols] + board[(rows-1)*cols] + board[cols+1] + board[(rows*cols)-1] + board[(2*cols)-1] + board[(rows-1)*cols +1];
			}
			else if(i == ((rows-1)*cols))
			{
				neighbors = board[i+1] + board[rows*cols -1] + [(rows-2)*cols] + board[0] + board[(rows-2)*cols +1] + board[cols-1] + board[1] + board[i-1];
			}
			else
			{
				neighbors = board[i+1] + board[i+cols-1] + board[i-cols] + board[i+cols] + board[i+cols+1] + board[i-1] + board[i-1+2*cols] + board[i-cols+1];
			}
		}
		else if((i+1)%cols == 0)
		{
			if(i == (cols-1))
			{
				neighbors = board[i-1] + board[i-cols+1] + board[rows*cols -1] + board[i+cols] + board[i+cols-1] + board[(rows-1)*cols] + board[rows*cols -2] + board[i+1];
			}
			else if(i == (rows*cols -1))
			{
				neighbors = board[i-1] + board[i-cols+1] + board[i-cols] + board[cols-1] + board[cols-2] + board[i+1-2*cols] + board[i-cols-1] + board[0];
			}
			else
			{
				neighbors = board[i-1] + board[i-cols+1] + board[i-cols] + board[i+cols] + board[i+cols-1] + board[i+1-2*cols] + board[i-cols-1] + board[i+1];
			}
		}
		else
		{
			neighbors = board[i+1] + board[i-1] + board[i-cols] + board[i+cols] + board[i+cols-1] + board[i+cols+1] + board[i-cols-1] + board[i-cols+1];
		}
		if((board[i] == 1) && (neighbors<2 || neighbors>3)) next[i] = 0;
		else if((board[i] == 0) && (neighbors == 3)) next[i] = 1;
		else next[i] = board[i];
	}
}

int continue_running(int *board, int *next, int* same, int* empty, int rows, int cols)
{
	int i;
	int same = 1;
	int sum = 0;
	for(i = 0; i < rows*cols; i++)
	{
		sum = sum + next[i];
		if(board[i] != next[i])
			same = 0;
	}
	if((sum > 0) && (!same)) return 1;
	else return 0;
}

int main(int argc, char* argv[])
{
	int *board, *next;
	int *d_board, *d_next;
	int M, N, gens, i;
	N = strtol(argv[1], NULL, 10);
	gens = strtol(argv[2], NULL, 10);
	
	M = N*N;
	board = (int *)malloc(M*sizeof(int));
	next = (int *)malloc(M*sizeof(int));

	cudaMalloc(&d_board, M*sizeof(int));
	cudaMalloc(&d_next, M*sizeof(int));

	for(i = 0; i < M; ++i )
	{
		srand(time(NULL)); 
		board[i] = rand()%2;
		next[i] = 0;
	}
	int *temp;
	for(i = 0; i < gens; i++)
	{
		cudaMemcpy(d_board, board, M*sizeof(int), cudaMemcpyHostToDevice );
		cudaMemcpy(d_next, next, M*sizeof(int), cudaMemcpyHostToDevice );
		
		run<<(N+(threads_per_block-1))/threads_per_block, threads_per_block>>(N, N, d_board, d_next);
		
		cudaMemcpy(board, d_board, M*sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(next, d_next, M*sizeof(int), cudaMemcpyDeviceToHost);
		if(continue_running(board, next, N, N))
		{
			temp = board;
			board = next;
			next = temp;
		}
		else break;
	}
	
	cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) 
    printf("Error: %s\n", cudaGetErrorString(err));

	free(board);
	free(next);

	cudaFree(d_board);
	cudaFree(d_next);

	return 0;
}
