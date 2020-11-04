#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>

#define TAG 1

enum direction
{
	up,
	down,
	left,
	right,
	upleft,
	upright,
	downleft,
	downright
};

void check(int** board, int** next, int* cont, int numprocs, int myrank, int x, int y) //tsekarei an einai adeios o pinakas h idios me prin
{
	int empty = 0, same = 1, i = 1, j = 1;
	int allempty = 0, allsame = 0;
	int *emptyprocs, *sameprocs, *allconts;
	emptyprocs = NULL;
	sameprocs = NULL;
	
//#pragma omp parallel for collapse(2)
	for(i = 1; i < x-1; i++)
	{
		for(j = 1; j < y-1; j++)
		{
			empty = empty + board[i][j];
			if(board[i][j] != next[i][j])
				same = 0;
		}
	}

	if(myrank == 0)
	{
		emptyprocs = malloc(sizeof(int) * numprocs);
		sameprocs = malloc(sizeof(int) * numprocs);
		allconts = malloc(sizeof(int)*numprocs);
	}
		MPI_Gather(&empty, 1, MPI_INT, emptyprocs, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Gather(&same, 1, MPI_INT, sameprocs, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(myrank == 0)
	{
		for(i = 0; i < numprocs; i++)
		{
			allempty = allempty + emptyprocs[i];
			allsame = allsame + sameprocs[i];
		}
		if((allempty == 0) || (allsame)) *cont = 0; //an einai ola dead h ola idia me ta prohgoumena
		else *cont = 1;
		for(i = 0; i < numprocs; i++)
		{
			allconts[i] = *cont;
		}	
	}
	MPI_Scatter(allconts, 1, MPI_INT, cont, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

void compute(int i, int j, int** board, int** next)
{
	int k = 0;
	k = board[i+1][j] + board[i-1][j] + board[i][j+1] + board[i][j-1] + 
		board[i+1][j+1] + board[i+1][j-1] + board[i-1][j+1] + board[i-1][j-1];
	if((k < 0) || (k > 8))
	{
		printf("Error! Calculation went wrong\n");
	}
	if((k == 3) && (board[i][j] == 0)) //3 geitones kai nekros (ara adeia thesi)
	{
		next[i][j] = 1;
	}
	else if((board[i][j] == 1) && ((k < 2) || (k > 3)))
	{
			next[i][j] = 0;
	}
	else next[i][j] = board[i][j];
}

void run(int** board, int** next, long int x, long int y, int myrank, int numprocs, int gens, MPI_Comm comm_cart)
{	
	int coordinates[2], neighbor_coordinates[2], source_neighbors[8], dest_neighbors[8];
	MPI_Request sends[8], recvs[8];
	MPI_Status status;
	int i, j;
	int xlimit = x - 2; //to orio tou x, dld mexri to shmeio pou exei dika tou stoixeia
	int ylimit = y - 2; //to orio tou y, dld mexri to shmeio pou exei dika tou stoixeia
	int runs = 0;
	int cont = 1;
	int** temp;
	//dimiourgia tipou dedomenwn sthlhs
	MPI_Datatype column;
	MPI_Type_vector(xlimit, 1, y, MPI_INT, &column);
	MPI_Type_commit(&column);
	MPI_Cart_coords(comm_cart, myrank, 2, coordinates);
	while(runs < gens)
	{
		MPI_Cart_shift(comm_cart, 1, 1, &source_neighbors[up], &dest_neighbors[up]);
		MPI_Cart_shift(comm_cart, 1, -1, &source_neighbors[down], &dest_neighbors[down]);
		MPI_Cart_shift(comm_cart, 0, -1, &source_neighbors[left], &dest_neighbors[left]);
		MPI_Cart_shift(comm_cart, 0, 1, &source_neighbors[right], &dest_neighbors[right]);
		

		neighbor_coordinates[0] = coordinates[0] + 1; 
		neighbor_coordinates[1] = coordinates[1] -1; 
		MPI_Cart_rank(comm_cart, neighbor_coordinates, &dest_neighbors[upleft]);
		source_neighbors[downright] = dest_neighbors[upleft]; 
		
		neighbor_coordinates[0] = coordinates[0] + 1; 
		neighbor_coordinates[1] = coordinates[1] + 1; 
		MPI_Cart_rank(comm_cart, neighbor_coordinates, &dest_neighbors[upright]);
		source_neighbors[downleft] = dest_neighbors[upright]; 
		
		neighbor_coordinates[0] = coordinates[0] - 1; 
		neighbor_coordinates[1] = coordinates[1] -1; 
		MPI_Cart_rank(comm_cart, neighbor_coordinates, &dest_neighbors[downleft]);
		source_neighbors[upright] = dest_neighbors[downleft]; 
		
		neighbor_coordinates[0] = coordinates[0] - 1; 
		neighbor_coordinates[1] = coordinates[1] +1; 
		MPI_Cart_rank(comm_cart, neighbor_coordinates, &dest_neighbors[downright]);
		source_neighbors[upleft] = dest_neighbors[downright]; 
		
		MPI_Isend(&board[1][1], ylimit, MPI_INT, dest_neighbors[up], TAG, comm_cart, &sends[up]);
		MPI_Isend(&board[xlimit][1], ylimit, MPI_INT, dest_neighbors[down], TAG, comm_cart, &sends[down]);
		MPI_Isend(&board[1][1], 1, column, dest_neighbors[left], TAG, comm_cart, &sends[left]);
		MPI_Isend(&board[1][ylimit], 1, column, dest_neighbors[right], TAG, comm_cart, &sends[right]);
		MPI_Isend(&board[1][1], 1, MPI_INT, dest_neighbors[upleft], TAG, comm_cart, &sends[upleft]);
		MPI_Isend(&board[1][ylimit], 1, MPI_INT, dest_neighbors[upright], TAG, comm_cart, &sends[upright]);
		MPI_Isend(&board[xlimit][1], 1, MPI_INT, dest_neighbors[downleft], TAG, comm_cart, &sends[downleft]);
		MPI_Isend(&board[xlimit][ylimit], 1, MPI_INT, dest_neighbors[downright], TAG, comm_cart, &sends[downright]);
					
		MPI_Irecv(&board[x-1][1], ylimit, MPI_INT, source_neighbors[down], MPI_ANY_TAG, comm_cart, &recvs[down]);
		MPI_Irecv(&board[0][1], ylimit, MPI_INT, source_neighbors[up], MPI_ANY_TAG, comm_cart, &recvs[up]);
		MPI_Irecv(&board[1][y-1], 1, column, source_neighbors[right], MPI_ANY_TAG, comm_cart, &recvs[right]);
		MPI_Irecv(&board[1][0], 1, column, source_neighbors[left], MPI_ANY_TAG, comm_cart, &recvs[left]);
		MPI_Irecv(&board[x-1][y-1], 1, MPI_INT, source_neighbors[downright], MPI_ANY_TAG, comm_cart, &recvs[downright]);
		MPI_Irecv(&board[x-1][0], 1, MPI_INT, source_neighbors[downleft], MPI_ANY_TAG, comm_cart, &recvs[downleft]);
		MPI_Irecv(&board[0][y-1], 1, MPI_INT, source_neighbors[upright], MPI_ANY_TAG, comm_cart, &recvs[upright]);
		MPI_Irecv(&board[0][0], 1, MPI_INT, source_neighbors[upleft], MPI_ANY_TAG, comm_cart, &recvs[upleft]);

#pragma omp parallel for collapse(2)
		for(i = 2; i < xlimit; i++)                    //ESWTERIKA
		{
			for(j = 2; j < ylimit; j++)
			{
				compute(i, j, board, next);
			}
		}

		MPI_Wait(&recvs[down], &status);					//DOWN
#pragma omp parallel for 
		for(j = 1; j < ylimit; j++)
				compute(xlimit, j, board, next);
		
		MPI_Wait(&recvs[up], &status);						//UP
//#pragma omp parallel for 
		for(j = 1; j < ylimit; j++);
				compute(1, j, board, next);
		
		MPI_Wait(&recvs[right], &status);					//RIGHT
#pragma omp parallel for
		for(i = 1; i < xlimit; i++)
			compute(i, ylimit, board, next);
		
		
		MPI_Wait(&recvs[left], &status);					//LEFT
#pragma omp parallel for
		for(i = 1; i < xlimit; i++)
			compute(i, 1, board, next);
		
		MPI_Wait(&recvs[downright], &status);
		compute(xlimit, ylimit, board, next);
		MPI_Wait(&recvs[downleft], &status);
		compute(xlimit, 1, board, next);
		MPI_Wait(&recvs[upright], &status);
		compute(1, ylimit, board, next);
		MPI_Wait(&recvs[upleft], &status);
		compute(1, 1, board, next);
		for (i = 0; i < 8; i++)
		{
			MPI_Wait(&sends[i], &status);
		}
		runs++;
		check(board,next,&cont, numprocs, myrank, x, y);
		if(cont == 0)
		{
			gens = runs;
		}
		temp = board;
		board = next;
		next = temp;
	}
	if (myrank == 0)
	{
		printf("Generations run: %d\n", runs);
	}
}

int main(int argc, char *argv[])
{
	long int N, M, gens;
	int *data, *nextdata;
	int **board, **next;
	int numprocs, myrank, i, j;
	int periods[2] = {1, 1};
	int dims[2];
	double start_t, end_t;
	long int x, y; //for the "table" every processor is in charge of
	omp_set_num_threads(4);
	int required = MPI_THREAD_SERIALIZED;
	int provided;
	MPI_Comm comm_cart;
	MPI_Init_thread(&argc, &argv, required, &provided);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	printf("Process: %d \n", myrank);
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
	
	x = M + 2; 																					//for left and right neighbors
	y = N + 2;      																				//for up and down neighbors
	
	if((numprocs==4)||(numprocs==9)||(numprocs==16)||(numprocs==25)||(numprocs==36))  						//Einai arithmoi poy exoyn akeraia riza, se auth th periptwsh oi dims einai sqrt(numprocs)
	{	
		dims[0] = (int) sqrt(numprocs); 
		dims[1] = dims[0];
	}
	else
	{
		if((numprocs%2)==0)
		{
			dims[0]=2;
			dims[1]=numprocs/2;
		}
		else
		{
			for(i=(numprocs/2) ; i>1 ; i--)
			{
				if((numprocs%i)==0)
				{
					dims[0]=i;
					dims[1]=numprocs/i;
					break;
				}
				
			}
			if(i==1)
			{
				printf("ERROR: number of procceses is prime, please enter a different number \n");
				return 1;
			}
		}
	}
	
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &comm_cart);
	
	
	//allocating space for each processor's "board" and the next board
	
	data = malloc(x*y*sizeof(int));
	board = malloc(x * sizeof(int));
	board[0] = data;
	
	for(i = 1; i < x; i++) 
	{
		board[i] = board[i-1] + y;
	}	
	nextdata = malloc(x*y*sizeof(int));
	next = malloc(x * sizeof(int));
	next[0] = nextdata;
	
	for(i = 1; i < x; i++) 
	{
		next[i] = next[i-1] + y;
	}
	
	srand(time(NULL)); //init. srand
	for(i = 1; i < x-1; i++)
	{
		for(j = 1; j < y-1; j++)
		{
			board[i][j] = rand()%2; //0 = dead, 1 = alive
		}
	}
	
	MPI_Barrier(comm_cart);
	start_t = MPI_Wtime();
	
	run(board, next, x, y, myrank, numprocs, gens, comm_cart);
	
	MPI_Barrier(comm_cart);
	end_t = MPI_Wtime();
	if (myrank == 0)
	{
		printf("Ending with %f time \n", end_t - start_t);
		printf("Destroying board... \n");
	}
	
	free(board[0]);
	free(board);
	free(next[0]);
	free(next);
	MPI_Finalize();
	return 0;
}
