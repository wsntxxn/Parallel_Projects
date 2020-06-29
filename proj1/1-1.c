#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#define ALLGATHER_TAG 0

void My_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, 
				  void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm){
	int world_rank;
	int world_size;
	MPI_Status status;
	MPI_Request request;

	MPI_Comm_size(comm, &world_size);
	MPI_Comm_rank(comm, &world_rank);

	int tsize;
    MPI_Type_size(recvtype, &tsize);

	if (world_rank)
		MPI_Isend(sendbuf, sendcount, sendtype, 0, ALLGATHER_TAG, comm, &request);

	if (!world_rank) {
		for (int i = 0; i < world_size; ++i) {
			if (!i) {
				memcpy(recvbuf, sendbuf, recvcount);
			}
			else {
				MPI_Recv(recvbuf + i * recvcount * tsize, recvcount, recvtype, i, ALLGATHER_TAG, comm, &status);
			}
		}
		
		for (int i = 1; i < world_size; ++i) {
			MPI_Send(recvbuf, world_size * recvcount, sendtype, i, ALLGATHER_TAG, comm);
		}
	}
	else
	{		
		MPI_Recv(recvbuf, world_size * sendcount, sendtype, 0, ALLGATHER_TAG, comm, &status);
	}
}

int main(int argc, char* argv[]) {
	int i;
	int world_rank;
	int world_size;
	int isend;
	int *recv_buffer;
	int *send_buffer;


	MPI_Status status;

	double time_elapsed;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	if (argc != 2) {
		if (!world_rank) {
			printf("Usage: %s <send number>\n", argv[0]);
			MPI_Abort(MPI_COMM_WORLD, -1);
		}
	}

	int number_send = atoi(argv[1]);
	send_buffer = malloc(sizeof(int) * number_send);
	recv_buffer = malloc(sizeof(int) * number_send * world_size);
	
	for (int i = 0 ; i < number_send; ++i)
		send_buffer[i] = i + number_send * world_rank;
	// isend = world_rank + 1;

	time_elapsed = -MPI_Wtime();
	My_Allgather(send_buffer, number_send, MPI_INT, recv_buffer, number_send, MPI_INT, MPI_COMM_WORLD);
	time_elapsed += MPI_Wtime();

	if (!world_rank) {
		/*
		for (int i = 0; i < world_size * number_send; ++i)
			printf("%d ", recv_buffer[i]);
		printf("\n");
		*/
		printf("Time elapsed of my implemented all gather: %f\n", time_elapsed);
	}

	send_buffer = malloc(sizeof(int) * number_send);
	recv_buffer = malloc(sizeof(int) * number_send * world_size);

	for (int i = 0 ; i < number_send; ++i)
		send_buffer[i] = i + number_send * world_rank;
	time_elapsed = -MPI_Wtime();
	MPI_Allgather(send_buffer, number_send, MPI_INT, recv_buffer, number_send, MPI_INT, MPI_COMM_WORLD);
	time_elapsed += MPI_Wtime();
	if (!world_rank) {
		/*
		for (int i = 0; i < world_size * number_send; ++i)
			printf("%d ", recv_buffer[i]);
		printf("\n");
		*/
		printf("Time elapsed of original implemented all gather: %f\n", time_elapsed);
	}

	MPI_Finalize();
	return 0;
}
