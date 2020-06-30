#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#define ALLGATHER_TAG 0

void My_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, 
				  void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm){
	int rank;
	int size;
	MPI_Status status;
	MPI_Request request;

	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);

	int tsize;
    MPI_Type_size(recvtype, &tsize);

	if (rank)
		MPI_Isend(sendbuf, sendcount, sendtype, 0, ALLGATHER_TAG, comm, &request);

	if (!rank) {
		for (int id = 0; id < size; ++id) {
			if (!id) {
				memcpy(recvbuf, sendbuf, recvcount);
			}
			else {
				MPI_Recv(recvbuf + id * recvcount * tsize, recvcount, recvtype, id, ALLGATHER_TAG, comm, &status);
			}
		}
	}
    MPI_Bcast(recvbuf, size * recvcount, sendtype, 0, comm);
}

int main(int argc, char* argv[]) {
	int i;
	int rank;
	int size;
	int isend;
	int *recv_buffer;
	int *send_buffer;


	MPI_Status status;

	double time_elapsed;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (argc != 2) {
		if (!rank) {
			printf("Usage: %s <send number>\n", argv[0]);
			MPI_Abort(MPI_COMM_WORLD, -1);
		}
	}

	int number_send = atoi(argv[1]);
	send_buffer = malloc(sizeof(int) * number_send);
	recv_buffer = malloc(sizeof(int) * number_send * size);
	
	for (int i = 0 ; i < number_send; ++i)
		send_buffer[i] = i + number_send * rank;

    MPI_Barrier(MPI_COMM_WORLD);
	time_elapsed = -MPI_Wtime();
	My_Allgather(send_buffer, number_send, MPI_INT, recv_buffer, number_send, MPI_INT, MPI_COMM_WORLD);
	time_elapsed += MPI_Wtime();

	if (!rank) {
		printf("Time elapsed of my implemented all gather: %f\n", time_elapsed);
	}

	send_buffer = malloc(sizeof(int) * number_send);
	recv_buffer = malloc(sizeof(int) * number_send * size);

	for (int i = 0 ; i < number_send; ++i)
		send_buffer[i] = i + number_send * rank;

    MPI_Barrier(MPI_COMM_WORLD);
	time_elapsed = -MPI_Wtime();
	MPI_Allgather(send_buffer, number_send, MPI_INT, recv_buffer, number_send, MPI_INT, MPI_COMM_WORLD);
	time_elapsed += MPI_Wtime();
	if (!rank) {
		printf("Time elapsed of original implemented all gather: %f\n", time_elapsed);
	}

	MPI_Finalize();
	return 0;
}
