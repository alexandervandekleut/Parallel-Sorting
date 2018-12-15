#include <vector>
#include <mpi.h>
#include <stdio.h>
#include <string>

typedef std::vector<int> list;

list quicksort(list &A, int q, int r){
	if (q<r){
		int x = A[q];
		int s = q;
		for (int i=q+1; i<r; i++){
			if (A[i] <= x){
				s = s+1;
				int temp = A[s];
				A[s] = A[i];
				A[i] = temp;
			}
		}
		int temp =  A[q];
		A[q] = A[s];
		A[s] = temp;
		quicksort(A, q, s);
		quicksort(A, s+1, r);
	}
	return A;
};

void print(list A){
	for (int i=0; i<A.size(); i++){
		printf("%i ", A[i]);
	} printf("\n");
}

int main(int argc, char *argv[]){

	MPI_Init(NULL, NULL); 

	int communicationSize = 0; // number of processes
  	int processRank = 0; // rank of current process

 	MPI_Comm_size(MPI_COMM_WORLD, &communicationSize); // set communicationSize with the appropriate number of processes
	MPI_Comm_rank(MPI_COMM_WORLD, &processRank); // set processRank with the appropriate process number

	list A = {22, 7, 13, 18, 2, 17, 1, 14, 20, 6, 10, 24, 15, 9, 21, 3, 16, 19, 23, 4, 11, 12, 5, 8}; // sample list to be sorted



	/* Part 1: Partition the original list into sublists to make good choices for partitions */
	list B = list(A.size()/communicationSize); // local partition of list in each processor

	MPI_Scatter(A.data(), B.size(), MPI_INT, B.data(), B.size(), MPI_INT, 0, MPI_COMM_WORLD); // split array into chunks

	B = quicksort(B, 0, B.size()); // locally sort chunk

	int samples = (argc == 2) ? std::stoi(argv[1]) : communicationSize; // default to one sample per processor. rounded down to next lowest multiple of the number of processors.

	int samplesPerBucket = samples/communicationSize; // get samples per processor

	list s = list(samplesPerBucket);
	for (int i=s.size()-1; i>=0; i--){ // go backwards to sample this way
		s[i] = B[i*B.size()/s.size()];
	}

	list S = list(samples); // list of global splitters to choose from

	MPI_Gather(s.data(), s.size(), MPI_INT, S.data(), s.size(), MPI_INT, 0, MPI_COMM_WORLD); // collect local splitters into global splitters

	if (processRank == 0){
		quicksort(S, 0, S.size()); // sort list of global splitters

		list G = list(communicationSize-1);
		for (int i=0; i<G.size(); i++){ // go forwards to sample the final global splitters
			G[i] = S[i*S.size()/G.size()];
		}

		std::vector<list> partitions = std::vector<list>(G.size() + 1); // a list of lists that are split based on the list of final global splitters
		
		for (int i=0; i<A.size(); i++){ // for each element in the original list
			bool placed = false; // default to not placed
			for (int j=0; j<G.size(); j++){ // for each final global splitter
				if (!placed){ // so long as we havent placed this element yet
					if (A[i] < G[j]){ // if it is smaller than this splitter
					partitions[j].push_back(A[i]); // place it
					placed = true; // and move on to the next element
					}
				}
			}
			if (!placed){ // if it was never placed, then it must go in the highest bucket
				partitions[partitions.size()-1].push_back(A[i]);
			}
		}

		for (int i=0; i<communicationSize; i++){ // each process corresponds to a list in partitions
			list toSend = partitions[i]; // for each process/partition
			int sendSize = toSend.size(); // get the size

			MPI_Send(&sendSize, 1, MPI_INT, i, 0, MPI_COMM_WORLD); // first send over the size of the array we are sending
			MPI_Send(toSend.data(), sendSize, MPI_INT, i, 1, MPI_COMM_WORLD); // then send over the array, now that we know how much to recieve
		}
	}

	int recvSize; // how much we are going to recieve gets sent first
	MPI_Recv(&recvSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	list recvBuffer = list(recvSize); // then make a buffer of the appropriate size
	MPI_Recv(recvBuffer.data(), recvSize, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // fill it 

	recvBuffer = quicksort(recvBuffer, 0 ,recvBuffer.size()); // local sort on partition of array

	MPI_Send(&recvSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD); // send back size to expect of sorted partition of array
	MPI_Send(recvBuffer.data(), recvBuffer.size(), MPI_INT, 0, 1, MPI_COMM_WORLD); // send back sorted partition of array

	if (processRank == 0){
		list result = list(); // final sorted list to concatenate
		
		for (int i=0; i<communicationSize; i++){ // collect sorted sublists from each process

			MPI_Recv(&recvSize, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // reuse old buffer from above
			recvBuffer = list(recvSize); // then make a buffer of the appropriate size, reuse buffer declaration from above
			MPI_Recv(recvBuffer.data(), recvSize, MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // fill it 
			result.insert(result.end(), recvBuffer.begin(), recvBuffer.end()); // concatenate recieved array
		}
		print(result);
	}

	MPI_Finalize();
}