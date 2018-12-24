#include <vector>
#include <mpi.h>
#include <stdio.h>
#include <string>
#include <math.h>

typedef std::vector<int> list;

void print(list A){
	for (int i=0; i<A.size(); i++){
		printf("%i ", A[i]);
	} printf("\n");
}


/* Sequential mergesort routine */
void merge(list &A, int p, int q, int r){
	list L(&A[p], &A[q]);
	list R(&A[q], &A[r]);

	int i=0;
	int j=0;

	for (int k=p; k < r; k++){
		if (i >= L.size()){
			A[k] = R[j];
			j++;
		} else if (j >= R.size()){
			A[k] = L[i];
			i++;
		} else {
			if (L[i] <= R[j]){
				A[k] = L[i];
				i++;
			} else {
				A[k] = R[j];
				j++;
			}
		}
	}
}

void mergesort(list &A, int p, int r){
	if ((r-p) > 1){
		int q = (p+r)/2; // [p, q)
		mergesort(A, p, q); // [q, r)
		mergesort(A, q, r);
		merge(A, p, q, r);
	}
}

/* Parallel mergesort routine */
list merge(list &L, list &R){

	list A(L.size() + R.size()); // to store merged lists

	int i=0; // index to loop through L
	int j=0; // index to loop through R

	for (int k=0; k < A.size(); k++){ // k loops through A
		if (i >= L.size()){ // if we are done with L just fill A with R
			A[k] = R[j];
			j++;
		} else if (j >= R.size()){ // if we are done with R just fill A with L
			A[k] = L[i];
			i++;
		} else { // otherwise choose the minimum from the front of L or R and add it to the end of A
			if (L[i] <= R[j]){
				A[k] = L[i];
				i++;
			} else {
				A[k] = R[j];
				j++;
			}
		}
	}

	return A;
}


int main(int argc, char *argv[]){ 

	MPI_Init(NULL, NULL); 

	int communicationSize = 0; // number of processes (should be a power of two)
  	int processRank = 0; // rank of current process

 	MPI_Comm_size(MPI_COMM_WORLD, &communicationSize); // set communicationSize with the appropriate number of processes
	MPI_Comm_rank(MPI_COMM_WORLD, &processRank); // set processRank with the appropriate process number

	list A = {22, 7, 13, 18, 2, 17, 1, 14, 20, 6, 10, 24, 15, 9, 21, 3, 16, 19, 23, 4, 11, 12, 5, 8}; // sample list to be sorted. 24 = 3*2^3 so we can use max n=8 processors to sort it.
	
	list local_A(A.size()/communicationSize); // for now we assume that this is an even division

	MPI_Scatter(A.data(), local_A.size(), MPI_INT, local_A.data(), local_A.size(), MPI_INT, 0, MPI_COMM_WORLD); // partitioning phase done all at once

	mergesort(local_A, 0, local_A.size()); // lists locally sorted

	/* Example of binary tree topology used to sort
	R = receive
	S = send

	P   1  2  3   
	0 | R1 R2 R4 
	1 | S0
	2 | R3 S0
	3 | S2
	4 | R5 R6 S0
	5 | S4
	6 | R7 S4
	7 | S6

	*/

	for (int d = 2; d <= communicationSize; d*=2){ // merging procedure. d/2 is the difference in ranks between communicating processes.
		if (processRank % d == 0){ // receive
			list local_B(local_A.size()); // temp buffer to hold data received from receiving partner
			MPI_Recv(local_B.data(), local_B.size(), MPI_INT, processRank+d/2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // receive into local B from sending partner
			local_A = merge(local_A,local_B); // merge local A and local B together and store result as new local A
		} else if ((processRank - d/2)%d == 0){ // send
			MPI_Send(local_A.data(), local_A.size(), MPI_INT, processRank-d/2, 0, MPI_COMM_WORLD); // send the local A to the receiving partner
		}
	}

	if (processRank == 0){ // process 0 is where the results eventually pool to using the above process.
		print(local_A);
	}

	MPI_Finalize();
}