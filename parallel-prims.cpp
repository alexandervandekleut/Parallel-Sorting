#include <vector>
#include <mpi.h>
#include <stdio.h>
#include <string>
#include <math.h>
#include <limits.h>

const int INF = INT_MAX;
typedef std::vector<int> list;
typedef std::vector<list> matrix;

void print(matrix M){
	for (list A : M){
		for (int i=0; i<A.size(); i++){
			printf("%i ", A[i]);
		} printf("\n");
	}
}

void print(list A){
	for (int i=0; i<A.size(); i++){
		printf("%i ", A[i]);
	} printf("\n");
}

void print(int p, list A){
	printf("Process %i: ", p);
	for (int i=0; i<A.size(); i++){
		printf("%i ", A[i]);
	} printf("\n");
}



/* Sequential prim's */
matrix prims(matrix &distances){
	int n = distances.size();

	matrix tree(n); // initialize empty nxn matrix. we set tree[i][j] = 1 if (i,j) \in tree at the end
	for (int i=0; i<n; i++){
		tree[i] = list(n);
		for (int j=0; j<n; j++){
			tree[i][j] = 0;
		}
	}

	list unexplored(n-1); // {1, 2, 3, ..., n-1}
	for (int i=0; i<n-1; i++){
		unexplored[i] = i+1;
	}

	list explored = {0}; // seed tree at 0 vertex

	while (unexplored.size() > 0){ // while there are still more unexplored vertices
		int min_dist = INF; // determine the shortest distance to a vertex from explore to unexplored
		int f; // the from vertex
		int t; // the to vertex
		int idx; // the index of the to vertex, used to remove it from the unexplored list and add it to the explored list
		for (int i=0; i<explored.size(); i++){
			int from = explored[i];
			for (int j=0; j<unexplored.size(); j++){
				int to = unexplored[j];
				if (distances[from][to] < min_dist){
					min_dist = distances[from][to];
					f = from;
					t = to;
					idx = j;
				}
			}
		}
		tree[f][t] = 1;
		tree[t][f] = 1;
  		explored.push_back(t); // add the vertex we found with the minimum distance to the explored list
  		unexplored.erase(unexplored.begin() + idx); // erase the vertex we found with the minimum distance from the unexplored list
  	}

  	return tree;
  }

  int main(int argc, char *argv[]){ 

  	MPI_Init(NULL, NULL); 

	int communicationSize = 0; // number of processes (should be a power of two)
  	int processRank = 0; // rank of current process

  	MPI_Comm_size(MPI_COMM_WORLD, &communicationSize); // set communicationSize with the appropriate number of processes
	MPI_Comm_rank(MPI_COMM_WORLD, &processRank); // set processRank with the appropriate process number

  	matrix distances =  {{INF,	5,		2,		9,		3,		6,		INF,	1}, // symmetric 8x8 matrix corresponding to distances between vertices
					  	{5,		INF,	3,		3,		INF,	9,		4,		5},
					  	{2,		3,		INF,	6,		2,		9,		5,		3},
					  	{9,		3,		6,		INF,	INF,	8,		8,		7},
					  	{3,		INF,	2,		INF,	INF,	3,		1,		4},
					  	{6,		9,		9,		INF,	3,		INF,	2,		5},
					  	{INF,	4,		5,		8,		1,		2,		INF,	8},
					  	{1,		5,		3,		7,		4,		5,		8,		INF}};

	int n = distances.size();

	matrix tree(n); // initialize empty nxn matrix. we set tree[i][j] = 1 if (i,j) \in tree at the end
	for (int i=0; i<n; i++){
		tree[i] = list(n);
		for (int j=0; j<n; j++){
			tree[i][j] = 0;
		}
	}

	list unexplored; // list of unexplored vertices (cyclic assignment) {1 + process rank, 1+comm size + process rank, 1+2*comm size + process rank, ...}
	for (int i=1; i<n; i++){
		if ((i%communicationSize) == processRank){
			unexplored.push_back(i);
		}
	}

	list explored = {0}; // seed tree at 0 vertex

	int flag = 1; // used to determine if any process still has unexplored vertices

	while (flag){ // while there are still more unexplored vertices

		int min_dist = INF; // used to determine the shortest distance of an edge from an explored vertex to an unexplored vertex
		int f; // the from vertex
		int t; // the to vertex
		int idx; // the index of the to vertex, used to remove it from the unexplored list and add it to the explored list
		int from; // temp from during loops
		int to; // temp to during loops


		for (int i=0; i<explored.size(); i++){ // for each vertex in explored
			from = explored[i];
			for (int j=0; j<unexplored.size(); j++){ // for each vertex in unexplored
				to = unexplored[j];
				if (distances[from][to] < min_dist){ // if this beats the current shortest edge connected connected to unconnected
					min_dist = distances[from][to]; // update everything
					f = from;
					t = to;
					idx = j;
				}
			}
		}

		int u = (unexplored.size() > 0) ? 1 : 0; // update u to determine if there is still some elements left to explore

		MPI_Send(&u, 1, MPI_INT, 0, 0, MPI_COMM_WORLD); // for simplicity, everyone sends and recieves each iteration even if doing no work
		MPI_Send(&f, 1, MPI_INT, 0, 1, MPI_COMM_WORLD); // the u flag tells process 0 whether or not the data incoming is relevant
		MPI_Send(&t, 1, MPI_INT, 0, 2, MPI_COMM_WORLD); // every process sends u, f, and t whether or not its being used to sends and recieves can match up

		if (processRank == 0){ // process 0 determines the best (f, t) pair
			int min_dist = INF; // assume worst possible distance
			for (int p=0; p<communicationSize; p++){ // for each process
				MPI_Recv(&u, 1, MPI_INT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // store their u flag, from and to
				MPI_Recv(&from, 1, MPI_INT, p, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(&to, 1, MPI_INT, p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				if (u){ // ONLY if u is set (the process actually sent usedful info)
					int d = distances[from][to]; // get the distance for these (from, to) pair
					if (d < min_dist){ // if it beats the best distance
						min_dist = d; // consider the from and to to be the new (f,t) pair
						f = from;
						t = to;
					}
				}
			}
		}

		MPI_Bcast(&f, 1, MPI_INT, 0, MPI_COMM_WORLD); // send out (f,t) to everyone so they all update their own tree and explored/unexplored list
		MPI_Bcast(&t, 1, MPI_INT, 0, MPI_COMM_WORLD);

		tree[f][t] = 1; // set the tree
		tree[t][f] = 1;
  		explored.push_back(t); // add the vertex we found with the minimum distance to the explored list
  		
  		if (u){ // only if we sent useful info could we possibly remove a vertex
  			int idx = -1; // start with impossible index. if when we are done its still -1 then we didnt find it 
  			for (int i=0; i<unexplored.size(); i++){ // get the index of the element to be removed from unexplored
  				if (unexplored[i] == t){
  					idx = i;
  				}
  			}
  			if (idx != -1){
  				unexplored.erase(unexplored.begin() + idx); // erase the vertex we found with the minimum distance from the unexplored list
  			}
  		}

  		MPI_Allreduce(&u, &flag, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD); // update everyones flags. if eveyones unexplored is empty then its 0

  	}

  	if (processRank == 0){
  		print(tree);
  		int total = 0;
  		for (int i=0; i<n; i++){
  			for (int j=0; j<i; j++){
  				if (tree[i][j] == 1){
  					total += distances[i][j];
  				}
  			}
  		}
  		printf("total: %i\n",total);
  	}

	MPI_Finalize();
}











