#include "iostream"
#include "stdlib.h"
#include "mpi.h"
#include "math.h"
using namespace std;
#define INF 10000

int main63(int argc, char** argv) {

	MPI_Init(&argc, &argv);
	int size;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	const int N = 4;
	int arr[N * N] = {0, 3, INF, 5,
					 2, 0, INF, 4,
					 INF, 1, 0, INF,
					 INF, INF, 2, 0};
 
	//int m = 4;
	int rowsPerProc = N / size;
	int rem = N % size;
	int** arr2D;
	int root;


	arr2D = new int* [N];


	for (int i = 0; i < N; ++i) {
		arr2D[i] = &arr[(i * N)];
	}

	int* recBuff = new int[N * (rowsPerProc + rem)];
	MPI_Status status;
	
	int* kRow = new int[N];
 	//broadcast to each proc their respective rows

	if (rank == 0) {
		for (int i = 1; i < size; ++i) {
			if (i != size - 1) {
				MPI_Send(&arr[i * N * rowsPerProc], N * rowsPerProc, MPI_INT, i, 0, MPI_COMM_WORLD);
			}
			else {
				MPI_Send(&arr[i * N * rowsPerProc],  N * (rowsPerProc + rem), MPI_INT, i, 0, MPI_COMM_WORLD);
			}
		}
		for (int i = 0; i < N * rowsPerProc; ++i) {
			//cout << arr[i] << endl;
			recBuff[i] = arr[i];
		}
	}
	else {
		
		if (rank == size - 1) {
			MPI_Recv(recBuff, N * (rowsPerProc + rem), MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		}
		else {
			MPI_Recv(recBuff, N * rowsPerProc, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		}
	
	}
 
	        /*if (rank == 0) {
	          	int k = 13;
	          	cout << "Proc#" << k / size << " is the owner of the " << k << "th" << " row" << endl;
	        }*/
	/*at this point, all processes have the rows in their receive buffer*/

	
	//if (rank == size - 1)
		//  for (int i = 0; i < N ; ++i) {
		//	cout << kRow[i] << endl;
		//   }


	/*if (rank == size - 1) {
		int i = 17;
		int tmpDiff = i - (rowsPerProc * rank);
		for (int j = tmpDiff * N; j < (tmpDiff * N) + N; ++j) {
			cout << recBuff[j] << endl;
		}
	}*/
	double start, end;
	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();

	int tmpDiff;
	for (int i = 0; i < N; ++i) {
 
	 
		root = int(i / rowsPerProc);
		if (root >= size)
			root = size - 1;
	 
		if (root == rank) {
			if (root == size - 1) {
				int rowStart = root * rowsPerProc;
				tmpDiff = i - rowStart;
				int z = 0;
				for (int j = tmpDiff * N; j < (tmpDiff * N) + N; ++j) {
					kRow[z++] = recBuff[j];
				}
			}
			else {
				int rowStart = root * rowsPerProc;
				tmpDiff = i - rowStart;
				int z = 0;
				for (int j = tmpDiff * N; j < (tmpDiff * N) + N; ++j) {
					kRow[z++] = recBuff[j];
				}
			}
		}
		MPI_Bcast(kRow, N, MPI_INT, root, MPI_COMM_WORLD);

		/*if (i == 1 && rank == 3) {
			for (int j = 0; j < N; ++j) {
				cout << kRow[j] << endl;
			}
		}*/
		int lim;
		if (rank == size - 1) {
			lim = (rowsPerProc + rem);
		}
		else {
			lim = rowsPerProc;
		}
		//if (rank == 0) {
		//	if (i == 0) {
		//		for (int j = 0; j < N; ++j) {
		//			cout << kRow[j] << "\t";
		//		}
		//	}
		//}
		for (int k = 0; k < lim; ++k) {
			for (int j = 0; j < N; ++j) {
				int min;
				if (recBuff[(k * N) + j] < (recBuff[(k * N) + i] + kRow[j])) {
					min = recBuff[(k * N) + j];
				}
				else 
					min = recBuff[(k * N) + i] + kRow[j];

					recBuff[(k * N) + j] = min;
					/*if (rank == 0) {
						cout << "index : " << kRow[j] << "  where k = " << i << " and j = " << j << endl;
					}*/
			}
		}

		/*for (int j = 0; j < N; j++) {
			int min;
			
			if (recBuff[j] < (recBuff[i] + kRow[j]))
				min = recBuff[j];
			else
				min = recBuff[i] + kRow[j];

			recBuff[j] = min;
		}*/
	 }
	MPI_Barrier(MPI_COMM_WORLD);
	end = MPI_Wtime();

	for (int i = 0; i < size; ++i) {
		if (rank == i) {
			if (rank != size - 1) {
				for (int j = 0; j < rowsPerProc; ++j) {
					//cout << "Row#" << (j + i);
					for (int k = 0; k < N; ++k) {
						cout << recBuff[(j * N) + k] << "\t";
					}
					cout << "\n";
				}
			}
			else {
				for (int j = 0; j < rowsPerProc + rem; ++j) {
					for (int k = 0; k < N; ++k) {
						cout << recBuff[(j * N) + k] << "\t";
					}
					cout << "\n";
				}
			}
		}
	}

	delete[] arr2D;
	delete[] recBuff;
	delete[] kRow;
	if (rank == 0) { /* use time on master node */
		cout << "Runtime = " << end - start << "\n";
	}
	MPI_Finalize();
	return 0;
}