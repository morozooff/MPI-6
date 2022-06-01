#include <complex>
#include <iostream>
#include <mpi.h>
#include <stdio.h>
#include <vector>
#include <Windows.h>

using namespace std;

const double PI = 3.1415;

//fast Fourie transform for multyplying
void fastFourie(vector<complex<double>>& coef, bool invert) {
	//array of coef
	int n = (int)coef.size();

	//separating to two polinoms 
	for (int i = 1, j = 0; i < n; i++) {
		int bit = n / 2;
		for (; j >= bit; bit /= 2)
			j -= bit;
		j += bit;
		//swapping
		if (i < j)
			swap(coef[i], coef[j]);
	}

	for (int len = 2; len <= n; len *= 2) {
		double angle = 2 * PI / len * (invert ? -1 : 1);
		complex<double> wlen(cos(angle), sin(angle));
		for (int i = 0; i < n; i += len) {
			complex<double> w(1);
			for (int j = 0; j < len / 2; ++j) {
				complex<double> u = coef[i + j], v = coef[i + j + len / 2] * w;
				coef[i + j] = u + v;
				coef[i + j + len / 2] = u - v;
				w *= wlen;
			}
		}
	}
	if (invert)
		for (int i = 0; i < n; ++i)
			coef[i] /= n;
}

//function for multiplying two polinoms
void multiply(vector<int> a, vector<int> b, vector<int>& c) {
	//{begin, begin+1,.....,end=begin+n}
	vector<complex<double>> fourierA = vector<complex<double>>(a.begin(), a.end());
	vector<complex<double>> fourierB = vector<complex<double>>(b.begin(), b.end());

	int n = 1;

	while (n < max(a.size(), b.size())) {
		n *= 2;
	}
	n *= 2;

	while (fourierA.size() != n) {
		fourierA.emplace(fourierA.begin(), complex<double>(0, 0));
	}
	while (fourierB.size() != n) {
		fourierB.emplace(fourierB.begin(), complex<double>(0, 0));
	}

	fastFourie(fourierA, false);
	fastFourie(fourierB, false);

	for (int i = 0; i < n; i++) {
		fourierA[i] *= fourierB[i];
	}

	fastFourie(fourierA, true);

	c.resize(n, 0);

	for (int i = 0; i < n - 1; i++) {
		c[i + 1] = int(fourierA[i].real() + 0.5);
	}

	for (int i = c.size() - 1; i > 0; i--) {
		if (c[i] >= 10) {
			c[i - 1] += c[i] / 10;
			c[i] = c[i] % 10;
		}
	}

}

vector<int> int_to_vector(const int& src) {
	string s = to_string(src);
	vector<int> res(s.begin(), s.end());
	for (int& re : res) {
		re -= '0';
	}
	return res;
}

int vector_to_int(vector<int> src) {
	int res = 0;

	reverse(src.begin(), src.end());

	int i = 0;
	for (const int& elem : src) {
		res += elem * pow(10, i);
		i++;
	}

	return res;
}

template<typename T>
//delete null in the beggining
void trimmer(vector<T>& src) {
	vector<T> res;
	bool leadingZero = true;

	for (const T& element : src) {
		if (element != 0)
			leadingZero = false;

		if (leadingZero)
			continue;

		res.push_back(element);

	}
	src = res;
}

int main(int* argc, char** argv) {

	//initialization
	MPI_Init(argc, &argv);
	int ProcNum, ProcRank;
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Status status;

	//num of outgoing arc from each node
	int indexes[] = { 6, 7, 8, 9, 10, 11, 12};
	//sequential list of graph arcs
	int edges[] = { 1,2,3,4,5,6, 0, 0, 0, 0, 0, 0 };

	//create new subcommunicator for graphs topology
	//let ProcNum = nnode
	MPI_Comm graph_comm_world;
	MPI_Graph_create(MPI_COMM_WORLD, ProcNum, indexes, edges, 1, &graph_comm_world);

	// amount of procces, wich checkProc have outgoing arc 
	int neighbours_count;
	MPI_Graph_neighbors_count(graph_comm_world, ProcRank, &neighbours_count);

	// getting rangs of neighbors
	int* neighbours = (int*)malloc(neighbours_count * sizeof(int));
	MPI_Graph_neighbors(graph_comm_world, ProcRank, neighbours_count, neighbours);

	//initializate new datatype for multiplying long int numb
	MPI_Datatype MPI_LONGINT;
	int long_int_size = 64;
	//creating this type
	MPI_Type_contiguous(long_int_size, MPI_INT, &MPI_LONGINT);
	MPI_Type_commit(&MPI_LONGINT);
	srand(time(NULL) + ProcRank);

	//for
	if (neighbours_count == 1) {
		//randomize value
		int a = rand() % 1000;
		int b = rand() % 1000;
		//out random value
		printf("(%d * %d)", a, b);
		if (neighbours_count = 5) {
			printf("\n     *");
		}

		vector<int> A_vector = int_to_vector(a);
		vector<int> B_vector = int_to_vector(b);
		vector<int> C_vector;

		multiply(A_vector, B_vector, C_vector);
		MPI_Send(C_vector.data(), 1, MPI_LONGINT, neighbours[0], 0, graph_comm_world);

		int size = C_vector.size();
		// send array size
		MPI_Send(&size, 1, MPI_INT, neighbours[0], 1, graph_comm_world); 

	}
	else if (neighbours_count == 6)
	{
		int* res1 = new int[long_int_size];
		int* res2 = new int[long_int_size];
		int* res3 = new int[long_int_size];
		int* res4 = new int[long_int_size];
		int* res5 = new int[long_int_size];
		int* res6 = new int[long_int_size];

		//accept res
		MPI_Recv(res1, 1, MPI_LONGINT, neighbours[0], 0, graph_comm_world, &status);
		MPI_Recv(res2, 1, MPI_LONGINT, neighbours[1], 0, graph_comm_world, &status);
		MPI_Recv(res3, 1, MPI_LONGINT, neighbours[2], 0, graph_comm_world, &status);
		MPI_Recv(res4, 1, MPI_LONGINT, neighbours[3], 0, graph_comm_world, &status);
		MPI_Recv(res5, 1, MPI_LONGINT, neighbours[4], 0, graph_comm_world, &status);
		MPI_Recv(res6, 1, MPI_LONGINT, neighbours[5], 0, graph_comm_world, &status);

		//accept size of res
		int size1, size2, size3, size4, size5, size6;
		MPI_Recv(&size1, 1, MPI_INT, neighbours[0], 1, graph_comm_world, &status);
		MPI_Recv(&size2, 1, MPI_INT, neighbours[1], 1, graph_comm_world, &status);
		MPI_Recv(&size3, 1, MPI_INT, neighbours[2], 1, graph_comm_world, &status);
		MPI_Recv(&size4, 1, MPI_INT, neighbours[3], 1, graph_comm_world, &status);
		MPI_Recv(&size5, 1, MPI_INT, neighbours[4], 1, graph_comm_world, &status);
		MPI_Recv(&size6, 1, MPI_INT, neighbours[5], 1, graph_comm_world, &status);

		//(array, arrays size)
		vector<int> A_vector(res1, res1 + size1);
		vector<int> B_vector(res2, res2 + size2);
		vector<int> C_vector(res3, res3 + size3);
		vector<int> D_vector(res4, res4 + size4);
		vector<int> E_vector(res5, res5 + size5);
		vector<int> F_vector(res6, res6 + size6);


		// vector for resultig value's after the myltiplying
		vector<int> G_vector;
		vector<int> H_vector;
		vector<int> I_vector;
		vector<int> J_vector;
		vector<int> K_vector;
		

		multiply(A_vector, B_vector, G_vector);
		multiply(G_vector, C_vector, H_vector);
		multiply(H_vector, D_vector, I_vector);
		multiply(I_vector, E_vector, J_vector);
		multiply(J_vector, F_vector, K_vector);


		//trim nulls
		trimmer(K_vector);
		//out result vector
		printf("Result of multiplying = ");
		for (int i = 0; i < K_vector.size(); i++)
		{
			printf("%d", K_vector[i]);
		}
		printf("\n_________________________________________________________");
		cout << endl;
	}


	//finish
	MPI_Finalize();
	return 0;
}
