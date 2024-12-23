#include "emp-tool/emp-tool.h"
#include "emp-zk/emp-zk.h"
#include <iostream>

#define ZK_F 12

using namespace emp;
using namespace std;

int port, party;
const int threads = 1;

int rows = 100000;
int cols = 10;

uint64_t comm(BoolIO<NetIO> *ios[threads])
{
	uint64_t c = 0;
	for (int i = 0; i < threads; ++i)
		c += ios[i]->counter;
	return c;
}

void test_Softmax(BoolIO<NetIO> *ios[threads], int party, IntFp *a_zks, IntFp *output, int rows, int cols)
{
	int size = rows * cols;

	vector<Integer> a_int(size);
	arith2bool<BoolIO<NetIO>>(a_int.data(), a_zks, size);

	// softmax operation(row major)
	vector<Float> rows_ft_temp(cols, (0.0, PUBLIC));
	vector<Integer> res_int(cols);
	vector<IntFp> res_zks(cols);

	for (int i = 0; i < rows; i++)
	{
		// calc exp(a[ix]) & sum(exp(a[ix]...))
		Float sum(0., PUBLIC);
		for (int j = 0; j < cols; j++)
		{
			rows_ft_temp[j] = Int62ToFloat(a_int[i * cols + j], ZK_F);
		}

		// get max
		Float input_max = rows_ft_temp[0];
		for (size_t j = 1; j < cols; j++)
		{
			input_max = If(input_max.less_than(rows_ft_temp[j]), rows_ft_temp[j], input_max);
		}

		for (int j = 0; j < cols; j++)
		{
			rows_ft_temp[j] = (rows_ft_temp[j] - input_max).exp();
			sum = sum + rows_ft_temp[j];
		}

		// calc output[x0] = exp(a[x0]) / sum
		for (int j = 0; j < cols; j++)
		{
			rows_ft_temp[j] = rows_ft_temp[j] / sum;

			res_int[j] = FloatToInt62(rows_ft_temp[j], ZK_F);
		}

		// Bool to Arithmetic
		bool2arith<BoolIO<NetIO>>(res_zks.data(), res_int.data(), cols);

		sync_zk_bool<BoolIO<NetIO>>();

		for (int j = 0; j < cols; j++)
		{
			output[i * cols + j] = res_zks[j];
		}
	}
}

int main(int argc, char **argv)
{
	parse_party_and_port(argv, &party, &port);
	BoolIO<NetIO> *ios[threads];
	for (int i = 0; i < threads; ++i)
		ios[i] = new BoolIO<NetIO>(new NetIO(party == ALICE ? nullptr : "127.0.0.1", port + i), party == ALICE);

	std::cout << std::endl
			  << "------------ Mystique softmax zero-knowledge proof test ------------" << std::endl
			  << std::endl;

	setup_zk_bool<BoolIO<NetIO>>(ios, threads, party);
	setup_zk_arith<BoolIO<NetIO>>(ios, threads, party, true);

	sync_zk_bool<BoolIO<NetIO>>();

	cout << "start data generation" << endl;

	srand(time(NULL));

	int size = rows * cols;

	uint64_t *a = new uint64_t[size];
	for (int i = 0; i < size; ++i)
	{
		a[i] = rand() % PR;
	}

	IntFp *a_zks = new IntFp[size];
	IntFp *output = new IntFp[size];
	batch_feed(a_zks, a, size);

	cout << "start test" << endl;

	uint64_t com1 = comm(ios);
	auto start = clock_start();

	test_Softmax(ios, party, a_zks, output, rows, cols);

	finalize_zk_bool<BoolIO<NetIO>>();
	finalize_zk_arith<BoolIO<NetIO>>();

	double time2 = time_from(start);
	cout << "time - Mystique softmax: " << time2 / 1000000 << " s\t " << party << endl;
	uint64_t com2 = comm(ios) - com1;
	std::cout << "communication (KB) - Mystique softmax: " << com2 / 1024.0 << std::endl;

	cout << "finish test" << endl;

	for (int i = 0; i < threads; ++i)
	{
		delete ios[i]->io;
		delete ios[i];
	}
	return 0;
}
