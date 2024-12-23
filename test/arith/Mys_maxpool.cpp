#include "emp-tool/emp-tool.h"
#include "emp-zk/emp-zk.h"
#include <iostream>

#define ZK_F 12
#define ZK_INT_LEN 62

using namespace emp;
using namespace std;

int port, party;
const int threads = 1;

int rows = 100000;
int cols = 4;

uint64_t comm(BoolIO<NetIO> *ios[threads])
{
	uint64_t c = 0;
	for (int i = 0; i < threads; ++i)
		c += ios[i]->counter;
	return c;
}

void test_maxpool(BoolIO<NetIO> *ios[threads], int party, IntFp *x, IntFp *y, int rows, int cols)
{
	int dim = rows * cols;

	// Step 1: A2F
	vector<Integer> x_bool(dim);
	arith2bool<BoolIO<NetIO>>(x_bool.data(), x, dim); // 该function内部有check阶段
	vector<Float> x_float(dim);
	for (int i = 0; i < dim; i++){
		x_float[i] = Int62ToFloat(x_bool[i], ZK_F);
	}

	// step 2: max
	Float *y_float = new Float[rows];
	Integer *y_bool = new Integer[rows];
	for (int i = 0; i < rows; i++){
		Float input_max = x_float[i * cols + 0];
		for (size_t j = 1; j < cols; j++){
			input_max = If(input_max.less_than(x_float[i * cols + j]), x_float[i * cols + j], input_max);
		}
		y_float[i] = input_max;
		y_bool[i] = FloatToInt62(y_float[i], ZK_F);
	}
  	
	bool2arith<BoolIO<NetIO>>(y, y_bool, rows); // 该function内部有check阶段
}

int main(int argc, char **argv)
{
	parse_party_and_port(argv, &party, &port);
	BoolIO<NetIO> *ios[threads];
	for (int i = 0; i < threads; ++i)
		ios[i] = new BoolIO<NetIO>(new NetIO(party == ALICE ? nullptr : "127.0.0.1", port + i), party == ALICE);

	std::cout << std::endl
			  << "------------ Mystique - maxpool zero-knowledge proof test ------------" << std::endl
			  << std::endl;

	setup_zk_bool<BoolIO<NetIO>>(ios, threads, party);
	setup_zk_arith<BoolIO<NetIO>>(ios, threads, party, true);

	sync_zk_bool<BoolIO<NetIO>>();

	cout << "start data generation" << endl;

	srand(time(NULL));

	uint64_t *witness = new uint64_t[rows * cols];
	memset(witness, 0, rows * cols * sizeof(uint64_t));

	IntFp *x = new IntFp[rows * cols];
	IntFp *y = new IntFp[rows];
	// for (int i = 0; i < rows * cols; i++){
	// 	if (party == ALICE){
	// 		witness[i] = rand() % PR;
	// 	} 
	// 	x[i] = IntFp(witness[i], ALICE);
	// }
	__uint128_t *randomness = new __uint128_t[rows * cols]; 
	PRG prg(fix_key);
	prg.random_block((block*)randomness, rows * cols);
    for (int i = 0; i < rows * cols; i++){
		if (party == ALICE){
			// witness[i] = randomness[i] % PR;
			witness[i] = 1;
		}
		x[i] = IntFp(witness[i], ALICE);
	}

	cout << "start test" << endl;

	uint64_t com1 = comm(ios);
	auto start = clock_start();

	test_maxpool(ios, party, x, y, rows, cols);

	finalize_zk_bool<BoolIO<NetIO>>();
	finalize_zk_arith<BoolIO<NetIO>>();

	double time1 = time_from(start);
	cout << "time - maxpool: " << time1 / 1000000 << " s\t " << party << endl;
	uint64_t com11 = comm(ios) - com1;
	std::cout << "communication - maxpool (KB): " << com11 / 1024.0 << std::endl;

	cout << "finish test" << endl;

	for (int i = 0; i < threads; ++i)
	{
		delete ios[i]->io;
		delete ios[i];
	}
	return 0;
}
