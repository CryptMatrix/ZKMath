#include "emp-tool/emp-tool.h"
#include "emp-zk/emp-zk.h"
#include <iostream>

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

int main(int argc, char **argv)
{
	parse_party_and_port(argv, &party, &port);
	BoolIO<NetIO> *ios[threads];
	for (int i = 0; i < threads; ++i)
		ios[i] = new BoolIO<NetIO>(new NetIO(party == ALICE ? nullptr : "127.0.0.1", port + i), party == ALICE);

	std::cout << std::endl
			  << "------------ ZKMax test ------------" << std::endl
			  << std::endl;

	setup_zk_bool<BoolIO<NetIO>>(ios, threads, party);
	setup_zk_arith<BoolIO<NetIO>>(ios, threads, party, true);

	sync_zk_bool<BoolIO<NetIO>>();

	uint64_t *witness = new uint64_t[rows * cols];
	memset(witness, 0, rows * cols * sizeof(uint64_t));

	startComputation(party);

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

	uint64_t com = comm(ios);
	auto start = clock_start();

	ZKMax(party, x, y, rows, cols); 

	endComputation(party);

	double time = time_from(start);
	cout << "Performance without LUT construction" << endl;
	cout << "time - ZKMax (s): " << time / 1000000 << " s\t " << party << endl;
	uint64_t com1 = comm(ios) - com;
	std::cout << "communication - ZKMax (KB): " << com1 / 1024.0 << std::endl;

	/****************************/
	/**** verify correctness ****/
	/****************************/
	// if (party == ALICE){
	// 	int64_t witness_int = 0;
	// 	uint64_t output_real = 0;
	// 	uint64_t output_prot = 0;
	// 	for (int i = 0; i < rows; i++){
	// 		int64_t max_y = witness[i * cols] > (PR-1)/2 ? witness[i * cols] - PR : witness[i * cols];
	// 		for (int j = 1; j < cols; j++){
	// 			witness_int = witness[i * cols + j] > (PR-1)/2 ? witness[i * cols + j] - PR : witness[i * cols + j];
	// 			// if (witness[i * cols + j] > (PR-1)/2 && witness_int < 0){
	// 			// 	cout << "this is a negative value" << endl;
	// 			// }
	// 			if (witness_int > max_y){
	// 				max_y = witness_int;
	// 			}
	// 		}
	// 		output_real = max_y < 0 ? PR + max_y : max_y;
	// 		output_prot = (uint64_t)HIGH64(y[i].value);
	// 		if (output_real != output_prot){
	// 			cout << "fault !!!" << endl;
	// 		}
	// 	}
	// }

	cout << "finish test" << endl;

	finalize_zk_bool<BoolIO<NetIO>>();
	finalize_zk_arith<BoolIO<NetIO>>();

	for (int i = 0; i < threads; i++)
	{
		delete ios[i]->io;
		delete ios[i];
	}
	return 0;
}
