#include "emp-tool/emp-tool.h"
#include "emp-zk/emp-zk.h"
#include <iostream>

#define ZK_F 12
#define ZK_INT_LEN 62

using namespace emp;
using namespace std;

int port, party;
const int threads = 1;

int size = 100000;

uint64_t comm(BoolIO<NetIO> *ios[threads])
{
	uint64_t c = 0;
	for (int i = 0; i < threads; ++i)
		c += ios[i]->counter;
	return c;
}

void test_relu(BoolIO<NetIO> *ios[threads], int party, IntFp *a_zks, IntFp *output, int size)
{
	// Step 1: calc Invert b result
	vector<Integer> a_Int(size);

	arith2bool<BoolIO<NetIO>>(a_Int.data(), a_zks, size); // 该function内部有check阶段

	Integer zero(ZK_INT_LEN, 0, PUBLIC);
  	Integer one(ZK_INT_LEN, 1, PUBLIC);

  	vector<Integer> relu(size);
  	Integer smallest_neg(ZK_INT_LEN, (uint64_t)((PR - 1) / 2) + 1, PUBLIC); // ??? this is a positive?
  	for (size_t i = 0; i < size; ++i)
    	relu[i] = a_Int[i].select(a_Int[i].geq(smallest_neg), zero);

	// Step 3: Bool to Arithmetic
	bool2arith<BoolIO<NetIO>>(output, relu.data(), size); // 该function内部有check阶段
}

int main(int argc, char **argv)
{
	parse_party_and_port(argv, &party, &port);
	BoolIO<NetIO> *ios[threads];
	for (int i = 0; i < threads; ++i)
		ios[i] = new BoolIO<NetIO>(new NetIO(party == ALICE ? nullptr : "127.0.0.1", port + i), party == ALICE);

	std::cout << std::endl
			  << "------------ Mystique - relu - zero-knowledge proof test ------------" << std::endl
			  << std::endl;

	setup_zk_bool<BoolIO<NetIO>>(ios, threads, party);
	setup_zk_arith<BoolIO<NetIO>>(ios, threads, party, true);

	sync_zk_bool<BoolIO<NetIO>>();

	cout << "start data generation" << endl;

	uint64_t *witness = new uint64_t[size];
	memset(witness, 0, size * sizeof(uint64_t));
	IntFp *x = new IntFp[size];
	IntFp *y = new IntFp[size];
	__uint128_t *randomness = new __uint128_t[size]; 
	PRG prg(fix_key);
	prg.random_block((block*)randomness, size);
    for (int i = 0; i < size; i++){
		if (party == ALICE){
			witness[i] = randomness[i] % PR;
		}
		x[i] = IntFp(witness[i], ALICE);
	}

	cout << "start test" << endl;

	uint64_t com1 = comm(ios);
	auto start = clock_start();

	test_relu(ios, party, x, y, size);

	finalize_zk_bool<BoolIO<NetIO>>();
	finalize_zk_arith<BoolIO<NetIO>>();

	double time1 = time_from(start);
	cout << "time - Mystique - relu: " << time1 / 1000000 << " s\t " << party << endl;
	uint64_t com11 = comm(ios) - com1;
	std::cout << "communication - Mystique - relu (KB): " << com11 / 1024.0 << std::endl;

	cout << "finish test" << endl;

	for (int i = 0; i < threads; ++i)
	{
		delete ios[i]->io;
		delete ios[i];
	}
	return 0;
}
