#include "emp-tool/emp-tool.h"
#include "emp-zk/emp-zk.h"
#include <iostream>

#define ZK_F 16
#define ZK_INT_LEN 62

using namespace emp;
using namespace std;

int port, party;
const int threads = 1;

int size = 1;

uint64_t comm(BoolIO<NetIO> *ios[threads])
{
	uint64_t c = 0;
	for (int i = 0; i < threads; ++i)
		c += ios[i]->counter;
	return c;
}


void test_Div(BoolIO<NetIO> *ios[threads], int party, IntFp *a, IntFp *b, IntFp *output, int size) 
{

	uint64_t com1 = comm(ios);
	auto start = clock_start();

	// Step 1: calc Invert b result
  	vector<IntFp> b_invert_zks(size);
	vector<Integer> a_Int(size);
	
	arith2bool<BoolIO<NetIO>>(a_Int.data(), a, size);

  	// Step 1-2: calc invert
  	Float ir_res(0., PUBLIC);
  	Float one(1.0, PUBLIC);

  	for (auto i = 0; i < size; i++) 
	{
    	// integer to float, and calc invert
    	ir_res = Int62ToFloat(a_Int[i], ZK_F);
    	ir_res = (one / ir_res);
    	a_Int[i] = FloatToInt62(ir_res, ZK_F);
  	}

  	// Step 1-3: Bool to Arithmetic
  	bool2arith<BoolIO<NetIO>>((IntFp*)b_invert_zks.data(), a_Int.data(), size);


	vector<IntFp> c_zks(size);
    // Step 2-1: a * b => c
	for (auto i = 0; i < size; i++) 
	{
      	c_zks[i] = a[i] * b_invert_zks[i];
  	}

  	// Step 2-2: result(c_zks) to bool
  	sync_zk_bool<BoolIO<NetIO>>();
  	vector<Integer> c_res(size, Integer(ZK_INT_LEN, 0, PUBLIC));
  	arith2bool<BoolIO<NetIO>>(c_res.data(), (IntFp*)c_zks.data(), size);

  	// Step 2-3: truncate
  	// Float ir_res(0., PUBLIC);
  	for (auto i = 0; i < size; i++) 
	{
    	ir_res = Int62ToFloat(c_res[i], ZK_F);
    	c_res[i] = FloatToInt62(ir_res, 0);
  	}

  	// Step 2-4: Bool to Arithmetic
  	bool2arith<BoolIO<NetIO>>(output, c_res.data(), size);

	double time2 = time_from(start);

	cout << "time: " << time2 / 1000 << " ms\t " << party << endl;

	uint64_t com2 = comm(ios) - com1;
	std::cout << "communication (KB): " << com2 / 1024.0 << std::endl;
}

int main(int argc, char **argv)
{
	parse_party_and_port(argv, &party, &port);
	BoolIO<NetIO> *ios[threads];
	for (int i = 0; i < threads; ++i)
		ios[i] = new BoolIO<NetIO>(new NetIO(party == ALICE ? nullptr : "127.0.0.1", port + i), party == ALICE);

	std::cout << std::endl
			  << "------------ Exp_exact zero-knowledge proof test ------------" << std::endl
			  << std::endl;

	setup_zk_bool<BoolIO<NetIO>>(ios, threads, party);
	setup_zk_arith<BoolIO<NetIO>>(ios, threads, party, true);

	sync_zk_bool<BoolIO<NetIO>>();

	cout << "start data generation" << endl;

	srand(time(NULL));



	uint64_t *a = new uint64_t[size];
	uint64_t *b = new uint64_t[size];
	for (int i = 0; i < size; ++i)
	{
		// cout << i << endl;
		a[i] = rand() % PR;
		b[i] = rand() % PR;
	}

	IntFp *a_zks = new IntFp[size];
	IntFp *b_zks = new IntFp[size];
	IntFp *output = new IntFp[size];
	batch_feed(a_zks, a, size);
	batch_feed(b_zks, b, size);

	cout << "start test" << endl;

	test_Div(ios, party, a_zks, b_zks, output, size);

	cout << "finish test" << endl;

	sync_zk_bool<BoolIO<NetIO>>();

	finalize_zk_bool<BoolIO<NetIO>>();
	finalize_zk_arith<BoolIO<NetIO>>();

	for (int i = 0; i < threads; ++i)
	{
		delete ios[i]->io;
		delete ios[i];
	}
	return 0;
}
