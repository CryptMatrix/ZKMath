#include "emp-tool/emp-tool.h"
#include "emp-zk/emp-zk.h"
#include <iostream>

#define ZK_F 12

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

void test_Div(BoolIO<NetIO> *ios[threads], int party, IntFp *a_zks, IntFp *output, int size)
{
	// Step 1: calc Invert b result
  	vector<IntFp> b_invert_zks(size);
	vector<Integer> a_Int(size);

	arith2bool<BoolIO<NetIO>>(a_Int.data(), a_zks, size); // 该function内部有check阶段

	// double time1 = time_from(start);
	// cout << "time - A2B: " << time1 / 1000 << " ms\t " << party << endl;

	Float ir_res(0., PUBLIC);
  	Float one(1.0, PUBLIC);

  	for (auto i = 0; i < size; i++) 
	{
    	// integer to float, and calc invert
    	ir_res = Int62ToFloat(a_Int[i], ZK_F);
    	ir_res = (one / ir_res);
    	a_Int[i] = FloatToInt62(ir_res, ZK_F);
  	}

	// double time2 = time_from(start) - time1;
	// cout << "time - Circuit eveluation: " << time2 / 1000 << " ms\t " << party << endl;

	// Step 3: Bool to Arithmetic
	bool2arith<BoolIO<NetIO>>(output, a_Int.data(), size); // 该function内部有check阶段

	// double time3 = time_from(start) - time1 - time2;
	// cout << "time - B2A: " << time3 / 1000 << " ms\t " << party << endl;
}

int main(int argc, char **argv)
{
	parse_party_and_port(argv, &party, &port);
	BoolIO<NetIO> *ios[threads];
	for (int i = 0; i < threads; ++i)
		ios[i] = new BoolIO<NetIO>(new NetIO(party == ALICE ? nullptr : "127.0.0.1", port + i), party == ALICE);

	std::cout << std::endl
			  << "------------ Mystique - Div zero-knowledge proof test ------------" << std::endl
			  << std::endl;

	setup_zk_bool<BoolIO<NetIO>>(ios, threads, party);
	setup_zk_arith<BoolIO<NetIO>>(ios, threads, party, true);

	sync_zk_bool<BoolIO<NetIO>>();

	cout << "start data generation" << endl;

	uint64_t *witness = new uint64_t[size];
	memset(witness, 0, size * sizeof(uint64_t));

	startComputation(party);

	uint64_t constant = (PR + 1)/2;
	IntFp *x = new IntFp[size];
	IntFp *y = new IntFp[size];
	// for (int i = 0; i < size; i++){
	// 	if (party == ALICE){
	// 		while (witness[i] == 0)
	// 		{
	// 			witness[i] = rand();
	// 		}
	// 		witness[i] = witness[i] % PR;
	// 	}
	// 	x[i] = IntFp(witness[i], ALICE);
	// }
	__uint128_t *randomness = new __uint128_t[size]; 
	PRG prg(fix_key);
	prg.random_block((block*)randomness, size);
    for (int i = 0; i < size; i++){
		if (party == ALICE){
			witness[i] = randomness[i] % constant;
			// if (witness[i] == 0){
			// 	witness[i] = witness[i] + 1;
			// }
			if (witness[i] < 1){
				witness[i] = witness[i] + 1;
			}
		}
		x[i] = IntFp(witness[i], ALICE);
	}

	cout << "start test" << endl;

	uint64_t com1 = comm(ios);
	auto start = clock_start();

	test_Div(ios, party, x, y, size);

	finalize_zk_bool<BoolIO<NetIO>>();
	finalize_zk_arith<BoolIO<NetIO>>();
	
	double time1 = time_from(start);
	cout << "time - Mystique - div: " << time1 / 1000000 << " s\t " << party << endl;
	uint64_t com11 = comm(ios) - com1;
	std::cout << "communication - Mystique - div (KB): " << com11 / 1024.0 << std::endl;

	uint64_t total_err_fixed = 0;
    uint64_t max_ULP_err_fixed = 0;
	if (party == ALICE){
		for (int i = 0; i < size; i++){
			uint64_t ori_y = (uint64_t)HIGH64(y[i].value);

			double witness_real = Field2Real(witness[i], ZK_F);
			double div = 1.0/witness_real;
			uint64_t div_field = Real2Field(div, ZK_F);
			uint64_t err_fixed = computeULPErr(ori_y, div_field);
      		if (err_fixed > 1)
      		{
        		cout << "ZKDiv ULP Error Fixed: " << ori_y << "," << div_field << ","
             		<< err_fixed << endl;
      		}
      		total_err_fixed += err_fixed;
      		max_ULP_err_fixed = std::max(max_ULP_err_fixed, err_fixed);
		}
		cout << "Average ULP error fixed: " << total_err_fixed / size << endl;
    	cout << "Total ULP error fixed: " << total_err_fixed << endl;
    	cout << "Max ULP error fixed: " << max_ULP_err_fixed << endl;
    	cout << "Number of tests fixed: " << size << endl;
	}

	cout << "finish test" << endl;

	for (int i = 0; i < threads; ++i)
	{
		delete ios[i]->io;
		delete ios[i];
	}
	return 0;
}
