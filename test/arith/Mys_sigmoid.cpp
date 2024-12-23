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

void test_sigmoid(BoolIO<NetIO> *ios[threads], int party, IntFp *x, IntFp *y, int size)
{
	Integer zero(ZK_INT_LEN, 0, PUBLIC);
  	Integer one(ZK_INT_LEN, 1, PUBLIC);
	Integer smallest_neg(ZK_INT_LEN, (uint64_t)((PR - 1) / 2) + 1, PUBLIC); 

	// Step 1: compute b
	vector<Integer> x_integer(size);
	IntFp *b_fp = new IntFp[size];
	arith2bool<BoolIO<NetIO>>(x_integer.data(), x, size); 
	vector<Integer> b_integer(size);
  	for (size_t i = 0; i < size; ++i){
		b_integer[i] = one.select(x_integer[i].geq(smallest_neg), zero);
	}
	bool2arith<BoolIO<NetIO>>(b_fp, b_integer.data(), size); 

	// step 2: compute x_bar
	IntFp *x_bar_fp = new IntFp[size];
	for (int i = 0; i < size; i++){
		x_bar_fp[i] = (b_fp[i] * 2 + (PR - 1)) * x[i];
	}

	// step 3: compute exp
	IntFp *z_fp = new IntFp[size];
	vector<Integer> z_integer(size);
	vector<Integer> x_bar_integer(size);
	arith2bool<BoolIO<NetIO>>(x_bar_integer.data(), x_bar_fp, size); 
	Float ir_res(0., PUBLIC);
	for (auto i = 0; i < size; i++)
	{
		ir_res = Int62ToFloat(x_bar_integer[i], ZK_F); 
		ir_res = ir_res.exp();				  
		z_integer[i] = FloatToInt62(ir_res, ZK_F); 
	}
	bool2arith<BoolIO<NetIO>>(z_fp, z_integer.data(), size); 

	// step 4: compute division
	IntFp *zaddone_fp = new IntFp[size];
	for (int i = 0; i < size; i++){
		zaddone_fp[i] = z_fp[i] + 1;
	}
	IntFp *d1_fp = new IntFp[size];
	vector<Integer> zaddone_integer(size);
	vector<Integer> d1_integer(size);
	Float one_float(1.0, PUBLIC);
	arith2bool<BoolIO<NetIO>>(zaddone_integer.data(), zaddone_fp, size); 
	for (auto i = 0; i < size; i++)
	{
		ir_res = Int62ToFloat(zaddone_integer[i], ZK_F); 
		ir_res = (one_float / ir_res);				  
		d1_integer[i] = FloatToInt62(ir_res, ZK_F); 
	}
	bool2arith<BoolIO<NetIO>>(d1_fp, d1_integer.data(), size); 

	// step 5: compute d2 and then truncation
	IntFp *d2_fp = new IntFp[size];
	for (int i = 0; i < size; i++){
		d2_fp[i] = z_fp[i] * d1_fp[i];
	}
	vector<Integer> d2_integer(size);
	arith2bool<BoolIO<NetIO>>(d2_integer.data(), d2_fp, size);
	for (size_t i = 0; i < size; i++){
		Float fl = std::move(Int62ToFloat(d2_integer[i], 2 * ZK_F));
		d2_integer[i] = FloatToInt62(fl, ZK_F);
	}

	bool2arith<BoolIO<NetIO>>(d2_fp, d2_integer.data(), size);

	// step 6: compute y
	for (size_t i = 0; i < size; i++){
		y[i] = b_fp[i] * d1_fp[i] + (b_fp[i].negate() + 1) * d2_fp[i];
	}
}

int main(int argc, char **argv)
{
	parse_party_and_port(argv, &party, &port);
	BoolIO<NetIO> *ios[threads];
	for (int i = 0; i < threads; ++i)
		ios[i] = new BoolIO<NetIO>(new NetIO(party == ALICE ? nullptr : "127.0.0.1", port + i), party == ALICE);

	std::cout << std::endl
			  << "------------ Mystique - sigmoid zero-knowledge proof test ------------" << std::endl
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
			// witness[i] = randomness[i] % PR;
			witness[i] = 1;
		}
		x[i] = IntFp(witness[i], ALICE);
	}

	cout << "start test" << endl;

	uint64_t com1 = comm(ios);
	auto start = clock_start();

	test_sigmoid(ios, party, x, y, size);

	finalize_zk_bool<BoolIO<NetIO>>();
	finalize_zk_arith<BoolIO<NetIO>>();
	
	double time1 = time_from(start);
	cout << "time - Mystique - sigmoid: " << time1 / 1000000 << " s\t " << party << endl;
	uint64_t com11 = comm(ios) - com1;
	std::cout << "communication - Mystique - sigmoid (KB): " << com11 / 1024.0 << std::endl;

	cout << "finish test" << endl;

	for (int i = 0; i < threads; ++i)
	{
		delete ios[i]->io;
		delete ios[i];
	}
	return 0;
}
