#include "emp-tool/emp-tool.h"
#include "emp-zk/emp-zk.h"
#include <iostream>

#define ZK_F 12

using namespace emp;
using namespace std;

int port, party;
const int threads = 1;

int size = 10000;

uint64_t comm(BoolIO<NetIO> *ios[threads])
{
	uint64_t c = 0;
	for (int i = 0; i < threads; ++i)
		c += ios[i]->counter;
	return c;
}

void test_Trunc(BoolIO<NetIO> *ios[threads], int party, IntFp *in, int size)
{
	vector<Integer> bc(size);
	arith2bool<BoolIO<NetIO>>(bc.data(), (IntFp *)in, size);
	for (size_t i = 0; i < bc.size(); i++){
		Float fl = std::move(Int62ToFloat(bc[i], 2 * ZK_F));
		bc[i] = FloatToInt62(fl, ZK_F);
	}
	bool2arith<BoolIO<NetIO>>(in, bc.data(), bc.size());
}

int main(int argc, char **argv)
{
	parse_party_and_port(argv, &party, &port);
	BoolIO<NetIO> *ios[threads];
	for (int i = 0; i < threads; ++i)
		ios[i] = new BoolIO<NetIO>(new NetIO(party == ALICE ? nullptr : "127.0.0.1", port + i), party == ALICE);

	std::cout << std::endl
			  << "------------ Truncation zero-knowledge proof test ------------" << std::endl
			  << std::endl;

	setup_zk_bool<BoolIO<NetIO>>(ios, threads, party);
	setup_zk_arith<BoolIO<NetIO>>(ios, threads, party, true);

	sync_zk_bool<BoolIO<NetIO>>();

	cout << "start data generation" << endl;

	srand(time(NULL));

	uint64_t *a = new uint64_t[size];
	for (int i = 0; i < size; ++i)
	{
		// cout << i << endl;
		a[i] = rand() % PR;
	}

	IntFp *a_zks = new IntFp[size];
	batch_feed(a_zks, a, size);

	cout << "start test" << endl;

	uint64_t com1 = comm(ios);
	auto start = clock_start();

	test_Trunc(ios, party, a_zks, size);

	finalize_zk_bool<BoolIO<NetIO>>();
	finalize_zk_arith<BoolIO<NetIO>>();

	double time2 = time_from(start);
	cout << "time: " << time2 / 1000 << " ms\t " << party << endl;
	uint64_t com2 = comm(ios) - com1;
	std::cout << "communication (KB): " << com2 / 1024.0 << std::endl;

	cout << "finish test" << endl;

	for (int i = 0; i < threads; ++i)
	{
		delete ios[i]->io;
		delete ios[i];
	}
	return 0;
}
