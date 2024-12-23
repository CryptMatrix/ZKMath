#include <iostream>
#include "emp-tool/emp-tool.h"
#include "emp-zk/emp-zk.h"
using namespace emp;
using namespace std;

int port, party;
const int threads = 1;
// int index_sz = 5, val_sz = 32;

int LUT_bitlength = 4;
int index_sz = LUT_bitlength, val_sz = 32;

uint64_t comm(BoolIO<NetIO> *ios[threads]) {
	uint64_t c = 0;
	for(int i = 0; i < threads; ++i)
		c += ios[i]->counter;
	return c;
}
uint64_t comm2(BoolIO<NetIO> *ios[threads]) {
	uint64_t c = 0;
	for(int i = 0; i < threads; ++i)
		c += ios[i]->io->counter;
	return c;
}

void LUT_Exp_Gen(ROZKRAM<BoolIO<NetIO>> *ram, int LUT_bitlength){
	// TODO: one-time step-generate LUT and it's MAC in F2K
	int LUT_size = (1 << LUT_bitlength);
    vector<Integer> value;
	for(int i = 0; i < LUT_size; ++i)
		value.push_back(Integer(val_sz, exp(i), ALICE)); 
	ram->init(value);
}

void OneLUT_Exp(ROZKRAM<BoolIO<NetIO>> *ram, Integer *index, Integer *res, int n){
	// TODO: look up and check (n times)
	for (int i = 0; i < n; i++){
		res[i] = ram->read(index[i]);  // 这个index其实是只有Prover知道的  这一点与ram中不同！！？？
	}
}

void LUT_NegExp(IntFp *index, int size){
	// TODO: 
	// 1) mult -1
	// 2) A2B
	// 3) digit decomposition
	// 4) LUT generation
	// 5) OneLUT_Exp
	// 6) Mult
	// 7) Truncation
}

void test(BoolIO<NetIO> *ios[threads], int party) {
	setup_zk_bool<BoolIO<NetIO>>(ios, threads, party);
	vector<Integer> data;
	int test_n = (1<<index_sz);
	for(int i = 0; i < test_n; ++i)
		data.push_back(Integer(val_sz, 2*i, ALICE));  // 第i个位置赋值为2i
	ROZKRAM<BoolIO<NetIO>> *ram = new ROZKRAM<BoolIO<NetIO>>(party, index_sz, val_sz);
	ram->init(data);
	for(int i = 0; i < test_n; ++i) {
		Integer res = ram->read(Integer(index_sz, i, PUBLIC));   // 读第i个位置，i是public的！！
		Bit eq = (res == Integer(val_sz, i*2, ALICE));           // 第i个位置对应值为2i
		if(!eq.reveal<bool>(PUBLIC)) { 
			cout <<i<<"something is wrong!!\n";
		}
	}
	ram->check();   // 上面ram->read()中当一个batch存满时会调用check()，这里为何还要单独调用一遍？  因为可能出现最后一次 batch_size没放满的情况
	delete ram;
	finalize_zk_bool<BoolIO<NetIO>>();
	cout <<"done\n";
}

int main(int argc, char** argv) {
	parse_party_and_port(argv, &party, &port);
	BoolIO<NetIO>* ios[threads];
	for(int i = 0; i < threads; ++i)
		ios[i] = new BoolIO<NetIO>(new NetIO(party == ALICE?nullptr:"127.0.0.1",port), party==ALICE);

	if (argc > 3)
		index_sz = atoi(argv[3]);

	test(ios, party);

	for(int i = 0; i < threads; ++i) {
		delete ios[i]->io;
		delete ios[i];
	}
	return 0;
}
