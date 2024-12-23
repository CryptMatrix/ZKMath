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

void test_ExpExact(BoolIO<NetIO> *ios[threads], int party, IntFp *a_zks, IntFp *output, int size)
{
	// Step 1:FP to Integer
	vector<Integer> a_Int(size);
	// vector<IntFp> res_Int(size);
	// sync_zk_bool<BoolIO<NetIO>>();

	arith2bool<BoolIO<NetIO>>(a_Int.data(), a_zks, size); // 该function内部有check阶段

	// double time1 = time_from(start);
	// cout << "time - A2B: " << time1 / 1000 << " ms\t " << party << endl;
	// uint64_t com11 = comm(ios) - com1;
	// std::cout << "communication - A2B (KB): " << com11 / 1024.0 << std::endl;

	// Step 2:calc exp
	Float ir_res(0., PUBLIC);
	for (auto i = 0; i < size; i++)
	{
		ir_res = Int62ToFloat(a_Int[i], ZK_F); // 该过程所有的操作符都已被rewrite，最终都会表示成circuit的形式，check在AND gate协议中给出
		ir_res = ir_res.exp();				   // 有一个全局变量统计AND门次数，一旦大于batchcheck size, 在调用AND gate协议时，将先check
		a_Int[i] = FloatToInt62(ir_res, ZK_F); // 同样表示成circuit的形式，check在AND gate协议中给出
	}

	// double time2 = time_from(start) - time1;
	// cout << "time - Circuit eveluation: " << time2 / 1000 << " ms\t " << party << endl;
	// uint64_t com2 = comm(ios) - com11;
	// std::cout << "communication - Circuit eveluation (KB): " << com2 / 1024.0 << std::endl;

	// Step 3: Bool to Arithmetic
	bool2arith<BoolIO<NetIO>>(output, a_Int.data(), size); // 该function内部有check阶段

	// double time3 = time_from(start) - time1 - time2;
	// cout << "time - B2A: " << time3 / 1000 << " ms\t " << party << endl;
	// uint64_t com3 = comm(ios) - com2;
	// std::cout << "communication - B2A (KB): " << com3 / 1024.0 << std::endl;
}

int main(int argc, char **argv)
{
	parse_party_and_port(argv, &party, &port);
	BoolIO<NetIO> *ios[threads];
	for (int i = 0; i < threads; ++i)
		ios[i] = new BoolIO<NetIO>(new NetIO(party == ALICE ? nullptr : "127.0.0.1", port + i), party == ALICE);

	std::cout << std::endl
			  << "------------ Mystique_Exp_exact zero-knowledge proof test ------------" << std::endl
			  << std::endl;

	setup_zk_bool<BoolIO<NetIO>>(ios, threads, party);
	setup_zk_arith<BoolIO<NetIO>>(ios, threads, party, true);

	sync_zk_bool<BoolIO<NetIO>>();

	cout << "start data generation" << endl;

	srand(time(NULL));

	uint64_t constant = (PR + 1)/2;
	uint64_t *a = new uint64_t[size];
	IntFp *a_zks = new IntFp[size];
	__uint128_t *randomness = new __uint128_t[size]; 
	PRG prg(fix_key);
	prg.random_block((block*)randomness, size);
    for (int i = 0; i < size; i++){
		if (party == ALICE){
			a[i] = randomness[i] % constant;
		}
		a_zks[i] = IntFp(a[i], ALICE);
	}
	IntFp *output = new IntFp[size];

	cout << "start test" << endl;

	uint64_t com1 = comm(ios);
	auto start = clock_start();

	test_ExpExact(ios, party, a_zks, output, size);

	finalize_zk_bool<BoolIO<NetIO>>();
	finalize_zk_arith<BoolIO<NetIO>>();

	double time1 = time_from(start);
	cout << "time - Mystique_exp: " << time1 / 1000000 << " s\t " << party << endl;
	uint64_t com11 = comm(ios) - com1;
	std::cout << "communication - Mystique_exp (KB): " << com11 / 1024.0 << std::endl;

	/****************************/
	/**** verify correctness ****/
	/****************************/
	uint64_t total_err_fixed = 0;
    uint64_t max_ULP_err_fixed = 0;
	if (party == ALICE){
		uint64_t last_digit_bitlen = EXP_N - (EXP_DIGIT_LEN) * (EXP_LUT_NUM - 1);
		for (int i = 0; i < size; i++){
			// 计算每个分块的结果
			uint64_t *value_field = new uint64_t[EXP_LUT_NUM];
			for (int j = 0; j < EXP_LUT_NUM - 1; j++){
				uint64_t digit = (a[i] >> (j * EXP_DIGIT_LEN)) & ((1ULL << EXP_DIGIT_LEN) - 1);
				double digit_real = double(digit) / double(1ULL << SCALE);
            	double value_real = exp(-1 * (digit_real * (1ULL << EXP_DIGIT_LEN * j)));
				value_field[j] = value_real * (1ULL << SCALE);
			}
			uint64_t digit = (a[i] >> ((EXP_LUT_NUM - 1) * EXP_DIGIT_LEN)) & ((1ULL << last_digit_bitlen) - 1);
			double digit_real = double(digit) / double(1ULL << SCALE);
            double value_real = exp(-1 * (digit_real * (1ULL << EXP_DIGIT_LEN * (EXP_LUT_NUM - 1))));
			value_field[EXP_LUT_NUM - 1] = value_real * (1ULL << SCALE);
			// 相乘
			uint64_t ori_y = value_field[0];
			for (int j = 1; j < EXP_LUT_NUM; j++){
				ori_y = ori_y * value_field[j] >> SCALE;
			}
			uint64_t nm_y = (uint64_t)HIGH64(output[i].value);
			if (ori_y != nm_y){
				cout << "exp fault !!!" << "real: " << ori_y << ", protocol: " << nm_y << endl;
			}

			/****************************/
			/****** verify ULPerror *****/
			/****************************/
			double witness_real = Field2Real(a[i], SCALE);
			double exp = std::exp(-1 * witness_real);
			uint64_t exp_field = Real2Field(exp, SCALE);
			uint64_t err_fixed = computeULPErr(ori_y, exp_field);
      		if (err_fixed > 1)
      		{
        		cout << "ULP Error Fixed: " << ori_y << "," << exp_field << ","
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
