#include "emp-tool/emp-tool.h"
#include "emp-zk/emp-zk.h"
#include <iostream>

#define ZK_F 12

using namespace emp;
using namespace std;

int port, party;
const int threads = 1;

int rows = 100000;
int cols = 16;

uint64_t comm(BoolIO<NetIO> *ios[threads])
{
	uint64_t c = 0;
	for (int i = 0; i < threads; ++i)
		c += ios[i]->counter;
	return c;
}

void test_layernorm(BoolIO<NetIO> *ios[threads], int party, IntFp *input, IntFp *output, IntFp *gamma, IntFp *beta, int rows, int cols)
{
	int size = rows * cols;

	int64_t cols_fix = floor(cols * (1ULL << ZK_F));
    uint64_t cols_field = cols_fix < 0 ? PR + cols_fix : cols_fix; 

	// step 1: 计算mu
	IntFp *sum = new IntFp[rows];
	IntFp *mu = new IntFp[rows];
	for (int i = 0; i < rows; i++){
		sum[i] = input[i * cols + 0];
		for (int j = 1; j < cols; j++){
			sum[i] = sum[i] + input[i * cols + j];
		}
		mu[i] = sum[i] * cols_field;
	}
	vector<Integer> mu_bool(rows);
	arith2bool<BoolIO<NetIO>>(mu_bool.data(), mu, rows);
	for (size_t i = 0; i < mu_bool.size(); i++){
		Float fl = std::move(Int62ToFloat(mu_bool[i], 2 * ZK_F));
		mu_bool[i] = FloatToInt62(fl, ZK_F);
	}
	bool2arith<BoolIO<NetIO>>(mu, mu_bool.data(), rows);
	cout << "end step 1" << endl;

	// step 2: 计算sigma
	IntFp *sigma = new IntFp[rows];
	IntFp *tmp = new IntFp[rows * cols];
	for (int i = 0; i < rows; i++){
		tmp[i * cols + 0] = input[i * cols + 0] + mu[i].negate();
		sum[i] = tmp[i * cols + 0] * tmp[i * cols + 0];
		for (int j = 1; j < cols; j++){
			tmp[i * cols + j] = input[i * cols + j] + mu[i].negate();
			sum[i] = sum[i] + (tmp[i * cols + j] * tmp[i * cols + j]);
		}
		sigma[i] = sum[i] * cols_field;
	}
	vector<Integer> sigma_bool(rows);
	arith2bool<BoolIO<NetIO>>(sigma_bool.data(), sigma, rows);
	for (size_t i = 0; i < sigma_bool.size(); i++){
		Float fl = std::move(Int62ToFloat(sigma_bool[i], 3 * ZK_F));
		sigma_bool[i] = FloatToInt62(fl, ZK_F);
	}
	// bool2arith<BoolIO<NetIO>>(sigma, sigma_bool.data(), rows);
	cout << "end step 2" << endl;

	// step 3: 计算rsqrt
	IntFp *t = new IntFp[rows];
	// arith2bool<BoolIO<NetIO>>(sigma_bool.data(), sigma, rows);
	Float ir_res(0., PUBLIC);
	Float one(1.0, PUBLIC);
	for (int i = 0; i < rows; i++){
		ir_res = Int62ToFloat(sigma_bool[i], ZK_F); 
		ir_res = ir_res.sqrt();				   
		ir_res = (one / ir_res);
		sigma_bool[i] = FloatToInt62(ir_res, ZK_F); 
	}
	bool2arith<BoolIO<NetIO>>(t, sigma_bool.data(), rows);
	cout << "end step 3" << endl;

	// step 4: 计算z
	IntFp *z = new IntFp[rows * cols];
	for (int i = 0; i < rows; i++){
		for (int j = 0; j < cols; j++){
			z[i * cols + j] = t[i] * tmp[i * cols + j];
		}
	}
	vector<Integer> z_bool(rows * cols);
	arith2bool<BoolIO<NetIO>>(z_bool.data(), z, rows * cols);
	for (size_t i = 0; i < z_bool.size(); i++){
		Float fl = std::move(Int62ToFloat(z_bool[i], 2 * ZK_F));
		z_bool[i] = FloatToInt62(fl, ZK_F);
	}
	bool2arith<BoolIO<NetIO>>(z, z_bool.data(), rows * cols);
	cout << "end step 4" << endl;

	// step 5: 计算output
	for (int i = 0; i < rows; i++){
		for (int j = 0; j < cols; j++){
			output[i * cols + j] = z[i * cols + j] * gamma[i] + beta[i];
		}
	}
	vector<Integer> output_bool(rows * cols);
	arith2bool<BoolIO<NetIO>>(output_bool.data(), output, rows * cols);
	for (size_t i = 0; i < output_bool.size(); i++){
		Float fl = std::move(Int62ToFloat(output_bool[i], 2 * ZK_F));
		output_bool[i] = FloatToInt62(fl, ZK_F);
	}
	bool2arith<BoolIO<NetIO>>(output, output_bool.data(), rows * cols);
	cout << "end step 5" << endl;
}

int main(int argc, char **argv)
{
	parse_party_and_port(argv, &party, &port);
	BoolIO<NetIO> *ios[threads];
	for (int i = 0; i < threads; ++i)
		ios[i] = new BoolIO<NetIO>(new NetIO(party == ALICE ? nullptr : "127.0.0.1", port + i), party == ALICE);

	std::cout << std::endl
			  << "------------ Mystique layernorm zero-knowledge proof test ------------" << std::endl
			  << std::endl;

	setup_zk_bool<BoolIO<NetIO>>(ios, threads, party);
	setup_zk_arith<BoolIO<NetIO>>(ios, threads, party, true);

	sync_zk_bool<BoolIO<NetIO>>();

	cout << "start data generation" << endl;

	srand(time(NULL));

	int size = rows * cols;

	uint64_t *a = new uint64_t[size];
	for (int i = 0; i < size; ++i){
		a[i] = 10;
	}
	IntFp *a_zks = new IntFp[size];
	batch_feed(a_zks, a, size);

	IntFp *gamma = new IntFp[rows];
	IntFp *beta = new IntFp[rows];
	uint64_t *b = new uint64_t[rows];
	for (int i = 0; i < rows; ++i){
		b[i] = 1;
	}
	batch_feed(gamma, b, rows);
	for (int i = 0; i < rows; ++i){
		b[i] = 1;
	}
	batch_feed(beta, b, rows);

	IntFp *output = new IntFp[size];

	cout << "start test" << endl;

	uint64_t com1 = comm(ios);
	auto start = clock_start();

	test_layernorm(ios, party, a_zks, output, gamma, beta, rows, cols);

	finalize_zk_bool<BoolIO<NetIO>>();
	finalize_zk_arith<BoolIO<NetIO>>();

	double time2 = time_from(start);
	cout << "time - Mystique layernorm: " << time2 / 1000000 << " s\t " << party << endl;
	uint64_t com2 = comm(ios) - com1;
	std::cout << "communication (KB) - Mystique layernorm: " << com2 / 1024.0 << std::endl;

	cout << "finish test" << endl;

	for (int i = 0; i < threads; ++i)
	{
		delete ios[i]->io;
		delete ios[i];
	}
	return 0;
}
