#include "emp-zk/emp-zk.h"
#include <iostream>
#include "emp-tool/emp-tool.h"
#if defined(__linux__)
	#include <sys/time.h>
	#include <sys/resource.h>
#elif defined(__APPLE__)
	#include <unistd.h>
	#include <sys/resource.h>
	#include <mach/mach.h>
#endif

using namespace emp;
using namespace std;

int port, party;
const int threads = 1;

uint64_t dim1 = 256;
uint64_t dim2 = 768;
uint64_t dim3 = 768; 

uint64_t comm(BoolIO<NetIO> *ios[threads]) {
	uint64_t c = 0;
	for(int i = 0; i < threads; ++i)
		c += ios[i]->counter;
	return c;
}

void test_circuit_zk(BoolIO<NetIO> *ios[threads], int party, uint64_t dim1, uint64_t dim2, uint64_t dim3, uint64_t *ar, uint64_t *br, uint64_t *cr) 
{
    setup_zk_arith<BoolIO<NetIO>>(ios, threads, party);

    IntFp *mat_a = new IntFp[dim1 * dim2];
	IntFp *mat_b = new IntFp[dim2 * dim3];
	IntFp *mat_c = new IntFp[dim1 * dim3];
	
	for(int i = 0; i < dim1 * dim2; ++i) {
		mat_a[i] = IntFp((uint64_t)i, ALICE);
	}
    for(int i = 0; i < dim2 * dim3; ++i) {
		mat_b[i] = IntFp((uint64_t)i, ALICE);
	}
    for(int i = 0; i < dim1 * dim3; ++i) {
		mat_c[i] = IntFp((uint64_t)0, PUBLIC);
	}

    uint64_t com1 = comm(ios);
    auto start = clock_start();

    for(int i = 0; i < dim1; ++i) {
		for(int j = 0; j < dim2; ++j) {
			for(int k = 0; k < dim3; ++k) {
				uint64_t tmp = (ar[i*dim2+j] * br[j*dim3+k])%pr;
				cr[i*dim3+k] = (cr[i*dim3+k] + tmp)%pr;
                IntFp tmp1 = mat_a[i*dim2+j] * mat_b[j*dim3+k];
                mat_c[i*dim3+k] = mat_c[i*dim3+k] + tmp1;
			}
		}
	}

	batch_reveal_check(mat_c, cr, dim1 * dim3);
	double timeuse = time_from(start);
	
	cout << "time: " << timeuse / 1000 << " ms\t " << party << endl;
	
    uint64_t com2 = comm(ios) - com1;
	std::cout << "communication (KB): " << com2 / 1024.0 << std::endl;

    finalize_zk_arith<BoolIO<NetIO>>();

	delete[] ar;
	delete[] br;
	delete[] cr;
	delete[] mat_a;
	delete[] mat_b;
	delete[] mat_c;


#if defined(__linux__)
	struct rusage rusage;
	if (!getrusage(RUSAGE_SELF, &rusage))
		std::cout << "[Linux]Peak resident set size: " << (size_t)rusage.ru_maxrss << std::endl;
	else std::cout << "[Linux]Query RSS failed" << std::endl;
#elif defined(__APPLE__)
	struct mach_task_basic_info info;
	mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
	if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &count) == KERN_SUCCESS)
		std::cout << "[Mac]Peak resident set size: " << (size_t)info.resident_size_max << std::endl;
	else std::cout << "[Mac]Query RSS failed" << std::endl;
#endif
}

int main(int argc, char** argv) {
	parse_party_and_port(argv, &party, &port);
	BoolIO<NetIO>* ios[threads];
	for(int i = 0; i < threads; ++i)
		ios[i] = new BoolIO<NetIO>(new NetIO(party == ALICE?nullptr:"127.0.0.1",port+i), party==ALICE);

	std::cout << std::endl << "FC zero-knowledge proof test ------------" << std::endl << std::endl;;

	int num = 0;
	if(argc < 3) {
		std::cout << "usage: bin/arith/matrix_mul_arith PARTY PORT DIMENSION" << std::endl;
		return -1;
	} else if (argc == 3) {
		num = 10;
	} else {
		num = atoi(argv[3]);
	}


    uint64_t *ar, *br, *cr;
	ar = new uint64_t[dim1 * dim2];
	br = new uint64_t[dim2 * dim3];
	cr = new uint64_t[dim1 * dim3];

    for(int i = 0; i < dim1 * dim2; ++i) {
		ar[i] = i;
	}
    for(int i = 0; i < dim2 * dim3; ++i) {
		br[i] = i;
	}
    for(int i = 0; i < dim1 * dim3; ++i) {
		cr[i] = 0;
	}
	

	test_circuit_zk(ios, party, dim1, dim2, dim3, ar, br, cr);

	for(int i = 0; i < threads; ++i) {
		delete ios[i]->io;
		delete ios[i];
	}
	return 0;
}
