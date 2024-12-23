#include "emp-tool/emp-tool.h"
#include "emp-zk/emp-zk.h"
#include <iostream>

#define ZK_F 16
#define LUT_size 8
#define Bit_width 32

using namespace emp;
using namespace std;

int port, party;
const int threads = 1;

int size = 4000;

uint64_t comm(BoolIO<NetIO> *ios[threads])
{
    uint64_t c = 0;
    for (int i = 0; i < threads; ++i)
        c += ios[i]->counter;
    return c;
}

// void LUT_Exp_Gen(ROZKRAM<BoolIO<NetIO>> *ram)
// { // 可以传入一个ram类型的vector分别对应不同的LUT
//     // TODO: one-time step-generate LUT and it's MAC in F2K
//     int digit_size = LUT_size;
//     int num_digits = ceil(double(Bit_width) / digit_size);
//     int last_digit_size = Bit_width - (num_digits - 1) * digit_size;

//     int item_num = (1 << LUT_size);

//     vector<Integer> value;

//     for (int j = 0; j < num_digits; j++)
//     {
//         for (int i = 0; i < item_num; ++i)
//             value.push_back(Integer(Bit_width, exp(i), ALICE));
//     }
//     ram->init(value);
// }


static uint64_t lookup_neg_exp(uint64_t val_in, uint32_t s_in)
{
    if (s_in < 0)
    {
        s_in *= -1;
        val_in *= (1ULL << (s_in));
        s_in = 0;
    }
    uint64_t res_val = exp(-1.0 * (val_in / double(1LL << s_in))) * (1LL << ZK_F);
    return res_val;
}

void digitDecomp(int party, IntFp *input, IntFp *input_digits, int size, int num_digits, int digit_bitlen) {
    LUTRangeIntFp *LUTdigitDcomp[num_digits];
    for (int i = 0; i < num_digits; i++){
        LUTdigitDcomp[i] = new LUTRangeIntFp(party);
        LUTdigitDcomp[i]->LUTRangeinit(1ULL<<digit_bitlen);
    }
    for (auto i = 0; i < size; i++)
    {
        for (auto j = 0; j < num_digits; j++)
        {
            input_digits[j + i * num_digits] = IntFp((input[i].value >> j * digit_bitlen) & (1ULL << digit_bitlen - 1), ALICE);
            LUTdigitDcomp[j]->LUTRangeread(input_digits[j + i * num_digits]); 
        }
    }
}

void test_negExp(BoolIO<NetIO> *ios[threads], int party, IntFp *input, IntFp *output, int size)
{
    int digit_bitlen = LUT_size;
    int num_digits = ceil(double(Bit_width) / digit_bitlen);
    // int last_digit_size = Bit_width - (num_digits - 1) * digit_size;

    // convert to positive
    for (auto i = 0; i < size; i++)
        input[i] = input[i] * (PR - 1);

    LUTIntFp *LUTnegExp[num_digits];
    for (int i = 0; i < num_digits; i++){
        LUTnegExp[i] = new LUTIntFp(party);
        vector<uint64_t> data;
        for (int j = 0; j < (1ULL<<digit_bitlen); j++){
            uint64_t result = lookup_neg_exp(j, ZK_F - i * digit_bitlen);
            data.push_back(result);
        }
        LUTnegExp[i]->LUTinit(data);
    }

    uint64_t com1 = comm(ios);
    auto start = clock_start();


    IntFp *input_digits = new IntFp[num_digits * size];

    // store all decompositions
    digitDecomp(party, input, input_digits, size, num_digits, digit_bitlen);

    for (auto i = 0; i < size; i++)
    {
        IntFp *output_digits = new IntFp[num_digits];
        for (auto j = 0; j < num_digits; j++)
        {
            output_digits[j] = IntFp(lookup_neg_exp(j, ZK_F - i * digit_bitlen), ALICE);
            LUTnegExp[j]->LUTread(input_digits[j + i * num_digits], output_digits[j]); 
        }

        for (int j = 1; j < num_digits; j *= 2)
        {
            for (int k = 0; k < num_digits and k + j < num_digits; k += 2 * j)
            {
                output_digits[k] = output_digits[k + j] * output_digits[k];
            }
        }
        output[i] = output_digits[0]; // TODO truncation, only one-time?
    }

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
    for (int i = 0; i < size; ++i)
    {
        // cout << i << endl;
        a[i] = rand() % PR;
    }

    IntFp *a_zks = new IntFp[size];
    IntFp *output = new IntFp[size];
    batch_feed(a_zks, a, size);

    cout << "start test" << endl;

    test_negExp(ios, party, a_zks, output, size);

    cout << "finish test" << endl;

    finalize_zk_bool<BoolIO<NetIO>>();
    finalize_zk_arith<BoolIO<NetIO>>();

    for (int i = 0; i < threads; ++i)
    {
        delete ios[i]->io;
        delete ios[i];
    }
    return 0;
}
