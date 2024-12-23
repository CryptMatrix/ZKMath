#include "emp-tool/emp-tool.h"
#include "emp-zk/emp-zk.h"
#include <iostream>

#define ZK_F 12

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

// only for positive inputs
// TODO add argument-scale
void test_Trunc_positive(BoolIO<NetIO> *ios[threads], int party, IntFp *in, IntFp *out, int size)
{
    vector<Integer> inB(size);

    uint64_t com1 = comm(ios);
    auto start = clock_start();

    arith2bool<BoolIO<NetIO>>(inB.data(), (IntFp *)in, size);

    for (size_t i = 0; i < inB.size(); i++)
    {
        inB[i] = inB[i] >> ZK_F; // logic shift
    }

    bool2arith<BoolIO<NetIO>>(out, inB.data(), inB.size());

    double time2 = time_from(start);

    cout << "time: " << time2 / 1000 << " ms\t " << party << endl;

    uint64_t com2 = comm(ios) - com1;
    std::cout << "communication (KB): " << com2 / 1024.0 << std::endl;

    int incorrect_cnt = 0;
    for (size_t i = 0; i < inB.size(); i++)
    {
        int64_t input = in[i].reveal();
        int64_t output = out[i].reveal();
        if (input >> ZK_F != output)
            incorrect_cnt++;
    }

    if (incorrect_cnt)
        std::cout << "incorrect boolean: " << incorrect_cnt << std::endl;
    std::cout << "end check boolean" << std::endl;
}

// general truncation
void test_Trunc(BoolIO<NetIO> *ios[threads], int party, IntFp *in, IntFp *out, int size)
{
    vector<Integer> inB(size);

    uint64_t com1 = comm(ios);
    auto start = clock_start();

    arith2bool<BoolIO<NetIO>>(inB.data(), (IntFp *)in, size);

    const Integer zero(62, 0, PUBLIC);
    const Integer prime(62, PR, PUBLIC);
    const Integer power2(62, 1ULL << (61 - ZK_F), PUBLIC);
    const Bit oneB(true, PUBLIC);

    for (size_t i = 0; i < inB.size(); i++)
    {
        // compute msb. Here approximate comparison is used.
        Bit msb = inB[i].bits[60];

        // convert to real value from Fp value
        Integer offset = If(msb, prime, zero);
        inB[i] = inB[i] - offset;

        // compute truncation
        inB[i].bits[61] ^= oneB;
        inB[i] = (inB[i] >> ZK_F) - power2 + offset;
    }

    bool2arith<BoolIO<NetIO>>(out, inB.data(), inB.size());

    double time2 = time_from(start);

    cout << "time: " << time2 / 1000 << " ms\t " << party << endl;

    uint64_t com2 = comm(ios) - com1;
    std::cout << "communication (KB): " << com2 / 1024.0 << std::endl;

    int incorrect_cnt = 0;
    for (size_t i = 0; i < inB.size(); i++)
    {
        int64_t input = in[i].reveal();
        int64_t output = out[i].reveal();
        if (input >> ZK_F != output)
        {
            incorrect_cnt++;
            // std::cout << input << std::endl;
        }
    }

    if (incorrect_cnt)
        std::cout << "incorrect boolean: " << incorrect_cnt << std::endl;
    std::cout << "end check boolean" << std::endl;
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

    size = 10000;

    srand(time(NULL));

    uint64_t *a = new uint64_t[size];
    for (int i = 0; i < size; ++i)
    {
        // cout << i << endl;
        a[i] = rand() % PR;
        // a[i] = rand() % (uint64_t)((PR + 1) / 2); // only positive values
    }

    IntFp *a_zks = new IntFp[size];
    batch_feed(a_zks, a, size);

    IntFp *b_zks = new IntFp[size];
    batch_feed(b_zks, a, size);

    cout << "start test" << endl;

    // test_Trunc_positive(ios, party, a_zks, b_zks, size);

    test_Trunc(ios, party, a_zks, b_zks, size);

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
