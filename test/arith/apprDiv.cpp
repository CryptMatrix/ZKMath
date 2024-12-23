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


// // A1 \in (1/4, 1)  b
static int64_t lookup_div_A1(int64_t index, int m)
{
  int32_t k = 1ULL << m;
  // p = 1.y1...ym
  double p = 1 + (double(index) / double(k));
  double A1 = 1.0 / (p * (p + 1.0 / double(k)));
  // TODO 提高精度？
  int32_t scale = m + 3;
  int64_t val = (A1 * (1LL << scale));
  // val = 1.xxxxx scale=m+3
  return val;
}

// A0 \in (1/2, 1)  a 
static int64_t lookup_div_A0(int64_t index, int m)
{
  int32_t k = 1ULL << m;
  double p = 1 + (double(index) / double(k));
  double z = (p * (p + (1.0 / double(k))));
  double A0 = ((1.0 / double(k * 2)) + sqrt(z)) / z;
  int32_t scale = 2 * m + 2;
  int64_t val = (A0 * (1LL << scale));
  //val = 0.1xxxxx       scale=2 * m + 2
  return val;
}

// input format = 1.xxxxx
void div_plaintext(int32_t size, int32_t bwA, int32_t bwB, int32_t bwC, int32_t sA, int32_t sB, int32_t sC, vector<int64_t> &A, vector<int64_t> &B, vector<int64_t> &C, bool compute_msnzb)
{
  int32_t m, iters;
  // m-bit 用于lookup table求近似。scale=12，m=5.
  // TODO set large m??
  m = (sC <= 18 ? ceil((sC - 2) / 2.0) : ceil((ceil(sC / 2.0) - 2) / 2.0));
  iters = (sC <= 18 ? 0 : 1);
  assert(m <= 8);
  assert(sB >= m);

  // in ZK, bwB表示输入数据的上界

  int64_t m_mask = (1LL << m) - 1;
  uint64_t s_min_m_mask = (1ULL << (sB - m)) - 1;
  for (int i = 0; i < size; i++)
  {
    int32_t tmp_sB;
    int32_t msnzb, s_adjust;
    int64_t adjust, tmpB;
    if (compute_msnzb)
    {
      msnzb = floor(log2(int64_t(B[i])));
      tmp_sB = bwB - 1;
      s_adjust = bwB - 1 - sB; // Note，这里是-sB
      adjust = (1LL << (bwB - 1 - msnzb));
      tmpB = (adjust * int64_t(B[i])); // TODO 这里需要调整回来
      s_min_m_mask = (1ULL << (tmp_sB - m)) - 1;
    }
    else
    {
      assert((B[i] >= (1LL << sB)) && (B[i] < (2LL << sB)));
      tmpB = int64_t(B[i]);
      tmp_sB = sB;
    }
    int64_t B_m = (tmpB >> (tmp_sB - m)) & m_mask;
    // scale_A0 = 2m+2, \in (1/2, 1) 开区间，0.1xxx
    // 原因是 scale = 2m+2
    int64_t A0 = lookup_div_A0(B_m, m);
    // scale_A1 = m+3, \in (1/4, 1) 开区间 0.x1xxx
    int64_t A1 = lookup_div_A1(B_m, m);
    int64_t Q = tmpB & s_min_m_mask;
    // reciprocal approximation: Y=1/B, scale = sC
    // scale A1=m+3 Q=tmp_sB A0=2m+2
    // A1*Q scale = m+3+tmp_sB
    int64_t Y = ((A0 << (tmp_sB + 1 - m)) - A1 * Q) >> (tmp_sB + m + 3 - sC);
    int64_t a0 = (int64_t(A[i]) * Y) >> sA;
    if (compute_msnzb)
    {
      // 本质上就是 a0 * 2^(sB-msnzb)
      a0 = (a0 * adjust) >> s_adjust;
    }
    // if (iters > 0)
    // {
    //   int64_t e = (tmpB * Y) >> sB;
    //   // e0 = 1 - e
    //   int64_t e0 = (1LL << sC) - e;
    //   int64_t a_prev = a0, e_prev = e0, a_curr, e_curr;
    //   for (int j = 0; j < iters - 1; j++)
    //   {
    //     e_curr = (e_prev * e_prev) >> sC;
    //     a_curr = (a_prev * ((1LL << sC) + e_prev)) >> sC;
    //     a_prev = a_curr;
    //     e_prev = e_curr;
    //   }
    //   a0 = (a_prev * ((1LL << sC) + e_prev)) >> sC;
    // }
    C[i] = (int64_t)((int32_t)(a0));
  }
}


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
