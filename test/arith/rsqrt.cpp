#include "emp-tool/emp-tool.h"
#include "emp-zk/emp-zk.h"
#include <iostream>

using namespace emp;
using namespace std;

#define MAX_THREADS 4
#define SQRT_LOOKUP_SCALE 2
#define KKOT_LIMIT 8

int port, party;
const int threads = 1;

int LUT_size = 5;

inline int64_t Saturate(int32_t inp) { return (int64_t)inp; }

static int64_t lookup_sqrt(int32_t index, int m, int exp_parity)    // index = e; m = g; exp_parity = B (in paper)
{
  int32_t k = 1ULL << m;
  double u = (1.0 + (double(index) / double(k))) * (1ULL << exp_parity);
  double Y = 1.0 / sqrt(u);
  int32_t scale = m + SQRT_LOOKUP_SCALE;
  int64_t val = (Y * (1ULL << scale));
  return val;
}

void cleartext_rsqrt(int64_t *A, int32_t I, int32_t J, int32_t shrA,
                    int32_t shrB, int32_t bwA, int32_t bwB, int64_t *B)
{
  int32_t sA = log2(shrA);
  int32_t sB = log2(shrB);
  int32_t m = ceil(sB / 2.0);   // m = g (in paper)
  int32_t iters = 1;

  assert(m <= KKOT_LIMIT);
  assert(sA >= m);

  int64_t m_mask = (1LL << m) - 1;
  int64_t error = 0;

  for (int i = 0; i < I * J; i++)
  {
// Normalizes A to shifted_A in [1, 2)
    assert(A[i] >= 0);

    int32_t msnzb = floor(log2(int64_t(A[i])));   // msnzb: 最高非零位的index (从0开始)

    int32_t adjust_amt = bwA - 2 - msnzb;    
    int64_t adjust = (1LL << adjust_amt);   // adjust = A (in paper) 

    int64_t sqrt_adjust_scale = floor((bwA - 1 - sA) / 2.0);
    int64_t sqrt_adjust = (1LL << int(floor((sA - msnzb + 1) / 2.0) + sqrt_adjust_scale));   // sqrt_adjust = C (in paper)

    int64_t exp_parity = ((msnzb - sA) & 1);    // exp_parity = B (in paper)
    int64_t shifted_A = int64_t(A[i]) * adjust;    // shifted_A = x' (in paper). i.e., x' = x after range reduction 

// Computes the initial approximation Y
    int64_t A_m = (int64_t(shifted_A) >> (bwA - 2 - m)) & m_mask;   // A_m = e (in paper)
    int64_t Y = lookup_sqrt(A_m, m, exp_parity);    // Y = w (in paper)

// Goldschmidt’s method for t iterations
    int64_t X = (exp_parity ? 2 * shifted_A : shifted_A);

    int64_t x_prev = Y;    // x_prev = a0 (in paper)
    int64_t b_prev = X;    // b_prev = q0 (in paper)
    if (bwA - 2 > sB)
    {
      b_prev = b_prev >> (bwA - 2 - sB);    
    }
    else
    {
      b_prev = b_prev << (sB + 2 - bwA);
    }

    int64_t Y_prev = Y;    // Y_prev = p0 (in paper)
    int64_t x_curr, b_curr, Y_curr;
    for (int j = 0; j < iters; j++)    // iters = 1
    {
      int64_t Y_sq = Y_prev * Y_prev;  // Y_sq = Yi (in paper)

      b_curr = (Y_sq * b_prev) >> (2 * m + 2 * SQRT_LOOKUP_SCALE);    // b_curr = qi (in paper)

      Y_curr = (1.5 * (1LL << (sB + 1))) - b_curr;   // Y_curr = pi (in paper)
      x_curr = x_prev * Y_curr;

      x_curr = x_curr >> (m + SQRT_LOOKUP_SCALE + 1);    // x_curr = ai (in paper)

      x_prev = x_curr;
      b_prev = b_curr;
      Y_prev = Y_curr;
    }

// Uses reciprocal square root of x' to compute the same for x
    int64_t res;
    res = x_curr * sqrt_adjust;
    res = res >> sqrt_adjust_scale;

    B[i] = Saturate(res);

    double flA = (A[i] / double(1LL << sA));
    double sqrt_A = 1.0 / sqrt(flA);
    int64_t actualB = sqrt_A * (1LL << sB);
    error += abs(actualB - B[i]);
  }
  if (party == 2)
  {
    std::cout << "Sqrt Average Error: " << error / (I * J) << std::endl;
  }
}


// void ZK_rsqrt(IntFp *A, int32_t I, int32_t J, 
//               int32_t sA, int32_t sB, int32_t bwA, int32_t bwB, IntFp *B)
// {
//   uint64_t m = ceil(sB / 2.0);   // m = g (in paper)
//   uint64_t iters = 1;

//   int64_t m_mask = (1LL << m) - 1;

//   for (int i = 0; i < I * J; i++){
//     if (party == ALICE){
//       uint64_t msnzb = floor(log2((uint64_t)HIGH64(A[i])));
//     }
//   }
// }