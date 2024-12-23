#ifndef RO_ZK_RAM_H__
#define RO_ZK_RAM_H__
#include "emp-zk/extensions/ram-zk/ostriple.h"

template <typename IO>
class ROZKRAM
{
public:
	double check1 = 0, check2 = 0, check3 = 0; // 时间开销： check1：read的开销（无check）；check2: 生成sorted_list，以及check是否升序的开销；check3: check是否互为permutation的开销
	int party;
	int index_sz, val_sz; // index_sz指的是一个index用多少个bit表示，并不是list中index的数量
	uint64_t step = 0;
	vector<uint64_t> clear_mem;
	vector<block> check_MAC;
	vector<__uint128_t> list;
	GaloisFieldPacking gp;
	IO *io;
	block Delta;
	F2kOSTriple<IO> *ostriple = nullptr;

	ROZKRAM(int _party, int index_sz, int val_sz) : party(_party), index_sz(index_sz), val_sz(val_sz)
	{
		ZKBoolCircExec<IO> *exec = (ZKBoolCircExec<IO> *)(CircuitExecution::circ_exec);
		io = exec->ostriple->io;
		Delta = exec->ostriple->delta;
		ostriple = new F2kOSTriple<IO>(party, exec->ostriple->threads, exec->ostriple->ios, exec->ostriple->ferret, exec->ostriple->pool);
	}

	~ROZKRAM()
	{
		delete ostriple;
	}

	void init(vector<Integer> &data)
	{
		// 1）list: 将index+原始值(不带MAC)放到队列中；
		// 2）check_MAC: 将index+原始值(带MAC)packing到F2k上再放到队列中
		uint64_t val = 0;
		block m;
		for (size_t i = 0; i < data.size(); ++i)
		{
			val = data[i].reveal<uint64_t>(ALICE);	// val是data中的witness
			list.push_back(pack((uint64_t)i, val)); // pack(1，2)：将i和val连接到一起（不是F2k上的意思，只是单纯的连在一起）   push_back：感觉就是放到表中的意思
			clear_mem.push_back(val);
			pack(m, Integer(index_sz + 1, (uint64_t)i, PUBLIC), data[i]); // pack(1，2，3)：将 index 和 data 一起 packing到F2k上
			check_MAC.push_back(m);										  // check_MAC 是为了 verify 两个队列 (list 和 sorted_list) 存在permutation关系； sorted_list将在check()函数中生成
			assert(data[i].size() == val_sz);
		}
	}

	Integer read(const Integer &index)
	{
		// 1）读出index位上对应的value；
		// 2）将该index+value(不带MAC)append到list；
		// 3）将该index+value(带MAC，F2k)append到check_MAC
		// 4）一个batch_size放满后，check
		auto start = clock_start();
		uint64_t clear_index = index.reveal<uint64_t>(ALICE);
		uint64_t tmp = 0;
		if (party == ALICE)
		{
			tmp = clear_mem[clear_index];			// tmp 就是 index 对应的 value
			list.push_back(pack(clear_index, tmp)); // 将该次查询append到原始list中
		}
		else
			list.push_back(0);
		Integer res = Integer(val_sz, tmp, ALICE);
		block m;
		pack(m, index, res);
		check_MAC.push_back(m); // 放到 check_MAC 中
		++step;
		check1 += time_from(start);
		if (step == clear_mem.size() * 8)
			check(); //
		return res;
	}

	void pack(block &mac, const Integer &index, const Integer &val)
	{
		// 1）将index和val合并pack到F2k
		block res1, res2;
		vector_inn_prdt_sum_red(&res1, (block *)(val.bits.data()), gp.base, val_sz); //  gp.base 应该是gp的基
		vector_inn_prdt_sum_red(&res2, (block *)(index.bits.data()), gp.base + val_sz, index_sz);
		mac = res1 ^ res2; // 将index和value依次packing到 一个mac 中
	}

	__uint128_t pack(uint64_t index, uint64_t value)
	{
		// 1）简单的将index和value连接到一起，明文操作
		return (((__uint128_t)index) << val_sz) | value; // "|"是或运算符，只要有1，结果就为1
	}

	void check()
	{
		// 1）生成sorted_list和sorted_MAC
		//    - 生成sorted_list（明文，不带MAC）
		//    - 生成sort_value：生成sorted_list中的value的MAC；sort_index：生成sorted_list中的index的MAC
		//    - 生成sorted_MAC（带MAC，F2k）：将index+value packing到F2k上再放到队列中
		// 2）第一步check: check sorted list是否正确升序排序 (使用sort_value和sort_index)
		// 3）第二部check: check list 和 sorted_list是否互为permutation (使用sorted_list和sorted_MAC)
		auto start = clock_start();
		vector<__uint128_t> sorted_list;
		sorted_list = vector<__uint128_t>(list.begin(), list.end());
		if (party == ALICE)
		{
			sort(sorted_list.begin(), sorted_list.end());
		}
		vector<Integer> sort_value, sort_index;
		vector<block> sorted_X, sorted_MAC;
		block m;
		__uint128_t val_mask = ((__uint128_t)1 << val_sz) - 1;
		for (size_t i = 0; i < list.size(); ++i)
		{
			__uint128_t item;
			if (party == ALICE)
				item = sorted_list[i];
			else
				item = 0;
			uint64_t i1 = item & val_mask;
			sort_value.push_back(Integer(val_sz, i1, ALICE)); // Integer(val_sz, i1, ALICE)：将sorted_list中的每个value生成MAC
			uint64_t i2 = LOW64(item >> val_sz);
			sort_index.push_back(Integer(index_sz + 1, i2, ALICE)); // Integer(index_sz+1, i2, ALICE)：将sorted_list中的每个index生成MAC
			pack(m, sort_index[i], sort_value[i]);					// 将 index 和 data 一起 packing到F2k上
			sorted_MAC.push_back(m);								// sorted_MAC
		}

		bool cheat = true;
		for (size_t i = 0; i < list.size() - 1; ++i)
		{ // 进行第一步的check （check sorted list是否正确升序排序）
			Bit eq = !(sort_index[i].geq(sort_index[i + 1])) | (sort_index[i].equal(sort_index[i + 1]) & sort_value[i].equal(sort_value[i + 1]));
			bool res = eq.reveal<bool>(PUBLIC);
			cheat = cheat and res; // 如果没有cheat， eq应=1，即res=1
		}
		if (!cheat)
			error("cheat!");
		check2 += time_from(start);

		start = clock_start();
		sync_zk_bool<IO>();
		check_set_euqality(sorted_list, sorted_MAC, list, check_MAC); // 进行第二步的check （check list 和 sorted list是否互为permutation）
		check3 += time_from(start);
		list.resize(clear_mem.size());
		check_MAC.resize(clear_mem.size());
		step = 0;
	}

	void vector_inn_prdt(block &xx, block &mm, vector<__uint128_t> &X, vector<block> &MAC, block r)
	{
		block x, m;
		size_t i = 1;
		block tmp = (block)X[0];
		ostriple->compute_add_const(xx, mm, tmp, MAC[0], r);
		while (i < list.size())
		{
			tmp = (block)X[i];
			ostriple->compute_add_const(x, m, tmp, MAC[i], r);
			ostriple->compute_mul(xx, mm, xx, mm, x, m);
			++i;
		}
	}

	void vector_inn_prdt_bch2(block &xx, block &mm, vector<__uint128_t> &X, vector<block> &MAC, block r)
	{
		block x[2], m[2], t[2];
		size_t i = 1;
		block tmp = (block)X[0];
		ostriple->compute_add_const(xx, mm, tmp, MAC[0], r);
		while (i < list.size() - 1)
		{
			for (int j = 0; j < 2; ++j)
			{
				t[j] = (block)X[i + j];
				ostriple->compute_add_const(x[j], m[j], t[j], MAC[i + j], r);
			}
			ostriple->compute_mul3(xx, mm, x[0], m[0],
								   x[1], m[1], xx, mm);
			i += 2;
		}
		while (i < list.size())
		{
			t[0] = (block)X[i];
			ostriple->compute_add_const(x[0], m[0], t[0], MAC[i], r);
			ostriple->compute_mul(xx, mm, xx, mm, x[0], m[0]);
			++i;
		}
	}

	void vector_inn_prdt_bch3(block &xx, block &mm, vector<__uint128_t> &X, vector<block> &MAC, block r)
	{
		block x[3], m[3], t[3];
		size_t i = 1;
		block tmp = (block)X[0];
		ostriple->compute_add_const(xx, mm, tmp, MAC[0], r);
		while (i < list.size() - 2)
		{
			for (int j = 0; j < 3; ++j)
			{
				t[j] = (block)X[i + j];
				ostriple->compute_add_const(x[j], m[j], t[j], MAC[i + j], r);
			}
			ostriple->compute_mul4(xx, mm, x[0], m[0],
								   x[1], m[1], x[2], m[2], xx, mm);
			i += 3;
		}
		while (i < list.size())
		{
			t[0] = (block)X[i];
			ostriple->compute_add_const(x[0], m[0], t[0], MAC[i], r);
			ostriple->compute_mul(xx, mm, xx, mm, x[0], m[0]);
			++i;
		}
	}

	void vector_inn_prdt_bch4(block &xx, block &mm, vector<__uint128_t> &X, vector<block> &MAC, block r)
	{
		// TODO: 向量X中所有元素（其MAC为MAC）+r后再依次相乘，得到最终结果xx和相应的MAC值mm
		// 1）计算X+r和其MAC
		// 2）每5个一组做乘法
		block x[4], m[4], t[4];
		size_t i = 1;
		block tmp = (block)X[0];
		// compute_add_const: (valb, macb, vala, maca, c) 计算 [a] + c = [b]  c为常数
		ostriple->compute_add_const(xx, mm, tmp, MAC[0], r); // F2kOSTriple<IO> *ostriple  [xx] = [tmp] + r
		while (i < list.size() - 3)
		{ // 依次将X中每四个元素一组：首先，每个元素加一个常数r；随后，所有-r后的元素做乘法
			for (int j = 0; j < 4; ++j)
			{
				t[j] = (block)X[i + j];
				ostriple->compute_add_const(x[j], m[j], t[j], MAC[i + j], r);
			}
			// compute_mul5：5个MACed values 做乘法，结果（value，MAC）放入第一二个参数（xx, mm）
			ostriple->compute_mul5(xx, mm, x[0], m[0],
								   x[1], m[1], x[2], m[2], x[3], m[3], xx, mm);
			i += 4; //
		}			// 上述整个循环完成后，得到的是：X-r后的向量中的逐元素乘法，输出（xx, mm）

		// 下述几行应该是忘记注释了，所以我给它注释掉了
		// while(i < list.size()) {
		// 	t[0] = (block)X[i];
		// 	ostriple->compute_add_const(x[0], m[0], t[0], MAC[i], r);
		// 	ostriple->compute_mul(xx, mm, xx, mm, x[0], m[0]);
		// 	++i;
		// }
	}

	// mult batch 4
	void check_set_euqality(vector<__uint128_t> &sorted_X, vector<block> &sorted_MAC, vector<__uint128_t> &check_X, vector<block> &check_MAC)
	{
		// TODO: check 两个 lists 是否互为 permutation
		// 1）将两个 lists 分别转换成 单变量多项式；输入变量r评估多项式
		// 2）check 两个多项式所得评估结果的相等性
		block r, val[2], mac[2];									   // block：128-bit的整数
		r = io->get_hash_block();									   // 从F2k中随机采样r
		vector_inn_prdt_bch4(val[0], mac[0], sorted_X, sorted_MAC, r); // 将 sorted list 转换成 多项式 后 输入r 得到计算结果
		vector_inn_prdt_bch4(val[1], mac[1], check_X, check_MAC, r);   // 将 list 转换成 多项式 后 输入r 得到计算结果

		// TODO comparison
		if (party == ALICE)
		{ // check
			io->send_data(mac, 2 * sizeof(block));
			io->flush();
		}
		else
		{
			block macrecv[2];
			io->recv_data(macrecv, 2 * sizeof(block));
			mac[0] ^= macrecv[0];
			mac[1] ^= macrecv[1];
			if (memcmp(mac, mac + 1, 16) != 0)
			{ // 16 * 8 = 128; int memcmp(const void *str1, const void *str2, size_t n)) 把存储区str1 和存储区str2 的前n个字节进行比较
				error("check set equality failed!\n");
			}
		}
	}

	void check_MAC_valid(block X, block MAC)
	{
		if (party == ALICE)
		{
			io->send_data(&X, 16);
			io->send_data(&MAC, 16);
		}
		else
		{
			block M = zero_block, x = zero_block;
			io->recv_data(&x, 16);
			io->recv_data(&M, 16);
			gfmul(x, Delta, &x);
			x = x ^ MAC;
			if (memcmp(&x, &M, 16) != 0)
			{
				error("check_MAC failed!\n");
			}
		}
	}
};
#endif // RO_ZK_RAM_H__
