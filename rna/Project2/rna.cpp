#include "rna.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>

using namespace std;

int rna::get_nucl(int ind) const
{
	if (ind < 0 || ind>this->nucl_num) {
		exit(0);
	}

	int p_1 = ind / (4 * sizeof(size_t));
	int p_2 = ind % (4 * sizeof(size_t));

	size_t tmp = this->rna_arr[p_1];
	tmp >>= 2 * (sizeof(size_t) * 4 - p_2);
	tmp &= 3;

	switch (tmp) {
	case(0):
		return A;
	case(1):
		return G;
	case(2):
		return C;
	case(3):
		return T;
	}
	exit(0);
}

void rna::insert(Nucleotide nucl, int ind)
{
	size_t letter;
	switch (nucl) {
	case(A):
		letter = 0;
	case(G):
		letter = 1;
	case(C):
		letter = 2;
	case(T):
		letter = 3;
	}

	if (ind = -1) {
		for (int i = 0; i < this->nucl_num; i++) {
			int p_1 = i / (4 * sizeof(size_t));
			int p_2 = i % (4 * sizeof(size_t));

			letter <<= (sizeof(size_t) * 8) - p_2 * 2 - 2;
			size_t tmp = this->rna_arr[p_1];
			tmp |= letter;
			this->rna_arr[p_1] = tmp;
		}
	}
	else {
		int bit_n = ind % 4 * sizeof(size_t);
		int pos = ind / (4 * sizeof(size_t));
		size_t shift = 8 * sizeof(size_t) - 2 - bit_n * 2;
		size_t mask = (size_t)3 << (shift);
		size_t tmp = this->rna_arr[pos];
		this->rna_arr[pos] = (this->rna_arr[pos] & (~mask)) | (letter << shift);
	}
}


rna::rna()
{
	nucl_num = 1;
	capacity = 1;
	rna_arr = new size_t[1];
}

int rna::count_capacity(int numb)
{	
	if (numb * 2 <= sizeof(size_t)) {
		return 1;
	}
	return (numb % 4 == 0 ? numb / 4 : numb / 4 + 1);
}

rna::rna(int numb)
{
	nucl_num = numb;
	capacity = count_capacity(numb);
	rna_arr = new size_t[capacity];
}

rna::rna(Nucleotide nucl, int numb)
{
	nucl_num = numb;
	capacity = count_capacity(numb);
	rna_arr = new size_t[capacity];
	insert(nucl, -1);
}
rna::~rna()
{
	delete[] rna_arr;
}

rna::reference rna::operator[] (int ind) 
{
	reference ref(ind, this);
	return ref;
}

rna::const_reference rna::operator[](int ind) const
{
	return const_reference(ind, this);
}

rna operator+(const rna & rna1, const rna & rna2)
{
	unsigned int last_number1 = (rna1.nucl_num) % 4;
	unsigned int last_number2 = (rna2.nucl_num) % 4;
	unsigned int extra_1 = sizeof(size_t) - last_number1;

	rna rna_new;
	rna_new.capacity = rna1.capacity + rna2.capacity;
	rna_new.nucl_num = rna1.nucl_num + rna2.nucl_num;
	rna_new.rna_arr = new size_t[rna_new.nucl_num / 4];

	memcpy(rna_new.rna_arr, rna1.rna_arr, sizeof(size_t)*rna1.capacity);
	memcpy(&rna_new.rna_arr[rna1.capacity], rna2.rna_arr, sizeof(size_t)*rna2.capacity);

	if (last_number1 != 0) {
		for (int i = rna1.capacity; i < rna_new.capacity; i++) {
			if (i >= 1) {																//is mul 8 necessary?
				rna_new.rna_arr[i] <<= extra_1*8;
			}
			rna_new.rna_arr[i] |= rna_new.rna_arr[i + 1] >> last_number1*8;
		}
	}
	return rna_new;
}

bool operator==(const rna & nucl_1, const rna & nucl_2)
{
	if (nucl_1.nucl_num != nucl_2.nucl_num) {
		return false;
	}
	for (int i = 0; i < nucl_1.capacity; i++) {
		if (nucl_1.rna_arr[i] != nucl_2.rna_arr[i]) {
			return false;
		}
	}
	return true;

}

bool operator!=(const rna & nucl_1, const rna & nucl_2)
{
	return !(nucl_1 == nucl_2);
}

rna operator!(const rna &rna1)
{
	rna rna2;
	rna2.nucl_num = rna1.nucl_num;
	rna2.capacity = rna1.capacity;
	rna2.rna_arr = new size_t[rna1.capacity];

	if (rna1.nucl_num% rna1.capacity == 0) {
		for (int i = 0; i < rna1.capacity; i++) {
			rna2.rna_arr[i] = ~rna1.rna_arr[i];
		}
	}
	else {
		for (int i = 0; i < rna1.capacity; i++) {
			rna2.rna_arr[i] = ~rna1.rna_arr[i];
		}
		rna2.rna_arr[rna2.capacity - 1] >>=sizeof(size_t)- rna2.nucl_num%rna2.capacity;   // 1->0
		rna2.rna_arr[rna2.capacity - 1] <<= sizeof(size_t) - rna2.nucl_num%rna2.capacity;	//0 ->num
	}
	return rna2;
}

bool is_complimentary(const rna & rna1, const rna & rna2)
{
	if (rna1.nucl_num != rna2.nucl_num) {
		return false;
	}
	if (rna1.capacity == rna2.capacity == 0) {
		return true;
	}

	if (!(rna1) == rna2) {
		return true;
	}
	return false;
}
