#include "rna.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;


void rna::count_capacity()
{
	if (capacity == 0) rna_arr = nullptr;
	if (ceil(1.0 * nucl_num / NUCL_SIZE_T) <= capacity) return;

	size_t old_capacity = capacity;
	capacity = ceil(2.0 * nucl_num / NUCL_SIZE_T);
	auto* tmp = new size_t[capacity];
	for (size_t i = 0; i < old_capacity; i++)
		tmp[i] = rna_arr[i];
	for (size_t i = old_capacity; i < capacity; i++)
		tmp[i] = fill_nucl(nucl);
	delete[] rna_arr;
	rna_arr = tmp;
}

size_t rna::fill_nucl(Nucleotide nucl)
{
	size_t nucleotide = 0;
	for (size_t i = 0; i < NUCL_SIZE_T; i++) {
		nucleotide <<= 2u;
		nucleotide += nucl;
	}
	return nucleotide;
}


void rna::add_nucl(Nucleotide nucl)
{
	nucl_num++;
	count_capacity();
	write_nucl(nucl,nucl_num - 1);
}

Nucleotide rna::get_nucl(size_t index) const
{
	if (index >= nucl_num) {
		return nucl;
	}
	size_t block_idx = index / NUCL_SIZE_T;
	size_t pair_idx = index % NUCL_SIZE_T;
	size_t block = rna_arr[block_idx];
	return (Nucleotide)(((block >> (2u * (NUCL_SIZE_T - pair_idx - 1u)))) & 3u);
}

void rna::write_nucl(Nucleotide nucl, size_t index)
{
	if (index >= nucl_num) {
		size_t from = nucl_num;
		nucl_num = index + 1;
		count_capacity();
		for (size_t i = from; i < index; i++)
			write_nucl(nucl, i);
	}

	size_t block_idx = index / NUCL_SIZE_T;
	size_t pair_idx = index % NUCL_SIZE_T;
	size_t block = rna_arr[block_idx];
	size_t mask = ((size_t)3) << (2u * (NUCL_SIZE_T - pair_idx - 1));
	size_t _nucl = fill_nucl(nucl) & mask;

	rna_arr[block_idx] = (block & (~mask)) + nucl;
}


rna::rna()
{
	nucl = A;
	nucl_num = 1;
	capacity = 1;
	rna_arr = new size_t[1];
}


rna::rna(Nucleotide _nucl, size_t numb)
{
	nucl = _nucl;
	nucl_num = numb;
	count_capacity();
}
rna::rna(const rna& rnk)
{
	nucl = rnk.nucl;
	nucl_num = rnk.nucl_num;
	capacity = rnk.capacity;
	rna_arr = new size_t[capacity];
	for (size_t i = 0; i < capacity; i++)
		rna_arr[i] = rnk.rna_arr[i];
}
rna::~rna()
{
	delete[] rna_arr;
}

rna& rna::operator=(const rna& _rna)
{
	if (&_rna == this) return *this;

	nucl = _rna.nucl;
	nucl_num = _rna.nucl_num;
	capacity = _rna.capacity;

	delete[] rna_arr;
	rna_arr = nullptr;

	rna_arr = new size_t[capacity];
	for (size_t i = 0; i < capacity; i++)
		rna_arr[i] = _rna.rna_arr[i];
	return *this;
}

rna::reference rna::operator[] (size_t ind)
{
	return reference (ind, *this);
}

rna::const_reference rna::operator[](size_t ind) const
{
	return const_reference(ind, *this);
}

rna& rna::operator+=(Nucleotide nucl)
{
	add_nucl(nucl);
	return *this;
}

rna operator+(const rna& r1, const rna& r2)
{
	rna rnk(r1);
	for (size_t i = 0; i < r2.nucl_num; i++)
		rnk.add_nucl(r2.get_nucl(i));
	return rnk;
}

bool operator==(const rna& nucl_1, const rna& nucl_2)
{
	if (nucl_1.nucl_num != nucl_2.nucl_num) {
		return false;
	}
	for (size_t i = 0; i < nucl_1.capacity; i++) {
		if (nucl_1.rna_arr[i] != nucl_2.rna_arr[i]) {
			return false;
		}
	}
	return true;

}

bool operator!=(const rna& nucl_1, const rna& nucl_2)
{
	return !(nucl_1 == nucl_2);
}

rna rna::operator!()
{
	rna rnk(*this);
	for (size_t i = 0; i < nucl_num; i++) {
		Nucleotide cur = get_nucl(i);
		rnk.write_nucl((Nucleotide)(3 - cur), i);
	}
	return rnk;
}

rna rna::split(size_t ind)
{
	rna rna1(nucl, 0);
	for (size_t i = 0; i < min(ind, nucl_num); i++)
		rna1 += get_nucl(i);
	return rna1;

}

bool is_complimentary(const rna& rna1, const rna& rna2)
{
	if (rna1.nucl_num != rna2.nucl_num) {
		return false;
	}
	for (size_t i = 0; i < rna1.nucl_num; i++) {
		if (rna1.get_nucl(i) + rna2.get_nucl(i) != 3) {
			return false;
		}
	}
	return true;
	
}

