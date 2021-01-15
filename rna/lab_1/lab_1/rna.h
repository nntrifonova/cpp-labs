#pragma once

#ifndef RNA_RNA_H
#define RNA_RNA_H

#include <iostream>
#define NUCL_SIZE_T (sizeof(size_t) * 4)

using namespace std;
enum Nucleotide { A, G, C, T };

class rna
{
private:
	Nucleotide nucl;
	size_t nucl_num;
	size_t  capacity;
	size_t* rna_arr;

	void count_capacity();
	size_t fill_nucl(Nucleotide nucl);

	void add_nucl(Nucleotide nucl);
	Nucleotide get_nucl(size_t index) const;
	void write_nucl(Nucleotide nucl, size_t index);

public:
	class const_reference;
	class reference {
	private:
		size_t index;
		rna* rna_array;
	public:
		reference(size_t n, rna &rna_ar) {
			index = n;
			rna_array = &rna_ar;
		}
		reference& operator=(Nucleotide nucl) {
			(*rna_array).write_nucl(nucl, index);
			return *this;
		}
		
	};
	class const_reference {
	private:
		int num;
		const rna* rna_array;
	public:
		const_reference(int n, const rna& rna_ar) {
			num = n;
			rna_array = &rna_ar;
		}

		operator Nucleotide() const {
			return (Nucleotide)(*rna_array).get_nucl(num);
		}
	};
	rna();
	rna(Nucleotide nucl, size_t nucl_num);
	rna(const rna& rnk);
	virtual ~rna();

	rna& operator=(const rna& _rna);
	reference operator [](size_t ind);
	const_reference operator [](size_t ind) const;
	rna& operator +=(Nucleotide nucl);
	friend rna operator+(const rna& r1, const rna& r2);
	friend bool operator == (const rna& nucl_1, const rna& nucl_2);
	friend bool operator != (const rna& nucl_1, const rna& nucl_2);
	rna operator ! ();
	friend bool is_complimentary(const rna& nucleotide1, const rna& nucleotide2);
	rna split(size_t ind);
};
#endif //RNA_RNA_H

