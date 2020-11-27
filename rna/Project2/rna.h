#pragma once

#ifndef RNA_RNA_H
#define RNA_RNA_H

enum Nucleotide { A, G, C, T };

class rna
{
private:
	int nucl_num;
	int  capacity;
	size_t* rna_arr;
	void insert(Nucleotide nucl, int ind);
	Nucleotide get_nucl(int ind);
	class reference {
	private:
		int num;
		rna* rna_array;
	public:
		reference(int n, rna* rna_ar) {
			num = n;
			rna_array = rna_ar;
		}
		void operator=(Nucleotide nucl) {
			return (*rna_array).insert(nucl, num);
		}
		operator Nucleotide() const {
			return (Nucleotide)(*rna_array).get_nucl(num);
		}
	};
	class reference_const {
	private:
		int ind;
		const rna* rna_array;
	public:
		reference_const(int index, rna* rna_arr) {
			ind = index;
			rna_array = rna_arr;
		}
		operator Nucleotide() const {
			return (Nucleotide)(*rna_array).get_nucl(ind);
		}

	};
public:
	rna();
	int count_capacity(int numb);
	rna(int numb);
	rna(Nucleotide nucl, int nucl_num);
	virtual ~rna();

	reference operator [](int ind);
	reference_const operator [](int ind) const;
	friend rna operator+(const rna &rna1, const rna &rna2);
	friend bool operator == (const rna& nucl_1, const rna& nucl_2);
	friend bool operator != (const rna& nucl_1, const rna& nucl_2);
	friend rna operator ! (const rna &rna1);
	friend bool is_complimentary(const rna& nucleotide1, const rna& nucleotide2);
	reference operator [](int index);
};
#endif //RNA_RNA_H
