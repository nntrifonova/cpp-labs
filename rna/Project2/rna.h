#pragma once

#ifndef RNA_RNA_H
#define RNA_RNA_H

class rna
{
public:

	enum Nucleotide { A, G, C, T };

private:
	int nucl_num;
	int  capacity;
	size_t* rna_arr;

public:
	void insert(Nucleotide nucl, int ind);
	int get_nucl(int ind) const;
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
		operator Nucleotide() {
			return (Nucleotide)(*rna_array).get_nucl(num);
		}
	};
	class const_reference {
	private:
		int num;
		const rna* rna_array;
	public:
		const_reference(int n, const rna* rna_ar) {
			num = n;
			rna_array = rna_ar;
		}

		operator Nucleotide() const {
			return (Nucleotide)(*rna_array).get_nucl(num);
		}
	};
	rna();
	int count_capacity(int numb);
	rna(int numb);
	rna(Nucleotide nucl, int nucl_num);
	virtual ~rna();

	reference operator [](int ind);
	const_reference operator [](int ind) const;
	friend rna operator+(const rna &rna1, const rna &rna2);
	friend bool operator == (const rna& nucl_1, const rna& nucl_2);
	friend bool operator != (const rna& nucl_1, const rna& nucl_2);
	friend rna operator ! (const rna &rna1);
	friend bool is_complimentary(const rna& nucleotide1, const rna& nucleotide2);
};
#endif //RNA_RNA_H
