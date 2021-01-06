#include "rna.h"


int main() {
	rna rna1(A, 10), rna2(T, 10);
	rna sum = rna1 + rna2;
	return 0;
}