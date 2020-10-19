#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "uthash.h"
#include <chrono>
#include <thread>
#include <mutex>
#include <omp.h>

typedef unsigned char byte;

#define SIGNATURE_LEN 64
#define N 4 

int DENSITY = 21;
int PARTITION_SIZE;

int inverse[256];
const char* alphabet = "CSTPAGNDEQHRKMILVFYW"; // unique animo acids

/* CSTPAGNDEQHRKMILVFYW = 01234567890123456789
   e.g. MKI = 0, -1, 1, 0, 0, 1, 1, -1, 0, 0, -1
   M = 13, K = 12, I = 14 = 13*20^2 + 12*20^1 + 1*20^0
   This number will represent the specific kmer where that is 0 between 20^WORDLEN
*/

void seed_random(char* term, int length);
short random_num(short max);
void Init();

int doc_sig[SIGNATURE_LEN];

int WORDLEN;
FILE* sig_file;

short* compute_new_term_sig(char* term, short* term_sig)
{
	seed_random(term, WORDLEN); // rand nums are seeded based on the term itself and WORDLEN

	int non_zero = SIGNATURE_LEN * DENSITY / 100; // set everything to 0


	// 10.5% +ve
	int positive = 0;
	while (positive < non_zero / 2)
	{
		short pos = random_num(SIGNATURE_LEN);
		if (term_sig[pos] == 0)
		{
			term_sig[pos] = 1;
			positive++;
		}
	}

	// 10.5% -ve
	int negative = 0;
	while (negative < non_zero / 2)
	{
		short pos = random_num(SIGNATURE_LEN);
		if (term_sig[pos] == 0)
		{
			term_sig[pos] = -1;
			negative++;
		}
	}
	return term_sig;
}

#pragma omp parallel
{
	hash_term* vocab = NULL;

	#pragma omp parallel for
		for (int i = 0; i < N; i++) {
			run_parahash(&vocab);
		}
	}
}

void run_parahash(hash_term** vocab)
{
	hash_term* entry;
	HASH_FIND(hh, *vocab, term, WORDLEN, entry);
	HASH_ADD(hh, *vocab, term, WORDLEN, entry);
}

/*
	hash find has the highest hotpath and cpu usage, one possible way to fix this is having to
	find functions. Where hash find -> returns null -> lock -> hash find -> then adds the entry.
	While running different hash tables concurrently on separate threads. To then combine the entries
	towards the end.
*/
short* find_sig(char* term, hash_term** vocab)
{
	hash_term* entry; // cache that remembers kmers used before

	HASH_FIND(hh, *vocab, term, WORDLEN, entry); // hotpath
	if (entry == NULL)
	{
		entry = (hash_term*)malloc(sizeof(hash_term));
		strncpy_s(entry->term, sizeof(entry->term), term, WORDLEN);
		memset(entry->sig, 0, sizeof(entry->sig));
		compute_new_term_sig(term, entry->sig);
		HASH_ADD(hh, *vocab, term, WORDLEN, entry);
	}
	return entry->sig;
}

void signature_add(char* term)
{
	short* term_sig = find_sig(term, );
	for (int i = 0; i < SIGNATURE_LEN; i++)
		doc_sig[i] += term_sig[i];
}

int doc = 0;

void compute_signature(char* sequence, int length)
{
	memset(doc_sig, 0, sizeof(doc_sig)); // sets total to all 0s

	// goes through each kmer that makes the parition and adds them to the above doc sig
	for (int i = 0; i < length - WORDLEN + 1; i++)
		signature_add(sequence + i); // Adds a complete row each loop iteration

	// save document number to sig file
	fwrite(&doc, sizeof(int), 1, sig_file);

	// flatten and output to sig file
	for (int i = 0; i < SIGNATURE_LEN; i += 8)
	{
		byte c = 0;
		for (int j = 0; j < 8; j++)
			c |= (doc_sig[i + j] > 0) << (7 - j);
		fwrite(&c, sizeof(byte), 1, sig_file);
	}
}

#define min(a,b) ((a) < (b) ? (a) : (b))

void partition(char* sequence, int length)
{
	int i = 0;
	do
	{
		compute_signature(sequence + i, min(PARTITION_SIZE, length - i)); // function computes sig for just one partition
		i += PARTITION_SIZE / 2; // increments start of the partition by half the partition size (half overlap between one partition and the next)
	} while (i + PARTITION_SIZE / 2 < length);
	doc++;
}

int power(int n, int e)
{
	int p = 1;
	for (int j = 0; j < e; j++)
		p *= n;
	return p;
}

int main(int argc, char* argv[])
{
	//const char* filename = "qut2.fasta";
	const char* filename = "qut3.fasta";

	WORDLEN = 3;
	PARTITION_SIZE = 16;
	int WORDS = power(20, WORDLEN); // only 20 letters of the aplphabet are used, hence base 20

	for (int i = 0; i < strlen(alphabet); i++)
		inverse[alphabet[i]] = i; // for each position in the aplhabet, the inverse is calculated to be position i

	auto start = std::chrono::high_resolution_clock::now();

	FILE* file;
	errno_t OK = fopen_s(&file, filename, "r");

	if (OK != 0)
	{
		fprintf(stderr, "Error: failed to open file %s\n", filename);
		return 1;
	}

	/* writing to file */
	char outfile[256];
	sprintf_s(outfile, 256, "%s.part%d_sigs%02d_%d", filename, PARTITION_SIZE, WORDLEN, SIGNATURE_LEN);
	fopen_s(&sig_file, outfile, "w");

	char buffer[10000];
	while (!feof(file))
	{
		fgets(buffer, 10000, file); // skip meta data line - every first or so line that contains irrelevant info  
		fgets(buffer, 10000, file);
		int n = (int)strlen(buffer) - 1; // reads the sig data and returns n characters 
		buffer[n] = 0;
		partition(buffer, n); // Written out here
	}
	fclose(file);

	fclose(sig_file);

	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end - start;

	printf("%s %f seconds\n", filename, duration.count());

	return 0;
}
