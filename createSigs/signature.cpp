#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "uthash.h"
#include <chrono>
#include <vector>
#include <string>
#include <numeric>

typedef unsigned char byte;

#define SIGNATURE_LEN 64

int DENSITY = 21;
int PARTITION_SIZE;

int inverse[256];
const char* alphabet = "CSTPAGNDEQHRKMILVFYW";

void seed_random(const char* term, int length);
short random_num(short max);
void Init();

//int doc_sig[SIGNATURE_LEN];

int WORDLEN;
FILE* sig_file;

struct signature
{
    int doc = 0;
    int doc_sig[SIGNATURE_LEN];
    uint8_t sig_size[SIGNATURE_LEN / 8];

    // pack the contents of the int array into the byte array 
    void Flatten()
    {
        for (int i = 0; i < SIGNATURE_LEN / 8; i++)
        {
            sig_size[i] = 0;
            for (int j = 0; j < 8; j++)
                sig_size[i] |= (doc_sig[8 * i + j] > 0) << (7 - j);
        }
    }

    // saves doc int and dumps byte array to sig file
    void Write()
    {
        fwrite(&doc, sizeof(int), 1, sig_file);
        fwrite(sig_size, sizeof(sig_size), 1, sig_file);
    }

};

typedef struct
{
    char term[100];
    short sig[SIGNATURE_LEN];
    UT_hash_handle hh;
} hash_term;

hash_term* vocab = NULL;

short* compute_new_term_sig(const char* term, short* term_sig)
{
    seed_random(term, WORDLEN);

    int non_zero = SIGNATURE_LEN * DENSITY / 100;

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

short* find_sig(const char* term)
{
    hash_term* entry;
    HASH_FIND(hh, vocab, term, WORDLEN, entry);
    if (entry == NULL)
    {
        entry = (hash_term*)malloc(sizeof(hash_term));
        strncpy_s(entry->term, sizeof(entry->term), term, WORDLEN);
        memset(entry->sig, 0, sizeof(entry->sig));
        compute_new_term_sig(term, entry->sig);
        HASH_ADD(hh, vocab, term, WORDLEN, entry);
    }

    return entry->sig;
}

signature s;

void signature_add(int* doc_sig, const char* term)
{
    short* term_sig = find_sig(term);
    for (int i = 0; i < SIGNATURE_LEN; i++)
        doc_sig[i] += term_sig[i];
}

void compute_signature(const char* sequence, int length, signature& s)
{
    //memset(s.doc_sig, 0, sizeof(s.doc_sig));
    std::fill(s.doc_sig, s.doc_sig + SIGNATURE_LEN, 0);
    for (int i = 0; i < length - WORDLEN + 1; i++)
        signature_add(s.doc_sig, sequence + i);
    s.Flatten(); // flatten and output to sig file
}

#define min(a,b) ((a) < (b) ? (a) : (b))

// Process arrays + accepts a vector of signatures.
void partition(std::vector<std::string> input, std::vector<std::vector<signature>> output, int length)
{
#pragma omp parallel
    {
        int m = (input.size() - 1) / (PARTITION_SIZE / 2); // signature blocks needed 
        output.resize(m); // resize output to hold m signature blocks
#pragma omp for
        for (int i = 0; i < m; i += PARTITION_SIZE / 2) {
            compute_signature((char*)&input + i, min(PARTITION_SIZE, length - i), s);
        }
    }
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
    int WORDS = power(20, WORDLEN);

    for (int i = 0; i < strlen(alphabet); i++)
        inverse[alphabet[i]] = i;

    auto start = std::chrono::high_resolution_clock::now();

    FILE* file;
    errno_t OK = fopen_s(&file, filename, "r");

    if (OK != 0)
    {
        fprintf(stderr, "Error: failed to open file %s\n", filename);
        return 1;
    }

    char outfile[256];
    sprintf_s(outfile, 256, "%s.part%d_sigs%02d_%d", filename, PARTITION_SIZE, WORDLEN, SIGNATURE_LEN);
    fopen_s(&sig_file, outfile, "w");

    std::vector<std::string> in(1000); // vector of 1000 input strings
    std::vector<std::vector<signature>> out(1000); // matrix of signature structs with 1000

    static char buffer[10000];

    int counter = 0;

    // reads input file in batches. To then later process in the partition function
    while (!feof(file))
    {
        fgets(buffer, 10000, file); // skip the meta data line
        fgets(buffer, 10000, file);
        int n = (int)strlen(buffer) - 1;
        buffer[n] = 0;
        in[counter++] = buffer;
        if (counter == 1000) {
            partition(in, out, counter);
            counter = 0;
        }
    }
    if (counter > 0)
        partition(in, out, counter);

    s.Write();

    fclose(file);
    fclose(sig_file);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    printf("%s %f seconds\n", filename, duration.count());

    return 0;
}
