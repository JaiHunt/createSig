#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "uthash.h"
#include <chrono>
#include <Windows.h>

typedef unsigned char byte;

#define SIGNATURE_LEN 64

int DENSITY = 21;
int PARTITION_SIZE;
int NUM_THREADS = 4;

int inverse[256];
const char* alphabet = "CSTPAGNDEQHRKMILVFYW"; // unique animo acids

void seed_random(char* term, int length);
short random_num(short max);
void Init();

int doc_sig[SIGNATURE_LEN];

int WORDLEN;
FILE* sig_file;

typedef struct
{
    char term[100];
    short sig[SIGNATURE_LEN];
    UT_hash_handle hh;
} hash_term;

hash_term* vocab = NULL;


short* compute_new_term_sig(char* term, short* term_sig)
{
    seed_random(term, WORDLEN); // rand nums are seeded based on the term itself and WORDLEN

    int non_zero = SIGNATURE_LEN * DENSITY / 100; // set everything to 0

    int positive = 0;
    while (positive < non_zero / 2) // 10.5% +ve
    {
        short pos = random_num(SIGNATURE_LEN); 
        if (term_sig[pos] == 0)
        {
            term_sig[pos] = 1;
            positive++;
        }
    }

    int negative = 0;
    while (negative < non_zero / 2) // 10.5% -ve
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

DWORD WINAPI ParallelHash(LPVOID lpThreadParameter) {

    return 0;
}

short* find_sig(char* term)
{
    hash_term* entry; // cache that remembers kmers used before

    HASH_FIND(hh, vocab, term, WORDLEN, entry);

    // create parallel hash tables
    if (entry == NULL)
    {
        // lock hash function
        HASH_FIND(hh, vocab, term, WORDLEN, entry);
        entry = (hash_term*)malloc(sizeof(hash_term));
        strncpy_s(entry->term, sizeof(entry->term), term, WORDLEN);
        memset(entry->sig, 0, sizeof(entry->sig));
        compute_new_term_sig(term, entry->sig);
        HASH_ADD(hh, vocab, term, WORDLEN, entry);
    }

    return entry->sig;
}


void signature_add(char* term)
{
    /* Creates Handle */
    HANDLE* handles = new HANDLE[NUM_THREADS];

    short* term_sig; 

    for (int i = 0; i < NUM_THREADS; i++) {
        handles[i] = CreateThread(NULL, 0, ParallelHash, (LPVOID)i, 0, NULL);
        term_sig = find_sig(term);
    }
    for (int i = 0; i < SIGNATURE_LEN; i++)
        doc_sig[i] += term_sig[i];
}

int doc = 0;

void compute_signature(char* sequence, int length)
{
    memset(doc_sig, 0, sizeof(doc_sig)); // sets total to all 0s

    for (int i = 0; i < length - WORDLEN + 1; i++) // goes through each kmer that makes the parition and adds them to the above doc sig
        signature_add(sequence + i); // Adds a complete row each loop iteration

    fwrite(&doc, sizeof(int), 1, sig_file); // save document number to sig file

    for (int i = 0; i < SIGNATURE_LEN; i += 8) // flatten and output to sig file
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
