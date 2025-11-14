#include "prefixTree.h"
#include <fstream>//work with files
#include <iostream>
#include <cstdlib>//stop the program
#include <cstring>//string functions 
using namespace std;

prefixTree::prefixTree() {}/*******/

prefixTree::prefixTree(const string& path, int, int):genomeFilePath(path)
{
    genomeLength = getLengthOfGenome();// get the length 
    genome = new char[genomeLength];// make array for the length
    saveGenome();
    root = new Node();
}

prefixTree::~prefixTree() {
    deleteTrie(root);
    delete[] genome;
    if (genomeKmerArray) {// not null
        for (int i = 0; i < rowCountForGenomeKmerArray; ++i)
            delete[] genomeKmerArray[i];
        delete[] genomeKmerArray;
    }
}


void prefixTree::buildTrieFromGenome(int kmerLen) {
    splitGenomeIntoKmer(kmerLen);
    createTrieFromGenomeKmers(kmerLen);
    std::cout << "Trie built from genome has " << trieNodeCount << " nodes\n";
}

void prefixTree::splitGenomeIntoKmer(int kmer) {
    //l=6 & kmer=3 +1=4 >>> the possible output
    rowCountForGenomeKmerArray = genomeLength - kmer + 1;
    genomeKmerArray = new char* [rowCountForGenomeKmerArray];//rw=4 >> [ 0, 1, 2, 3 ]
    // 2D array >> row figure
    for (int r = 0; r < rowCountForGenomeKmerArray; ++r) {
        genomeKmerArray[r] = new char[kmer];//kmer=3 [ 0,1,2]
        for (int i = 0; i < kmer - 1; ++i) { //loop till fill the 3 places
            genomeKmerArray[r][i] = genome[r + i];
            genomeKmerArray[r][kmer - 1] = '$';
        }
    }
}

void prefixTree::createTrieFromGenomeKmers(int kmerLen) {
    trieNodeCount = 0;// had one root then add the children 
    for (int r = 0; r < rowCountForGenomeKmerArray; ++r) {
        Node* cur = root;
        for (int i = 0; i < kmerLen; ++i) {
            char c = genomeKmerArray[r][i];
            Node** slot = nullptr;
            switch (c) { // this for search when u get the result just break and update the slot pointer
            case 'A': slot = &cur->A; break;// slot is node** which point to the place 
            case 'C': slot = &cur->C; break;
            case 'G': slot = &cur->G; break;
            case 'T': slot = &cur->T; break;
            case '$': slot = &cur->endMarker; break;
            }
            if (!*slot) {// the value of cur->A for example if equall null +1 count 
                //as the start point
                *slot = new Node();// make new node and +1
                ++trieNodeCount;
            }
            cur = *slot;// move cur to child node
        }
       
        cur->terminalCount += 1;// for the repeatation
    }
}

// ——— exact search (≤0 mismatches) ——————————————————————————————————————

int prefixTree::exactSearchTrie(const char* query) const {
    Node* cur = root;

    // Walk each base by index until you hit '\0' or '$'
    for (int i = 0; query[i] != '\0' && query[i] != '$'; ++i) {
        char c = query[i];
        switch (c) {//“AC$”
        case 'A': cur = cur->A; break;
        case 'C': cur = cur->C; break;//i=1: c='C' → cur = A-node->C
        case 'G': cur = cur->G; break;
        case 'T': cur = cur->T; break;
        default:  return 0;        // illegal character
        }
        if (!cur) return 0;            //IF null return 0
    }
    //RETURN count
    return cur->endMarker->terminalCount;
}


// ——— fuzzy search (≤maxMismatch mismatches) ————————————————————————————

int prefixTree::fuzzySearchTrie(const char* query, int maxMismatch) const {
    return fuzzySearchRec(root, query, 0, maxMismatch);
}

int prefixTree::fuzzySearchRec(Node* node,
    const char* q,
    int mismatches,
    int maxMismatch) const
{
    // base case
    if (*q == '$') {//completly identical
        if (!node->endMarker || mismatches > maxMismatch) {
            return 0;
        }
        return node->endMarker->terminalCount;
    }

    int total = 0;
    // try each possible base
    const char bases[4] = { 'A','C','G','T' };
    for (int b = 0; b < 4; ++b) {
        char letter = bases[b];
        Node* child = nullptr;
        switch (letter) {
        case 'A': child = node->A; break;
        case 'C': child = node->C; break;
        case 'G': child = node->G; break;
        case 'T': child = node->T; break;
        }
        if (!child) continue;//jump back to the top
        int newMis = mismatches + (letter != *q);//1 IF they are different and 0 if same
        if (newMis <= maxMismatch) {
            total += fuzzySearchRec(child, q + 1, newMis, maxMismatch);
        }
    }
    return total;
}

// ——— FASTA I/O & cleanup ——————————————————————————————————————————

int prefixTree::getLengthOfGenome() {
    ifstream in(genomeFilePath);
    if (!in) { cerr << "Cannot open " << genomeFilePath << "\n"; exit(1); }
    string line, seq;
    //boolean context is true
    //keep reading lines from the file until there are no more 
    while ( getline(in, line))
        if (!line.empty() && line[0] != '>')
            seq += line;//seq = "ACGTAC" + "TTGG"
    return int(seq.size());
}

void prefixTree::saveGenome() {
    ifstream in(genomeFilePath);
    string line, seq;
    // skip header
    if (getline(in, line) && line.size()/*ensures it’s not empty*/ && line[0] == '>') {}
    while (getline(in, line))
        if (!line.empty())
            seq += line;
    for (int i = 0; i < genomeLength; ++i)
        genome[i] = seq[i];
}

void prefixTree::deleteTrie(Node* n) {
    if (!n) return;
    deleteTrie(n->A);
    deleteTrie(n->C);
    deleteTrie(n->G);
    deleteTrie(n->T);
    deleteTrie(n->endMarker);
    delete n;
}
