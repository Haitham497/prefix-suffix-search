#ifndef PREFIXTREE_H
#define PREFIXTREE_H

#include <string>
using namespace std;
struct Node {
    Node* A = nullptr;
    Node* C = nullptr;
    Node* G = nullptr;
    Node* T = nullptr;
    Node* endMarker = nullptr;
    int   terminalCount = 0;// count for how many the repetion of the kmer
};

class prefixTree {
public:
    prefixTree();
    prefixTree(const string& path, int /*readCount*/, int /*readLen*/);
    ~prefixTree();

    // build a trie over *all* genome k-mers of length = kmerLen
    void buildTrieFromGenome(int kmerLen);

    // exact search (≤0 mismatches)
    int  exactSearchTrie(const char* query) const;

    // fuzzy search (≤ maxMismatch mismatches; default = 1)
    int  fuzzySearchTrie(const char* query, int maxMismatch = 1) const;
    int getNodeCount() const { return trieNodeCount; }

private:
    string genomeFilePath;
    char* genome = nullptr;
    int genomeLength = 0;

    char** genomeKmerArray = nullptr;
    int rowCountForGenomeKmerArray = 0;

    Node* root = nullptr;
    int trieNodeCount = 0;

    int getLengthOfGenome();
    void saveGenome();
    void splitGenomeIntoKmer(int kmer);
    void createTrieFromGenomeKmers(int kmerLen);
    void deleteTrie(Node* n);

    
    int  fuzzySearchRec(Node* node,
        const char* q,
        int mismatches, int maxMismatch) const;// check and count mismatch if it reach the max stop
};

#endif // PREFIXTREE_H
