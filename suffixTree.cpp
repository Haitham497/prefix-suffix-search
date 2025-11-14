#include "suffixTree.h"
#include <stdexcept>

suffixTree::suffixTree(string inputGenome)
{

    nodecount = 0;
    filepath = inputGenome;
    genomesize = grablength();
    genome = new char[genomesize];
    savinggenome();

    root = new Node();
    root->A = nullptr;
    root->C = nullptr;
    root->G = nullptr;
    root->T = nullptr;
    root->$ = nullptr;
}

suffixTree::~suffixTree()
{
    if (root == nullptr) return;

    stack<Node*> freestack;
    freestack.push(root);

    while (!freestack.empty())
    {
        Node* curr = freestack.top();
        freestack.pop();

        if (curr == nullptr) continue;

        if (curr->A != nullptr) freestack.push(curr->A);
        if (curr->C != nullptr) freestack.push(curr->C);
        if (curr->G != nullptr) freestack.push(curr->G);
        if (curr->T != nullptr) freestack.push(curr->T);
        if (curr->$ != nullptr) freestack.push(curr->$);

        delete curr;
    }

    if (genome != nullptr) {
        delete[] genome;
        genome = nullptr;
    }
}

int suffixTree::grablength()
{
    string header;
    string sequence;
    string line;
    ifstream file(filepath);
    if (!file.is_open()) {
        throw runtime_error("Failed to open genome file: " + filepath);
    }
    if (!getline(file, header)) {
        file.close();
        throw runtime_error("Empty genome file");
    }
    while (getline(file, line)) {
        if (line.empty() || line[0] == '>') continue;
        sequence += line;
    }

    file.close();

    if (sequence.empty()) {
        throw runtime_error("No sequence data found in file");
    }

    return sequence.length(); 
}

void suffixTree::savinggenome()
{
    ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open genome file: " + filepath);
    }

    string header;
    string sequence;
    string line;

    if (!getline(file, header)) {
        file.close();
        throw runtime_error("Empty genome file");
    }
    while (getline(file, line)) {
        if (line.empty() || line[0] == '>') continue;
        sequence += line;
    }

    file.close();

    if (sequence.empty()) {
        throw runtime_error("No sequence data found in file");
    }

    if (sequence.length() + 1 != genomesize) {
        genomesize = sequence.length() + 1;
        delete[] genome;
        genome = new char[genomesize];
    }

    for (int i = 0; i < genomesize - 1; ++i) {
        genome[i] = sequence[i];
    }

    genome[genomesize - 1] = '$';
}


void suffixTree::createsuffixtree()
{
    for (int i = 0; i < genomesize; i++)
    {
        Node* curr = root;
        for (int j = i; j < genomesize; j++)
        {
            char currchar = genome[j];
            Node** child = nullptr;

            switch (currchar)
            {
            case 'A': child = &curr->A; break;
            case 'C': child = &curr->C; break;
            case 'G': child = &curr->G; break;
            case 'T': child = &curr->T; break;
            case '$': child = &curr->$; break;
            default:
                throw runtime_error("Invalid character in genome at position " +
                    to_string(j) + ": '" + string(1, currchar) + "'");
            }

            if (*child == nullptr)
            {
                *child = new Node();
                (*child)->A = nullptr;
                (*child)->C = nullptr;
                (*child)->G = nullptr;
                (*child)->T = nullptr;
                (*child)->$ = nullptr;
                nodecount++;
            }
            curr = *child;
        }
    }
    cout << "Suffix tree created with " << nodecount << " nodes." << endl;
}

bool suffixTree::searchsuffix(char* inputQuery, int inputQueryLength)
{
    if (inputQuery == nullptr || inputQueryLength <= 0 || root == nullptr) {
        return false;
    }

    Node* currentNode = root;
    int remainlength = inputQueryLength;
    char* currentquery = inputQuery;
    int matched = 0;

    while (matched < inputQueryLength && remainlength > 0)
    {
        char currentchar = currentquery[0];
        Node* nextNode = nullptr;

        switch (currentchar)
        {
        case 'A': nextNode = currentNode->A; break;
        case 'C': nextNode = currentNode->C; break;
        case 'G': nextNode = currentNode->G; break;
        case 'T': nextNode = currentNode->T; break;
        case '$': nextNode = currentNode->$; break;
        default:
            cerr << "Invalid character in query at position " << matched
                << ": '" << currentchar << "'" << endl;
            return false;
        }

        if (nextNode == nullptr)
        {
            return false; 
        }

        matched++;
        remainlength--;
        currentquery++;
        currentNode = nextNode;
    }

    return (matched == inputQueryLength);
}


char** suffixTree::generate_random_kmers_from_genome(int readCount, int readLength)
{
    if (genome == nullptr || readCount <= 0 || readLength <= 0 ||
        readLength >= genomesize - 1) {
        return nullptr;
    }

    char** fragments = new char* [readCount];

    for (int i = 0; i < readCount; i++) {
        fragments[i] = new char[readLength + 1]; 

        int startingposition = rand() % (genomesize - readLength - 1);

        for (int j = 0; j < readLength; j++) {
            fragments[i][j] = genome[startingposition + j];
        }

        fragments[i][readLength] = '\0';
    }

    return fragments;
}

int suffixTree::search_suffix_for_kmers(int readCount, int readLength)
{
    char** fragments = generate_random_kmers_from_genome(readCount, readLength);

    if (fragments == nullptr) {
        return 0;
    }

    int matchCount = 0;

    for (int i = 0; i < readCount; i++) {
        if (searchsuffix(fragments[i], readLength)) {
            matchCount++;
        }

        delete[] fragments[i];
    }

    delete[] fragments;

    return matchCount;
}