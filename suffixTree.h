#pragma once
#include <cmath>
#include <ctime>
#include <stack>
#include <chrono>
#include <random>
#include <string>
#include <cstring>
#include <fstream>
#include <iostream>

using namespace std;

struct Node
{
	Node* A; Node* C; Node* G; Node* T; Node* $;
};

class suffixTree
{
private:
	Node* root;
	char* genome;
	int genomesize, nodecount;
	ifstream file;
	string filepath;

	void savinggenome();
	int grablength();

public:
	suffixTree();
	suffixTree(string inputGenome);
	~suffixTree();
	void createsuffixtree();
	bool searchsuffix(char* inputQuery, int inputQueryLength);
	char** generate_random_kmers_from_genome(int readCount, int readLength);
	int search_suffix_for_kmers(int readCount, int readLength);
};