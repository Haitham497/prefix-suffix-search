#include "suffixTree.h"
#include <iostream>
#include <chrono>

int main() {
    try {
        suffixTree fileTree(R"(C:\\Users\\ZC\\Downloads\\genomic2.txt)");
        cout << "Genome loaded from file and tree initialized.\n";

        auto start_build = chrono::high_resolution_clock::now();
        fileTree.createsuffixtree();
        auto end_build = chrono::high_resolution_clock::now();
        cout << "Suffix tree built.\n";
        cout << "Time to build suffix tree: "
            << chrono::duration<double>(end_build - start_build).count()
            << " seconds\n";

        int readCount = 5;
        int readLength = 36;

        char** reads = fileTree.generate_random_kmers_from_genome(readCount, readLength);
        cout << "Selected " << readCount << " random reads of length " << readLength << ":\n";

        std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < readCount; ++i) {
            cout << reads[i] << "\n";
            bool found = fileTree.searchsuffix(reads[i], readLength);
            cout << "Read " << i + 1 << (found ? " found" : " not found") << " in suffix tree.\n";
            delete[] reads[i];
        }
        delete[] reads;
        std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = std::chrono::duration<double>(end - start);
        std::cout << "Time to search all reads: " << elapsed.count() << " seconds\n";


    }
    catch (const exception& ex) {
        cerr << "Exception: " << ex.what() << endl;
    }

    return 0;
}
