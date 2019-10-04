#include <iostream>
#include <chrono>
#include <ctime>
#include "PBasedALG.h"

using namespace std;
using namespace std::chrono;

int main(int argc, char *argv[]){
    if(argc != 4){
        cerr << "[!] Invalid syntax. Use: " << argv[0] << " path/file.dat seed mode" << endl;
        return 1;
    }

    chrono::high_resolution_clock::time_point tBefore, tAfter;
    duration<double> duration5, duration6, duration7;

    PBasedQAP::PBasedALG myMemeticALG(argv[1], atoi(argv[2]), 10);

    // AM with pLS = 1
    tBefore = chrono::high_resolution_clock::now();
    PBasedQAP::QAP_solution S5 = myMemeticALG.AM(10, 1.0);
    tAfter = chrono::high_resolution_clock::now();
    duration5 = duration_cast<duration<double> >(tAfter - tBefore);

    // AM with pLS = 0.1
    tBefore = chrono::high_resolution_clock::now();
    PBasedQAP::QAP_solution S6 = myMemeticALG.AM(10, 0.1);
    tAfter = chrono::high_resolution_clock::now();
    duration6 = duration_cast<duration<double> >(tAfter - tBefore);

    // AM best 0.1
    tBefore = chrono::high_resolution_clock::now();
    PBasedQAP::QAP_solution S7 = myMemeticALG.AM_best(10, 0.1);
    tAfter = chrono::high_resolution_clock::now();
    duration7 = duration_cast<duration<double> >(tAfter - tBefore);

    if(atoi(argv[3]) == 1){
        cout << "[+] AM pLS=1 solution for " << argv[1] << endl << "(";
        for(unsigned int i=0; i<S5.S.size(); ++i)
            cout << S5.S[i] << ", ";

        cout << ")" << endl;
        cout << fixed << "Solution cost: " << S5.cost << endl; 
        cout << fixed << "Elapsed Time: " << duration5.count() << " s" << endl << endl;

        cout << "[+] AM pLS=0.1 solution for " << argv[1] << endl << "(";
        for(unsigned int i=0; i<S6.S.size(); ++i)
            cout << S6.S[i] << ", ";

        cout << ")" << endl;
        cout << fixed << "Solution cost: " << S6.cost << endl; 
        cout << fixed << "Elapsed Time: " << duration6.count() << " s" << endl << endl;

        cout << "[+] AM BEST 0.1 solution for " << argv[1] << endl << "(";
        for(unsigned int i=0; i<S7.S.size(); ++i)
            cout << S7.S[i] << ", ";

        cout << ")" << endl;
        cout << fixed << "Solution cost: " << S7.cost << endl; 
        cout << fixed << "Elapsed Time: " << duration7.count() << " s" << endl << endl;

    } else {
        cout << argv[1] << endl;
        cout << fixed << S5.cost << "\t" << duration5.count()  << "\t" 
             << S6.cost << "\t" << duration6.count()  << "\t"
             << S7.cost << "\t" << duration7.count()  << endl;
    }
    
}