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
    duration<double> duration1, duration2, duration3, duration4;

    PBasedQAP::PBasedALG myGeneticALG(argv[1], atoi(argv[2]), 50);

    // AGG with position cross operator
    tBefore = chrono::high_resolution_clock::now();
    PBasedQAP::QAP_solution S1 = myGeneticALG.AGG_pos();
    tAfter = chrono::high_resolution_clock::now();
    duration1 = duration_cast<duration<double> >(tAfter - tBefore);

    // AGG with OX cross operator
    tBefore = chrono::high_resolution_clock::now();
    PBasedQAP::QAP_solution S2 = myGeneticALG.AGG_OX();
    tAfter = chrono::high_resolution_clock::now();
    duration2 = duration_cast<duration<double> >(tAfter - tBefore);

    // AGE with position cross operator
    tBefore = chrono::high_resolution_clock::now();
    PBasedQAP::QAP_solution S3 = myGeneticALG.AGE_pos();
    tAfter = chrono::high_resolution_clock::now();
    duration3 = duration_cast<duration<double> >(tAfter - tBefore);

    // AGE with OX cross operator
    tBefore = chrono::high_resolution_clock::now();
    PBasedQAP::QAP_solution S4 = myGeneticALG.AGE_OX();
    tAfter = chrono::high_resolution_clock::now();
    duration4 = duration_cast<duration<double> >(tAfter - tBefore);

    if(atoi(argv[3]) == 1){
        cout << "[+] AGG position solution for " << argv[1] << endl << "(";
        for(unsigned int i=0; i<S1.S.size(); ++i)
            cout << S1.S[i] << ", ";

        cout << ")" << endl;
        cout << fixed << "Solution cost: " << S1.cost << endl; 
        cout << fixed << "Elapsed Time: " << duration1.count() << " s" << endl << endl;

        cout << "[+] AGG OX solution for " << argv[1] << endl << "(";
        for(unsigned int i=0; i<S2.S.size(); ++i)
            cout << S2.S[i] << ", ";

        cout << ")" << endl;
        cout << fixed << "Solution cost: " << S2.cost << endl; 
        cout << fixed << "Elapsed Time: " << duration2.count() << " s" << endl << endl;

        cout << "[+] AGE position solution for " << argv[1] << endl << "(";
        for(unsigned int i=0; i<S3.S.size(); ++i)
            cout << S3.S[i] << ", ";

        cout << ")" << endl;
        cout << fixed << "Solution cost: " << S3.cost << endl; 
        cout << fixed << "Elapsed Time: " << duration3.count() << " s" << endl << endl;

        cout << "[+] AGE OX solution for " << argv[1] << endl << "(";
        for(unsigned int i=0; i<S4.S.size(); ++i)
            cout << S4.S[i] << ", ";

        cout << ")" << endl;
        cout << fixed << "Solution cost: " << S4.cost << endl; 
        cout << fixed << "Elapsed Time: " << duration4.count() << " s" << endl << endl;

    } else {
        cout << argv[1] << endl;
        cout << fixed << S1.cost << "\t" << duration1.count()  << "\t" 
             << S2.cost << "\t" << duration2.count()  << "\t"
             << S3.cost << "\t" << duration3.count()  << "\t"
             << S4.cost << "\t" << duration4.count()  << endl;
    }
    
}