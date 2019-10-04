///////////////////////////////////////////////////////////////////////////////////
// Project: Metaheurísticas: Práctica 1, Búsqueda Local y Greedy
// Choosen problem: QAP
// Author: Álvaro Fernández García
// Date: March, 2018
// Compilation: g++ -o QAP -I. *.cpp -std=c++11
///////////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <climits>
#include <cstdlib>
#include <chrono>
#include <ctime>
#include "random.h"

using namespace std;
using namespace std::chrono;

///////////////////////////////////////////////////////////////////////////////////

// Define a short name for a float matrix:
typedef vector< vector<float> > _fmatrix;

// Create a QAP solutión: Defines a permutation where S[i] = j,
// means that i unit is located on j location:
struct QAP_solution{
    vector<int> S;
    float cost;
};

///////////////////////////////////////////////////////////////////////////////////

// Load .dat file:
void load(const string path, _fmatrix &F, _fmatrix &D){
    ifstream input;
    int n;
    float val;

    input.open(path.c_str());

    if(input.is_open()){
        input >> n;  //Read problem size
        F.resize(n);  //Resize rows

        for(int i=0; i<n; ++i){
            F[i].resize(n);  //Resize columns
            for(int j=0; j<n; ++j){
                input >> val;
                F[i][j] = val;
            }
        }
    
        D.resize(n);  //Resize rows

        for(int i=0; i<n; ++i){
            D[i].resize(n);  //Resize columns
            for(int j=0; j<n; ++j){
                input >> val;
                D[i][j] = val;
            }
        }
    } else{
        cerr << "[!] Error opening file " << path << endl;
        exit(1);
    }
}

///////////////////////////////////////////////////////////////////////////////////

// Computes QAP cost solution:
void QAP_cost(QAP_solution &S, const _fmatrix &F, const _fmatrix &D){
    float cost = 0.0;
    int n = S.S.size();

    for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
            cost += F[i][j] * D[S.S[i]][S.S[j]];

    S.cost = cost; 
}

///////////////////////////////////////////////////////////////////////////////////

// A Greedy based algorithm for QAP problem:
QAP_solution GQAP(const _fmatrix &F, const _fmatrix &D){
    int n = F[0].size();  // Get problem size
    vector<float> _fp, _dp;  //Vectors who will contain potential
    float sumf, sumd, min, max;
    int posMax, posMin;
    QAP_solution out;
    out.S.resize(n);

    // Compute potential:
    for(int i=0; i<n; ++i){
        sumf = sumd = 0.0;
        for(int j=0; j<n; ++j){
            sumf += F[i][j];
            sumd += D[i][j];
        }
        _fp.push_back(sumf);
        _dp.push_back(sumd);
    }

    // Start Greedy:
    for(int i=0; i<n; ++i ){
        max = numeric_limits<float>::min();
        min = numeric_limits<float>::max();
        for(int j=0; j<n; ++j){
            if(_fp[j] > max){
                max = _fp[j];
                posMax = j;
            }
            if(_dp[j] < min){
                min = _dp[j];
                posMin = j;
            }
        }
        out.S[posMax] = posMin;
        //"Remove" candidates: assing a value that ensures it never gonna be chosen again
        _fp[posMax] = numeric_limits<float>::min();
        _dp[posMin] = numeric_limits<float>::max();
    }

    QAP_cost(out, F, D);
    return out;
}

///////////////////////////////////////////////////////////////////////////////////
// Local Search functions.
///////////////////////////////////////////////////////////////////////////////////

// Generate a random solution for QAP problem:
void GenerateRandom(QAP_solution &S){
    int n = S.S.size();
    vector<int> aux;
    int lim = n-1;
    int index;

    for(int i=0; i<n; ++i)
        aux.push_back(i);

    for(int i=0; i<n; ++i){
        index = Randint(0, lim);
        S.S[i] = aux[index];
        auto it = find(aux.begin(), aux.end(), aux[index]);
        aux.erase(it);
        --lim;
    }
}


// Factorice cost:
float DeltaCost(const QAP_solution& pi, int r, int s, const _fmatrix &F, const _fmatrix &D){
    float variation = 0.0;
    int n = pi.S.size();

    for(int k=0; k<n; ++k){
        if(k != r && k != s){
            variation += F[r][k] * (D[pi.S[s]][pi.S[k]] - D[pi.S[r]][pi.S[k]]) + 
                         F[s][k] * (D[pi.S[r]][pi.S[k]] - D[pi.S[s]][pi.S[k]]) +
                         F[k][r] * (D[pi.S[k]][pi.S[s]] - D[pi.S[k]][pi.S[r]]) + 
                         F[k][s] * (D[pi.S[k]][pi.S[r]] - D[pi.S[k]][pi.S[s]]);
        }
    }

    return variation;
}


// Update the current solution:
void ApplyMove(QAP_solution& S, int i, int j, float var){
    int aux;

    aux = S.S[i];
    S.S[i] = S.S[j];
    S.S[j] = aux;

    S.cost += var;
}


// Local Search implementation:
QAP_solution LSQAP(const _fmatrix &F, const _fmatrix &D){
    int n = F[0].size();
    QAP_solution out;
    out.S.resize(n);
    bool improveFlag;
    bool _continue = true;
     float var;
    
    //Generate random initial solution
    GenerateRandom(out);
    QAP_cost(out, F, D);

    // Define the 'Don't look Bits' array:
    vector<bool> DLB;
    for(int i=0; i<n; ++i)
        DLB.push_back(false);

    // Only 50000 iterations or while there is improve
    for(int k=0; k<50000 && _continue; ++k){
        _continue = false;

        for(int i=0; i<n; ++i){
        
            if(DLB[i] == false){
                improveFlag = false;
        
                for(int j=0; j<n; ++j){
                    var = DeltaCost(out, i, j, F, D);  
        
                    if(var < 0){  //If cost of new neighbour improves
                        ApplyMove(out, i, j, var);
                        DLB[i] = false;
                        DLB[j] = false;
                        improveFlag = true;
                        _continue = true;
                        break;  //First Better
                    }
                }
                if(!improveFlag)  // if there is no improve
                    DLB[i] = true;
                else
                    break;  //First Better
            }
        }
    }

    return out;

}

///////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]){
    if(argc != 4){
        cerr << "[!] Invalid syntax. Use: " << argv[0] << " path/file.dat seed mode" << endl;
        return 1;
    }

    Set_random(atoi(argv[2]));  //Set random Seed.

    chrono::high_resolution_clock::time_point tBefore, tAfter;
    duration<double> durationG, durationLS;

    _fmatrix F, D;
    load(argv[1], F, D);
 
    tBefore = chrono::high_resolution_clock::now();
    QAP_solution S = GQAP(F, D);
    tAfter = chrono::high_resolution_clock::now();
    durationG = duration_cast<duration<double> >(tAfter - tBefore);

    tBefore = chrono::high_resolution_clock::now();
    QAP_solution S2 = LSQAP(F, D);
    tAfter = chrono::high_resolution_clock::now();
    durationLS = duration_cast<duration<double> >(tAfter - tBefore);

    if(atoi(argv[3]) == 1){
        cout << "[+] Greedy solution for " << argv[1] << endl << "(";
        for(int i=0; i<S.S.size(); ++i)
            cout << S.S[i] << ", ";

        cout << ")" << endl;
        cout << fixed << "Solution cost: " << S.cost << endl; 
        cout << fixed << "Elapsed Time: " << durationG.count() << " s" << endl << endl;

        cout << "[+] Local Search solution for " << argv[1] << endl << "(";
        for(int i=0; i<S2.S.size(); ++i)
            cout << S2.S[i] << ", ";

        cout << ")" << endl;
        cout << fixed << "Solution cost: " << S2.cost << endl;
        cout << fixed << "Elapsed Time: " << durationLS.count() << " s" << endl << endl;
    } else {
        cout << argv[1] << endl;
        cout << fixed << S.cost << " " << durationG.count() << " " << S2.cost << " " << durationLS.count() << endl;
    }
}