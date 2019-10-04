///////////////////////////////////////////////////////////////////////////////////
// Project: Metaheurísticas: Práctica 3, Técnicas basadas en trayectorias:
// Choosen problem: QAP
// Author: Álvaro Fernández García
// Date: May, 2018
// Compilation: g++ -o QAP -I. *.cpp -std=c++11 -O3 -march=native
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
#include <cmath>
#include "random.h"

using namespace std;
using namespace std::chrono;


///////////////////////////////////////////////////////////////////////////////////
// Algorithms parameters:
const float ALPHA = 0.3;             //For quality threshold in Randomized Greedy.
const int MAX_EVAL_LS = 50000;       //LS max evaluations
const int MAX_ITER_GRASP = 25;       //GRASP max iterations
const int MAX_ITER_BMB = 25;         //BMB max iterations
const double MU = 0.3;               //For T_0 in ES
const double FI = 0.3;               //For T_0 in ES
const float K = 1;                   //For metropoli equation in ES
const int MAX_ITER_ILS = 24;         //ILS max iterations
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


// Generate a random solution for QAP problem:
void GenerateRandom(QAP_solution &S, const _fmatrix &F, const _fmatrix &D){
    int n = F[0].size(); //Get problem size
    vector<int> aux;
    int lim = n-1;
    int index;

    for(int i=0; i<n; ++i)
        aux.push_back(i);

    S.S.resize(n);

    for(int i=0; i<n; ++i){
        index = Randint(0, lim);
        S.S[i] = aux[index];
        auto it = find(aux.begin(), aux.end(), aux[index]);
        aux.erase(it);
        --lim;
    }

    QAP_cost(S, F, D);
}


///////////////////////////////////////////////////////////////////////////////////
// Local Search functions.
///////////////////////////////////////////////////////////////////////////////////


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
QAP_solution LS(QAP_solution S, const _fmatrix &F, const _fmatrix &D){
    int n = S.S.size();
    bool improveFlag;
    bool _continue = true;
    float var;
    int evaluations = 0;

    // Define the 'Don't look Bits' array:
    vector<bool> DLB;
    for(int i=0; i<n; ++i)
        DLB.push_back(false);

    // Only 50000 evaluations or while there is improve
    do{
        _continue = false;

        for(int i=0; i<n && evaluations < MAX_EVAL_LS; ++i){
        
            if(DLB[i] == false){
                improveFlag = false;
        
                for(int j=0; j<n && evaluations < MAX_EVAL_LS; ++j){
                    var = DeltaCost(S, i, j, F, D); 
                    evaluations++; 
        
                    if(var < 0){  //If cost of new neighbour improves
                        ApplyMove(S, i, j, var);
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
    } while(evaluations < MAX_EVAL_LS && _continue);

    return S;

}


///////////////////////////////////////////////////////////////////////////////////
// Functions for GRASP
///////////////////////////////////////////////////////////////////////////////////


// Randomized Greedy for GRASP
QAP_solution RandomizedGreedy(const _fmatrix &F, const _fmatrix &D){
    int n = F[0].size();  // Get problem size
    vector<float> _fp, _dp;  //Vectors who will contain potential
    float sumf, sumd, minf, maxf, mind, maxd, sum;
    float thresholdF, thresholdD;
    vector< pair<float, pair<int,int> > > C;
    
    QAP_solution out;
    out.S.resize(n);

    //***************************************************
    // STAGE 1:
    //***************************************************

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

    // Find max and min values:
    minf = *min_element(_fp.begin(), _fp.end());
    maxf = *max_element(_fp.begin(), _fp.end());
    mind = *min_element(_dp.begin(), _dp.end());
    maxd = *max_element(_dp.begin(), _dp.end());

    //Compute threshold:
    thresholdF = maxf - ALPHA * (maxf - minf);
    thresholdD = mind + ALPHA * (maxd - mind);

    //Create LCR:
    vector<int> LCRf, LCRd;

    for(int i=0; i<n; ++i){
        if(_fp[i] >= thresholdF)
            LCRf.push_back(i);
        if(_dp[i] <= thresholdD)
            LCRd.push_back(i);
    }

    //Generate random numbers:
    int ub1, ub2, loc1, loc2;
    
    ub1 = Randint(0, LCRf.size()-1);
    loc1 = Randint(0, LCRd.size()-1);
    out.S[ub1] = loc1;

    do{
        ub2 = Randint(0, LCRf.size()-1);
    } while(ub1 == ub2);

    do{
        loc2 = Randint(0, LCRd.size()-1);
    } while(loc1 == loc2);

    out.S[ub2] = loc2;

    //***************************************************
    // STAGE 2:
    //***************************************************

    vector<pair<int,int> > S, LC, LCR;

    //Temporal S:
    S.push_back(pair<int,int>(ub1, loc1));
    S.push_back(pair<int,int>(ub2, loc2));

    //Build initial candidate list:
    for(int i=0; i<n; ++i){
        if(i != ub1 && i != ub2){
            for(int j=0; j<n; ++j)
                if(j != loc1 && j != loc2)
                    LC.push_back(pair<int,int>(i,j));

        }
    }

    //Start greedy selection:
    for(int z=0; z<n-2; ++z){

        //Clear temporal lists:
        C.clear();
        LCR.clear();

        //Compute C matrix
        for(auto LCit=LC.begin(); LCit!=LC.end(); ++LCit){
            sum = 0.0;
            for(auto Sit=S.begin(); Sit!=S.end(); ++Sit)
                sum += F[LCit->first][Sit->first] * D[LCit->second][Sit->second];

            C.push_back(pair<float, pair<int,int> >(sum, pair<int,int>(LCit->first,LCit->second)));
        }

        //Find max and min element:
        float max = numeric_limits<float>::min();
        float min = numeric_limits<float>::max();

        for(unsigned int i=0; i<C.size(); ++i){
                if(C[i].first > max)
                    max = C[i].first;
                if(C[i].first < min)
                    min = C[i].first;
        }

        //Compute threshold:
        float mu = min + ALPHA * (max - min);

        //Construct LCR:
        LCR.clear();
        for(unsigned int i=0; i<C.size(); ++i)
            if(C[i].first <= mu)
                LCR.push_back(pair<int,int>(C[i].second.first, C[i].second.second));
            

        //Select random element:
        int element = Randint(0, LCR.size()-1);
        pair<int,int> selected = LCR[element];

        //Add to S:
        S.push_back(selected);
        out.S[selected.first] = selected.second;
        
        //Remove from candidates:
        auto it = LC.begin();
        while(it!=LC.end()){
            if(it->first == selected.first || it->second == selected.second)
                LC.erase(it);
            else
                ++it;
        }        
    }

    //Cost:
    QAP_cost(out, F, D);

    return out;
}


//GRASP for QAP:
QAP_solution GRASP(const _fmatrix &F, const _fmatrix &D){
    QAP_solution best, SGreedy, S;
    best.cost = numeric_limits<float>::max();

    for(int i=0; i<MAX_ITER_GRASP; ++i){
        SGreedy = RandomizedGreedy(F, D);
        S = LS(SGreedy, F, D);

        if(S.cost < best.cost)
            best = S;
    }

    return best;
}


///////////////////////////////////////////////////////////////////////////////////
// BMB:
///////////////////////////////////////////////////////////////////////////////////


// BMB for QAP:
QAP_solution BMB(const _fmatrix &F, const _fmatrix &D){
    QAP_solution best, SRandom, S;
    best.cost = numeric_limits<float>::max();

    for(int i=0; i<MAX_ITER_BMB; ++i){
        GenerateRandom(SRandom, F, D);
        S = LS(SRandom, F, D);

        if(S.cost < best.cost)
            best = S;
    }

    return best;
}


///////////////////////////////////////////////////////////////////////////////////
// Simulated annealing
///////////////////////////////////////////////////////////////////////////////////


//Neighbour operator for simulated annealing:
QAP_solution NeighbourOP(QAP_solution S, const _fmatrix &F, const _fmatrix &D){
    int i, j, aux, n = S.S.size();
    QAP_solution out;

    i = Randint(0, n-1);
    do{
        j = Randint(0, n-1);
    }while(i==j);

    aux = S.S[i];
    S.S[i] = S.S[j];
    S.S[j] = aux;

    out = S;
    QAP_cost(out, F, D);
    return out;
}


// ES for QAP:
QAP_solution ES(QAP_solution S, const _fmatrix &F, const _fmatrix &D, bool cauchy=true){
    int n = F[0].size();
    QAP_solution SN, best;
    float T, T0, deltaF, BETA;
    int currentSuccesses, currentNeighbours;
    
    //Compute L(T), T0 and BETA:
    const int MAX_NEIGHBOURS = 10 * n;
    const int MAX_SUCCESSES = round(0.1 * MAX_NEIGHBOURS);

    best = S;

    float TF = 0.001;
    T0 = (MU * S.cost) / (-log(FI));  //T_0
    
    if(T0 < TF)
        TF = T0 * 0.001;

    if(cauchy)
        BETA = (T0 - TF) / ( (50000/MAX_NEIGHBOURS) * T0 * TF ); 
    
    T = T0;

    do{
        currentNeighbours = 0;
        currentSuccesses = 0;

        while(currentSuccesses < MAX_SUCCESSES && currentNeighbours < MAX_NEIGHBOURS){
            SN = NeighbourOP(S, F, D);
            currentNeighbours++;
            deltaF = SN.cost - S.cost;

            if( (deltaF < 0) || (Rand() < exp((-deltaF) / (K * T))) ){
                S = SN;
                currentSuccesses++;
                if(S.cost < best.cost)
                    best = S;
            }
        }

        if(currentSuccesses == 0)
            break;

        //Annealing
        if(cauchy)
            T = T / (1 + BETA * T);
        else
            T = 0.99 * T;

    } while(T > TF);

    return best;
}


///////////////////////////////////////////////////////////////////////////////////
// Iterative Local Search
///////////////////////////////////////////////////////////////////////////////////


//Mutation operator for ILS:
QAP_solution Mutate(QAP_solution S, const _fmatrix &F, const _fmatrix &D){
    int n = S.S.size();
    QAP_solution out;
    vector<int> sublist;

    int index = Randint(0, (3/4)*n);
    int offset = n/4;

    for(int i=index; i<=index+offset; ++i)
        sublist.push_back(S.S[i]);

    random_shuffle(sublist.begin(), sublist.end());

    int k = 0;
    for(int i=index; i<=index+offset; ++i)
        S.S[i] = sublist[k++];

    QAP_cost(S, F, D);
    out = S;
    return out;
}


// ILS for QAP:
QAP_solution ILS(const _fmatrix &F, const _fmatrix &D){
    QAP_solution S, best;

    GenerateRandom(S, F, D);
    S = LS(S, F, D);  //Optimize
    best = S;

    for(int i=0; i<MAX_ITER_ILS; ++i){
        S = Mutate(best, F, D);
        S = LS(S, F, D);
        if(S.cost < best.cost)
            best = S;
    }
    return best;
}


// ILS-ES for QAP:
QAP_solution ILSES(const _fmatrix &F, const _fmatrix &D){
    QAP_solution S, best;

    GenerateRandom(S, F, D);
    S = ES(S, F, D);  //Optimize
    best = S;

    for(int i=0; i<MAX_ITER_ILS; ++i){
        S = Mutate(best, F, D);
        S = ES(S, F, D);
        if(S.cost < best.cost)
            best = S;
    }
    return best;
}

///////////////////////////////////////////////////////////////////////////////////


int main(int argc, char *argv[]){
    if(argc != 4){
        cerr << "[!] Invalid syntax. Use: " << argv[0] << " path/file.dat seed mode" << endl;
        return 1;
    }

    Set_random(atoi(argv[2]));  //Set random Seed.
    srand(unsigned(atoi(argv[2])));

    chrono::high_resolution_clock::time_point tBefore, tAfter;
    duration<double> durationGRASP, durationBMB, durationES1, durationES2, durationILS, durationILSES;

    _fmatrix F, D;
    load(argv[1], F, D);
 
    //GRASP
    tBefore = chrono::high_resolution_clock::now();
    QAP_solution S = GRASP(F, D);
    tAfter = chrono::high_resolution_clock::now();
    durationGRASP = duration_cast<duration<double> >(tAfter - tBefore);

    //BMB
    tBefore = chrono::high_resolution_clock::now();
    QAP_solution S2 = BMB(F, D);
    tAfter = chrono::high_resolution_clock::now();
    durationBMB = duration_cast<duration<double> >(tAfter - tBefore);

    //ES Cauchy:
    QAP_solution aux;
    GenerateRandom(aux, F, D);
    tBefore = chrono::high_resolution_clock::now();
    QAP_solution S3 = ES(aux, F, D);
    tAfter = chrono::high_resolution_clock::now();
    durationES1 = duration_cast<duration<double> >(tAfter - tBefore);

    //ES proportional:
    aux;
    GenerateRandom(aux, F, D);
    tBefore = chrono::high_resolution_clock::now();
    QAP_solution S5 = ES(aux, F, D, false);
    tAfter = chrono::high_resolution_clock::now();
    durationES2 = duration_cast<duration<double> >(tAfter - tBefore);

    //ILS
    tBefore = chrono::high_resolution_clock::now();
    QAP_solution S4 = ILS(F, D);
    tAfter = chrono::high_resolution_clock::now();
    durationILS = duration_cast<duration<double> >(tAfter - tBefore);

    //ILS-ES
    tBefore = chrono::high_resolution_clock::now();
    QAP_solution S6 = ILSES(F, D);
    tAfter = chrono::high_resolution_clock::now();
    durationILSES = duration_cast<duration<double> >(tAfter - tBefore);

    if(atoi(argv[3]) == 1){
        cout << "[+] GRASP solution for " << argv[1] << endl << "(";
        for(int i=0; i<S.S.size(); ++i)
            cout << S.S[i] << ", ";

        cout << ")" << endl;
        cout << fixed << "Solution cost: " << S.cost << endl; 
        cout << fixed << "Elapsed Time: " << durationGRASP.count() << " s" << endl << endl;

        cout << "[+] BMB solution for " << argv[1] << endl << "(";
        for(int i=0; i<S2.S.size(); ++i)
            cout << S2.S[i] << ", ";

        cout << ")" << endl;
        cout << fixed << "Solution cost: " << S2.cost << endl; 
        cout << fixed << "Elapsed Time: " << durationBMB.count() << " s" << endl << endl;

        cout << "[+] ES CAUCHY solution for " << argv[1] << endl << "(";
        for(int i=0; i<S3.S.size(); ++i)
            cout << S3.S[i] << ", ";

        cout << ")" << endl;
        cout << fixed << "Solution cost: " << S3.cost << endl; 
        cout << fixed << "Elapsed Time: " << durationES1.count() << " s" << endl << endl;

        cout << "[+] ES PROPORTIONAL solution for " << argv[1] << endl << "(";
        for(int i=0; i<S5.S.size(); ++i)
            cout << S5.S[i] << ", ";

        cout << ")" << endl;
        cout << fixed << "Solution cost: " << S5.cost << endl; 
        cout << fixed << "Elapsed Time: " << durationES2.count() << " s" << endl << endl;

        cout << "[+] ILS solution for " << argv[1] << endl << "(";
        for(int i=0; i<S4.S.size(); ++i)
            cout << S4.S[i] << ", ";

        cout << ")" << endl;
        cout << fixed << "Solution cost: " << S4.cost << endl; 
        cout << fixed << "Elapsed Time: " << durationILS.count() << " s" << endl << endl;

        cout << "[+] ILS-ES solution for " << argv[1] << endl << "(";
        for(int i=0; i<S6.S.size(); ++i)
            cout << S6.S[i] << ", ";

        cout << ")" << endl;
        cout << fixed << "Solution cost: " << S6.cost << endl; 
        cout << fixed << "Elapsed Time: " << durationILSES.count() << " s" << endl << endl;

    } else {
        cout << argv[1] << endl;
        cout << fixed << S.cost << " " << durationGRASP.count() << " "
                      << S2.cost << " " << durationBMB.count() << " " 
                      << S3.cost << " " << durationES1.count() << " "
                      << S5.cost << " " << durationES2.count() << " "
                      << S4.cost << " " << durationILS.count() << " "
                      << S6.cost << " " << durationILSES.count() << endl;
    }
}