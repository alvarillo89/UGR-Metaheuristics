#include "PBasedALG.h"

namespace PBasedQAP{

    void PBasedALG::load(const std::string path){
        std::ifstream input;
        float val;

        input.open(path.c_str());

        if(input.is_open()){
            input >> _sizeProblem;  //Read problem size
            _F.resize(_sizeProblem);  //Resize rows

            for(int i=0; i<_sizeProblem; ++i){
                _F[i].resize(_sizeProblem);  //Resize columns
                for(int j=0; j<_sizeProblem; ++j){
                    input >> val;
                    _F[i][j] = val;
                }
            }
        
            _D.resize(_sizeProblem);  //Resize rows

            for(int i=0; i<_sizeProblem; ++i){
                _D[i].resize(_sizeProblem);  //Resize columns
                for(int j=0; j<_sizeProblem; ++j){
                    input >> val;
                    _D[i][j] = val;
                }
            }
        } else{
            std::cerr << "[!] Error opening file " << path << std::endl;
            exit(1);
        }
    }


    void PBasedALG::EvaluateSolution(PBasedQAP::QAP_solution &S){
        float cost = 0.0;

        for(int i=0; i<_sizeProblem; ++i)
            for(int j=0; j<_sizeProblem; ++j)
                cost += _F[i][j] * _D[S.S[i]][S.S[j]];

        S.cost = cost; 
        S.needReevaluation = false;
        ++_currentEvaluations;
    }


    void PBasedALG::EvaluatePopulation(){
        for(int i=0; i<_POPULATION_SIZE; ++i)
            if(_population[i].needReevaluation)
                PBasedALG::EvaluateSolution(_population[i]);
    }


    void PBasedALG::GenerateRandomInitialPopulation(){
        PBasedQAP::QAP_solution S;
        std::vector<int> aux;
        int lim; 
        int index;

        S.S.resize(_sizeProblem);

        for(int j=0; j<_POPULATION_SIZE; ++j){
            aux.clear();
            lim = _sizeProblem-1;

            for(int i=0; i<_sizeProblem; ++i)
                aux.push_back(i);

            for(int i=0; i<_sizeProblem; ++i){
                index = Randint(0, lim);
                S.S[i] = aux[index];
                std::vector<int>::iterator it = std::find(aux.begin(), aux.end(), aux[index]);
                aux.erase(it);
                --lim;
            }

            PBasedALG::EvaluateSolution(S);
            _population.push_back(S);
        }
        _currentEvaluations = 0;
    }


    PBasedALG::PBasedALG(std::string path, int seed, int pSize) : _POPULATION_SIZE(pSize){
        Set_random(seed);
        std::srand(unsigned(seed));  // For random_shuffle
        this->PBasedALG::load(path);
    }


    void PBasedALG::Mutate(PBasedQAP::QAP_solution &S, int pos){
        int other, aux;

        do{
            other = Randint(0, _sizeProblem-1);
        }while(other == pos);

        aux = S.S[pos];
        S.S[pos] = S.S[other];
        S.S[other] = aux;

        S.needReevaluation = true;
    }


    std::vector<PBasedQAP::QAP_solution> PBasedALG::Select(int size){
        std::vector<PBasedQAP::QAP_solution> newPopulation;
        int firstC, secondC;

        for(int i=0; i<size; ++i){
            firstC = Randint(0, _POPULATION_SIZE-1);
            secondC = Randint(0, _POPULATION_SIZE-1);

            if(_population[firstC].cost < _population[secondC].cost)
                newPopulation.push_back(_population[firstC]);
            else
                newPopulation.push_back(_population[secondC]);
        }

        return newPopulation;
    }


    void PBasedALG::CrossPosition(PBasedQAP::QAP_solution &S1, PBasedQAP::QAP_solution &S2){
        PBasedQAP::QAP_solution H1, H2;
        std::vector<int> rest, restCopy;
        std::vector<bool> freePos;

        freePos.resize(_sizeProblem);
        H1.S.resize(_sizeProblem);
        H1.needReevaluation = true;
        H2.S.resize(_sizeProblem);
        H2.needReevaluation = true;

        for(int i=0; i<_sizeProblem; ++i){
            if(S1.S[i] == S2.S[i]){
                H1.S[i] = S1.S[i];
                H2.S[i] = S2.S[i];
                freePos[i] = false;
            }
            else{
                freePos[i] = true;
                rest.push_back(S1.S[i]);
            }
        }

        restCopy = rest;
        random_shuffle(rest.begin(), rest.end());
        random_shuffle(restCopy.begin(), restCopy.end());

        int k = 0;
        for(int i=0; i<_sizeProblem; ++i){
            if(freePos[i]){
                H1.S[i] = rest[k];
                H2.S[i] = restCopy[k];
                ++k;
            }
        }
    
        S1 = H1;
        S2 = H2;
    }


    void PBasedALG::CrossOX(PBasedQAP::QAP_solution &S1, PBasedQAP::QAP_solution &S2){
        PBasedQAP::QAP_solution H1, H2;
        int low, top, k;
        std::vector<bool> freePos(_sizeProblem, true);
        std::vector<int> tmp, rest1, rest2;
        
        low = Randint(1, _sizeProblem/2);
        top = Randint((_sizeProblem/2)+1, _sizeProblem-2);

        freePos.resize(_sizeProblem);
        H1.S.resize(_sizeProblem);
        H1.needReevaluation = true;
        H2.S.resize(_sizeProblem);
        H2.needReevaluation = true;

        for(int i=low; i<=top; ++i)
            freePos[i] = false;

        for(int i=0; i<_sizeProblem; ++i){
            if(!freePos[i]){
                H1.S[i] = S1.S[i];
                H2.S[i] = S2.S[i];
            } else{
                rest1.push_back(S1.S[i]);
                rest2.push_back(S2.S[i]);
            }
        }

        for(int i=0; i<_sizeProblem; ++i){
            if(std::find(rest1.begin(), rest1.end(), S2.S[i]) != rest1.end())
                tmp.push_back(S2.S[i]);
        }

        k = 0;
        for(int i=top+1; i<low || i>top; i=(i+1)%_sizeProblem)
            H1.S[i] = tmp[k++];

        tmp.clear();
        for(int i=0; i<_sizeProblem; ++i){
            if(std::find(rest2.begin(), rest2.end(), S1.S[i]) != rest2.end())
                tmp.push_back(S1.S[i]);
        }

        k = 0;
        for(int i=top+1; i<low || i>top; i=(i+1)%_sizeProblem)
            H2.S[i] = tmp[k++];
        
        S1 = H1;
        S2 = H2;
    }


    void PBasedALG::findBestSolution(){
        PBasedQAP::QAP_solution tmp;
        tmp.cost = std::numeric_limits<float>::max();
        
        for(int i=0; i<_POPULATION_SIZE; ++i)
            if(_population[i].cost < tmp.cost)
                tmp = _population[i];

        this->_bestSolution = tmp;
    }


    PBasedQAP::QAP_solution PBasedALG::AGG_pos(){
        std::vector<PBasedQAP::QAP_solution> intermediatePopulation;
        int expectedCrosses = ceilf((_POPULATION_SIZE/2.0) * _CROSS_P_G); 
        int expectedMutations = ceilf(_POPULATION_SIZE * _sizeProblem * _MUT_P_G);
        int k, chr, gen;

        _population.clear();

        this->PBasedALG::GenerateRandomInitialPopulation();
        this->PBasedALG::findBestSolution();  //Elitism

        while(_currentEvaluations < _MAX_EVALUATIONS){
            intermediatePopulation = PBasedALG::Select(_POPULATION_SIZE);

            //Cross:
            k = 0;
            for(int i=0; i<expectedCrosses; ++i){
                PBasedALG::CrossPosition(intermediatePopulation[k], intermediatePopulation[k+1]);
                k += 2;
            }

            //Mutate:
            for(int i=0; i<expectedMutations; ++i){
                chr = Randint(0, _POPULATION_SIZE-1);
                gen = Randint(0, _sizeProblem-1);
                PBasedALG::Mutate(intermediatePopulation[chr], gen);
            }

            //Replace:
            _population = intermediatePopulation;
            intermediatePopulation.clear();

            //Evaluate:
            this->PBasedALG::EvaluatePopulation();

            //Elitism:
            _population[_POPULATION_SIZE-1] = _bestSolution;
            PBasedALG::findBestSolution();
        }
        return _bestSolution;
    }


    PBasedQAP::QAP_solution PBasedALG::AGG_OX(){
        std::vector<PBasedQAP::QAP_solution> intermediatePopulation;
        int expectedCrosses = ceilf((_POPULATION_SIZE/2.0) * _CROSS_P_G); 
        int expectedMutations = ceilf(_POPULATION_SIZE * _sizeProblem * _MUT_P_G);
        int k, chr, gen;

        _population.clear();

        this->PBasedALG::GenerateRandomInitialPopulation();
        this->PBasedALG::findBestSolution();  //Elitism

        while(_currentEvaluations < _MAX_EVALUATIONS){
            intermediatePopulation = PBasedALG::Select(_POPULATION_SIZE);

            //Cross:
            k = 0;
            for(int i=0; i<expectedCrosses; ++i){
                PBasedALG::CrossOX(intermediatePopulation[k], intermediatePopulation[k+1]);
                k += 2;
            }

            //Mutate:
            for(int i=0; i<expectedMutations; ++i){
                chr = Randint(0, _POPULATION_SIZE-1);
                gen = Randint(0, _sizeProblem-1);
                PBasedALG::Mutate(intermediatePopulation[chr], gen);
            }

            //Replace:
            _population = intermediatePopulation;
            intermediatePopulation.clear();

            //Evaluate:
            this->PBasedALG::EvaluatePopulation();

            //Elitism:
            _population[_POPULATION_SIZE-1] = _bestSolution;
            PBasedALG::findBestSolution();
        }
        return _bestSolution;
    }


    PBasedQAP::QAP_solution PBasedALG::AGE_pos(){
        std::vector<PBasedQAP::QAP_solution> intermediatePopulation;
        int expectedMutations = ceilf(2 * _sizeProblem * _MUT_P_G);
        int chr, gen;

        _population.clear();
        this->PBasedALG::GenerateRandomInitialPopulation();

        while(_currentEvaluations < _MAX_EVALUATIONS){
            // Select two parents:
            intermediatePopulation = PBasedALG::Select(2);

            // Always cross
            PBasedALG::CrossPosition(intermediatePopulation[0], intermediatePopulation[1]);

            // Mutate:
            for(int i=0; i<expectedMutations; ++i){
                chr = Randint(0, 1);
                gen = Randint(0, _sizeProblem-1);
                PBasedALG::Mutate(intermediatePopulation[chr], gen);
            }

            //Evaluate:
            PBasedALG::EvaluateSolution(intermediatePopulation[0]);
            PBasedALG::EvaluateSolution(intermediatePopulation[1]);

            //Compete:
            _population.push_back(intermediatePopulation[0]);
            _population.push_back(intermediatePopulation[1]);
            std::sort(_population.begin(), _population.end());
            _population.pop_back();
            _population.pop_back();
            intermediatePopulation.clear();
        }

        PBasedALG::findBestSolution();
        return _bestSolution;    
    }


    PBasedQAP::QAP_solution PBasedALG::AGE_OX(){
        std::vector<PBasedQAP::QAP_solution> intermediatePopulation;
        int expectedMutations = ceilf(2 * _sizeProblem * _MUT_P_G);
        int chr, gen;

        _population.clear();
        this->PBasedALG::GenerateRandomInitialPopulation();

        while(_currentEvaluations < _MAX_EVALUATIONS){
            // Select two parents:
            intermediatePopulation = PBasedALG::Select(2);

            // Always cross
            PBasedALG::CrossOX(intermediatePopulation[0], intermediatePopulation[1]);

            // Mutate:
            for(int i=0; i<expectedMutations; ++i){
                chr = Randint(0, 1);
                gen = Randint(0, _sizeProblem-1);
                PBasedALG::Mutate(intermediatePopulation[chr], gen);
            }

            //Evaluate:
            PBasedALG::EvaluateSolution(intermediatePopulation[0]);
            PBasedALG::EvaluateSolution(intermediatePopulation[1]);

            //Compete:
            _population.push_back(intermediatePopulation[0]);
            _population.push_back(intermediatePopulation[1]);
            std::sort(_population.begin(), _population.end());
            _population.pop_back();
            _population.pop_back();
            intermediatePopulation.clear();
        }

        PBasedALG::findBestSolution();
        return _bestSolution;   
    }


    float PBasedALG::DeltaCost(const PBasedQAP::QAP_solution& pi, int r, int s){
        float variation = 0.0;
        int n = pi.S.size();

        for(int k=0; k<n; ++k){
            if(k != r && k != s){
                variation += _F[r][k] * (_D[pi.S[s]][pi.S[k]] - _D[pi.S[r]][pi.S[k]]) + 
                             _F[s][k] * (_D[pi.S[r]][pi.S[k]] - _D[pi.S[s]][pi.S[k]]) +
                             _F[k][r] * (_D[pi.S[k]][pi.S[s]] - _D[pi.S[k]][pi.S[r]]) + 
                             _F[k][s] * (_D[pi.S[k]][pi.S[r]] - _D[pi.S[k]][pi.S[s]]);
            }
        }
        _currentEvaluations++;
        return variation;
    }


    void PBasedALG::ApplyMove(PBasedQAP::QAP_solution& S, int i, int j, float var){
        int aux;

        aux = S.S[i];
        S.S[i] = S.S[j];
        S.S[j] = aux;

        S.cost += var;
    }


    void PBasedALG::LS(PBasedQAP::QAP_solution &S){
        bool improveFlag;
        bool _continue = true;
        float var;
        int localEvaluations = 0;
        
        // Define the 'Don't look Bits' array:
        std::vector<bool> DLB;
        for(int i=0; i<_sizeProblem; ++i)
            DLB.push_back(false);

        // Only 400 evaluations or while there is improve
        do{
            _continue = false;

            for(int i=0; i<_sizeProblem && localEvaluations < 400; ++i){
            
                if(DLB[i] == false){
                    improveFlag = false;
            
                    for(int j=0; j<_sizeProblem && localEvaluations < 400; ++j){
                        var = DeltaCost(S, i, j);
                        localEvaluations++;
            
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
        } while(localEvaluations < 400 && _continue);
    }


    PBasedQAP::QAP_solution PBasedALG::AM(int iters, float pLS){
        std::vector<PBasedQAP::QAP_solution> intermediatePopulation;
        int expectedCrosses = ceilf((_POPULATION_SIZE/2.0) * _CROSS_P_G); 
        int expectedMutations = ceilf(_POPULATION_SIZE * _sizeProblem * _MUT_P_G);
        int k, chr, gen;
        int timeout = 0;

        _population.clear();

        this->PBasedALG::GenerateRandomInitialPopulation();

        //Optimize population:
        for(int i=0; i<_POPULATION_SIZE; ++i){
            if(Randfloat(0,1) <= pLS)
                PBasedALG::LS(_population[i]);
        }

        PBasedALG::findBestSolution();

        while(_currentEvaluations < _MAX_EVALUATIONS){
            intermediatePopulation = PBasedALG::Select(_POPULATION_SIZE);

            //Cross:
            k = 0;
            for(int i=0; i<expectedCrosses; ++i){
                PBasedALG::CrossPosition(intermediatePopulation[k], intermediatePopulation[k+1]);
                k += 2;
            }

            //Mutate:
            for(int i=0; i<expectedMutations; ++i){
                chr = Randint(0, _POPULATION_SIZE-1);
                gen = Randint(0, _sizeProblem-1);
                PBasedALG::Mutate(intermediatePopulation[chr], gen);
            }
            
            //Replace:
            _population = intermediatePopulation;
            intermediatePopulation.clear();

            //Evaluate:
            this->PBasedALG::EvaluatePopulation();

            //Optimize population:
            if(timeout == iters){
                for(int i=0; i<_POPULATION_SIZE; ++i){
                    if(Randfloat(0,1) <= pLS)
                        PBasedALG::LS(_population[i]);
                }
                timeout = 0;
            }

            //Elitism:
            _population[_POPULATION_SIZE-1] = _bestSolution;
            PBasedALG::findBestSolution();

            timeout++;
        }
        return _bestSolution;
    }


    PBasedQAP::QAP_solution PBasedALG::AM_best(int iters, float percentage){
        std::vector<PBasedQAP::QAP_solution> intermediatePopulation;
        int expectedCrosses = ceilf((_POPULATION_SIZE/2.0) * _CROSS_P_G); 
        int expectedMutations = ceilf(_POPULATION_SIZE * _sizeProblem * _MUT_P_G);
        int numberOfOptimizations = ceilf(_POPULATION_SIZE * percentage);
        int k, chr, gen;
        int timeout = 0;

        _population.clear();

        this->PBasedALG::GenerateRandomInitialPopulation();

        //Optimize population:
        std::sort(_population.begin(), _population.end());
        for(int i=0; i<numberOfOptimizations; ++i){
            PBasedALG::LS(_population[i]);
        }

        PBasedALG::findBestSolution();

        while(_currentEvaluations < _MAX_EVALUATIONS){
            intermediatePopulation = PBasedALG::Select(_POPULATION_SIZE);

            //Cross:
            k = 0;
            for(int i=0; i<expectedCrosses; ++i){
                PBasedALG::CrossPosition(intermediatePopulation[k], intermediatePopulation[k+1]);
                k += 2;
            }

            //Mutate:
            for(int i=0; i<expectedMutations; ++i){
                chr = Randint(0, _POPULATION_SIZE-1);
                gen = Randint(0, _sizeProblem-1);
                PBasedALG::Mutate(intermediatePopulation[chr], gen);
            }
            
            //Replace:
            _population = intermediatePopulation;
            intermediatePopulation.clear();

            //Evaluate:
            this->PBasedALG::EvaluatePopulation();

            //Optimize population:
            if(timeout == iters){
                std::sort(_population.begin(), _population.end());
                for(int i=0; i<numberOfOptimizations; ++i){
                    PBasedALG::LS(_population[i]);
                }
                timeout = 0;
            }

            //Elitism:
            _population[_POPULATION_SIZE-1] = _bestSolution;
            PBasedALG::findBestSolution();

            timeout++;
        }
        return _bestSolution;
    }
}