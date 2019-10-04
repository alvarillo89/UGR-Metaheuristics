#ifndef _PBASED_ALG_H
#define _PBASED_ALG_H


#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <limits>
#include <cmath>
#include "random.h"

namespace PBasedQAP{

    /** @struct QAP solutión: 
     * @brief Defines a permutation where S[i] = j, 
     * means that i unit is located on j location. */
    struct QAP_solution{
        std::vector<int> S;     /** Permutation */
        float cost;             /** Permutation Cost */
        bool needReevaluation;  /** If true its necessary to compute cost again */

        inline bool operator<(const QAP_solution &other) const{
            return (this->cost < other.cost);
        }
    };

    /** @typedef _fmatrix 
     * @brief Defines a short name for a float matrix */
    typedef std::vector< std::vector<float> > _fmatrix;  

    /** @class Population Based Algorithm
     * @brief  Population Based Algorithm class */
    class PBasedALG{
        public:
            
        /** @brief Constructor 
         * @param path to .dat file 
         * @param seed for random number generator */
        PBasedALG(std::string path, int seed, int pSize);

        /** @brief Start Genetic Algorithm  Generational
         * First with position Cross operator
         * Second with OX Cross operator
         * @return The best solution found by the Genetic Algorithm */
        PBasedQAP::QAP_solution AGG_pos();
        PBasedQAP::QAP_solution AGG_OX();

        /** @brief Start Genetic Algorithm  Stationary
         * First with position Cross operator
         * Second with OX Cross operator
         * @return The best solution found by the Genetic Algorithm */
        PBasedQAP::QAP_solution AGE_pos();
        PBasedQAP::QAP_solution AGE_OX();

        /** @brief Start Memetic Algorithm
         * First apply local search to a solution with probability
         * Second apply LS to a percentage of the best solutions in population 
         * @return The best solution found by the Memetic Algorithm */
        PBasedQAP::QAP_solution AM(int iters, float pLS);
        PBasedQAP::QAP_solution AM_best(int iters, float percentage);

        //private:
            // Atributes:
            // Define neccesary constants: 
            const int _POPULATION_SIZE;             /**< population's solutions */
            const float _CROSS_P_G = 0.7;           /**< Generational Cross probability */
            const float _MUT_P_G = 0.001;           /**< Generational Mutation probability */
            const int _MAX_EVALUATIONS = 50000;     /**< Max number of evaluations to do */

            int _currentEvaluations;
            int _sizeProblem;
            std::vector<PBasedQAP::QAP_solution> _population; /**< Population */
            PBasedQAP::_fmatrix _F;
            PBasedQAP::_fmatrix _D;

            PBasedQAP::QAP_solution _bestSolution;  /**< Elitism */
            
            // Private member functions:

            /** @brief Load .dat file 
             * @param path to .dat file*/
            void load(const std::string path);

            /** @brief Generate initial population */ 
            void GenerateRandomInitialPopulation();

            /** @brief Evaluate cost of a solution
             * @param S The QAP solution to evaluate */
            void EvaluateSolution(PBasedQAP::QAP_solution &S);

            /** @brief Evaluate cost of entire population */
            void EvaluatePopulation();

            /** @brief Computes the new best solution in population */
            void findBestSolution();

            /** @brief Selection operator (binary tournament) */
            std::vector<PBasedQAP::QAP_solution> Select(int size);

            /** @brief Mutation operator */
            void Mutate(PBasedQAP::QAP_solution &S, int pos);

            /** @brief Position Cross operator */
            void CrossPosition(PBasedQAP::QAP_solution &S1, PBasedQAP::QAP_solution &S2);

            /** @brief OX Cross Operator */
            void CrossOX(PBasedQAP::QAP_solution &S1, PBasedQAP::QAP_solution &S2);

            /** @brief factorized cost */
            float DeltaCost(const PBasedQAP::QAP_solution& pi, int r, int s);

            /** @brief Neighbour operator for local search */
            void ApplyMove(PBasedQAP::QAP_solution& S, int i, int j, float var);

            /** @brief Local Search for Memetic Algorithm */
            void LS(PBasedQAP::QAP_solution &S);
    };
}

#endif