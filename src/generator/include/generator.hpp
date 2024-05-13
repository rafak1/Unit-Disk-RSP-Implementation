#pragma once
#include "../grid/include/point.hpp"
#include <vector>
#include <random>

template <typename T>
class dataset{
    public:
        std::vector<Point<T>> points;
        int lambda;
        int s_i;
        int t_i;
};


template <typename T>
class generator{
    public:
        /**
         * @brief Generates random dataset for RSPP.
         * 
         * @param n number of points
         * @param min minimal coordinate
         * @param max maximal coordinate
         * @param lambda maximal lambda
         * @return dataset<T> generated dataset
         */
        dataset<T> generate_dataset(int n, T min, T max, int lambda);
};


template <typename T>
dataset<T> generator<T>::generate_dataset(int n, T min, T max, int lambda){
    dataset<T> d;
    
    std::uniform_real_distribution<T> unif(min, max);
    std::default_random_engine re;

    for(int i = 0; i < n; i++){
        d.points.push_back(Point<T>{unif(re), unif(re)});
    }

    d.lambda = (rand() % lambda) + 1;

    d.s_i = rand() % n;
    do{
        d.t_i = rand() % n;
    }while(d.s_i == d.t_i);
    return d;
}