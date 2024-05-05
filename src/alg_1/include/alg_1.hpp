#pragma once
#include <../grid/include/point.hpp>
#include <../cs/include/cs.hpp>
#include <iostream>

template <typename T>
class alg_1 {
    public:
        /**
         * @brief Solves the Reverse Shortest Path Problem by binary search on real numbers using the CS algorithm as a decision algorithm.
         * 
         * @param points Center of the circles
         * @param lambda At most length of the path between s_i and t_i
         * @param s_i Index of the source point
         * @param t_i Index of the target point
         * @param precision Number of iterations
         * @return T The minimal radius of the circles, so that the shortest path between s_i and t_i is at most lambda
         */
        T solve(std::vector<Point<T>> points, int lambda, T s_i, T t_i, int precision = 1000);
};



template <typename T>
T alg_1<T>::solve(std::vector<Point<T>> points, int lambda, T s_i, T t_i, int precision){
    cs<T> cs;
    // Binary search on real numbers
    T r = 1;
    T last_r = 0;
    for (int i = 0; i < precision; i++){
        std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAAAAA r: "<<r<<std::endl;
        if (cs.solve(points, r, s_i, t_i, lambda)){
            T temp = r;
            r = (r + last_r) / 2;
            last_r = temp;
        }else{
            last_r = r;
            r *= 2;
        }
    }

    return r;
}
