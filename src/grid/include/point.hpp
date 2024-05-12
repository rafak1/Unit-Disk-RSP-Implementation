#pragma once
#include <cmath>

template <typename T>
class Point {
    public:
        T x, y;
        int grid_x, grid_y, grid_cell, grid_index, original_index;

        /**
         * @brief if i == 0 return x else return y
         * 
         * @param i 
         * @return int x or y
         */
        T get(int i) const {
            if (i == 0) return x;
            return y;
        }

        T euclidean_distance(Point<T> &p){
            return sqrt(pow(x - p.x, 2) + pow(y - p.y, 2));
        }
};