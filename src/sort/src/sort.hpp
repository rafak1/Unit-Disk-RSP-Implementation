#pragma once
#include <limits>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>


template <typename T>
class c_point{
    public:
        /**
         *  0 -> blue point (independet of r)
         *  1 -> vertex of an envelope
         */
        short type;

        std::pair<T, T> point_1;
        std::pair<T, T> point_2;
};

template <typename T>
class c_sort{
    public:
        void sort(std::vector<c_point<T>> &points)

    private:
        T get_root(c_point<T> &a, c_point<T> &b);
};

template <typename T>
class comparision{
    public:
        c_point<T> a;
        c_point<T> b;
}


template <typename T>
void c_sort<T>::sort(std::vector<c_point<T>> &points){
    std::sort(points.begin(), points.end(), [](c_point<T> &a, c_point<T> &b){
        if(a.point_1.first == b.point_1.first){
            return a.point_1.second < b.point_1.second;
        }
        return a.point_1.first < b.point_1.first;
    });
}


template <typename T>
T c_sort<T>::get_root(c_point<T> &a, c_point<T> &b){
    short sum = a.type + b.type;
    if(sum == 0){
        return -1;
    }else if (sum == 1){
        return -2;
    }else if (a.type == 1 && b.type == 1){
        if(a.point.first == b.point.first){
            return -3;
        }
        T a = (a.point.second + b.point.second) / (a.point.first - b.point.first);
        T mid_x = (a.point.first + b.point.first) / 2;
        T mid_y = (a.point.second + b.point.second) / 2;


        

    }
}