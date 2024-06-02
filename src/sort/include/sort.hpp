#pragma once
#include <limits>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>
#include "../grid/include/point.hpp"
#include "../cs/include/cs.hpp"

template <typename T>
class c_point{
    public:
        /**
         *  0 -> blue point (independet of r)
         *  1 -> vertex of an envelope
         */
        short type;

        int instance;
        int instance_n;
        bool is_vertical;
        bool is_reverse;

        Point<T> point_1;
        Point<T> point_2;

        c_point(short type, int instance, int instance_n, bool is_vertical, bool is_reverse, std::pair<T, T> point_1, std::pair<T, T> point_2):
            type(type), instance(instance), instance_n(instance_n), is_vertical(is_vertical), is_reverse(is_reverse){
            this->point_1 = Point<T>{point_1.first, point_1.second};
            this->point_2 = Point<T>{point_2.first, point_2.second};
        }
};

template <typename T>
class c_sort{
    public:
        T interval_l;
        T interval_r;
        c_sort() = default;
        c_sort(T l, T r): interval_l(l), interval_r(r){};
        std::pair<T,T> sort(std::vector<c_point<T>> &points, std::vector<int> & og_by_x, std::vector<int> & og_by_y, std::vector<Point<T>> &og_points, int s_i, int t_i, int lambda, int x_or_y);
        const T eps = std::numeric_limits<T>::epsilon();
    private:
        std::pair<T,bool> get_root(c_point<T> &a, c_point<T> &b, int x_or_y);
        T get_r(){
            return interval_l + (interval_r - interval_l) / 2;
        }
};

template <typename T>
class Comparision{
    public:
        c_point<T> a;
        c_point<T> b;
};

/**
 * @brief Does a parametric sdort of the points
 * 
 * @tparam T 
 * @param points c_points to be sorted
 * @param og_by_x values to run the cs algorithm:
 * @param og_by_y 
 * @param og_points 
 * @param s_i 
 * @param t_i 
 * @param lambda 
 * @param x_or_y determines if the sort is by x or y, see point.hpp
 */
template <typename T>
std::pair<T,T> c_sort<T>::sort(
    std::vector<c_point<T>> &points,
    std::vector<int> & og_by_x, 
    std::vector<int> & og_by_y, 
    std::vector<Point<T>> &og_points, 
    int s_i, int t_i, int lambda,
    int x_or_y){

    interval_l = 0;
    interval_r = std::numeric_limits<T>::max();

    int other = x_or_y == 0 ? 1 : 0;

    cs<T> cs;

    //std::cout<<"SORTING"<<std::endl;

    std::sort(points.begin(), points.end(), [&](c_point<T> &a, c_point<T> &b){

        std::pair<T, bool> root_pair = get_root(a, b, x_or_y);

        std::cout<<"ROOT: "<<root_pair.first<<" "<<root_pair.second<<std::endl;
        std::cout<<"A: "<<a.point_1.x<<" "<<a.point_1.y<<" "<<a.point_2.x<<" "<<a.point_2.y<<" type- "<<a.type<<std::endl;
        std::cout<<"B: "<<b.point_1.x<<" "<<b.point_1.y<<" "<<b.point_2.x<<" "<<b.point_2.y<<" type- "<<b.type<<std::endl;

        T root = root_pair.first;
        bool direction = root_pair.second;

        if(root == -1){
            // two blue points
            if(a.point_1.get(x_or_y) == b.point_1.get(x_or_y)){
                return a.point_1.get(other) < b.point_1.get(other);
            }
            return a.point_1.get(x_or_y) < b.point_1.get(x_or_y);
        }else if(root == -2){
             // two envelope vertices in the same instance
            T r = get_r();

            Point<T> p1 = cs.intersection_point(a.point_1, a.point_2, r, a.is_reverse, a.is_vertical);
            Point<T> p2 = cs.intersection_point(b.point_1, b.point_2, r, b.is_reverse, b.is_vertical);

            if(p1.get(x_or_y) == p2.get(x_or_y)){
                return p1.get(other) < p2.get(other);
            }
            return p1.get(x_or_y) < p2.get(x_or_y);
        }else if(root == -3){
            // one envelope vertex and one blue point, vertical line from vertex points(the same y)
            if(a.type == 1){
                T mid_x = (a.point_1.get(x_or_y) + a.point_2.get(x_or_y)) / 2;
                if(mid_x == b.point_1.get(x_or_y)){
                    return a.point_1.get(other) < b.point_1.get(other);
                }
                return mid_x < b.point_1.get(x_or_y);
            }else{
                T mid_x = (b.point_1.get(x_or_y) + b.point_2.get(x_or_y)) / 2;
                if(a.point_1.get(x_or_y) == mid_x){
                    return a.point_1.get(other) < b.point_1.get(other);
                }
                return a.point_1.get(x_or_y) < mid_x;
            }
        }else if(root == -4){
            // Two envelope vertices, not having points of interest with the same x (also when delta < 0)
            T mid_x1 = (a.point_1.get(x_or_y) + a.point_2.get(x_or_y)) / 2;
            T mid_x2 = (b.point_1.get(x_or_y) + b.point_2.get(x_or_y)) / 2;
            if(mid_x1 == mid_x2){
                T mid_y1 = (a.point_1.get(other) + a.point_2.get(other)) / 2;
                T mid_y2 = (b.point_1.get(other) + b.point_2.get(other)) / 2;
                return mid_y1 < mid_y2;
            }
            return mid_x1 < mid_x2;
        }else{
            bool res = cs.solve(og_points, og_by_x, og_by_y, root, s_i, t_i, lambda);

            if(res){
                interval_r = std::min(interval_r, root);
            }else{
                interval_l = std::max(interval_l, root);
            }

            if(direction){
                return res;
            }else{
                return !res;
            }
        }
    });

    std::cout<<"SORTED"<<std::endl;
    for(auto p : points){
        std::cout<<p.point_1.x<<" "<<p.point_1.y<<" "<<p.point_2.x<<" "<<p.point_2.y<<" type- "<<p.type<<std::endl;
    }

    return {interval_l, interval_r};
}


template <typename T>
std::pair<T,bool> c_sort<T>::get_root(c_point<T> &a_point, c_point<T> &b_point, int x_or_y){
    cs<T> cs;

    int other = x_or_y == 0 ? 1 : 0;

    short sum = a_point.type + b_point.type;
    if(sum == 0){
        return {-1, false};
    }else if (sum == 1){
        // one envelope vertex and one blue point
        c_point<T> vertex = a_point.type == 1 ? a_point : b_point;
        c_point<T> point = a_point.type == 0 ? a_point : b_point;


        if(vertex.point_1.get(other) == vertex.point_2.get(other)){
            return {-3, false};
        }

        T a = - (vertex.point_1.get(x_or_y) - vertex.point_2.get(x_or_y)) / (vertex.point_1.get(other) - vertex.point_2.get(other));
        T mid_x = (vertex.point_1.get(x_or_y) + vertex.point_2.get(x_or_y)) / 2;
        T mid_y = (vertex.point_1.get(other) + vertex.point_2.get(other)) / 2;
        T b = mid_y - a * mid_x;

        T d = sqrt(pow(point.point_1.get(x_or_y) - point.point_2.get(x_or_y), 2) + pow(point.point_1.get(other) - point.point_2.get(other), 2));

        Point<T> v = Point<T>{point.point_1.get(x_or_y), a * point.point_1.get(x_or_y) + b};

        T h = sqrt(pow(v.get(x_or_y) - mid_x, 2) + pow(v.get(other) - mid_y, 2));

        T r = sqrt(pow(d/2, 2) + pow(h, 2));

        bool direction;

        Point<T> inter = cs.intersection_point(vertex.point_1, vertex.point_2, r, vertex.is_reverse, vertex.is_vertical);
        
        if(abs(inter.get(x_or_y) - v.get(x_or_y)) > eps || abs(inter.get(other) - v.get(other)) > eps){
            return {-3, false};
        }

        direction = v.get(x_or_y) < mid_x;


        if(a_point.type == 1)
            direction = !direction;

        return {r, direction};

    }else if (a_point.instance != b_point.instance || a_point.instance_n != b_point.instance_n){
        // two envelope vertices in different instances
        T d1 = sqrt(pow(a_point.point_1.get(x_or_y) - a_point.point_2.get(x_or_y), 2) + pow(a_point.point_1.get(other) - a_point.point_2.get(other), 2));
        T d2 = sqrt(pow(b_point.point_1.get(x_or_y) - b_point.point_2.get(x_or_y), 2) + pow(b_point.point_1.get(other) - b_point.point_2.get(other), 2));
        T a = - (a_point.point_1.get(other) - a_point.point_2.get(other)) / d1;
        T b = - (b_point.point_1.get(other) - b_point.point_2.get(other)) / d2;
        T c =  (b_point.point_1.get(x_or_y) - b_point.point_2.get(x_or_y)) / 2  -  (a_point.point_1.get(x_or_y) - a_point.point_2.get(x_or_y)) / 2;
        T k1 = pow(d1/2 , 2);
        T k2 = pow(d2/2 , 2);
        T A = pow(a, 2);
        T B = pow(b, 2);
        T C = pow(c, 2);
        
        //quadratic equation
        T M = pow(A, 2) + pow(B, 2) - 2 * A * B;
        T N = 2*A*B*k2 + 2*A*B*k1 - 2*pow(A, 2)*k1 - 2*A*C - 2*pow(B, 2)*k2 - 2*B*C;
        T O = -k1*k2 + pow(A, 2)*pow(k1,2) + 2*A*B*k1*k2 + 2*A*C*k1 + pow(B, 2)*pow(k2,2) + 2*B*C*k2 + pow(C, 2);

        T delta = pow(N, 2) - 4 * M * O;
        if (delta < 0){
            return {-4, false};
        }

        T r1 = (-N + sqrt(delta)) / (2 * M);

        T r2 = (-N - sqrt(delta)) / (2 * M);

        if(r1 > 0){
            r1 = sqrt(r1);

            Point<T> inter1_r1 = cs.intersection_point(a_point.point_1, a_point.point_2, r1, a_point.is_reverse, a_point.is_vertical);
            Point<T> inter2_r1 = cs.intersection_point(b_point.point_1, b_point.point_2, r1, b_point.is_reverse, b_point.is_vertical);

            if(abs(inter1_r1.get(x_or_y) - inter2_r1.get(x_or_y)) < eps){
                Point<T> r1_moved1 = cs.intersection_point(a_point.point_1, a_point.point_2, r1 + 1, a_point.is_reverse, a_point.is_vertical);
                Point<T> r1_moved2 = cs.intersection_point(b_point.point_1, b_point.point_2, r1 + 1, b_point.is_reverse, b_point.is_vertical);

                return {r1, r1_moved1.get(x_or_y) < r1_moved2.get(x_or_y)};
            }
        }
        if( r2 > 0){
            r2 = sqrt(r2);

            Point<T> inter1_r2 = cs.intersection_point(a_point.point_1, a_point.point_2, r2, a_point.is_reverse, a_point.is_vertical);
            Point<T> inter2_r2 = cs.intersection_point(b_point.point_1, b_point.point_2, r2, b_point.is_reverse, b_point.is_vertical);

            if(abs(inter1_r2.get(x_or_y) - inter2_r2.get(x_or_y)) < eps){
                Point<T> r2_moved1 = cs.intersection_point(a_point.point_1, a_point.point_2, r2 + 1, a_point.is_reverse, a_point.is_vertical);
                Point<T> r2_moved2 = cs.intersection_point(b_point.point_1, b_point.point_2, r2 + 1, b_point.is_reverse, b_point.is_vertical);

                return {r2, r2_moved1.get(x_or_y) < r2_moved2.get(x_or_y)};
            }
        }
        return {-4, false};
    }else{
        //two envelope vertices in the same instance
        return {-2, false};
    }
}