#pragma once
#include <limits>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <stdlib.h>
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
        bool is_marked = false;

        Point<T> point_1;
        Point<T> point_2;

        c_point(short type, int instance, int instance_n, bool is_vertical, bool is_reverse, std::pair<T, T> point_1, std::pair<T, T> point_2):
            type(type), instance(instance), instance_n(instance_n), is_vertical(is_vertical), is_reverse(is_reverse){
            this->point_1 = Point<T>{point_1.first, point_1.second};
            this->point_2 = Point<T>{point_2.first, point_2.second};
        }
        c_point() = default;
};

template <typename T>
class node{
    public:
        c_point<T> point;
        node<T> *left = nullptr;
        node<T> *right = nullptr;
        int index;

        void delete_tree(){
            if(left != nullptr){
                left->delete_tree();
                delete left;
            }
            if(right != nullptr){
                right->delete_tree();
                delete right;
            }
        }
};

template <typename T>
class Comparision{
    public:
        c_point<T> a;
        c_point<T> b;
        bool result;
        T root;
        bool direction;
        int weight;
        int instance;
        node<T> *c_node = nullptr;
};

template <typename T>
class Instnace{
    public:
        std::vector<c_point<T>> points;
        T initial_weight;
};

template <typename T>
class c_sort{
    public:
        T interval_l;
        T interval_r;
        c_sort() = default;
        c_sort(T l, T r): interval_l(l), interval_r(r){};
        std::pair<T,T> brut_sort(std::vector<c_point<T>> &points, std::vector<int> & og_by_x, std::vector<int> & og_by_y, std::vector<Point<T>> &og_points, int s_i, int t_i, int lambda, int x_or_y);
        std::pair<T,T> sort(std::vector<c_point<T>> &points, std::vector<int> & og_by_x, std::vector<int> & og_by_y, std::vector<Point<T>> &og_points, int s_i, int t_i, int lambda, int x_or_y);
        const T eps = 1e-3;
    private:
        std::pair<T,bool> get_root(c_point<T> &a, c_point<T> &b, int x_or_y);
        T get_r(){
            return interval_l + (interval_r - interval_l) / 2;
        }
        std::pair<bool, T> resolve(c_point<T> &a, c_point<T> &b, int x_or_y);
        Comparision <T> weighted_median(std::vector<Comparision<T>> comparisions, T total_weight);
        node<T>* create_bst(std::vector<c_point<T>> &pivots, int l, int r);
        void move_comp(    bool first, Comparision<T>* comp,std::vector<Comparision<T>*> &comparisions, std::vector<Comparision<T>*> &active_comparisions,std::vector<std::vector<std::vector<c_point<T>>>> &buckets,  int instance_i, int x_or_y);
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

    T w = 1;
    cs<T> cs;

    int n = points.size();

    std::vector<Instnace<T>> instances;
    instances.push_back(Instnace<T>{points, w});

   //std::cout<<"SORTING"<<std::endl;

    while(true){
        std::vector< std::vector<Comparision<T>*>> comparisions;
        std::vector<Comparision<T>*> active_comparisions;
        
        comparisions.resize(instances.size(), std::vector<Comparision<T>*>());

        int instance_i = 0;
        for(auto &instance : instances){
            int sqrt_n = sqrt(instance.points.size());
            std::vector<c_point<T>> pivots;

            //std::cout<<"instance "<<instance_i<<std::endl;
            
            //randomly select sqrt(n) points
            std::random_shuffle(instance.points.begin(), instance.points.end());
            for(int i = 0; i < sqrt_n; i++){
                instance.points[i].is_marked = true;
                pivots.push_back(instance.points[i]);
            }

            //std::cout<<"pivots "<<pivots.size()<<std::endl;

            //create comparisions for sorting the marked points
            for(int i = 0; i < sqrt_n; i++){
                for(int j = 0; j < sqrt_n; j++){
                    Comparision<T> *comp = new Comparision<T>;
                    comp->a = instance.points[i];
                    comp->b = instance.points[j];
                    comp->weight = instance.initial_weight;
                    comparisions[instance_i].push_back(comp);
                }
            }

            //std::cout<<"comparisions "<<comparisions[instance_i].size()<<std::endl;

            for(int i=0;i<comparisions[instance_i].size();i++){
                std::pair<bool, T> comp_res = resolve(comparisions[instance_i][i]->a, comparisions[instance_i][i]->b, x_or_y);
                if(comp_res.second == -1){
                    comparisions[instance_i][i]->result = comp_res.first;
                }else{
                    comparisions[instance_i][i]->root = comp_res.second;
                    active_comparisions.push_back(comparisions[instance_i][i]);
                }
            }

            //std::cout<<"active_comparisions "<<active_comparisions.size()<<std::endl;
            
            instance_i++;
        }

        //resolve the comparisions
        while(!active_comparisions.empty()){

            //std::cout<<"active_comparisions "<<active_comparisions.size()<<std::endl;

            T total_weight = 0;
            for(auto c : active_comparisions){
                total_weight += c->weight;
            }
            
            //select a weghted median comparision
            std::vector<Comparision<T>> active_comparisions_copy;
            for(auto c : active_comparisions){
                active_comparisions_copy.push_back(*c);
            }
            Comparision<T> comp = weighted_median(active_comparisions_copy, total_weight);

            //std::cout<<"median "<<comp.root<<std::endl;

            //resolve the comparision
            std::pair<bool, T> comp_res = resolve(comp.a, comp.b, x_or_y);

            bool res = cs.solve(og_points, og_by_x, og_by_y, comp.root, s_i, t_i, lambda);

            if(res){
                interval_r = std::min(interval_r, comp.root);
            }else{
                interval_l = std::max(interval_l, comp.root);
            }

            //std::cout<<"res "<<res<<" "<<comp.root<<std::endl;

            std::vector<Comparision<T> *> new_comparisions;
            for(auto c : active_comparisions){
                if(res){
                    if(c->root < comp.root)
                        new_comparisions.push_back(c);
                    else{
                        std::pair<bool, T> comp_res_int = resolve(c->a, c->b, x_or_y);
                        c->result = comp_res_int.first;
                    }

                }else{
                    if(c->root > comp.root)
                        new_comparisions.push_back(c);
                    else{
                        std::pair<bool, T> comp_res_int = resolve(c->a, c->b, x_or_y);
                        c->result = comp_res_int.first;
                    }
                }
            }

            active_comparisions = new_comparisions;
        }


        //bubble sort but now with solved comparisions
        for(int inst = 0; inst < instances.size(); inst++){
            int sqrt_n = sqrt(instances[inst].points.size());
            for(int i = 0; i < sqrt_n; i++){
                for(int j = 0; j < sqrt_n; j++){
                    int ind = i * sqrt_n + j;
                    if(ind < comparisions[inst].size() && comparisions[inst][ind]->result){
                        std::swap(instances[inst].points[i], instances[inst].points[j]);
                    }
                }
            }
        }

        //std::cout<<"deleting"<<std::endl;

        for(int i=0;i<instances.size();i++){
            for(int j=0;j<comparisions[i].size();j++){
                delete comparisions[i][j];
            }
            comparisions[i].clear();
        }
        active_comparisions.clear();

        //std::cout<<"bst"<<std::endl;

        //create a bst tree with the sorted points
        std::vector<node<T>*> roots;
        for(auto &instance : instances){
            int sqrt_n = sqrt(instance.points.size());
            node<T> *root = create_bst(instance.points, 0, sqrt_n - 1);
            roots.push_back(root);
        }

        //std::cout<<"bst done"<<std::endl;

        
        std::vector<std::vector<std::vector<c_point<T>>>> buckets;
        buckets.resize(instances.size(), std::vector<std::vector<c_point<T>>>());
        for(int i = 0; i < instances.size(); i++){
            int sqrt_n = sqrt(instances[i].points.size());
            buckets[i].resize(sqrt_n+2, std::vector<c_point<T>>());
        }

        //route the rest of the points through the tree
        instance_i = 0;
        for(auto &instance : instances){
            for(int i = 0; i < instance.points.size(); i++){
                if(instance.points[i].is_marked){
                    continue;
                }
                Comparision<T> *comp = new Comparision<T>;
                comp->a = instance.points[i];
                comp->b = roots[instance_i]->point;
                comp->c_node = roots[instance_i];
                comp->weight = instance.initial_weight / ( 2 * pow((long double)instance.points.size(), 2));
                comp->instance = instance_i;

                move_comp(true, comp, comparisions[instance_i], active_comparisions, buckets, instance_i, x_or_y);
            }
            instance_i++;
        }

        //std::cout<<"routing start done"<<std::endl;
        
        while(!active_comparisions.empty()){

            //std::cout<<"active_comparisions 2 -> "<<active_comparisions.size()<<std::endl;

            
            T total_weight = 0;
            for(auto c : active_comparisions){
                total_weight += c->weight;
            }

            //select a weghted median comparision
            std::vector<Comparision<T>> active_comparisions_copy;
            for(auto c : active_comparisions){
                active_comparisions_copy.push_back(*c);
            }
            Comparision<T> comp = weighted_median(active_comparisions_copy, total_weight);

           //std::cout<<"median "<<comp.root<<std::endl;

            //resolve the comparision
            std::pair<bool, T> comp_res = resolve(comp.a, comp.b, x_or_y);

            T r = comp.root;

            bool res = cs.solve(og_points, og_by_x, og_by_y, r, s_i, t_i, lambda);

            if(res){
                interval_r = std::min(interval_r, r);
            }else{
                interval_l = std::max(interval_l, r);
            }

            //std::cout<<"res "<<res<<" "<<r<<std::endl;



            std::vector<Comparision<T> *> new_comparisions;
            for(auto c : active_comparisions){

                if(c->root != c->root || c->root < 0)
                    continue;
                
                //std::cout<<"active_comparisions 2 -> "<<c->root<<std::endl;

                if(res){
                    if(c->root < r){
                        new_comparisions.push_back(c);
                    }else{
                        c->result = true;
                        if(comp_res.first){
                            if(c->c_node->right != nullptr){
                                Comparision<T> *comp = new Comparision<T>;
                                comp->a = c->a;
                                comp->b = c->c_node->right->point;
                                comp->c_node = c->c_node->right;
                                comp->weight = c->weight / 2;
                                comp->instance = c->instance;
                                
                                move_comp(true, comp, comparisions[c->instance], new_comparisions, buckets, c->instance, x_or_y);
                            }else{
                                buckets[c->instance][c->c_node->index].push_back(c->a);
                            }
                        }else{
                            if(c->c_node->left != nullptr){
                                Comparision<T>* comp = new Comparision<T>;
                                comp->a = c->a;
                                comp->b = c->c_node->left->point;
                                comp->c_node = c->c_node->left;
                                comp->weight = c->weight / 2;
                                comp->instance = c->instance;
                                
                                move_comp(false, comp, comparisions[c->instance], new_comparisions, buckets, c->instance, x_or_y);
                            }else{
                                int index = c->c_node->index - 1;
                                index = index < 0 ? 0 : index;
                                buckets[c->instance][index].push_back(c->a);
                            }
                        }
                    }
                }else{
                    if(c->root > r){
                        new_comparisions.push_back(c);
                    }else{
                        std::pair<bool, T> comp_res = resolve(c->a, c->b, x_or_y);
                        if(comp_res.first){
                            if(c->c_node->left != nullptr){
                                Comparision<T> *comp = new Comparision<T>;
                                comp->a = c->a;
                                comp->b = c->c_node->left->point;
                                comp->c_node = c->c_node->left;
                                comp->weight = c->weight / 2;
                                comp->instance = c->instance;

                                move_comp(false, comp, comparisions[c->instance], new_comparisions, buckets, c->instance, x_or_y);

                            }else{
                                int index = c->c_node->index - 1;
                                index = index < 0 ? 0 : index;
                                buckets[c->instance][index].push_back(c->a);
                            }
                        }else{
                            if(c->c_node->right != nullptr){
                                Comparision<T> *comp = new Comparision<T>;
                                comp->a = c->a;
                                comp->b = c->c_node->right->point;
                                comp->c_node = c->c_node->right;
                                comp->weight = c->weight / 2;
                                comp->instance = c->instance;

                                move_comp(true, comp, comparisions[c->instance], new_comparisions, buckets, c->instance, x_or_y);

                            }else{
                                buckets[c->instance][c->c_node->index].push_back(c->a);
                            }
                        }
                    }
                }
            }
            active_comparisions = new_comparisions;
        }

        //std::cout<<"routing done"<<std::endl;

        //recurse into subproblems
        std::vector<Instnace<T>> new_instances;
        for(int i = 0; i < instances.size(); i++){
            int sqrt_n = sqrt(instances[i].points.size());
            for(int j = 0; j <= sqrt_n; j++){
                if(buckets[i][j].size() > 1){
                    Instnace<T> new_instance;
                    new_instance.points = buckets[i][j];
                    new_instance.initial_weight = instances[i].initial_weight / ( 4 * pow((long double)new_instance.points.size(), 4.5));
                    new_instances.push_back(new_instance);
                }
            }
        }

        instances = new_instances;

        //std::cout<<"new instances "<<instances.size()<<std::endl;


        //clean up
        for(int i = 0; i < n; i++){
            if(points[i].is_marked){
                points[i].is_marked = false;
            }
        }

        for(int i = 0; i < comparisions.size(); i++){
            for(int j = 0;j < comparisions[i].size(); j++){
                delete comparisions[i][j];
            }
        }
        //std::cout<<"a\n";


        for(int i = 0; i < instances.size(); i++){
            int sqrt_n = sqrt(instances[i].points.size());
            for(int j = 0; j < sqrt_n; j++){
                instances[i].points[j].is_marked = false;
            }
        }

        for(auto root : roots){
            if(root != nullptr){
                root->delete_tree();
                delete root;
            }
        }
        roots.clear();

        if(instances.size() == 0){
            break;
        }
    }


    return {interval_l, interval_r};
}


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
std::pair<T,T> c_sort<T>::brut_sort(
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

        std::pair<bool, T> comp_res = resolve(a, b, x_or_y);

        if(comp_res.second == -1){
            return comp_res.first;
        }

        bool res = cs.solve(og_points, og_by_x, og_by_y, comp_res.second, s_i, t_i, lambda);

        if(res){
            interval_r = std::min(interval_r, comp_res.second);
        }else{
            interval_l = std::max(interval_l, comp_res.second);
        }

        if(comp_res.first){
            return res;
        }else{
            return !res;
        }
    });

    /*std::cout<<"SORTED"<<std::endl;
    for(auto p : points){
        std::cout<<p.point_1.x<<" "<<p.point_1.y<<" "<<p.point_2.x<<" "<<p.point_2.y<<" type- "<<p.type<<std::endl;
    }*/

    return {interval_l, interval_r};
}


template <typename T>
node<T>* c_sort<T>::create_bst(std::vector<c_point<T>> &pivots, int l, int r){
    if(l > r){
        return nullptr;
    }

    int m = l + (r - l) / 2;

    if (m >= pivots.size() || m < 0){
        return nullptr;
    }   

    node<T> *root = new node<T>;
    root->point = pivots[m];
    root->index = m;
    root->left = create_bst(pivots, l, m - 1);
    root->right = create_bst(pivots, m + 1, r);

    return root;
}

template <typename T>
std::pair<bool, T> c_sort<T>::resolve(c_point<T> &a, c_point<T> &b, int x_or_y){

    int other = x_or_y == 0 ? 1 : 0; 

    cs<T> cs;

    std::pair<T, bool> root_pair = get_root(a, b, x_or_y);

    /*std::cout<<"ROOT: "<<root_pair.first<<" "<<root_pair.second<<std::endl;
    std::cout<<"A: "<<a.point_1.x<<" "<<a.point_1.y<<" "<<a.point_2.x<<" "<<a.point_2.y<<" type- "<<a.type<<std::endl;
    std::cout<<"B: "<<b.point_1.x<<" "<<b.point_1.y<<" "<<b.point_2.x<<" "<<b.point_2.y<<" type- "<<b.type<<std::endl;*/

    T root = root_pair.first;
    bool direction = root_pair.second;

    if(root == -1){
        // two blue points
        if(a.point_1.get(x_or_y) == b.point_1.get(x_or_y)){
            return {a.point_1.get(other) < b.point_1.get(other), -1};
        }
        return {a.point_1.get(x_or_y) < b.point_1.get(x_or_y), -1};
    }else if(root == -2){
            // two envelope vertices in the same instance
        T r = get_r();

        Point<T> p1 = cs.intersection_point(a.point_1, a.point_2, r, a.is_reverse, a.is_vertical);
        Point<T> p2 = cs.intersection_point(b.point_1, b.point_2, r, b.is_reverse, b.is_vertical);

        if(p1.get(x_or_y) == p2.get(x_or_y)){
            return {p1.get(other) < p2.get(other), -1};
        }
        return {p1.get(x_or_y) < p2.get(x_or_y), -1};
    }else if(root == -3){
        // one envelope vertex and one blue point, vertical line from vertex points(the same y)
        if(a.type == 1){
            T mid_x = (a.point_1.get(x_or_y) + a.point_2.get(x_or_y)) / 2;
            if(mid_x == b.point_1.get(x_or_y)){
                return {a.point_1.get(other) < b.point_1.get(other), -1};
            }
            return {mid_x < b.point_1.get(x_or_y), -1};
        }else{
            T mid_x = (b.point_1.get(x_or_y) + b.point_2.get(x_or_y)) / 2;
            if(a.point_1.get(x_or_y) == mid_x){
                return {a.point_1.get(other) < b.point_1.get(other), -1};
            }
            return {a.point_1.get(x_or_y) < mid_x, -1};
        }
    }else if(root == -4){
        // Two envelope vertices, not having points of interest with the same x (also when delta < 0)
        T mid_x1 = (a.point_1.get(x_or_y) + a.point_2.get(x_or_y)) / 2;
        T mid_x2 = (b.point_1.get(x_or_y) + b.point_2.get(x_or_y)) / 2;
        if(mid_x1 == mid_x2){
            T mid_y1 = (a.point_1.get(other) + a.point_2.get(other)) / 2;
            T mid_y2 = (b.point_1.get(other) + b.point_2.get(other)) / 2;
            return {mid_y1 < mid_y2, -1};
        }
        return {mid_x1 < mid_x2, -1};
    }else{
        return {direction, root};
    }
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

        //std::cout<<"a "<<a<<" b "<<b<<std::endl;

        T d = sqrt(pow(vertex.point_1.get(x_or_y) - vertex.point_2.get(x_or_y), 2) + pow(vertex.point_1.get(other) - vertex.point_2.get(other), 2));

        Point<T> v = Point<T>{point.point_1.get(x_or_y), a * point.point_1.get(x_or_y) + b};

        T h = sqrt(pow(v.get(x_or_y) - mid_x, 2) + pow(v.get(other) - mid_y, 2));

        //std::cout<<"d "<<d<<" h "<<h<<std::endl;

        T r = sqrt(pow(d/2, 2) + pow(h, 2));

        bool direction;

        Point<T> inter = cs.intersection_point(vertex.point_1, vertex.point_2, r, vertex.is_reverse, vertex.is_vertical);

        //std::cout<<"inter "<<inter.x<<" "<<inter.y<<" v "<<v.x<<" "<<v.y<<" | "<<abs(inter.get(x_or_y) - v.get(x_or_y))<<" "<<abs(inter.get(other) - v.get(other))<<std::endl;
        
        if(abs(inter.get(x_or_y) - v.get(x_or_y)) > eps || abs(inter.get(other) - v.get(other)) > eps){
            return {-3, false};
        }

        Point<T> point_1_moved = cs.intersection_point(vertex.point_1, vertex.point_2, r + 1, vertex.is_reverse, vertex.is_vertical);

        //std::cout<<"dir "<<v.get(x_or_y)<<" "<<point_1_moved.get(x_or_y)<<" "<<inter.get(x_or_y)<<std::endl;

        direction = v.get(x_or_y) < point_1_moved.get(x_or_y);


        if(a_point.type == 0)
            direction = !direction;

        return {r, direction};

    }else if (a_point.instance != b_point.instance || a_point.instance_n != b_point.instance_n){
        // two envelope vertices in different instances
        T d1 = sqrt(pow(a_point.point_1.get(x_or_y) - a_point.point_2.get(x_or_y), 2) + pow(a_point.point_1.get(other) - a_point.point_2.get(other), 2));
        T d2 = sqrt(pow(b_point.point_1.get(x_or_y) - b_point.point_2.get(x_or_y), 2) + pow(b_point.point_1.get(other) - b_point.point_2.get(other), 2));
        T a = - (a_point.point_1.get(other) - a_point.point_2.get(other)) / d1;
        T b = (b_point.point_1.get(other) - b_point.point_2.get(other)) / d2;
        T c =  (b_point.point_1.get(x_or_y) + b_point.point_2.get(x_or_y)) / 2  -  (a_point.point_1.get(x_or_y) + a_point.point_2.get(x_or_y)) / 2;
        T k1 = pow(d1/2 , 2);
        T k2 = pow(d2/2 , 2);
        T A = pow(a, 2);
        T B = pow(b, 2);
        T C = pow(c, 2);

        //std::cout<<"a "<<a<<" b "<<b<<" c "<<c<<" k1 "<<k1<<" k2 "<<k2<<" A "<<A<<" B "<<B<<" C "<<C<<" "<<d1<<" "<<d2<<std::endl;
        
        //quadratic equation
        T M = pow(A, 2) + pow(B, 2) - 2 * A * B;
        T N = 2*A*B*k2 + 2*A*B*k1 - 2*pow(A, 2)*k1 - 2*A*C - 2*pow(B, 2)*k2 - 2*B*C;
        T O = -k1*k2 + pow(A, 2)*pow(k1,2) + 2*A*B*k1*k2 + 2*A*C*k1 + pow(B, 2)*pow(k2,2) + 2*B*C*k2 + pow(C, 2);

        T delta = pow(N, 2) - 4 * M * O;

        //std::cout<<"M "<<M<<" N "<<N<<" O "<<O<<" delta "<<delta<<std::endl;

        if (delta < 0){
            return {-4, false};
        }

        T r1 = (-N + sqrt(delta)) / (2 * M);

        T r2 = (-N - sqrt(delta)) / (2 * M);

        //std::cout<<"r1 "<<r1<<" r2 "<<r2<<std::endl;

        if(r1 > 0){
            r1 = sqrt(r1);

            Point<T> inter1_r1 = cs.intersection_point(a_point.point_1, a_point.point_2, r1, a_point.is_reverse, a_point.is_vertical);
            Point<T> inter2_r1 = cs.intersection_point(b_point.point_1, b_point.point_2, r1, b_point.is_reverse, b_point.is_vertical);

            //std::cout<<"inter1_r1 "<<inter1_r1.x<<" "<<inter1_r1.y<<" inter2_r1 "<<inter2_r1.x<<" "<<inter2_r1.y<<std::endl;

            if(abs(inter1_r1.get(x_or_y) - inter2_r1.get(x_or_y)) < eps){
                Point<T> r1_moved1 = cs.intersection_point(a_point.point_1, a_point.point_2, r1 + 1, a_point.is_reverse, a_point.is_vertical);
                Point<T> r1_moved2 = cs.intersection_point(b_point.point_1, b_point.point_2, r1 + 1, b_point.is_reverse, b_point.is_vertical);

                return {r1, r1_moved1.get(x_or_y) > r1_moved2.get(x_or_y)};
            }
        }
        if( r2 > 0){
            r2 = sqrt(r2);

            Point<T> inter1_r2 = cs.intersection_point(a_point.point_1, a_point.point_2, r2, a_point.is_reverse, a_point.is_vertical);
            Point<T> inter2_r2 = cs.intersection_point(b_point.point_1, b_point.point_2, r2, b_point.is_reverse, b_point.is_vertical);

            //std::cout<<"inter1_r2 "<<inter1_r2.x<<" "<<inter1_r2.y<<" inter2_r2 "<<inter2_r2.x<<" "<<inter2_r2.y<<std::endl;

            if(abs(inter1_r2.get(x_or_y) - inter2_r2.get(x_or_y)) < eps){
                Point<T> r2_moved1 = cs.intersection_point(a_point.point_1, a_point.point_2, r2 + 1, a_point.is_reverse, a_point.is_vertical);
                Point<T> r2_moved2 = cs.intersection_point(b_point.point_1, b_point.point_2, r2 + 1, b_point.is_reverse, b_point.is_vertical);

                return {r2, r2_moved1.get(x_or_y) > r2_moved2.get(x_or_y)};
            }
        }
        return {-4, false};
    }else{
        //two envelope vertices in the same instance
        return {-2, false};
    }
}

template< typename T>
void c_sort<T>::move_comp( 
    bool first,
    Comparision<T>* comp,
    std::vector<Comparision<T>*> &comparisions, 
    std::vector<Comparision<T>*> &active_comparisions,
    std::vector<std::vector<std::vector<c_point<T>>>> &buckets, 
    int instance_i, int x_or_y){

    while(true){    
        std::pair<bool, T> comp_res = resolve(comp->a, comp->b, x_or_y);

        //std::cout<<"comp_res "<<comp_res.first<<" "<<comp_res.second<<std::endl;

        if(comp_res.second == -1){
            if((first && comp_res.first) || (!first && !comp_res.first)){
                if(comp->c_node->right != nullptr){
                    comp->b = comp->c_node->right->point;
                    comp->c_node = comp->c_node->right;
                    comp->weight /= 2;
                }else{
                    buckets[instance_i][comp->c_node->index].push_back(comp->a);
                    break;
                }
            }else{
                if(comp->c_node->left != nullptr){
                    comp->b = comp->c_node->left->point;
                    comp->c_node = comp->c_node->left;
                    comp->weight /= 2;
                }else{
                    int index = comp->c_node->index - 1;
                    index = index < 0 ? 0 : index;
                    buckets[instance_i][index].push_back(comp->a);
                    break;
                }
            }
        }else{
            comp->root = comp_res.second;
            //std::cout<<"root "<<comp->root<<std::endl;
            comparisions.push_back(comp);
            active_comparisions.push_back((comparisions.back()));
            break;
        }
    }
}

template< typename T>
Comparision <T> c_sort<T>::weighted_median(std::vector<Comparision<T>> comparisions, T total_weight){
    if(comparisions.size() == 1){
        return comparisions[0];
    }else if(comparisions.size() == 2){
        if (comparisions[0].weight > comparisions[1].weight){
            return comparisions[0];
        }else{
            return comparisions[1];
        }
    }
 
    //find the lower median
    std::nth_element(comparisions.begin(), comparisions.begin() + comparisions.size() / 2, comparisions.end(), [&](Comparision<T>& a, Comparision<T>& b){
        return a.root < b.root;
    });

    T median = comparisions[comparisions.size() / 2].weight;
    T median_root = comparisions[comparisions.size() / 2].root;


    //partiton around the median
    std::partition(comparisions.begin(), comparisions.end(), [&](Comparision<T> a){
        return a.root < median_root;
    });

    T sum_low = std::accumulate(comparisions.begin(), comparisions.begin() + comparisions.size() / 2, 0, [&](T sum, Comparision<T> a){
        return sum + a.weight;
    });

    T sum_high = std::accumulate(comparisions.begin() + comparisions.size() / 2, comparisions.end(), 0, [&](T sum, Comparision<T> a){
        return sum + a.weight;
    });


    if(sum_low < total_weight / 2 && sum_high < total_weight / 2){
        return comparisions[comparisions.size() / 2];
    }else{
        if(sum_low >= total_weight / 2){
            comparisions[comparisions.size() / 2].weight += sum_low;
            std::vector<Comparision<T>> new_comparisions(comparisions.begin(), comparisions.begin() + comparisions.size() / 2);
            return weighted_median(new_comparisions, total_weight);
        }else{
            comparisions[comparisions.size() / 2].weight += sum_high;
            std::vector<Comparision<T>> new_comparisions(comparisions.begin() + comparisions.size() / 2, comparisions.end());
            return weighted_median(new_comparisions, total_weight);
        }
    }
}