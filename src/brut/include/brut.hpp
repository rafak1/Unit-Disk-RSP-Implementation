#pragma once
#include "../grid/include/point.hpp"
#include <vector>
#include <algorithm>
#include <cmath>
#include <queue>
#include <iostream>


template <typename T>
class brut{
    public:
        /**
         * @brief Solves the Reverse Shortest Path Problem by brute force.
         * 
         * @param points centers of the circles
         * @param lambda at most length of the path between s_i and t_i
         * @param s_i index of the source point
         * @param t_i index of the target point
         * @return T minimal radius of the circles, so that the shortest path between s_i and t_i is at most lambda
         */
        T solve(std::vector<Point<T>> points, int lambda, T s_i, T t_i);

    private:
        T decision_algorithm(std::vector<Point<T>> points, int lambda, T r, int s_i, int t_i);
};



template <typename T>
T brut<T>::solve(std::vector<Point<T>> points, int lambda, T s_i, T t_i){
    
    //calculate all possible radii
    std::vector<T> possible_radii;
    for(int i = 0; i < points.size(); i++){
        for(int j = 0; j < points.size(); j++){
            if(i == j) continue;
            T dist = sqrt(pow(points[i].x - points[j].x, 2) + pow(points[i].y - points[j].y, 2));
            possible_radii.push_back(dist);
        }
    }

    //sort radii
    std::sort(possible_radii.begin(), possible_radii.end());


    //binsearch
    int l = 0;
    int r = possible_radii.size();
    while(l < r){
        int m = (l + r) / 2;
        if(decision_algorithm(points, lambda, possible_radii[m], s_i, t_i)){
            r = m;
        }else{
            l = m + 1;
        }
    }

    return possible_radii[l];
}


template <typename T>
T brut<T>::decision_algorithm(std::vector<Point<T>> points, int lambda, T r, int s_i, int t_i){

    //construct graph
    std::vector<std::vector<T>> graph(points.size(), std::vector<T>());

    for(int i = 0; i < points.size(); i++){
        for(int j = 0; j < points.size(); j++){
            if(i == j) continue;
            T dist = points[i].euclidean_distance(points[j]);
            if(dist <= r){
                graph[i].push_back(j);
            }
        }
    }

    //bfs
    std::vector<T> dist(points.size(), -1);
    dist[s_i] = 0;
    std::queue<int> queue;
    queue.push(s_i);
    while(!queue.empty()){
        int u = queue.front();
        queue.pop();
        for(int i = 0; i < graph[u].size(); i++){
            int v = graph[u][i];
            if(dist[v] == -1){
                dist[v] = dist[u] + 1;
                queue.push(v);
            }
        }
    }

    return dist[t_i] != -1 && dist[t_i] <= lambda;
}