#include "grid/include/point.hpp"
#include "alg_1/include/alg_1.hpp"
#include <iostream>
#include <chrono>

int PRECISION = 20;

int main(){
    int n, lambda, s_i, t_i;
    std::cin>>n>>lambda>>s_i>>t_i;
    
    std::vector<Point<long double>> points;
    for(int i = 0; i < n; i++){
        long double x, y;
        std::cin>>x>>y;
        points.push_back(Point<long double>{x, y});
    }

    alg_1<long double> alg;

    auto t1 = std::chrono::high_resolution_clock::now();
    long double res = alg.solve(points, lambda, s_i, t_i, PRECISION);
    auto t2 = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();

    std::cout<<"r^* = "<<res<<" ; calculated in: "<<duration<<" ms; with precision: "<<PRECISION<<std::endl;
}