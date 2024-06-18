#include "grid/include/point.hpp"
#include "brut/include/brut.hpp"
#include <iostream>

int main(){
    int n, lambda, s_i, t_i;
    std::cin>>n>>lambda>>s_i>>t_i;
    
    std::vector<Point<long double>> points;
    for(int i = 0; i < n; i++){
        long double x, y;
        std::cin>>x>>y;
        points.push_back(Point<long double>{x, y});
    }

    brut<long double> brut;
    long double res = brut.solve(points, lambda, s_i, t_i);

    std::cout<<res<<std::endl;
}