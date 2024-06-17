#include "generator/include/generator.hpp"
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>

int main(){
    int n, lambda;
    std::cin>>n>>lambda;

    generator<long double> gen;
    dataset<long double> d = gen.generate_dataset(n, -100, 100, lambda);

    std::cout<<n<<" "<<d.lambda<<" "<<d.s_i<<" "<<d.t_i<<std::endl;
    for(auto p : d.points){
        std::cout << p.x << " " << p.y << std::endl;
    }
}