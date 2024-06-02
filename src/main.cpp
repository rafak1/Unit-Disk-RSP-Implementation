#include "grid/include/grid.hpp"
#include "grid/include/point.hpp"
#include "alg_1/include/alg_1.hpp"
#include "alg_2/include/alg_2.hpp"
#include "brut/include/brut.hpp"
#include "generator/include/generator.hpp"
#include "cs/include/cs.hpp"
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>

#include <iostream>
int main(){
    std::vector<Point<long double>> points = {
        Point<long double>{0, 0}, 
        Point<long double>{0.1, 0.1}, 
        Point<long double>{-1, 1}, 
        Point<long double>{1, 1}, 
        Point<long double>{2, 2}, 
        Point<long double>{3, 3},
        Point<long double>{2, 1}, 
        Point<long double>{0, 1},
        Point<long double>{0.6, 0.4},
        };

    alg_2<long double> alg;
    long double res = alg.solve(points, 10, 0, 4);

    std::cout<<" ANSWER: "<<res<<std::endl;
    return 0;

   /*std::vector<Point<long double>> points = {
        Point<long double>{0.101366, -5.86374}, 
        Point<long double>{4.93276, -0.217757}, 
        Point<long double>{8.54065, -2.63342}, 
        Point<long double>{8.17979, 6.83221}, 
        Point<long double>{2.86956, 1.10816}, 
        Point<long double>{6.68214, -3.50696},
        Point<long double>{5.03714, -0.0729278}, 
        Point<long double>{-9.57615, -3.76574},
        Point<long double>{-7.16822, 0.758335},
        Point<long double>{1.16247, -3.33508}
        };*/
    /*std::vector<Point<long double>> points = {
        Point<long double>{0.306317, 6.34287},
        Point<long double>{4.61501, -7.1088},
        Point<long double>{-4.88258, -4.20902},
        Point<long double>{6.41196, 6.62846},
        Point<long double>{-2.35905, -4.86236},
        Point<long double>{-3.6755, 0.371716},
        Point<long double>{4.60791, 9.62049},
        Point<long double>{-1.98389, 8.65976},
        Point<long double>{8.20154, 2.77801},
        Point<long double>{9.33924, -9.72999}
        };*/
    /*std::vector<Point<long double>> points = {
        Point<long double>{-5.72779, 8.3437},
        Point<long double>{-9.66477, 0.909673},
        Point<long double>{9.71026, -5.19153},
        Point<long double>{-7.05702, -2.56154},
        Point<long double>{8.44864, 2.11344},
        Point<long double>{-7.17252, 2.2892},
        Point<long double>{-4.67754, 7.80791},
        Point<long double>{0.159981, 0.749154},
        Point<long double>{-1.80044, -1.34415},
        Point<long double>{0.199616, -3.2837}
        };*/
    /*std::vector<Point<long double>> points = {
        Point<long double>{8.9393, -9.39539},
        Point<long double>{-4.19049, -7.81817},
        Point<long double>{9.21043, -8.70884},
        Point<long double>{9.17925, -6.77818},
        Point<long double>{5.53497, -6.48097}
        };

    cs<long double> cs;
    long double res3 = cs.solve(points, 10, 2, 1, 7);
    std::cout<<" ANSWER: "<<res3<<std::endl;
    return 0;
    brut<long double> brut2;
    long double resa = brut2.solve(points, 4,2, 0);

    alg_1<long double> alg2;
    long double resa2 = alg2.solve(points, 4,2, 0, 10);
    std::cout<<" ANSWER: brut: "<<resa<<" algo: "<<resa2<<std::endl;

    return 0;*/

    while(true){

        generator<long double> gen;
        dataset<long double> d = gen.generate_dataset(10, -100, 100, 10);


        brut<long double> brut;
        long double res = brut.solve(d.points, d.lambda, d.s_i, d.t_i);

        alg_1<long double> alg;
        long double res2 = alg.solve(d.points, d.lambda, d.s_i, d.t_i, 100);

        if(res != res2){
            std::cout<<" ANSWER: brut: "<<res<<" algo: "<<res2<<std::endl;
            for(auto p : d.points){
                std::cout << p.x << " " << p.y << std::endl;
            }
            std::cout << d.lambda << std::endl;
            std::cout << d.s_i << std::endl;
            std::cout << d.t_i << std::endl;
            break;
        }
    }
    
    //std::cout<<" ANSWER: brut: "<<res<<" algo: "<<res2<<std::endl;

    return 0;
}




