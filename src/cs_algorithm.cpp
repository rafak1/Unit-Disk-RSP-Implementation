
#include "point.hpp"
#include "grid.hpp"
#include <vector>
#include <algorithm>



// need description
long double cs_solve(std::vector<Point<long double>> &input_points, long double r, int s_i, int t_i){

    std::vector<Point<long double>> sorted_by_x = input_points;
    std::sort(sorted_by_x.begin(), sorted_by_x.end(), [](Point<long double> a, Point<long double> b) { return a.x < b.x; });

    std::vector<Point<long double>> sorted_by_y = input_points;
    std::sort(sorted_by_y.begin(), sorted_by_y.end(), [](Point<long double> a, Point<long double> b) { return a.y < b.y; });

    Grid<long double> grid(sorted_by_x, sorted_by_y, input_points[s_i], r);


}




