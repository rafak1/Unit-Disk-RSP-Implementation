#include "grid/include/grid.hpp"
#include "grid/include/point.hpp"
#include <vector>
#include <algorithm>
#include <numeric>
#include "cs/include/cs.hpp"

#include <iostream>
int main(){
    std::vector<Point<long double>> points = {
        Point<long double>{0, 0}, 
        Point<long double>{-1, 1}, 
        Point<long double>{1, 1}, 
        Point<long double>{2, 2}, 
        Point<long double>{3, 3},
        Point<long double>{2, 1}, 
        Point<long double>{0, 1}
        };
    /*
    Vertical lines: -1.41421 -0.707107 0 0.707107 1.41421 2.12132 2.82843 3.53553
    Horizontal lines: -0.707107 0 0.707107 1.41421 2.12132 2.82843 3.53553
    Cells:
    Cell: 0 2
    Points: -1 1
    Neighbors: 1 2
    Cell: 2 1
    Points: 0 0
    Neighbors: 4 3 2 0
    Cell: 2 2
    Points: 0 1
    Neighbors: 4 5 3 1 0
    Cell: 3 2
    Points: 1 1
    Neighbors: 4 5 1 2
    Cell: 4 2
    Points: 2 1
    Neighbors: 5 3 1 2
    Cell: 4 3
    Points: 2 2
    Neighbors: 4 3 2
    Cell: 6 5
    Points: 3 3
    Neighbors:*/

    cs<long double> cs;
    cs.solve(points, 1, 0, 4);
    return 0;
}




