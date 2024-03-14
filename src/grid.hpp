#pragma once
#include <vector>
#include "point.hpp"

template <typename T>
class Line{
    public:
        T coordinate;
};

template <typename T>
class Cell{
    public:
        int line_x, line_y;
        std::vector<Point<T>> points;
};


template <typename T>
class Grid {
    public:
        std::vector<Line<T>> vertical_lines;
        std::vector<Line<T>> horizontal_lines;
        std::vector<Cell<T>> cells;


        /**
         * @brief Construct a new Grid object on a set on given points. If the points are not sorted, or if the points are not the same in both vectors, the behavior is undefined.
         *
         * @param sorted_by_x points sorted by x coordinate
         * @param sorted_by_y points sorted by y coordinate
         * @param s_i the index of the point that will be the starting point of the algorithm
         * @param r the radius of the circle
         */
        Grid(std::vector<Point<T>> &sorted_by_x, std::vector<Point<T>> &sorted_by_y, Point<T> &s, T r);





    private:
        const T grid_size;
        void prepare_one_way_grid(std::vector<Point<T>> points, Point<T> &s, T r);
        void assign_grid_to_points(std::vector<Point<T>> &sorted_by_x, std::vector<Point<T>> &sorted_by_y);
        void compute_cells(std::vector<Point<T>> &sorted_by_x, std::vector<Point<T>> &sorted_by_y);

};