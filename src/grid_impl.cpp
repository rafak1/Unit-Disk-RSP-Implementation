#include "grid.hpp"
#include "point.hpp"
#include <vector>
#include <list>
#include <cmath>


template <typename T>
Grid<T>::Grid(std::vector<Point<T>> &sorted_by_x, std::vector<Point<T>> &sorted_by_y, Point<T> &s, T r){
    
    this->grid_size = r/sqrt((T)2);

    this->vertical_lines = prepare_one_way_grid(sorted_by_x, s, r);
    this->horizontal_lines = prepare_one_way_grid(sorted_by_y, s, r);

    this->assign_grid_to_points(sorted_by_x, sorted_by_y);

    this->compute_cells(sorted_by_x, sorted_by_y);
}

template <typename T>
void Grid<T>::compute_cells(std::vector<Point<T>> &sorted_by_x, std::vector<Point<T>> &sorted_by_y){
    std::vector<std::pair<T,T>> points;
    points.reserve(sorted_by_x.size());

    //get all cells that have a point in them
    for(int i = 0; i < sorted_by_x.size(); i++){
        points.push_back(std::pair<T,T>{sorted_by_x[i].grid_x, sorted_by_x[i].grid_y});
    }

    //radix sort the cells to find duplicates

    

}


template <typename T>
void Grid<T>::assign_grid_to_points(std::vector<Point<T>> &sorted_by_x, std::vector<Point<T>> &sorted_by_y){

    //assign vertical grid coordinates to points
    int current_line = 0;
    for(int i = 0; i < points.size(); i++){
        while(points[i].x >= this->vertical_lines[current_line+1].coordinate){
            current_line++;
        }
        points[i].grid_x = current_line;
    }

    //assign horizontal grid coordinates to points
    current_line = 0;
    for(int i = 0; i < points.size(); i++){
        while(points[i].y >= this->horizontal_lines[current_line+1].coordinate){
            current_line++;
        }
        points[i].grid_y = current_line;
    }
}




template <typename T>
void Grid<T>::prepare_one_way_grid(std::vector<Point<T>> points, Point<T> &s, T r){    

    int p1_start = 0;
    int p2_end = 0;

    //binsearch TODO
    //find points to the left, and to the right of s
    for(int i = 0; i < sorted_by_x.size(); i++){
        if(sorted_by_x[i].x > s.x){
            p1_start = i;
            p2_end = i-1;
            break;
        }
    }

    int i1 = p1_start;

    //exclude 'right' points that are not connected to s
    while(i1 < sorted_by_x.size()-1 && sorted_by_x[i1+1].x - sorted_by_x[i1].x <= r){
        i1++;
    }

    int i2 = p2_end;

    //exclude 'left' points that are not connected to s
    while(i2 >= 1 && sorted_by_x[i2].x - sorted_by_x[i2-1].x <= r){
        i2--;
    }

    std::vector<Line<T>> lines;
    lines.reserve(2*points.size());

    std::vector<Line<T>> left_lines;
    long double current_line = s.x;
    left_lines.reserve(2*(p2_end-i2+1));

    //create lines to the left of s
    for(int j = i2; j <= p2_end; j++){
        current_line -= this->grid_size;
        left_lines.push_back(Line<T>{current_line});
    }

    //reverse left_lines into lines
    for(int j = left_lines.size()-1; j >= 0; j--){
        lines.push_back(left_lines[j]);
    }

    lines.push_back(Line<T>{s.x});

    current_line = lines[0].coordinate;

    //add lines to the right of s
    for(int j = p1_start; j <= i; j++){
        current_line += this->grid_size;
        lines.push_back(Line<T>{current_line});
    }
    
    return lines;
}

