#pragma once
#include <vector>
#include "point.hpp"
#include "point.hpp"
#include <vector>
#include <list>
#include <cmath>
#include <iostream>



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
        std::vector<Point<T>> points_by_y;
        std::vector<int> neighbors;
        bool visited = false;

        /**
         * @brief prints important information about the cell on the standard output
         * 
         */
        void print();
};


template <typename T>
class Grid {
    public:
        std::vector<Line<T>> vertical_lines;
        std::vector<Line<T>> horizontal_lines;
        std::vector<Cell<T>> cells;
        T grid_size;

        /**
         * @brief Construct a new Grid object on a set on given points. If the points are not sorted, or if the points are not the same in both vectors, the behavior is undefined.
         *
         * @param sorted_by_x points sorted by x coordinate
         * @param sorted_by_y points sorted by y coordinate
         * @param s_i the index of the point that will be the starting point of the algorithm
         * @param r the radius of the circle
         */
        Grid(std::vector<int> &sorted_by_x, std::vector<int> &sorted_by_y,std::vector<Point<T>>& points, Point<T> &s, T r);

        /**
         * @brief prints important information about the grid on the standard output
         * 
         */
        void print();



    private:
        std::vector<Line<T>> prepare_one_way_grid(std::vector<int>& i_points, std::vector<Point<T>>& points, Point<T> &s, T r, int x_or_y);
        void assign_grid_to_points(std::vector<int> &sorted_by_x, std::vector<int> &sorted_by_y, std::vector<Point<T>>& points);
        void compute_cells(std::vector<int> &sorted_by_x, std::vector<int> &sorted_by_y, std::vector<Point<T>>& points);
        void add_neighbors(std::vector<Cell<T>> &cells, int max_distance, int i, int j, int curr_line_y, int curr_line_x);

};


template <typename T>
Grid<T>::Grid(std::vector<int> &sorted_by_x, std::vector<int> &sorted_by_y, std::vector<Point<T>>& points, Point<T> &s, T r){

    this->grid_size = r/sqrt((T)2);
    std::cout<<"building grid with grid size: "<<this->grid_size<<std::endl;

    this->vertical_lines = prepare_one_way_grid(sorted_by_x, points, s, r, 0);
    for(int i = 0; i < this->vertical_lines.size(); i++){
        std::cout<<this->vertical_lines[i].coordinate<<" ";
    }
    std::cout<<std::endl;

    this->horizontal_lines = prepare_one_way_grid(sorted_by_y, points, s, r, 1);
    for(int i = 0; i < this->horizontal_lines.size(); i++){
        std::cout<<this->horizontal_lines[i].coordinate<<" ";
    }
    std::cout<<std::endl;

    this->assign_grid_to_points(sorted_by_x, sorted_by_y, points);
    for(int i = 0; i < sorted_by_x.size(); i++){
        std::cout<<points[sorted_by_x[i]].x<<" "<<points[sorted_by_x[i]].y<<" => "<<points[sorted_by_x[i]].grid_x<<" "<<points[sorted_by_x[i]].grid_y<<" | ";
    }
    std::cout<<std::endl;


    this->compute_cells(sorted_by_x, sorted_by_y, points);

    std::cout<<"end"<<std::endl;
}

//TODO redundant argument
template <typename T>
void Grid<T>::compute_cells(std::vector<int> &sorted_by_x, std::vector<int> &sorted_by_y, std::vector<Point<T>>& points){
    std::vector<Point<T>> points_cells_by_x;
    points_cells_by_x.reserve(points.size());

    //get all cells that have a point in them
    for(int i = 0; i < points.size(); i++){
        if(points[sorted_by_x[i]].grid_x == -1 || points[sorted_by_x[i]].grid_y == -1)
            continue;

        points[sorted_by_x[i]].original_index = sorted_by_x[i];
        points_cells_by_x.push_back(points[sorted_by_x[i]]);
    }


    //radix sort the cells to find duplicates
    radix_sort(points_cells_by_x, std::max(this->vertical_lines.size()+1, this->horizontal_lines.size()+1));

    for(int i = 0; i < points_cells_by_x.size(); i++){
        std::cout<<points_cells_by_x[i].x<<" "<<points_cells_by_x[i].y<<" "<<points_cells_by_x[i].grid_x<<" "<<points_cells_by_x[i].grid_y<<std::endl;
    }

    //remove duplicates and create lists of points for each cell
    for(int i = 0; i < points_cells_by_x.size(); i++){
        if(i == 0 || (points_cells_by_x[i].grid_x != points_cells_by_x[i-1].grid_x || points_cells_by_x[i].grid_y != points_cells_by_x[i-1].grid_y)){
            this->cells.push_back(Cell<T>{points_cells_by_x[i].grid_x, points_cells_by_x[i].grid_y});
        }
        this->cells.back().points.push_back(points_cells_by_x[i]);
    }

    //assign cell index to points
    for(int i = 0; i < this->cells.size(); i++){
        for(int j = 0; j < this->cells[i].points.size(); j++){
            this->cells[i].points[j].grid_cell = i;
            points[this->cells[i].points[j].original_index].grid_cell = i;
            this->cells[i].points[j].grid_index = j;
            points[this->cells[i].points[j].original_index].grid_index = j;
        }
    }

    //add points to cells sorted by y
    for(int i = 0; i < sorted_by_y.size(); i++){
        this->cells[points[sorted_by_y[i]].grid_cell].points_by_y.push_back(points[sorted_by_y[i]]);
        //std::cout<<"point "<<points[sorted_by_y[i]].x<<" "<<points[sorted_by_y[i]].y<<"was added to cell "<<points[sorted_by_y[i]].grid_cell<<std::endl;
    }


    //calculate neighbors for each cell
    int neighbor_ll = 0, neighbor_l = 0, neighbor_r = 0, neighbor_rr = 0;

    for(int i = 0; i < this->cells.size(); i++){
        int curr_line_x = this->cells[i].line_x;
        int curr_line_y = this->cells[i].line_y;
        
        //move neighbor_rr and add neighbors on that line
        if(curr_line_x+2 < this->vertical_lines.size()){
            move_to_next_neighbor(this->cells, neighbor_rr, curr_line_x, curr_line_y, 2);
            add_neighbors(this->cells, 1,neighbor_rr, i, curr_line_y, curr_line_x+2);
        }

        //move neighbor_r and add neighbors on that line
        if(curr_line_x+1 < this->vertical_lines.size()){
            move_to_next_neighbor(this->cells, neighbor_r, curr_line_x, curr_line_y, 1);
            add_neighbors(this->cells, 2,neighbor_r, i, curr_line_y, curr_line_x+1);
        }

        //add neighbors on current line
        add_neighbors(this->cells, 2,i, i, curr_line_y, curr_line_x);

        //move neighbor_l and add neighbors on that line
        if(curr_line_x-1 >= 0){
            move_to_next_neighbor(this->cells, neighbor_l, curr_line_x, curr_line_y, -1);
            add_neighbors(this->cells, 2,neighbor_l, i, curr_line_y, curr_line_x-1);
        }

        //move neighbor_ll and add neighbors on that line
        if(curr_line_x-2 >= 0){
            move_to_next_neighbor(this->cells, neighbor_ll, curr_line_x, curr_line_y, -2);
            add_neighbors(this->cells, 1,neighbor_ll, i, curr_line_y, curr_line_x-2);
        }        

        //std::cout<<i<<" "<<curr_line_x<<" "<<curr_line_y<<" "<<neighbor_ll<<" "<<neighbor_l<<" "<<i<<" "<<neighbor_r<<" "<<neighbor_rr<<std::endl;
    }
}


template <typename T>
void Grid<T>::assign_grid_to_points(std::vector<int> &sorted_by_x, std::vector<int> &sorted_by_y, std::vector<Point<T>>& points){

    //assign vertical grid coordinates to points
    int current_line = 0;
    for(int i = 0; i < sorted_by_x.size(); i++){
        while(current_line < this->vertical_lines.size()-1 && points[sorted_by_x[i]].x >= this->vertical_lines[current_line+1].coordinate){
            current_line++;
        }
        if(points[sorted_by_x[i]].x >= this->vertical_lines.back().coordinate){ //TODO +grid_size?
            points[sorted_by_x[i]].grid_x = -1;
        }else{
            points[sorted_by_x[i]].grid_x = current_line;
        }
    }

    //assign horizontal grid coordinates to points
    current_line = 0;
    //std::cout<<horizontal_lines[0].coordinate<<std::endl;
    for(int i = 0; i < sorted_by_y.size(); i++){
        while(current_line < this->horizontal_lines.size()-1 && points[sorted_by_y[i]].y >= this->horizontal_lines[current_line+1].coordinate){
            current_line++;
        }
        if(points[sorted_by_y[i]].y >= this->horizontal_lines.back().coordinate){
            points[sorted_by_y[i]].grid_y = -1;
        }else{
            points[sorted_by_y[i]].grid_y = current_line;
        }
    }
}




template <typename T>
std::vector<Line<T>> Grid<T>::prepare_one_way_grid(std::vector<int>& i_points, std::vector<Point<T>>& points, Point<T> &s, T r, int x_or_y){

    int p1_start = 0;
    int p2_end = 0;

   /*for(int i=0;i<points.size();i++){
        std::cout<<points[i_points[i]].x<<" "<<points[i_points[i]].y<< " | ";
    }
    std::cout<<std::endl;*/

    //binsearch TODO
    //find points to the left, and to the right of s
    for(int i = 0; i < points.size(); i++){
        if(points[i_points[i]].x == s.x && points[i_points[i]].y == s.y){
            p1_start = i;
            p2_end = i;
            break;
        }
    }

    //std::cout<<p1_start<<" "<<p2_end<<" starting points"<<std::endl;

    int i1 = p1_start;

    //exclude 'right' points that are not connected to s
    while(i1 < points.size()-1 && points[i_points[i1+1]].get(x_or_y) - points[i_points[i1]].get(x_or_y) <= r){
        i1++;
    }

    int i2 = p2_end;

    //exclude 'left' points that are not connected to s
    while(i2 >= 1 && points[i_points[i2]].get(x_or_y) - points[i_points[i2-1]].get(x_or_y) <= r){
        i2--;
    }

    //std::cout<<i1<<" "<<i2<<" points after"<<std::endl;

    std::vector<Line<T>> lines;
    lines.reserve(2*points.size());

    std::vector<Line<T>> left_lines;
    long double current_line = s.get(x_or_y);
    left_lines.reserve(2*(p2_end-i2+1));

    //create lines to the left of s
    while(current_line >= points[i_points[i2]].get(x_or_y)){
        current_line -= this->grid_size;
        left_lines.push_back(Line<T>{current_line});
    }

    //reverse left_lines into lines
    for(int j = left_lines.size()-1; j >= 0; j--){
        lines.push_back(left_lines[j]);
    }

    //add line on s
    lines.push_back(Line<T>{s.get(x_or_y)});

    current_line = s.get(x_or_y);

    //add lines to the right of s
    while(current_line <= points[i_points[i1]].get(x_or_y)){
        current_line += this->grid_size;
        lines.push_back(Line<T>{current_line});
    }

    for(int x = 0; x < lines.size(); x++){
        std::cout<<lines[x].coordinate<<" ";
    }
    std::cout<<std::endl;

    return lines;
}

template <typename T>
void radix_sort(std::vector<Point<T>> &points_cells, int n){
    bucket_sort(points_cells, n, 1);
    bucket_sort(points_cells, n, 0);
}

template <typename T>
void bucket_sort(std::vector<Point<T>> &points_cells, int n, int x_or_y){
    std::vector<std::list<Point<T>>> buckets(n , std::list<Point<T>>());

    for(int i = 0; i < points_cells.size(); i++){
        int index = x_or_y == 0 ? points_cells[i].grid_x : points_cells[i].grid_y;
        buckets[index].push_back(points_cells[i]);
    }

    int index = 0;
    for(int i = 0; i < n; i++){
        for(auto it = buckets[i].begin(); it != buckets[i].end(); it++){
            points_cells[index] = *it;
            index++;
        }
    }
}


template <typename T>
void move_to_next_neighbor(std::vector<Cell<T>> &cells, int &neighbor, int curr_line_x, int curr_line_y, int diff){
    //std::cout<<"move "<<neighbor<<" "<<curr_line_x<<" "<<curr_line_y<<" "<<diff<<std::endl;
    while(true){
        if(neighbor < cells.size() && cells[neighbor].line_x < curr_line_x+diff){
            neighbor++;
        }else if(cells[neighbor].line_x == curr_line_x+diff){
            if(neighbor + 1 >= cells.size() || cells[neighbor+1].line_x != curr_line_x+diff){
                break;
            }else{
                int curr_diff = abs(cells[neighbor].line_y - curr_line_y);
                int next_diff = abs(cells[neighbor+1].line_y - curr_line_y);
                if(next_diff < curr_diff){
                    neighbor++;
                }else{
                    break;
                }
            }
        }else{
            break;
        }
    }
}

template <typename T>
inline void Grid<T>::add_neighbors(std::vector<Cell<T>> &cells, int r, int suspect,int i, int curr_line_y, int curr_line_x){
    //std::cout<<"add neighbors: "<<i<<" "<<suspect<<" "<<curr_line_x<<" "<<curr_line_y<<" -r-"<<r<<std::endl;
    for(int j = -r; j <= r; j++){
        if(suspect + j >= 0 && suspect + j < this->cells.size() && suspect + j != i){
            //std::cout<<this->cells[suspect+j].line_x<<" "<<this->cells[suspect+j].line_y<<std::endl;
            if(this->cells[suspect+j].line_x == curr_line_x && abs(this->cells[suspect+j].line_y - curr_line_y) <= r){
                this->cells[i].neighbors.push_back(suspect+j);
                //std::cout<<"added: "<<suspect+j<<std::endl;
            }
        }
    }
}


template <typename T>
void Grid<T>::print(){
    std::cout<<"Vertical lines: ";
    for(int i = 0; i < this->vertical_lines.size(); i++){
        std::cout<<this->vertical_lines[i].coordinate<<" ";
    }
    std::cout<<std::endl;

    std::cout<<"Horizontal lines: ";
    for(int i = 0; i < this->horizontal_lines.size(); i++){
        std::cout<<this->horizontal_lines[i].coordinate<<" ";
    }
    std::cout<<std::endl;

    std::cout<<"Cells: "<<std::endl;
    for(int i = 0; i < this->cells.size(); i++){
        this->cells[i].print();
    }
}

template <typename T>
void Cell<T>::print(){
    std::cout<<"Cell: "<<this->line_x<<" "<<this->line_y<<std::endl;
    std::cout<<"Points: ";
    for(int i = 0; i < this->points.size(); i++){
        std::cout<<this->points[i].x<<" "<<this->points[i].y<<" ";
    }
    std::cout<<std::endl;

    std::cout<<"Points by y: ";
    for(int i = 0; i < this->points_by_y.size(); i++){
        std::cout<<this->points_by_y[i].x<<" "<<this->points_by_y[i].y<<" ";
    }
    std::cout<<std::endl;

    std::cout<<"Neighbors: ";
    for(int i = 0; i < this->neighbors.size(); i++){
        std::cout<<this->neighbors[i]<<" ";
    }
    std::cout<<std::endl;
}