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
    std::vector<Point<T>> points_cells;
    points_cells.reserve(sorted_by_x.size());

    //get all cells that have a point in them
    for(int i = 0; i < sorted_by_x.size(); i++){
        points_cells.push_back(sorted_by_x[i]);
    }

    //radix sort the cells to find duplicates
    radix_sort(points_cells, this->vertical_lines.size()+1);

    //remove duplicates and create lists of points for each cell
    for(int i = 0; i < points_cells.size(); i++){
        if(i == 0 || (points_cells[i].grid_x != points_cells[i-1].grid_x || points_cells[i].grid_y != points_cells[i-1].grid_y)){
            this->cells.push_back(Cell<T>{points_cells[i].grid_x, points_cells[i].grid_y});
        }
        this->cells.back().points.push_back(points_cells[i]);
    }

    //assign cell index to points
    for(int i = 0; i < this->cells.size(); i++){
        for(int j = 0; j < this->cells[i].points.size(); j++){
            this->cells[i].points[j].grid_cell = i;
        }
    }


    //calculate neighbors for each cell
    int neighbor_ll = 0, neighbor_l = 0, neighbor_r = 0, neighbor_rr = 0;

    for(int i = 0; i < this->cells.size(); i++){
        int curr_line_x = this->cells[i].line_x;
        int curr_line_y = this->cells[i].line_y;
        
        //move neighbor_rr and add neighbors on that line
        if(curr_line_x+2 < this->vertical_lines.size()){
            move_to_next_neighbot(this->cells, neighbor_rr, curr_line_x, curr_line_y, 2);
            add_neighbors(this->cells, 1,neighbor_rr, i, curr_line_y, curr_line_x);
        }

        //move neighbor_r and add neighbors on that line
        if(curr_line_x+1 < this->vertical_lines.size()){
            move_to_next_neighbot(this->cells, neighbor_r, curr_line_x, curr_line_y, 1);
            add_neighbors(this->cells, 2,neighbor_r, i, curr_line_y, curr_line_x);
        }

        //add neighbors on current line
        add_neighbors(this->cells, 2,i, i, curr_line_y, curr_line_x);

        //move neighbor_l and add neighbors on that line
        if(curr_line_x-1 >= 0){
            move_to_next_neighbot(this->cells, neighbor_l, curr_line_x, curr_line_y, -1);
            add_neighbors(this->cells, 2,neighbor_l, i, curr_line_y, curr_line_x);
        }

        //move neighbor_ll and add neighbors on that line
        if(curr_line_x-2 >= 0){
            move_to_next_neighbot(this->cells, neighbor_ll, curr_line_x, curr_line_y, -2);
            add_neighbors(this->cells, 1,neighbor_ll, i, curr_line_y, curr_line_x);
        }        
    }
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

template <typename T>
void radix_sort(std::vector<Point<T>> &points_cells, int n){

    bucket_sort(points_cells, n, 0);
    bucket_sort(points_cells, n, 1);
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
void move_to_next_neighbot(std::vector<Cell<T>> &cells, int &neighbor, int curr_line_x, int curr_line_y, int diff){
    while(true){
        if(cells[neighbor].line_x < curr_line_x+diff){
            neighbor++;
        }else if(cells[neighbor].line_x == curr_line_x+diff){
            if(neighbor + 1 >= cells.size()){
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
void move_to_next_neighbot(std::vector<Cell<T>> &cells, int &neighbor, int curr_line_x, int curr_line_y, int diff){
    while(true){
        if(cells[neighbor].line_x < curr_line_x+diff){
            neighbor++;
        }else if(cells[neighbor].line_x == curr_line_x+diff){
            if(neighbor + 1 >= cells.size()){
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
inline void add_neighbors(std::vector<Cell<T>> &cells, int r, int suspect,int i, int curr_line_y, int curr_line_x){
    for(int j = r; j <= r; j++){
        if(suspect + j >= 0 && suspect + j < this->cells.size() && suspect + j != i){
            if(this->cells[suspect+j] == curr_line_x && (this->cells[suspect+j].line_y - curr_line_y) <= r){
                this->cells[i].neighbors.push_back(suspect+j);
            }
        }
    }
}