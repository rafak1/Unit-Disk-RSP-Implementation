#pragma once
#include "../grid/include/point.hpp"
#include "../grid/include/grid.hpp"
#include <vector>
#include <algorithm>
#include <numeric>
#include <stack>

template <typename T>
class cs{
    public:
        bool solve(std::vector<Point<T>> input_points, T r, int s_i, int t_i, int lambda);
        bool solve(std::vector<Point<T>> input_points, std::vector<int> sorted_by_x, std::vector<int> sorted_by_y, T r, int s_i, int t_i, int lambda);
        std::vector<Point<T>> graham_scan(std::vector<Point<T>> &points, T r, T line, bool is_reverse, bool is_vertical);

    private:
        void prepare_distance(Grid<T> &grid, std::vector<std::vector<int>>& distance);
        bool run_bfs(Grid<T> &grid, Point<T> &start, Point<T> &end, T r, int lambda);
        bool do_step(Grid<T> &grid, std::vector<Point<T>> &S, std::vector<std::vector<int>>& distance, T r, Point<T> &end, int step);
        void solve_subproblem(std::vector<Point<T>>& C_points, std::vector<Point<T>>& C_prim_points, T r, std::vector<Point<T>>& result, bool is_reverse, bool is_vertical, T line, std::vector<std::vector<int>>& distance);

        bool do_intersect(Point<T> &p1, Point<T> &p2, T r, T line, bool is_reverse, bool is_vertical);
        Point<T> intersection_point(Point<T> &p1, Point<T> &p2, T r, bool is_reverse, bool is_vertical);
};

template <typename T>
bool cs<T>::solve(std::vector<Point<T>> input_points, T r, int s_i, int t_i, int lambda){
    std::vector<int> sorted_by_x(input_points.size());
    std::iota (std::begin(sorted_by_x), std::end(sorted_by_x), 0);
    std::sort(sorted_by_x.begin(), sorted_by_x.end(), [input_points](int a, int b) { 
        return input_points[a].x < input_points[b].x;
    });

    std::vector<int> sorted_by_y(input_points.size());
    std::iota (std::begin(sorted_by_y), std::end(sorted_by_y), 0);
    std::sort(sorted_by_y.begin(), sorted_by_y.end(), [input_points](int a, int b) { 
        return input_points[a].y < input_points[b].y;
    });

    Grid<T> grid(sorted_by_x, sorted_by_y,input_points, input_points[s_i], r);
    //grid.print();

    return run_bfs(grid, input_points[s_i], input_points[t_i], r, lambda);
}

template <typename T>
bool cs<T>::solve(std::vector<Point<T>> input_points, std::vector<int> sorted_by_x, std::vector<int> sorted_by_y, T r, int s_i, int t_i, int lambda){
    if( r <= 0){
        return false;
    }
    Grid<T> grid(sorted_by_x, sorted_by_y,input_points, input_points[s_i], r);
    //grid.print();

    return run_bfs(grid, input_points[s_i], input_points[t_i], r, lambda);
}

template <typename T>
bool cs<T>::run_bfs(Grid<T> &grid, Point<T> &start, Point<T> &end, T r, int lambda){
    std::vector<std::vector<int>> distance;
    prepare_distance(grid, distance);

    std::vector<Point<T>> S;

    distance[start.grid_cell][start.grid_index] = 0;

    S.push_back(start);

    for(int i = 0; i < lambda; i++){
        //std::cout<<"step\n";
        if(do_step(grid, S, distance, r, end, i)){
            //std::cout<<"found\n";
            return true;
        }else if(S.size() == 0){
            return false;
        }
    }
    return false;
}


template <typename T>
void cs<T>::prepare_distance(Grid<T> &grid, std::vector<std::vector<int>>& distance){
    distance.resize(grid.cells.size());
    for(int i = 0; i < grid.cells.size(); i++){
        distance[i].resize(grid.cells[i].points.size());
        for(int j = 0; j < grid.cells[i].points.size(); j++){
            distance[i][j] = -1;
        }
    }

    for(int i=0; i<grid.cells.size(); i++){
        grid.cells[i].visited = false;
    }
}

template <typename T>
bool cs<T>::do_step(Grid<T> &grid, std::vector<Point<T>> &S, std::vector<std::vector<int>>& distance, T r, Point<T> &end, int step){
    std::vector<Point<T>> new_S;
    std::vector<int> step_cells;

    //And save each cell that had a point in S_{i-1}, and add all other points from that cell to S_{i}
    for(int i = 0; i < S.size(); i++){
        Point<T> p = S[i];
        //std::cout<<"S point "<<p.x<<" "<<p.y<<"\n";

        bool should_add_cell = false;

        for(int i = 0; i < grid.cells[p.grid_cell].points.size(); i++){

            if(distance[p.grid_cell][i] == -1){
                distance[p.grid_cell][i] = step+1;
                new_S.push_back(grid.cells[p.grid_cell].points[i]);
                //std::cout<<"pushing (same cell)"<<grid.cells[p.grid_cell].points[i].x<<" "<<grid.cells[p.grid_cell].points[i].y<<"\n";
                should_add_cell = true;
            }else if(distance[p.grid_cell][i] == step-1 || distance[p.grid_cell][i] == step){ //second step conducted in this cell
                distance[p.grid_cell][i] = -step-2;
                should_add_cell = true;
            }
        }

        
        if(should_add_cell){
            //std::cout<<"adding\n";
            step_cells.push_back(p.grid_cell);
        }
    }

    /*std::cout<<"step cells\n";
    for(int i = 0; i < step_cells.size(); i++){
        std::cout<<step_cells[i]<<" | ";
    }
    std::cout<<"\n";*/

    //solve subproblem 1 (Decide whether the point should go into S_{i} or not)
    for(int i = 0; i < step_cells.size();i++){

        auto C = grid.cells[step_cells[i]];
        std::vector<Point<T>> C_points;
        std::vector<Point<T>> C_points_by_y;
        for(int k = 0; k < C.points.size(); k++){
            //std::cout<<"C point "<<C.points[k].x<<" "<<C.points[k].y<<" "<<distance[C.points[k].grid_cell][C.points[k].grid_index]<<" "<<step<<"\n";
            if(distance[C.points[k].grid_cell][C.points[k].grid_index] == -step-2){
                C_points.push_back(C.points[k]);
            }
            if(distance[C.points_by_y[k].grid_cell][C.points_by_y[k].grid_index] == -step-2){
                C_points_by_y.push_back(C.points_by_y[k]);
            }
        }

        /*for(int k=0;k<C_points.size();k++){
            std::cout<<"C_points "<<C_points[k].x<<" "<<C_points[k].y<<" | ";
        }
        std::cout<<"\n";

        for(int k=0;k<C_points_by_y.size();k++){
            std::cout<<"C_points_by_y "<<C_points_by_y[k].x<<" "<<C_points_by_y[k].y<<" | ";
        }
        std::cout<<"\n";*/

        //std::cout<<"step cell "<<step_cells[i]<<"\n";
        for(int j = 0; j < grid.cells[step_cells[i]].neighbors.size(); j++){
            
            auto neighbor = grid.cells[grid.cells[step_cells[i]].neighbors[j]];
            //std::cout<<"for neigbor "<<grid.cells[step_cells[i]].neighbors[j]<<" "<< grid.horizontal_lines[C.line_y].coordinate<<" "<< grid.vertical_lines[C.line_x].coordinate<<" "<<grid.grid_size<< "\n";

            //check the relative position of the neighbor
            if(neighbor.line_y > C.line_y){
                //standard
                solve_subproblem(C_points, neighbor.points, r, new_S, false,false, (grid.horizontal_lines[C.line_y].coordinate) + grid.grid_size, distance);
            }else if (neighbor.line_y == C.line_y){
                //vertical
                if(neighbor.line_x > C.line_x){
                    //right
                    solve_subproblem(C_points_by_y, neighbor.points_by_y, r, new_S, false,true, (grid.vertical_lines[C.line_x].coordinate) + grid.grid_size, distance);
                }else{
                    //left
                    std::vector<Point<T>> C_temp = C_points_by_y;
                    std::reverse(C_temp.begin(), C_temp.end());
                    std::vector<Point<T>> neighbor_temp = neighbor.points_by_y;
                    std::reverse(neighbor_temp.begin(), neighbor_temp.end());
                    solve_subproblem(C_temp, neighbor_temp, r, new_S, true,true, (grid.vertical_lines[C.line_x].coordinate), distance);
                }
            }else{
                //reverse
                std::vector<Point<T>> C_temp = C_points;
                std::reverse(C_temp.begin(), C_temp.end());
                std::vector<Point<T>> neighbor_temp = neighbor.points;
                std::reverse(neighbor_temp.begin(), neighbor_temp.end());

                solve_subproblem(C_temp, neighbor_temp, r, new_S, true,false, (grid.horizontal_lines[C.line_y].coordinate), distance);
            }
            
        }
    }
    
    //set distances of new_S
    for(int i = 0; i< new_S.size(); i++){
        distance[new_S[i].grid_cell][new_S[i].grid_index] = step+1;
    }

    //check whether we hit the end
    for(int i = 0; i< new_S.size(); i++){
        if(new_S[i].x == end.x && new_S[i].y == end.y){
            return true;
        }
    }


    /*std::cout<<"new_S\n";
    for(auto p : new_S){
        std::cout<<p.x<<" "<<p.y<<" | ";
    }
    std::cout<<"\n";*/

    S.clear();
    S = new_S;

    return false;
}


template <typename T>
void cs<T>::solve_subproblem(
    std::vector<Point<T>>& C_points, 
    std::vector<Point<T>>& C_prim_points, 
    T r, 
    std::vector<Point<T>>& result, 
    bool is_reverse,
    bool is_vertical,
    T line,
    std::vector<std::vector<int>>& distance
    ){
    std::vector<Point<T>> U;

    //calculate U
    std::vector<Point<T>> graham = graham_scan(C_points, r, line, is_reverse, is_vertical);
    
    /*std::cout<<"graham for line "<< line <<" reverse= "<< is_reverse << " vert= " << is_vertical <<"\n";
    for(auto p : graham){
        std::cout<<p.x<<" "<<p.y<<" | ";
    }
    std::cout<<"\n";*/


    for(int i = 0; i < graham.size()-1; i++){
        U.push_back(intersection_point(graham[i], graham[(i+1)], r, is_reverse, is_vertical));
    }

    //scan to determine wether the point should go into S_{i} or not

    int current_blue = 0; //index of C_prim.points
    int current_envelope = 0; //index of U and graham

    while(true){

        if(!is_vertical){
            while( current_envelope < U.size() &&
                ((!is_reverse && C_prim_points[current_blue].x > U[current_envelope].x) ||
                (is_reverse && C_prim_points[current_blue].x < U[current_envelope].x))){
                current_envelope++;
            }
        }else{
            while( current_envelope < U.size() &&
                ((!is_reverse && C_prim_points[current_blue].y > U[current_envelope].y) ||
                (is_reverse && C_prim_points[current_blue].y < U[current_envelope].y))){
                current_envelope++;
            }
        }

        //std::cout<<"point "<<C_prim_points[current_blue].x<<" "<<C_prim_points[current_blue].y<<" for envelope "<<current_envelope<<" "<<graham[current_envelope].x<<" "<<graham[current_envelope].y<<"\n";
        if(C_prim_points[current_blue].euclidean_distance(graham[current_envelope]) <= r && 
            distance[C_prim_points[current_blue].grid_cell][C_prim_points[current_blue].grid_index] == -1){
            
            result.push_back(C_prim_points[current_blue]);
            //std::cout<<"pushing "<<C_prim_points[current_blue].x<<" "<<C_prim_points[current_blue].y<<"\n";
        }

        current_blue++;
        if(current_blue == C_prim_points.size()){
            break;
        }
    }
}

template <typename T>
Point<T> next_to_top(std::stack<Point<T>> &S){
    Point<T> p = S.top();
    S.pop();
    Point<T> res = S.top();
    S.push(p);
    return res;
}

template <typename T>
std::vector<Point<T>> cs<T>::graham_scan(std::vector<Point<T>> &points, T r, T line, bool is_reverse, bool is_vertical){
    //points on input are sorted by x

    std::stack<Point<T>> S;

    //take care of points with the same x (only the top matters)
    int start = 0;
    while(start < points.size() && points[start].x == points[0].x){
        start++;
    }
    start = std::max(start - 1,0);

    S.push(points[start]);

    for(int i = start+1; i < points.size(); i++){

        int j = i;
        while(i < points.size() && points[i].x == points[i-1].x){
            i++;
        }
        i = std::max(i - 1,j);

        //std::cout<<"i "<<i<<std::endl;
        while(true){
            if(!do_intersect(S.top(), points[i], r, line, is_reverse, is_vertical)){
                Point<T> p = S.top();
                S.pop();
                if((!is_vertical && ((!is_reverse && p.y > points[i].y ) || (is_reverse && p.y < points[i].y))) ||
                    (is_vertical && ((!is_reverse && p.x > points[i].x) || (is_reverse && p.x < points[i].x)))){
                    S.push(p);
                }else{
                    S.push(points[i]);
                }
                break;
            }else{
                if(S.size() <= 1){
                    //std::cout<<"pushing "<<points[i].x<<" "<<points[i].y<<std::endl;
                    S.push(points[i]);
                    break;
                }else{
                    Point<T> next_next = next_to_top(S);
                    Point<T> v = intersection_point(S.top(), next_next, r, is_reverse, is_vertical);
                    Point<T> w = intersection_point(S.top(), points[i], r, is_reverse, is_vertical);
                    //std::cout<<"v "<<v.x<<" "<<v.y<<" w "<<w.x<<" "<<w.y<<std::endl;
                    if((!is_vertical && ((!is_reverse && w.x > v.x ) || (is_reverse && w.x < v.x))) ||
                        (is_vertical && ((!is_reverse && w.y > v.y) || (is_reverse && w.y < v.y)))){
                        S.push(points[i]);
                        //std::cout<<"pushing "<<points[i].x<<" "<<points[i].y<<std::endl;
                        break;
                    }else{
                        //std::cout<<"popping "<<S.top().x<<" "<<S.top().y<<std::endl;
                        S.pop();    
                    }
                }
            }
        }
    }

    std::vector<Point<T>> res;
    while(S.size() > 0){
        res.push_back(S.top());
        S.pop();
    }

    std::reverse(res.begin(), res.end());

    return res;
}

template <typename T>
bool cs<T>::do_intersect(Point<T> &p1, Point<T> &p2, T r, T line, bool is_reverse, bool is_vertical){
    Point<T> p = intersection_point(p1, p2, r, is_reverse, is_vertical);
    return (!is_vertical && ((!is_reverse && p.y >= line) || (is_reverse && p.y <= line))) ||
            (is_vertical && ((!is_reverse && p.x >= line) || (is_reverse && p.x <= line)));
} 

template <typename T>
Point<T> cs<T>::intersection_point(Point<T> &p1, Point<T> &p2, T r, bool is_reverse, bool is_vertical){
    T d = p1.euclidean_distance(p2) ;
    T h = sqrt(r*r - (d/2)*(d/2));

    T x1 = (p1.x + p2.x)/2 + h*(p1.y - p2.y)/d;
    T y1 = (p1.y + p2.y)/2 - h*(p1.x - p2.x)/d;

    T x2 = (p1.x + p2.x)/2 - h*(p1.y - p2.y)/d;
    T y2 = (p1.y + p2.y)/2 + h*(p1.x - p2.x)/d;

    //std::cout<<"intersection of "<<p1.x<<" "<<p1.y<<" and "<<p2.x<<" "<<p2.y<<" is "<<x1<<" "<<y1<<" or "<<x2<<" "<<y2<<"\n";

    if(!is_vertical){
        if((!is_reverse && y1 > y2) || (is_reverse && y1 < y2)){
            return Point<T>{x1, y1};
        }else{
            return Point<T>{x2, y2};
        }
    }else{
        if((!is_reverse && x1 > x2) || (is_reverse && x1 < x2)){
            return Point<T>{x1, y1};
        }else{
            return Point<T>{x2, y2};
        }
    }
}