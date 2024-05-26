#pragma once
#include <../grid/include/point.hpp>
#include <../cs/include/cs.hpp>
#include "../grid/include/grid.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>
#include <limits>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel                  K;
typedef CGAL::Delaunay_triangulation_2<K>                                    DT;
typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT>                 AT;
typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT> AP;
typedef CGAL::Voronoi_diagram_2<DT,AT,AP>                                    VD;
typedef K::Point_2                                                           Point_2;

template <typename T>
class alg_2 {
    public:
        /**
         * @brief Solves the Reverse Shortest Path Problem by binary search on real numbers using the CS algorithm as a decision algorithm.
         * 
         * @param points Center of the circles
         * @param lambda At most length of the path between s_i and t_i
         * @param s_i Index of the source point
         * @param t_i Index of the target point
         * @return T The minimal radius of the circles, so that the shortest path between s_i and t_i is at most lambda
         */
        T solve(std::vector<Point<T>> points, int lambda, int s_i, int t_i);
        
        alg_2() : my_cs() {}

    private:
        T interval_l, interval_r;
        Grid<T> grid;
        std::vector<Point<T>> input_points;
        std::vector<int> sorted_by_x;
        std::vector<int> sorted_by_y;
        cs<T> my_cs;
        int s_i, t_i, lambda;

        void shrink(T new_l, T new_r);
        void build_grid(int s_i, int t_i, int lambda);
        T get_value(int i, int j, int n, int s_pos, int size, int x_or_y, std::vector<int>& sorted);
        T matrix_search(int n, int s_pos, int s_i, int t_i, int lambda, bool is_vertical, bool look_for_one, std::vector<int>& sorted);
        void run_bfs(Point<T> &start, Point<T> &end, int lambda);
        void prepare_distance(std::vector<std::vector<int>>& distance);
        void pre_step(std::vector<Point<T>> &S, std::vector<std::vector<int>>& distance, Point<T> &end, int step);
        bool do_step(std::vector<Point<T>> &S, std::vector<std::vector<int>>& distance, Point<T> &end, int step);
};



template <typename T>
T alg_2<T>::solve(std::vector<Point<T>> input_points, int lambda, int s_i, int t_i){
    
    // set starting interval
    interval_l = 0;
    interval_r = std::numeric_limits<T>::max();
    this->s_i = s_i;
    this->t_i = t_i;
    this->lambda = lambda;

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

    this->input_points = input_points;
    this->sorted_by_x = sorted_by_x;
    this->sorted_by_y = sorted_by_y;

    // build the grid, and shrink the interval accordingly
    build_grid(s_i, t_i, lambda);

    Point<T> start = input_points[s_i];
    Point<T> end = input_points[t_i];

    run_bfs(start, end, lambda);

    return 0;
}

template <typename T>
void alg_2<T>::shrink(T new_l, T new_r){
    if(new_l != -1)
        interval_l = std::max(interval_l, new_l);

    if(new_r != -1)
        interval_r = std::min(interval_r, new_r);
}


template <typename T>
void alg_2<T>::run_bfs(Point<T> &start, Point<T> &end, int lambda){
    std::vector<std::vector<int>> distance;
    prepare_distance(grid, distance);

    std::vector<Point<T>> S;

    distance[start.grid_cell][start.grid_index] = 0;

    S.push_back(start);

    for(int i = 0; i < lambda; i++){
        pre_step(grid, S, distance, end, i);
        do_step(grid, S, distance, r, end, i);
    }
}


template <typename T>
void pre_step(std::vector<Point<T>> &S, std::vector<std::vector<int>>& distance, Point<T> &end, int step){
     std::vector<int> step_cells;

    //Save each cell that had a point in S_{i-1}, and add all other points from that cell to S_{i}
    for(int i = 0; i < S.size(); i++){
        Point<T> p = S[i];
        bool should_add_cell = false;
        for(int i = 0; i < grid.cells[p.grid_cell].points.size(); i++){
            if(distance[p.grid_cell][i] == -1){
                distance[p.grid_cell][i] = step+1;
                new_S.push_back(grid.cells[p.grid_cell].points[i]);
                should_add_cell = true;
            }else if(distance[p.grid_cell][i] == step-1 || distance[p.grid_cell][i] == step){ //second step conducted in this cell
                distance[p.grid_cell][i] = -step-2;
                should_add_cell = true;
            }
        }
        if(should_add_cell){
            step_cells.push_back(p.grid_cell);
        }

        //compute set Q of possible radii
        std::vector<T> Q;

        for(int i = 0; i < step_cells.size();i++){

            auto C = grid.cells[step_cells[i]];
            std::vector<Point<T>> C_points;
            std::vector<Point<T>> C_points_by_y;
            for(int k = 0; k < C.points.size(); k++){
                if(distance[C.points[k].grid_cell][C.points[k].grid_index] == -step-2){
                    C_points.push_back(C.points[k]);
                }
                if(distance[C.points_by_y[k].grid_cell][C.points_by_y[k].grid_index] == -step-2){
                    C_points_by_y.push_back(C.points_by_y[k]);
                }
            }

            //calculate the Voronoi diagram of the points in C
            std::vector<Point_2> v_points;
            for(int i = 0; i < C_points.size(); i++){
                v_points.push_back(Point_2(C_points[i].x, C_points[i].y));
            }

            VD vd;
            vd.insert(v_points.begin(), v_points.end());

            //for each vertex of the Voronoi diagram, calculate the distance to the closest point in C
            for(auto it = vd.vertices_begin(); it != vd.vertices_end(); it++){
                Point_2 p = it->point();

                auto dt_face = p->dual()->face();

                Point p1 = dt_face->vertex(0)->point();

                T dist = sqrt(CGAL::squared_distance(p, p1));
                
                Q.push_back(dist);
            }
        }

        //sort Q
        std::sort(Q.begin(), Q.end());

        //find the smallest radius that is feasible using binsearch
        auto it = std::lower_bound(Q.begin(), Q.end(), false, [this](T r){
            return this->cs.solve(this->input_points, this->sorted_by_x, this->sorted_by_y,r, this->s_i, this->t_i, this->lambda);
        });

        T res_l = -1;
        if(it != Q.end()){
            res_l = *it;
        }
        T res_r = -1;
        if(it+1 != Q.end()){
            res_r = *(it+1);
        }
        shrink(res_l, res_r);
    }
}

template <typename T>
bool alg_2<T>::do_step(std::vector<Point<T>> &S, std::vector<std::vector<int>>& distance, Point<T> &end, int step){
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

    S.clear();
    S = new_S;

    return false;
}

template <typename T>
void alg_2<T>::prepare_distance(std::vector<std::vector<int>>& distance){
    distance.resize(this->grid.cells.size());
    for(int i = 0; i < this->grid.cells.size(); i++){
        distance[i].resize(this->grid.cells[i].points.size());
        for(int j = 0; j < this->grid.cells[i].points.size(); j++){
            distance[i][j] = -1;
        }
    }

    for(int i=0; i<this->grid.cells.size(); i++){
        this->grid.cells[i].visited = false;
    }
}


















// BUILDIN THE GRID 

template <typename T>
void alg_2<T>::build_grid(int s_i, int t_i, int lambda){
    int s_pos = 0;

    Point<T> s = input_points[s_i];

    //find points to the left, and to the right of s
    for(int i = 0; i < input_points.size(); i++){
        if(input_points[sorted_by_x[i]].x == s.x && input_points[sorted_by_x[i]].y == s.y){
            s_pos = i;
            break;
        }
    }
    
    int m = input_points.size() - s_pos - 1;
    

    T r_res = matrix_search(m,s_pos, s_i, t_i, lambda, false, true, sorted_by_x);

    T l_res = matrix_search(m,s_pos, s_i, t_i, lambda, false, false, sorted_by_x);
    
    shrink(l_res, r_res);

    std::vector<int> left_by_x;
    for(int i = s_pos; i >= 0; i--){
        left_by_x.push_back(sorted_by_x[i]);
    }

    r_res = matrix_search(left_by_x.size(),0, s_i, t_i, lambda, false, true, left_by_x);

    l_res = matrix_search(left_by_x.size(),0, s_i, t_i, lambda, false, false, left_by_x);

    shrink(l_res, r_res);
    
    std::cout<<"RESULT after x :"<<interval_l<<" "<<interval_r<<std::endl;

    for(int i = 0; i < input_points.size(); i++){
        if(input_points[sorted_by_y[i]].x == s.x && input_points[sorted_by_y[i]].y == s.y){
            s_pos = i;
            break;
        }
    }
    m = input_points.size() - s_pos - 1;

    r_res = matrix_search(m,s_pos, s_i, t_i, lambda, true, true, sorted_by_y);

    l_res = matrix_search(m,s_pos, s_i, t_i, lambda, true, false, sorted_by_y);
    
    shrink(l_res, r_res);

    std::vector<int> left_by_y = {};
    for(int i = s_pos; i >= 0; i--){
        left_by_y.push_back(sorted_by_y[i]);
    }

    r_res = matrix_search(left_by_y.size(),0, s_i, t_i, lambda, true, true, left_by_y);

    l_res = matrix_search(left_by_y.size(),0, s_i, t_i, lambda, true, false, left_by_y);

    shrink(l_res, r_res);
    
    std::cout<<"RESULT after y :"<<interval_l<<" "<<interval_r<<std::endl;

    std::cout<<"building the grid for r = "<<interval_l<<std::endl;
    grid = Grid<T>(sorted_by_x, sorted_by_y, input_points, s, interval_l);
}

// MATRIX SEARCH ALGORITHM O(n log n)

template <typename T> 
class s_cell{
    public:
        int bot_x, bot_y;
        int top_x, top_y;
};

template <typename T> 
T alg_2<T>::get_value(int i, int j, int n, int s_pos, int size, int x_or_y, std::vector<int>& sorted){
    int i_start = size - n;
    if( i < i_start || j >= 2*n)
        return 0;
    T diff = this->input_points[sorted[s_pos + i - i_start]].get(x_or_y) - this->input_points[sorted[s_pos]].get(x_or_y);
    if(diff < 0)
        diff = -diff;
    return sqrt(2) * (diff / (j + 1));
}

template <typename T> 
T alg_2<T>::matrix_search(int n, int s_pos, int s_i, int t_i, int lambda, bool is_vertical, bool look_for_one, std::vector<int>& sorted){
    std::vector<s_cell<T>> active_cells;
    int x_or_y = is_vertical ? 1 : 0;

    if(n <= 1){
        return -1;
    }


    //input matrix is of size n x 2n
    //padding so that the matrix is a power of 2
    int padding = 1;
    while(padding < 2*n){
        padding *= 2;
    }

    /*std::cout<<"pad: "<<padding<<" for m="<<n<<" spos="<<s_pos<<std::endl;

    //print matrix
    for(int j = 0; j < padding; j++){
        for(int i=0;i<padding;i++){
            std::cout<<get_value(i,j,n,s_pos, padding, x_or_y, sorted)<<" ";
        }
        std::cout<<std::endl;
    }*/
    

    active_cells.push_back(s_cell<T>{0,padding-1,padding-1,0});


    while(true){
        std::vector<std::pair<T, int>> bot_values;
        std::vector<std::pair<T, int>> top_values;

        //get values of the top and bottom values of the active cells
        for(int i = 0; i < active_cells.size(); i++){
            bot_values.push_back(std::make_pair(get_value(active_cells[i].bot_x, active_cells[i].bot_y, n, s_pos, padding, x_or_y, sorted), i));
            top_values.push_back(std::make_pair(get_value(active_cells[i].top_x, active_cells[i].top_y, n, s_pos, padding, x_or_y, sorted), i));
        }


        //get the lower median of the bottom values
        nth_element(bot_values.begin(), bot_values.begin() + (int) ceil(((double) (bot_values.size()+1))/2) - 1, bot_values.end(), [](std::pair<T, int> a, std::pair<T, int> b){
            return a.first < b.first;
        });
        std::pair<T, int> x_s = bot_values[ceil(((double) (bot_values.size()+1))/2) - 1];

        //get the upper median of the top values
        nth_element(top_values.begin(), top_values.begin() + (top_values.size()+1)/2 - 1, top_values.end(), [](std::pair<T, int> a, std::pair<T, int> b){
            return a.first < b.first;
        });
        std::pair<T, int> x_l = top_values[(top_values.size()+1)/2 - 1];

        //std::cout<<"x_s: "<<x_s.first<<" x_l: "<<x_l.first<<std::endl;

        bool x_s_feasible = my_cs.solve(this->input_points, this->sorted_by_x, this->sorted_by_y, x_s.first, s_i, t_i, lambda);
        
        //get rid of the unimportant cells in regard to x_s
        std::vector<s_cell<T>> new_cells;
        bool one_equal = false;
        for(int i = 0; i < active_cells.size(); i++){
            T s_value = get_value(active_cells[i].bot_x, active_cells[i].bot_y, n, s_pos, padding, x_or_y, sorted);
            T l_value = get_value(active_cells[i].top_x, active_cells[i].top_y, n, s_pos, padding, x_or_y, sorted);
            if(look_for_one){
                if(x_s_feasible){
                    if(s_value < x_s.first)
                        new_cells.push_back(s_cell<T>{active_cells[i].bot_x, active_cells[i].bot_y, active_cells[i].top_x, active_cells[i].top_y});
                    else if(s_value == x_s.first && !one_equal){
                        one_equal = true;
                        new_cells.push_back(s_cell<T>{active_cells[i].bot_x, active_cells[i].bot_y, active_cells[i].top_x, active_cells[i].top_y});
                    }
                }else{
                    if(l_value > x_s.first)
                        new_cells.push_back(s_cell<T>{active_cells[i].bot_x, active_cells[i].bot_y, active_cells[i].top_x, active_cells[i].top_y});
                }
            }else{
                if(x_s_feasible){
                    if(s_value < x_s.first)
                        new_cells.push_back(s_cell<T>{active_cells[i].bot_x, active_cells[i].bot_y, active_cells[i].top_x, active_cells[i].top_y});
                }else{
                    if(l_value > x_s.first){
                        new_cells.push_back(s_cell<T>{active_cells[i].bot_x, active_cells[i].bot_y, active_cells[i].top_x, active_cells[i].top_y});
                    }else if(l_value == x_s.first && !one_equal){
                        one_equal = true;
                        new_cells.push_back(s_cell<T>{active_cells[i].bot_x, active_cells[i].bot_y, active_cells[i].top_x, active_cells[i].top_y});
                    }
                }
            }
        }
        
        /*std::cout<<"cells after x_s: "<<x_s.first<<" feas?: "<<x_s_feasible<<'\n';
        for(int i = 0; i < new_cells.size(); i++){
            std::cout<<"("<<new_cells[i].bot_x<<","<<new_cells[i].bot_y<<")-("<<new_cells[i].top_x<<","<<new_cells[i].top_y<<") ";
        }
        std::cout<<std::endl;*/

        bool x_l_feasible = my_cs.solve(this->input_points, this->sorted_by_x, this->sorted_by_y, x_l.first, s_i, t_i, lambda);

        //get rid of the unimportant cells in regard to x_l
        std::vector<s_cell<T>> new_new_cells;
        one_equal = false;
        for(int i = 0; i < new_cells.size(); i++){
            T s_value = get_value(new_cells[i].bot_x, new_cells[i].bot_y, n, s_pos, padding, x_or_y, sorted);
            T l_value = get_value(new_cells[i].top_x, new_cells[i].top_y, n, s_pos, padding, x_or_y, sorted);
            if(look_for_one){
                if(x_l_feasible){
                    if(s_value < x_l.first)
                        new_new_cells.push_back(s_cell<T>{new_cells[i].bot_x, new_cells[i].bot_y, new_cells[i].top_x, new_cells[i].top_y});
                    else if(s_value == x_l.first && !one_equal){
                        one_equal = true;
                        new_new_cells.push_back(s_cell<T>{new_cells[i].bot_x, new_cells[i].bot_y, new_cells[i].top_x, new_cells[i].top_y});
                    }
                }else{
                    if(l_value > x_l.first)
                        new_new_cells.push_back(s_cell<T>{new_cells[i].bot_x, new_cells[i].bot_y, new_cells[i].top_x, new_cells[i].top_y});
                }
            }else{
                if(x_l_feasible){
                    if(s_value < x_l.first)
                        new_new_cells.push_back(s_cell<T>{new_cells[i].bot_x, new_cells[i].bot_y, new_cells[i].top_x, new_cells[i].top_y});
                }else{
                    if(l_value > x_l.first)
                        new_new_cells.push_back(s_cell<T>{new_cells[i].bot_x, new_cells[i].bot_y, new_cells[i].top_x, new_cells[i].top_y});
                    else if(l_value == x_l.first && !one_equal){
                        one_equal = true;
                        new_new_cells.push_back(s_cell<T>{new_cells[i].bot_x, new_cells[i].bot_y, new_cells[i].top_x, new_cells[i].top_y});
                    }
                }
            }
        }

        /*std::cout<<"cells after x_l: "<<x_l.first<<" feas?: "<<x_l_feasible<<'\n';
        for(int i = 0; i < new_new_cells.size(); i++){
            std::cout<<"("<<new_new_cells[i].bot_x<<","<<new_new_cells[i].bot_y<<")-("<<new_new_cells[i].top_x<<","<<new_new_cells[i].top_y<<") ";
        }
        std::cout<<std::endl;*/

        //split every active cells into four subcells
        std::vector<s_cell<T>> new_active_cells;

        if(new_new_cells.size() == 1 && new_new_cells[0].bot_x == new_new_cells[0].top_x && new_new_cells[0].bot_y == new_new_cells[0].top_y){
            return get_value(new_new_cells[0].bot_x, new_new_cells[0].bot_y, n, s_pos, padding, x_or_y, sorted);
        }

        for(int i = 0; i < new_new_cells.size(); i++){
            s_cell<T> cell = new_new_cells[i];

            //is of size 1
            if(cell.bot_x == cell.top_x && cell.bot_y == cell.top_y){
                new_active_cells.push_back(cell);
                continue;
            }

            int size = cell.top_x - cell.bot_x + 1;

            //split into 4 subcells
            new_active_cells.push_back(s_cell<T>{cell.bot_x, cell.bot_y, cell.bot_x + size/2 - 1, cell.bot_y - size/2 + 1});
            new_active_cells.push_back(s_cell<T>{cell.bot_x, cell.bot_y - size/2, cell.bot_x + size/2 - 1, cell.top_y});
            new_active_cells.push_back(s_cell<T>{cell.bot_x + size/2, cell.bot_y, cell.top_x, cell.bot_y - size/2 + 1});
            new_active_cells.push_back(s_cell<T>{cell.bot_x + size/2, cell.bot_y - size/2, cell.top_x, cell.top_y});
        }

        active_cells = new_active_cells;

        /*std::cout<<"new active cells:";
        for(int i = 0; i < active_cells.size(); i++){
            std::cout<<"("<<active_cells[i].bot_x<<","<<active_cells[i].bot_y<<","<<get_value(active_cells[i].bot_x, active_cells[i].bot_y, n, s_pos, padding, x_or_y, sorted)<<")-("<<active_cells[i].top_x<<","<<active_cells[i].top_y<<","<<get_value(active_cells[i].top_x, active_cells[i].top_y, n, s_pos, padding, x_or_y, sorted)<<") ";
        }
        std::cout<<std::endl;*/
    }
}