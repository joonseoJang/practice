#include "level_set_extractor.h"
#include <float.h>
#include <stdexcept>

gridPoint2d::gridPoint2d(){
    x_ = y_ = 0;
    t_ = DBL_MAX;
    state_ = visitState::FAR;
}

gridPoint2d::gridPoint2d(double x, double y, double t=DBL_MAX){
    x_ = x;
    y_ = y;
    t_ = DBL_MAX;
    state_ = visitState::FAR;
}

void gridPoint2d::setX(double x){
    x_ = x;
}

void gridPoint2d::setY(double y){
    y_ = y;
}

bool gridPoint2d::setT(double t){
    if(state_ == visitState::FROZEN){
        return false;
    }
    t_ = t;
    return true;
}

bool gridPoint2d::setF(double f){
    if(f<0) return false;
    else{
        f_ = f;
        return true;
    }
}

void gridPoint2d::setState(visitState state){
    state_ = state;
}

void gridPoint2d::setPoint(double x, double y){
    x_ = x;
    y_ = y;
}

double gridPoint2d::x(){
    return x_;
}

double gridPoint2d::y(){
    return y_;
}

double gridPoint2d::t(){
    return t_;
}

double gridPoint2d::f(){
    return f_;
}

visitState gridPoint2d::state(){
    return state_;
}

minHeap::minHeap(){
}

minHeap::~minHeap(){
    binary_tree_.clear();
}

void minHeap::pushData(gridPoint2d* pt){
    binary_tree_.push_back(pt);
    int new_index = binary_tree_.size()-1;
    int parent_index = new_index/2;
    while(new_index>1 && *binary_tree_.at(new_index) < *binary_tree_.at(parent_index)){
        std::iter_swap(binary_tree_.begin() + new_index, binary_tree_.begin() + parent_index);
        new_index = parent_index;
        parent_index = new_index/2;
    }
}

void minHeap::deleteData(double i){
    if(i == 0) return;
    std::iter_swap(binary_tree_.begin() + i, binary_tree_.end() - 1);
    binary_tree_.resize(binary_tree_.size() - 1);
    int parent = i;
    int child = parent*2;
    if(child + 1 < binary_tree_.size()){
        child = (*binary_tree_.at(child) < *binary_tree_.at(child+1))? child : child + 1; 
    }

    while(child< binary_tree_.size() && *binary_tree_.at(child)<*binary_tree_.at(parent)){
        std::iter_swap(binary_tree_.begin() + parent, binary_tree_.begin() + child);
        parent = child;
        child = parent*2;
        if(child + 1 < binary_tree_.size()){
            child = (*binary_tree_.at(child) < *binary_tree_.at(child+1))? child : child + 1; 
        }
    }
}

void minHeap::sortData(int i=0){
    if(i>=binary_tree_.size()){
        return;
    }
    else if(i!=0){
        int parent = i;
        int child = parent*2;
        if(child + 1 < binary_tree_.size()){
            child = (*binary_tree_.at(child) < *binary_tree_.at(child+1))? child : child + 1; 
        }

        while(child< binary_tree_.size() && *binary_tree_.at(child)<*binary_tree_.at(parent)){
            std::iter_swap(binary_tree_.begin() + parent, binary_tree_.begin() + child);
            parent = child;
            child = parent*2;
            if(child + 1 < binary_tree_.size()){
                child = (*binary_tree_.at(child) < *binary_tree_.at(child+1))? child : child + 1; 
            }
        } 
        return;
    }
    else{
        for(i = 1; i<binary_tree_.size(); i++){
            int parent = i;
            int child = parent*2;
            if(child + 1 < binary_tree_.size()){
                child = (*binary_tree_.at(child) < *binary_tree_.at(child+1))? child : child + 1; 
            }

            while(child< binary_tree_.size() && *binary_tree_.at(child)<*binary_tree_.at(parent)){
                std::iter_swap(binary_tree_.begin() + parent, binary_tree_.begin() + child);
                parent = child;
                child = parent*2;
                if(child + 1 < binary_tree_.size()){
                    child = (*binary_tree_.at(child) < *binary_tree_.at(child+1))? child : child + 1; 
                }
            } 
        }
        return;
    }
}

gridPoint2d* minHeap::at(int i){
    return binary_tree_.at(i);
}

std::vector <gridPoint2d*>& minHeap::heap(){
    return binary_tree_;
}

levelSetExtractor::levelSetExtractor(){
    initialized = false;
}

levelSetExtractor::~levelSetExtractor(){
    gridMap_.clear();
    narrow_band_.heap().clear();
}

void levelSetExtractor::initialize(std::vector<std::vector<double>> polygon, 
                                          double step_x, 
                                          double step_y, 
                                          double alpha, 
                                          double beta=0){
    if(polygon.front().at(0)!=polygon.back().at(0) || polygon.front().at(1)!=polygon.back().at(1)){
        polygon.push_back(polygon.front());
    }

    step_x_ = step_x;
    step_y_ = step_y;
    alpha_ = alpha;
    beta_ = beta;
    int point_num = 0;
    double min_x = DBL_MAX, min_y = DBL_MAX;
    double max_x = -DBL_MAX, max_y = -DBL_MAX;
    for(std::vector<std::vector<double>>::iterator it = polygon.begin(); it!= polygon.end(); it++){
        if(it->at(0)>max_x) max_x = it->at(0);
        if(it->at(0)<min_x) min_x = it->at(0);
        if(it->at(0)>max_y) max_y = it->at(1);
        if(it->at(0)<min_y) min_y = it->at(1);
        polygon_.push_back({it->at(0), it->at(1)});
    }
    origin_x_index_ = ceil(min_x / step_x);
    origin_y_index_ = ceil(min_y / step_y);
    end_x_index_ = floor(max_x / step_x);
    end_y_index_ = floor(max_y / step_y);
    gridMap_.resize(end_x_index_ - origin_x_index_);
    for(std::vector<gridPoint2d*> v : gridMap_){
        v.resize(end_y_index_ - origin_y_index_);
    }
    
    for(int j = origin_y_index_; j <= end_y_index_; j++){
        double height = j*step_y;
        std::vector<double> x_points;
        for(std::vector<std::vector<double>>::iterator it = polygon.begin(); it!= polygon.end()-1; it++){
            if((it->at(1)-height)*((it+1)->at(1)-height)<0){
                double temp_x = it->at(0) + ((it+1)->at(0)-it->at(0))*(height - it->at(1))/((it+1)->at(1)-it->at(1));
                x_points.push_back(temp_x);
            }
        }
        if(x_points.size()%2==1) x_points.erase(x_points.end()-1);
        std::sort(x_points.begin(), x_points.end());
        for(std::vector<double>::iterator it = x_points.begin(); it!=x_points.end(); it+=2){
            for(int i = ceil(*it/step_x); i<=floor(*(it+1)/step_x); i++){
                double width = i*step_x;
                gridPoint2d pt(width, height);
                points_.push_back(pt);
                gridMap_.at(i-origin_x_index_).at(j-origin_y_index_) = &points_.back();
            }
        }
    }
    initialized = true;
}

bool levelSetExtractor::extractCenterLine(std::vector<std::vector<double>>& center_line){
    if(!initialized){
        std::cout<<"levelSetExtractor::extractCenterLine: you have to initialize extractor first"<<std::endl;
    }
    //TODO: starting point index should not be hardcoded
    current_point_ = &points_.at(0);
    current_point_->setT(0);
    current_point_->setState(visitState::FROZEN);
    setFastSpeed(*current_point_);
    neighbor2band();
    
    while(narrow_band_.heap().size()>1){
        neighbor2band();
        recomputeBand();
        current_point_ = narrow_band_.at(1);
        current_point_->setState(visitState::FROZEN);
        narrow_band_.deleteData(1);
    }
    center_line.clear();
    double dx = 0;
    double dy = 0;
    while(current_point_->t()>0){
        center_line.push_back({current_point_->x()+dx, current_point_->y()+dy});
        std::vector<double> grad;
        if(!getGradientT(*current_point_, grad, dx, dy)){
            std::cout<<"levelSetExtractor::extractCenterLine: Cannot get gradient of runge-kutta start point ("<<current_point_->x()<<", "<<current_point_->y()<<")"<<std::endl;
            return false;
        }
        double k_x = -step_x_*grad.at(0)/hypot(grad.at(0), grad.at(1));
        double k_y = -step_y_*grad.at(1)/hypot(grad.at(0), grad.at(1));
        if(!getGradientT(*current_point_, grad, k_x/2, k_y/2)){
            std::cout<<"levelSetExtractor::extractCenterLine: Cannot get gradient of runge-kutta middle  point ("<<current_point_->x() + k_x/2<<", "<<current_point_->y() + k_y/2<<")"<<std::endl;
            return false;
        }
        double next_x = current_point_->x() - step_x_*grad.at(0)/hypot(grad.at(0), grad.at(1));
        double next_y = current_point_->y() - step_y_*grad.at(1)/hypot(grad.at(0), grad.at(1));
        if(!getGridMapPoint(current_point_, round(next_x / step_x_), round(next_y / step_y_), visitState::FROZEN)){
            std::cout<<"levelSetExtractor::extractCenterLine: Cannot get gradient of runge-kutta end point ("<<next_x<<", "<<next_y<<")"<<std::endl;
            return false;
        }
        dx = next_x - current_point_->x();
        dy = next_y - current_point_->y();
    }
    return true;
}

void levelSetExtractor::neighbor2band(){
    int i = current_point_->x()/step_x_ - origin_x_index_;
    int j = current_point_->y()/step_y_ - origin_y_index_;
    for(int* v : neighbor_vector_){
        int n_i = i + *v;
        int n_j = j + *(v+1);
        gridPoint2d* point_ptr;
        if(getGridMapPoint(point_ptr, n_i, n_j, visitState::FAR)){
                setFastSpeed(*point_ptr);
                narrow_band_.pushData(point_ptr);
        }
    }
}

void levelSetExtractor::recomputeBand(){
    for(int pt_index = 1; pt_index<narrow_band_.heap().size(); pt_index++){
        gridPoint2d* pt = narrow_band_.at(pt_index);
        std::vector<double> coeff = {-pow(pt->f(), -2), 0, 0};
        int i = pt->x()/step_x_ - origin_x_index_;
        int j = pt->y()/step_y_ - origin_y_index_;
        for(int axis=0; axis<2; axis++){
            double val1 = DBL_MAX;
            double val2 = DBL_MAX;
            for(int dir=0; dir<2; dir++){
                int n_i = i + neighbor_vector_[2*axis + dir][0];
                int n_j = j + neighbor_vector_[2*axis + dir][1];
                gridPoint2d* temp_pt;
                if(getGridMapPoint(temp_pt, n_i, n_j, visitState::FROZEN)){
                    if(temp_pt->t()<val1){
                        val1 = temp_pt->t();
                        int n_2i = i + 2*neighbor_vector_[2*axis + dir][0];
                        int n_2j = j + 2*neighbor_vector_[2*axis + dir][1];
                        if(getGridMapPoint(temp_pt, n_2i, n_2j, visitState::FROZEN)){
                            if(temp_pt->t()<val1){
                                val2 = temp_pt->t();
                            }
                            else{
                                val2 = DBL_MAX;
                            }
                        }
                    }
                }
            }
            if(val2 != DBL_MAX){
                double a = (axis == 0)? 3/(2*step_x_) : 3/(2*step_y_);
                double K = (4*val1 - val2)/3;
                coeff.at(2) += a*a;
                coeff.at(1) -= 2*a*K;
                coeff.at(0) += pow(a*K,2);
            }
            else{
                double a = (axis == 0)? 1/step_x_ : 1/(step_y_);
                coeff.at(2) += a*a;
                coeff.at(1) -= 2*a*val1;
                coeff.at(0) += pow(a*val1, 2);
            }
        }
        std::vector<double> sol;
        if(int no_sol = quadEqSolver(coeff, sol)){
            if(no_sol == 2){
                pt->setT(std::min(pt->t(), std::max(sol.at(0), sol.at(1)) ) );
                narrow_band_.sortData(pt_index);
            }
            else{
                pt->setT(std::min(pt->t(), sol.at(0)));
                narrow_band_.sortData(pt_index);
            }
            continue;
        }
        else{
            continue;
        }
    }
}

bool levelSetExtractor::setModerateSpeed(gridPoint2d& pt){
    double min_val = DBL_MAX;
    for(std::vector<std::vector<double>>::iterator it = polygon_.begin(); it != polygon_.end()-1; it++){
        double line_x = (it+1)->at(0) - it->at(0);
        double line_y = (it+1)->at(1) - it->at(1);
        double point_x = pt.x() - it->at(0);
        double point_y = pt.y() - it->at(1);
        double dot = line_x*point_x + line_y*point_y;
        if(dot<0 || dot> line_x*line_x + line_y*line_y) continue;
        double dist = abs((line_x*point_y - line_y*point_x)/hypot(line_x, line_y));
        if(dist<min_val) min_val = dist;
    }
    return pt.setF(pow(min_val, beta_));
}

bool levelSetExtractor::setFastSpeed(gridPoint2d& pt){
    double min_val = DBL_MAX;
    for(std::vector<std::vector<double>>::iterator it = polygon_.begin(); it != polygon_.end()-1; it++){
        double line_x = (it+1)->at(0) - it->at(0);
        double line_y = (it+1)->at(1) - it->at(1);
        double point_x = pt.x() - it->at(0);
        double point_y = pt.y() - it->at(1);
        double dot = line_x*point_x + line_y*point_y;
        if(dot<0 || dot> line_x*line_x + line_y*line_y) continue;
        double dist = abs((line_x*point_y - line_y*point_x)/hypot(line_x, line_y));
        if(dist<min_val) min_val = dist;
    }
    return pt.setF(exp(alpha_*min_val));
}

bool levelSetExtractor::getGridMapPoint(gridPoint2d*& pt, int i, int j){
    try{
        if(gridMap_.at(i).at(j) == nullptr){
             pt = nullptr;
             return false;
        }
        else{
            pt = gridMap_.at(i).at(j);
            return true;
        }
    }
    catch(const std::out_of_range& oor){
        pt = nullptr;
        return false;
    }
}

bool levelSetExtractor::getGridMapPoint(gridPoint2d*& pt, int i, int j, visitState state){
    try{
        if(gridMap_.at(i).at(j) == nullptr){
             pt = nullptr;
             return false;
        }
        else if(gridMap_.at(i).at(j)->state() != state) return false;
        else{
            pt = gridMap_.at(i).at(j);
            return true;
        }
    }
    catch(const std::out_of_range& oor){
        pt = nullptr;
        return false;
    }
}

int levelSetExtractor::quadEqSolver(std::vector<double> coeff, std::vector<double>& sol){
    if(coeff.size()<3) return 0;
    sol.resize(2);
    if(coeff.at(2) == 0){
        sol.at(0) = -coeff.at(0)/coeff.at(1);
        return 1;
    }
    double checker = coeff.at(1)*coeff.at(1) - 4*coeff.at(2)*coeff.at(0);
    if(checker>0){
        sol.at(0) = (-coeff.at(1)+sqrt(checker))/(2*coeff.at(2));
        sol.at(1) = (-coeff.at(1)-sqrt(checker))/(2*coeff.at(2));
        return 1;
    }
    else if(checker == 0){
        sol.at(0) = -coeff.at(1)/(2*coeff.at(2));
    }
    else{
        return 0;
    }
}

bool levelSetExtractor::getGradientT(gridPoint2d& pt, std::vector<double>& grad, double dx = 0, double dy = 0){
    grad.resize(2);
    grad.at(0) = grad.at(1) = 0;
    int i = pt.x()/step_x_ - origin_x_index_;
    int j = pt.y()/step_y_ - origin_y_index_;
    {
    //x
    //1. get near points
    std::vector<gridPoint2d*> x_ptrs(5);
    int k = 0;
    for(std::vector<gridPoint2d*>::iterator it = x_ptrs.begin(); it!= x_ptrs.end();){
        if(k==2){
            *it = &pt;
            k++;
            continue;
        }
        else if(!getGridMapPoint(*it, i+k-2, j, visitState::FROZEN)){
            it = x_ptrs.erase(it);
            k++;
        }
        else{
            it++;
            k++;
        }
    }
    if(x_ptrs.size()<2) return false;
    //2. get \delta^nf_i 
    k = 1;
    std::vector<double> delta_f;
    std::vector<double> temp_delta;
    for(int l=0; l<x_ptrs.size()-k; l++){
        temp_delta.push_back((x_ptrs.at(l+k)->t() - x_ptrs.at(l)->t())/(x_ptrs.at(l+k)->x() - x_ptrs.at(l)->x()));
    }
    delta_f.push_back(temp_delta.at(0));
    k++;
    while(temp_delta.size()>1){
        for(int l = 0; l < x_ptrs.size()-k; l++){
            temp_delta.at(l) = (temp_delta.at(l+1) - temp_delta.at(l))/(x_ptrs.at(l+k)->x() - x_ptrs.at(l)->x());
        }
        delta_f.push_back(temp_delta.at(0));
        temp_delta.resize(temp_delta.size()-1);
        k++;
    }
    //3. f_n(x) = \delta_0 + \delta_1 (x-x_i) + \delta_2 (x-x_i)(x-x_{i+1}) ...
    //T_x = f'(x)
    //TODO when dx = 0 -> divide by zero
    for(int l=0; l<delta_f.size(); l++){
        double x_ = 1;
        double x_diff = 0;
        bool div_by_zero = false;
        for(int m=0; m<l+1; m++){
            if(pt.x() + dx - x_ptrs.at(m)->x() == 0){
                div_by_zero = true;
                continue;
            }
            else x_ *= pt.x() + dx - x_ptrs.at(m)->x();
        }
        if(div_by_zero){
            x_diff = x_;
        }
        else{
            for(int m=0; m<l+1; m++){
                x_diff += x_/(pt.x() + dx - x_ptrs.at(m)->x());
            }
        }
        grad.at(0) += delta_f.at(l) * x_diff;
    }
    } 
    
    {
    //y
    //1. get near points
    std::vector<gridPoint2d*> y_ptrs(5);
    int k = 0;
    for(std::vector<gridPoint2d*>::iterator it = y_ptrs.begin(); it!= y_ptrs.end();){
        if(k==2){
            *it = &pt;
            k++;
            continue;
        }
        else if(!getGridMapPoint(*it, i, j+k-2, visitState::FROZEN)){
            it = y_ptrs.erase(it);
            k++;
        }
        else{
            it++;
            k++;
        }
    }
    if(y_ptrs.size()<2) return false;
    //2. get \delta^nf_i 
    k = 1;
    std::vector<double> delta_f;
    std::vector<double> temp_delta;
    for(int l=0; l<y_ptrs.size()-k; l++){
        temp_delta.push_back((y_ptrs.at(l+k)->t() - y_ptrs.at(l)->t())/(y_ptrs.at(l+k)->y() - y_ptrs.at(l)->y()));
    }
    delta_f.push_back(temp_delta.at(0));
    k++;
    while(temp_delta.size()>1){
        for(int l = 0; l < y_ptrs.size()-k; l++){
            temp_delta.at(l) = (temp_delta.at(l+1) - temp_delta.at(l))/(y_ptrs.at(l+k)->y() - y_ptrs.at(l)->y());
        }
        delta_f.push_back(temp_delta.at(0));
        temp_delta.resize(temp_delta.size()-1);
        k++;
    }
    //3. f_n(y) = \delta_0 + \delta_1 (y-y_i) + \delta_2 (y-y_i)(y-y_{i+1}) ...
    //T_y = f'(y)
    //TODO when dy = 0 -> divide by zero
    for(int l=0; l<delta_f.size(); l++){
        double y_ = 1;
        double y_diff = 0;
        bool div_by_zero = false;
        for(int m=0; m<l+1; m++){
            if(pt.y() + dy - y_ptrs.at(m)->y() == 0){
                div_by_zero = true;
                continue;
            }
            else y_ *= pt.y() + dy - y_ptrs.at(m)->y();
        }
        if(div_by_zero){
            y_diff = y_;
        }
        else{
            for(int m=0; m<l+1; m++){
                y_diff += y_/(pt.y() + dy - y_ptrs.at(m)->y());
            }
        }
        grad.at(1) += delta_f.at(l) * y_diff;
    }
    }
    return true;
}

/*
int main(){
    std::cout<<"hello world";
    
    std::vector<std::vector<double>> polygon = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    levelSetExtractor lse;
    lse.initialize(polygon, 0.1, 0.1, 15, 0.5);
    
}
*/