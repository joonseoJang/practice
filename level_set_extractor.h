#ifndef LEVEL_SET_EXTRACTOR_H_
#define LEVEL_SET_EXTRACTOR_H_
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>

enum class visitState{
    FAR,
    NARROW_BAND,
    FROZEN
};

class gridPoint2d{
    public:
        // default constructor
        gridPoint2d();
        
        // constructor 
        gridPoint2d(double x, double y, double t);

        // set functions
        void setX(double x);
        void setY(double y);
        bool setT(double t);
        bool setF(double f);
        void setState(visitState state);
        void setPoint(double x, double y);
        
        //call functions
        double x();
        double y();
        double t();
        double f();
        visitState state();

        bool operator < (gridPoint2d& p){
            return this->t() < p.t();
        }

        bool operator > (gridPoint2d& p){
            return this->t() > p.t();
        }
    
    private:
        visitState state_;
        double x_;
        double y_;
        double t_;
        double f_;
};

class minHeap{
    public:
        minHeap();

        ~minHeap();

        void pushData(gridPoint2d* pt);

        void deleteData(double i);

        void sortData(int i);

        gridPoint2d* at(int i);

        std::vector <gridPoint2d*>& heap();
    
    private:
        std::vector <gridPoint2d*> binary_tree_;
        
};

class levelSetExtractor{
    public:
        levelSetExtractor();

        ~levelSetExtractor();

        void initialize(std::vector<std::vector<double>> polygon,
                        double step_x, 
                        double step_y, 
                        double alpha, 
                        double beta);

        bool extractCenterLine(std::vector<std::vector<double>>& center_line);

        void neighbor2band();
        
        void recomputeBand();

        bool setModerateSpeed(gridPoint2d& pt);

        bool setFastSpeed(gridPoint2d& pt);

        bool getGridMapPoint(gridPoint2d*& pt, int i, int j);
        
        bool getGridMapPoint(gridPoint2d*& pt, int i, int j, visitState state);
        
        int quadEqSolver(std::vector<double> coeff, std::vector<double>& sol);

        // get gradient value of gridPoint2d(pt.x + dx, pt.y + dy)
        bool getGradientT(gridPoint2d& pt, std::vector<double>& grad, double dx, double dy);

    private: 
        //all points
        std::vector<std::vector<double>> polygon_;
        std::vector<gridPoint2d> points_;
        std::vector<std::vector<gridPoint2d*>> gridMap_;
        int origin_x_index_, origin_y_index_, end_x_index_, end_y_index_;
        double step_x_, step_y_, alpha_, beta_;

        int neighbor_vector_[4][2] = 
        {
            {-1, 0},
            {1, 0},
            {0, -1},
            {0, 1}
        };
        
        //narrow band points
        minHeap narrow_band_;
        
        // current points
        gridPoint2d* current_point_;
        int current_point_index_;

        bool initialized;
        
        
};
#endif