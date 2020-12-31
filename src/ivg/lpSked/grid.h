#ifndef GRID_H
#define GRID_H
#include <vector>
#include <functional>
#include <iostream>
#include "geometry.h"
#include "ivg_const.h"
#include "date.h"
#include "definitions.h"


namespace lps{

class Wedge{
public:

    // e elevation, a azimuth, both degree
    Wedge(Point center, double e1, double e2, double a1, double a2) : center_(center),e1_(e1), e2_(e2), a1_(a1), a2_(a2){
        surfaceArea_ = abs((a2_-a1_)*ivg::d2rad*(sin(e2_*ivg::d2rad)-sin(e1_*ivg::d2rad)));
    }

    Wedge() {

    }

    bool contains(const Point& p) const;

    Path toPath(double precision = 2) const;

    double e1() const {return e1_;}
    double e2() const {return e2_;}
    double em() const {
//        double d = abs(a1_ - a2_);
//        if(d > 360.0-1e-3 && d < 360+1e-3)
//            return 90.0;
        return (e1_+e2_)/2;};
    double am() const {return (a1_+a2_)/2;}
    double a1() const {return a1_;}
    double a2() const {return a2_;}
    Point  center() const {return center_;}
    
    double surfaceArea() const{ return surfaceArea_;};
    
    void print() const;
    
private:
    Point  center_;
    double e1_;
    double e2_;
    double a1_;
    double a2_;
    double surfaceArea_;
};

Wedge quadrant(int i, const Wedge& wedge);
Wedge cakePiece( int i, int np, const Wedge& wedge);

std::vector<Wedge> createEqualAreaPartition(  std::vector<unsigned>& k );

template<class Object>
class DataPoint{

public:
    DataPoint(Point location, Object object) : location_(location), object_(object){

    }

    const Point& location()const {return location_;}
    const Object& object() const {return object_;}

private:
    Point  location_;
    Object object_;
};

template<class Object, class Shape>
class Node{

public:
    Node(){}
    
    // recursive constructor
    Node(const std::vector<DataPoint<Object>> & points, const Shape & boundingRect, int level=0, int max_level = std::numeric_limits<int>::max()) :
        boundingRect_(boundingRect), level_(level){
        
        // add all points within boundingRect to point list
        for( const DataPoint<Object>& p : points ){
            if( boundingRect.contains(p.location()) ){
                points_.push_back(p);
            }
        }
        
       // call constructor recursive for all 4 quadrants
        if(  level < max_level){
            for(int i=0; i < 4; ++i){
                Shape shape = quadrant(i, boundingRect_);
                children_.push_back( Node(points_, shape, level+1, max_level) );
            }

        }
        
        includesPointLevel1_ = level > 0 && points_.size() > 0;
              
    }
    
    // non recursive constructor for equal area size approach
    Node(const std::vector<DataPoint<Object>> & points, const Shape & boundingRect, std::vector< std::vector<unsigned> > levelPartition) :
        boundingRect_(boundingRect){
        includesPointLevel1_ = false;
        level_ = 0;
        for( const DataPoint<Object>& p : points ){
            if( boundingRect.contains(p.location()) ){
                points_.push_back(p);
            }
        }
                        
        for( unsigned level = 0; level < levelPartition.size(); ++level ){
            std::vector<Wedge> cells = lps::createEqualAreaPartition(levelPartition[level]);
            for(Wedge& w : cells ){
                Shape shape = w;
                children_.push_back( Node(points_, shape, level+1, level+1) );
            }
        }
              
    }
    
        
    void visitBound(std::function<void (int level, const Shape& boundingRect)> visitor) const{
        visitor(level_, boundingRect_);
        for(const Node<Object,Shape> & child : children_){
            child.visitBound(visitor);
        }
    }
    
    void visitBoundAndBool(std::function<void (int level, const Shape& boundingRect, bool includesPointLevel1)> visitor) const{
        visitor(level_, boundingRect_, includesPointLevel1_);
        for(const Node<Object,Shape> & child : children_){
            child.visitBoundAndBool(visitor);
        }
    }

    void visitPoints(std::function<void (const std::vector<DataPoint<Object>>&,int level)> visitor)const{
        visitor(points_, level_);
        for(const Node<Object,Shape> & child : children_){
            child.visitPoints(visitor);
        }
    }
    
    void visitPointsAndBool(std::function<void (const std::vector<DataPoint<Object>>&, bool includesPointLevel1, int level)> visitor)const{
        visitor(points_, includesPointLevel1_, level_);
        for(const Node<Object,Shape> & child : children_){
            child.visitPointsAndBool(visitor);
        }
    }
    
    void visitBool(std::function<void (int level, bool includesPointLevel1)> visitor)const{
        visitor(level_, includesPointLevel1_);
        for(const Node<Object,Shape> & child : children_){
            child.visitBool(visitor);
        }
    }
    
    void visitPointsAndBound(std::function<void (const std::vector<DataPoint<Object>>&, int level, const Shape& boundingRect)> visitor)const{
        visitor(points_,level_, boundingRect_);
        for(const Node<Object,Shape> & child : children_){
            child.visitPointsAndBound(visitor);
        }
    }
    
    void visitAll(std::function<void (const std::vector<DataPoint<Object>>&, int level, const Shape& boundingRect,  bool includesPointLevel1)> visitor)const{
        visitor(points_, level_, boundingRect_, includesPointLevel1_);
        for(const Node<Object,Shape> & child : children_){
            child.visitAll(visitor);
        }
    }

    const Shape & boundingRect() const {return boundingRect_;}
    int level() const {return level_;}
    std::vector<Node<Object,Shape>> children() const {return children_;}

private:
    Shape boundingRect_;
    std::vector<DataPoint<Object>> points_;
    std::vector<Node> children_;
    int level_;
    bool includesPointLevel1_; // true if level> 0 and at least one point is included

};


class TemporalGrid{
public:
    TemporalGrid(){};
    
    TemporalGrid( lps::Seconds sessionDuration, lps::Seconds resolution, lps::Seconds shift );
    
    TemporalGrid( lps::Seconds sessionDuration, lps::Seconds resolution, lps::Seconds shift , ivg::Date start);
    
    // Assume TemporalGrid and other Sampling have the same start epoch
    // Start is in the temporal grid. End not
    std::pair<unsigned, unsigned> getIndexOtherSampling( unsigned t, lps::Seconds intervalLength ) const;
    
    // returns the indices of the temporal cells including the point at time t
    std::vector<unsigned> getTemporalIndices( lps::Seconds t ) const;
    
    std::vector<unsigned> getTemporalIndices( ivg::Date d ) const;
    
    std::pair<lps::Seconds, lps::Seconds> getIntervalBorder( int i  ) const;
    
    unsigned get_number_of_intervals() const {return n_intervals; };
    
    std::vector<lps::Seconds> get_intervalBeginnings() const { return intervalBeginnings_; };
    
    std::vector<ivg::Date> get_intervalBeginningsDate() const;
    
    std::vector<lps::Seconds> get_intervalCenters() const;
    
    std::vector<ivg::Date> get_intervalCentersDate() const;
    
private:
    std::vector<lps::Seconds> intervalBeginnings_;
    lps::Seconds shift_;
    lps::Seconds resolution_;
    unsigned n_intervals;
    ivg::Date start_;
};

} //namespace


#endif // GRID_H
