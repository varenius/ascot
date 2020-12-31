#include "grid.h"
#include "matrix.h"


namespace lps{

bool Wedge::contains(const Point &p) const{
    
    Point v = p;
    
    while(v.x() < 0.0){
        v.rx() += 360.0;
    }
    while(v.x() > 360.0){
        v.rx() -= 360.0;
    } 
    
    if( v.x() == 0 )
        v.rx() += 1e-12;
                
    return v.x() <= a2_ && v.x() > a1_ && v.y() <= e2_ && v.y() > e1_;

}

Path Wedge::toPath(double precision)const{
    Path path;
  
    // check if wedge is full circle
    if(abs(a1_-a2_)<1e-3 || abs(abs(a1_-a2_)-360.0) <1e-3 ){
        bool first = true;
        for(double deg = a1_; deg <= a2_; deg+=precision){
            double rad = deg * M_PI/180;
            if(first){
                path.moveTo( (90-e1_)*std::sin(rad)+center_.x(), (90-e1_)*-std::cos(rad)+center_.y());
                first = false;
            }else{
                path.lineTo( (90-e1_)*std::sin(rad)+center_.x(), (90-e1_)*-std::cos(rad)+center_.y());
            }
        }
        path.closeSubpath();
        return path;
    } else {
    
        bool first = true;
        for(double deg = a1_; deg < a2_; deg+=precision){
            double rad = deg * M_PI/180;
            if(first){
                path.moveTo( (90-e1_)*std::sin(rad)+center_.x(), (90-e1_)*-std::cos(rad)+center_.y());
                first = false;
            }else{
                path.lineTo( (90-e1_)*std::sin(rad)+center_.x(), (90-e1_)*-std::cos(rad)+center_.y());
            }
        }
        path.lineTo( (90-e1_)*std::sin(a2_* M_PI/180)+center_.x(), (90-e1_)*-std::cos(a2_* M_PI/180)+center_.y());

        for(double deg = a2_; deg > a1_; deg-=precision){
            double rad = deg * M_PI/180;
            if(first){
                path.moveTo( (90-e2_)*std::sin(rad)+center_.x(), (90-e2_)*-std::cos(rad)+center_.y());
                first = false;
            }else{
                path.lineTo( (90-e2_)*std::sin(rad)+center_.x(), (90-e2_)*-std::cos(rad)+center_.y());
            }
        }
        path.lineTo( (90-e2_)*std::sin(a1_* M_PI/180)+center_.x(), (90-e2_)*-std::cos(a1_* M_PI/180)+center_.y());
        path.closeSubpath();
        return path;
    }
}

void Wedge::print() const{
    std::cout << "Elevation: " << this->e1_ << " - " << this->e2_ << " Azimuth: " << this->a1_ << " - " << this->a2_ << std::endl;
}

Wedge quadrant(int i, const Wedge &wedge){
    if(i ==0){
        return Wedge(wedge.center(),wedge.em(),wedge.e2(),wedge.a1(),wedge.am());
    }
    if(i ==1){
        return Wedge(wedge.center(),wedge.em(),wedge.e2(),wedge.am(),wedge.a2());
    }
    if(i ==2){
        return Wedge(wedge.center(),wedge.e1(),wedge.em(),wedge.am(),wedge.a2());
    }
    return Wedge(wedge.center(),wedge.e1(),wedge.em(),wedge.a1(),wedge.am());
}

Wedge cakePiece( int i, int np, const Wedge& wedge ){
    
    double az = 360.0/np;
    return Wedge(wedge.center(), wedge.e1(), wedge.e2(), i*az , (i+1)*az );
}

std::vector<Wedge> createEqualAreaPartition( std::vector<unsigned>& k ){
    
    unsigned n_rings = k.size();
    
    std::vector<double> cumk(n_rings);
    cumk[0] = k[0];
    for(unsigned i = 1; i < n_rings; ++i ){
        cumk[i] = cumk[i-1] + k[i];
    }
        
    ivg::Matrix r (n_rings, 1, 0.0);
    r(0) = 1.0;
    for(unsigned i = 1; i < n_rings; ++i ){
        r(i) = sqrt( (cumk[i]/cumk[i-1]) * std::pow(r(i-1), 2) );
    }
   
    r = r/(r(n_rings-1))*sqrt(2);
        
    std::vector<Wedge> result;
    
    for(unsigned i = 0; i < n_rings; ++i ){
        double epsilon2 = 90.0;
        if(i>0){
            epsilon2 = 90.0 - 2*asin(r(i-1)/2)*ivg::rad2d;
        }
        double epsilon1 = 90.0 - 2*asin(r(i)/2)*ivg::rad2d;
        double az = 360.0/k[i];
        for(unsigned j = 0; j < k[i]; ++j){
            double az1 = j*az;
            double az2 = (j+1)*az;
//            if(i%2==0){
//                az1 += az/2;
//                az2 += az/2;
//            }
                
            result.push_back(  Wedge(lps::Point(0,0), epsilon1, epsilon2, az1 , az2) );
        }
    }
    
    return result;
    
}


TemporalGrid::TemporalGrid( lps::Seconds sessionDuration, lps::Seconds resolution, lps::Seconds shift ){
    shift_ = shift;
    resolution_ = resolution;
    for( lps::Seconds t = 0; t + resolution <= sessionDuration + 1E-3; t+= shift){
        intervalBeginnings_.push_back(t);
    }
    n_intervals = intervalBeginnings_.size();
    
}

TemporalGrid::TemporalGrid( lps::Seconds sessionDuration, lps::Seconds resolution, lps::Seconds shift , ivg::Date start) : TemporalGrid( sessionDuration, resolution, shift ){
    start_ = start;
}
std::pair<unsigned, unsigned> TemporalGrid::getIndexOtherSampling( unsigned t, lps::Seconds intervalLength ) const{
    
    if(t > n_intervals){
        std::cerr << "Out of bounds in TemporalGrid::getIndexOtherSampling";
    }
    if(intervalLength > resolution_){
        std::cerr << "Interval length has to be smaller than resolution of temporal grid";
    }
    
    unsigned first = ceil( intervalBeginnings_[t]/intervalLength );
    unsigned second = ceil( (intervalBeginnings_[t]+resolution_)/intervalLength )-1;
    
    std::pair<unsigned, unsigned> pair = make_pair(first, second);
    return pair;
       
}

std::vector<unsigned> TemporalGrid::getTemporalIndices( lps::Seconds t ) const{
        
    
    std::vector<unsigned> indices;
    for( unsigned i = 0; i < n_intervals; ++i ){
        if(t >= intervalBeginnings_[i]  &&  t < intervalBeginnings_[i]+resolution_  ){
            indices.push_back(i);
        } else if( t < intervalBeginnings_[i] ) {
            break;
        }
    }
                    
//    std::cerr << t << " in :";
//    for(unsigned& i: indices ){
//        std::cerr << "\t" << intervalBeginnings_[i] << "-" << intervalBeginnings_[i]+resolution_ << std::endl;
//    }
    
    return indices;
}

 std::vector<unsigned> TemporalGrid::getTemporalIndices( ivg::Date d ) const{
    lps::Seconds t = ( d.get_double_mjd() - start_.get_double_mjd() )*24*3600;
    return getTemporalIndices(t);
 }
 
 std::pair<lps::Seconds, lps::Seconds> TemporalGrid::getIntervalBorder( int i  ) const{
    if(i > n_intervals){
        std::cerr << "Out of bounds in TemporalGrid::getIntervalBorder";
    }
    std::pair<lps::Seconds, lps::Seconds> pair = make_pair(intervalBeginnings_[i], intervalBeginnings_[i]+resolution_);
    return pair;
 }
 
 std::vector<ivg::Date> TemporalGrid::get_intervalBeginningsDate() const{
    std::vector<ivg::Date> dates(intervalBeginnings_.size(), start_);
    for( unsigned i = 0; i < intervalBeginnings_.size(); ++i ){
        dates[i].add_secs( intervalBeginnings_[i] );
    }
    return dates;
 }
 
std::vector<lps::Seconds> TemporalGrid::get_intervalCenters() const { 
    std::vector<lps::Seconds> ret = intervalBeginnings_;
    std::transform(ret.begin(), ret.end(), ret.begin(), bind2nd(std::plus<double>(), resolution_/2.0) );
    return ret;
};

std::vector<ivg::Date> TemporalGrid::get_intervalCentersDate() const { 
    std::vector<lps::Seconds> secs = get_intervalCenters();
    std::vector<ivg::Date> dates(secs.size());
    for( unsigned i = 0; i < secs.size(); ++i ){
        ivg::Date start = start_;
        start.add_secs( secs[i] );
        dates[i] = start;
    }
    return dates;
};

} //namespace
