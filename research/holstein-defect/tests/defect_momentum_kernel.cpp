#include "../src/defect_momentum.hpp"
#include <iostream>
void close(double a,double b,double tol=1e-11){if(!std::isfinite(a)||!std::isfinite(b)||std::abs(a-b)>tol)throw std::runtime_error("Open-endpoint identity failed");}
int main(){
    using namespace defect_momentum;
    Model m{1,1.2,.7,1.3,-3.,-1};std::vector<int> x;
    State empty{5,{},{}},one{5,{{1,1}},{}},pair{5,{{1,1},{2,1},{4,-1}}, {}};
    auto a=action(empty,m,2,x),b=action(one,m,2,x),c=action(pair,m,2,x);
    if(!a.valid||!b.valid||!c.valid)throw std::runtime_error("Valid open path rejected");
    close(b.residence,4);close(b.endpoint,0);close(c.residence,2);
    close(add_single_ratio(empty,a,b,m)+remove_single_ratio(one,b,a,m),0);
    close(add_single_ratio(empty,a,b,m),std::log(10.)+4*m.U);
    close(add_pair_ratio(one,b,c,m)+remove_pair_ratio(pair,c,b,m),0);
    State arc=one;arc.arcs.push_back({2,3});auto d=action(arc,m,2,x);
    close(defect::arc_add_ratio(one,b,d,arc.arcs[0],m)+defect::arc_remove_ratio(arc,d,b,arc.arcs[0],m),0);
    arc.arcs[0]={.5,3};if(action(arc,m,2,x).valid)throw std::runtime_error("Nonlocal phonon allowed");
    State beyond{5,{{1,-1},{2,-1}}, {}};
    if(action(beyond,m,2,x).valid)throw std::runtime_error("Endpoint outside window allowed");
    State excursion{5,{{1,-1},{2,-1},{3,1}}, {}};
    if(!action(excursion,m,2,x).valid)throw std::runtime_error("Intermediate path spuriously confined");
    Sampler sampler(m,6,2,128,256,3712);sampler.state=one;sampler.cached=action(one,m,2,x);
    double observed=0;const unsigned n=100000;
    for(unsigned i=0;i<n;++i){sampler.translate_origin();observed+=sampler.model.origin==-1;}
    // origins -2,-1,0,1: residence 0,4,1,0.
    double expected=std::exp(4*m.U)/(2+std::exp(4*m.U)+std::exp(m.U));
    close(observed/n,expected,.004);
    for(unsigned i=0;i<10000;++i){sampler.step();if(!sampler.cached.valid)throw std::runtime_error("Invalid accepted state");}
    std::cout<<"Open-endpoint balance, nonlocal sectors, heat bath and infinite-path checks passed\n";
}
