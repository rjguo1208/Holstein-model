#include "../src/defect.hpp"
#include <iostream>
void close(double a,double b,double tol=1e-10){if(!std::isfinite(a)||!std::isfinite(b)||std::abs(a-b)>tol)throw std::runtime_error("Defect kernel identity failed");}
int main(){
    using namespace defect;
    Model m{1,1.2,.7,1.3,-3.,0};
    State empty{5,{},{}},one{5,{{1,1},{4,-1}},{{2,3}}};
    auto a=action(one,m);if(!a.valid)throw std::runtime_error("Rejected valid path");
    close(a.residence,2);close(a.phonon_time,1);close(a.midpoint,1);close(a.nph_mid,1);
    close(a.value,1.2-2*1.3+15);
    State bad=one;bad.arcs={{.5,2}};if(action(bad,m).valid)throw std::runtime_error("Accepted nonlocal phonon line");
    State hop{5,{{1,1},{4,-1}},{}};
    close(hop_add_ratio(empty,action(empty,m),action(hop,m),m)+hop_remove_ratio(hop,action(hop,m),action(empty,m),m),0);
    close(hop_add_ratio(empty,action(empty,m),action(hop,m),m),2*std::log(5.)-3*m.U);
    close(arc_add_ratio(hop,action(hop,m),a,one.arcs[0],m)+arc_remove_ratio(one,a,action(hop,m),one.arcs[0],m),0);
    // Crossing contractions are allowed, including multiple phonons at one site.
    State crossed{5,{},{{1,3},{2,4}}};if(!action(crossed,m).valid)throw std::runtime_error("Lost crossing topology");
    for(double r:{.6,1.8}){
        auto scaled=one;scaled.tau*=r;for(auto& h:scaled.hops)h.time*=r;for(auto& l:scaled.arcs){l.u*=r;l.v*=r;}
        close(action(scaled,m).value,r*a.value);close(action(scaled,m).midpoint,a.midpoint);
    }
    // Independent quadrature checks the event Jacobian and energy derivative.
    for(const auto& s:{empty,hop,one,crossed}){
        const auto ac=action(s,m);const double rate=ac.value/s.tau,T=8.;const unsigned N=vertices(s);
        std::vector<double> p(16),energy(16);conditional(s,ac,m,T,p,energy);
        auto integral=[&](double lo,double hi,bool with_energy){
            double sum=0;const unsigned intervals=4000;
            for(unsigned j=0;j<=intervals;++j){
                const double t=lo+(hi-lo)*j/intervals;
                double v=(t==0 && N>0)?0:std::pow(t/T,N)*std::exp(-rate*t);
                if(with_energy)v*=t>0?rate+m.mu-N/t:rate+m.mu;
                sum+=(j==0||j==intervals?1:j%2?4:2)*v;
            }return sum*(hi-lo)/(3*intervals);
        };
        const double norm=integral(0,T,false);
        for(unsigned b=0;b<16;++b){close(p[b],integral(.5*b,.5*(b+1),false)/norm,2e-9);close(energy[b],integral(.5*b,.5*(b+1),true)/norm,2e-9);}
    }
    // Proposals preserve a closed, spatially consistent path over a live chain.
    Sampler sampler(m,10,32,64,117);
    for(unsigned i=0;i<10000;++i){sampler.step();if(!action(sampler.state,m).valid)throw std::runtime_error("Sampler left diagram space");}
    std::cout<<"Passed defect residence, locality, crossed diagrams, reversible insertion ratios, time Jacobian, conditional energy quadrature, live-path invariants.\n";
}
