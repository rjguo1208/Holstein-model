#include "../src/diagmc.hpp"
#include "../src/conditional.hpp"
#include <iostream>

void close(double a,double b,double tol=1e-10) {
    if (std::abs(a-b)>tol) throw std::runtime_error("Kernel identity failed");
}
int main() {
    using namespace diagmc;
    Model m{1,1.3,.8,.4,-3};
    auto eps=[&](double p) { return -2*m.t*std::cos(p); };
    State nested{5,{{1,4,.7},{2,3,-.4}}};
    State crossed{5,{{1,3,.7},{2,4,-.4}}};
    close(action(nested,m).kinetic,2*eps(m.k)+2*eps(m.k-.7)+eps(m.k-.3));
    close(action(crossed,m).kinetic,2*eps(m.k)+eps(m.k-.7)+eps(m.k-.3)+eps(m.k+.4));
    close(action(crossed,m).phonon,4*m.omega);
    auto permuted=crossed;std::reverse(permuted.arcs.begin(),permuted.arcs.end());
    close(action(permuted,m).value,action(crossed,m).value);
    State empty{5,{}};Arc l{1,3,.7};State one{5,{l}};
    const double add=log_add_ratio(empty,action(empty,m),action(one,m),l,m);
    const double expected=2*std::log(m.g)+std::log(5.)
        +std::log(-std::expm1(-4*m.omega))-std::log(m.omega)
        -2*(eps(m.k-.7)-eps(m.k));
    close(add,expected);
    close(add+log_remove_ratio(one,action(one,m),action(empty,m),l,m),0);
    for(double factor : {.6,1.8}) {
        auto scaled=crossed;scaled.tau*=factor;
        for(auto& arc:scaled.arcs) {arc.u*=factor;arc.v*=factor;}
        close(action(scaled,m).value,factor*action(crossed,m).value);
    }
    Sampler sampler(m,10,96,1);sampler.state=crossed;sampler.cached=action(crossed,m);
    if (sampler.crossing_pairs()!=1) throw std::runtime_error("Missed crossed topology");
    const double h=1e-5;
    {
        // Independent finite differences of a diagram's complete weight.
        auto plus=m,minus=m;const double step=1e-4;
        plus.k+=step;minus.k-=step;
        const auto a=action(crossed,m);
        const double wp=std::exp(a.value-action(crossed,plus).value);
        const double wm=std::exp(a.value-action(crossed,minus).value);
        close((wp-wm)/(2*step),-a.current,2e-6);
        close((wp-2+wm)/(step*step),a.current*a.current-a.curvature,2e-6);
    }
    auto target=[&](double tau) {
        auto s=crossed;const double r=tau/s.tau;s.tau=tau;
        for(auto& arc:s.arcs) {arc.u*=r;arc.v*=r;}
        return 2*s.arcs.size()*std::log(tau)-action(s,m).value;
    };
    close((target(5+h)-target(5-h))/(2*h),m.mu-sampler.energy_estimator(),1e-8);
    for (double x : {.001,.1,1.,5.,50.}) {
        close(gamma_p(1,x),-std::expm1(-x),1e-13);
        close(gamma_p(3,x),1-std::exp(-x)*(1+x+x*x/2),1e-13);
    }
    for (unsigned n : {0u,1u,5u,15u,40u}) {
        for (double rate : {.2,1.,4.}) {
            std::vector<double> histogram(16,0);
            conditional_histogram(n,rate,12.,histogram);
            double sum=0;
            for(double p:histogram) { if(p<0) throw std::runtime_error("Negative bin probability"); sum+=p; }
            close(sum,1,1e-12);
            // Independent Simpson integration of the unnormalized density.
            auto integrate=[&](double lo,double hi) {
                double result=0; const unsigned intervals=2000;
                for(unsigned j=0;j<=intervals;++j) {
                    const double tau=lo+(hi-lo)*j/intervals;
                    const double value=(tau==0 && n>0) ? 0 : std::pow(tau/12.,2*n)*std::exp(-rate*tau);
                    result+=(j==0 || j==intervals ? 1 : j%2 ? 4 : 2)*value;
                }
                return result*(hi-lo)/(3*intervals);
            };
            const double norm=integrate(0,12);
            for(unsigned j=0;j<16;++j) close(histogram[j],integrate(.75*j,.75*(j+1))/norm,2e-8);
        }
    }
    std::cout<<"Passed: diagram kinematics, phonon action, permutation invariance, detailed balance, stretch Jacobian, energy derivative, crossed topology.\n";
    std::cout<<"Passed: conditional imaginary-time histogram versus analytic gamma functions and independent quadrature.\n";
}
