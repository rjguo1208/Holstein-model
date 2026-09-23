#pragma once
#include "diagmc.hpp"
#include <limits>

namespace diagmc {
// Regularized lower incomplete gamma P(s,x), for positive integer s.
// Lower-tail series and upper-tail continued fraction avoid cancellation.
inline double gamma_p(unsigned s, double x) {
    if (!(s>0 && x>=0)) throw std::runtime_error("Invalid gamma arguments");
    if (x==0) return 0;
    const double prefactor=std::exp(s*std::log(x)-x-std::lgamma(static_cast<double>(s)));
    if (x<s+1.) {
        double term=1./s, sum=term;
        for (unsigned j=1;j<10000;++j) {
            term*=x/(s+j); sum+=term;
            if (std::abs(term)<std::abs(sum)*2e-15) return prefactor*sum;
        }
    } else {
        constexpr double tiny=1e-300;
        double b=x+1-s,c=1/tiny,d=1/b,h=d;
        for (unsigned j=1;j<10000;++j) {
            const double an=-static_cast<double>(j)*(static_cast<double>(j)-s);
            b+=2; d=an*d+b; if (std::abs(d)<tiny) d=tiny;
            c=b+an/c; if (std::abs(c)<tiny) c=tiny;
            d=1/d; const double delta=d*c; h*=delta;
            if (std::abs(delta-1)<2e-15) return 1-prefactor*h;
        }
    }
    throw std::runtime_error("Incomplete gamma failed to converge");
}

// In coordinates x_i=u_i/tau,y_i=v_i/tau, the conditional density is
// p(tau|n,x,y,q) proportional to tau^(2n)*exp(-rate*tau), 0<tau<T.
// This changes only the measurement; the sampled diagrams are unchanged.
inline void conditional_histogram(unsigned n,double rate,double tau_max,
                                  std::vector<double>& histogram) {
    if (!(rate>0)) throw std::runtime_error("Conditional estimator requires positive action rate");
    const unsigned shape=2*n+1;
    const double norm=gamma_p(shape,rate*tau_max);
    if (!(norm>std::numeric_limits<double>::min()))
        throw std::runtime_error("Conditional normalization underflow");
    double previous=0;
    for (unsigned j=0;j<histogram.size();++j) {
        const double next=(j+1==histogram.size()) ? 1. : gamma_p(shape,rate*tau_max*(j+1)/histogram.size())/norm;
        if (next<previous-1e-12) throw std::runtime_error("Nonmonotonic gamma CDF");
        histogram[j]+=std::max(0.,next-previous); previous=next;
    }
}
}
