#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

namespace diagmc {
constexpr double pi = 3.141592653589793238462643383279502884;
struct Model { double t=1, omega=1, g=1, k=0, mu=-2.6; };
struct Arc { double u, v, q; };
struct State { double tau=1; std::vector<Arc> arcs; };
struct Action { double kinetic=0, phonon=0, value=0, current=0, curvature=0; };

inline Action action(const State& s, const Model& m) {
    std::vector<std::pair<double,double>> events;
    events.reserve(2*s.arcs.size());
    Action a;
    for (const auto& l : s.arcs) {
        if (!(0 <= l.u && l.u < l.v && l.v <= s.tau))
            throw std::runtime_error("Invalid phonon endpoints");
        events.emplace_back(l.u,l.q);
        events.emplace_back(l.v,-l.q);
        a.phonon += m.omega*(l.v-l.u);
    }
    std::sort(events.begin(),events.end());
    double previous=0, qtotal=0;
    auto segment = [&](double length) {
        const double p=m.k-qtotal;
        a.kinetic += -2*m.t*std::cos(p)*length;
        a.current += 2*m.t*std::sin(p)*length;
        a.curvature += 2*m.t*std::cos(p)*length;
    };
    for (const auto& e : events) {
        segment(e.first-previous);
        qtotal += e.second;
        previous=e.first;
    }
    segment(s.tau-previous);
    a.value=a.kinetic+a.phonon-m.mu*s.tau;
    return a;
}

inline double log_add_ratio(const State& old, const Action& a,
        const Action& b, const Arc& l, const Model& m) {
    const double norm=-std::expm1(-m.omega*(old.tau-l.u));
    return 2*std::log(m.g)-(b.value-a.value)+std::log(old.tau)
        +std::log(norm)-std::log(m.omega)+m.omega*(l.v-l.u)
        -std::log(old.arcs.size()+1.0);
}
inline double log_remove_ratio(const State& old, const Action& a,
        const Action& b, const Arc& l, const Model& m) {
    const double norm=-std::expm1(-m.omega*(old.tau-l.u));
    return -2*std::log(m.g)-(b.value-a.value)-std::log(old.tau)
        -std::log(norm)+std::log(m.omega)-m.omega*(l.v-l.u)
        +std::log(static_cast<double>(old.arcs.size()));
}

class Sampler {
public:
    Model model;
    State state;
    Action cached;
    double tau_max;
    unsigned max_order;
    std::mt19937_64 rng;
    std::array<std::uint64_t,5> attempted{},accepted{};
    std::uint64_t cap_attempts=0;
    Sampler(Model m,double max_tau,unsigned cap,std::uint64_t seed)
        :model(m),tau_max(max_tau),max_order(cap),rng(seed) {
        state.tau=tau_max/2;
        cached=action(state,model);
    }
    double uniform() {
        double u;
        do { u=std::generate_canonical<double,53>(rng); } while (u==0);
        return u;
    }
    unsigned index(unsigned n) { return std::uniform_int_distribution<unsigned>(0,n-1)(rng); }
    bool accept(double ratio) { return ratio >= 0 || std::log(uniform()) < ratio; }
    void reset_counters() { attempted.fill(0); accepted.fill(0); cap_attempts=0; }
    void step() {
        // Fixed proposal probabilities, including null moves at n=0/cap.
        const double choice=uniform();
        const unsigned type=choice<.30 ? 0 : choice<.60 ? 1 : choice<.75 ? 2 : choice<.85 ? 3 : 4;
        ++attempted[type];
        State trial=state;
        const unsigned n=state.arcs.size();
        double ratio=0;
        if (type==0) {
            if (model.g==0) return;
            if (n==max_order) { ++cap_attempts; return; }
            const double u=state.tau*uniform();
            const double norm=-std::expm1(-model.omega*(state.tau-u));
            const double length=-std::log1p(-uniform()*norm)/model.omega;
            const Arc l{u,u+length,2*pi*uniform()-pi};
            if (!(l.u<l.v && l.v<state.tau)) return;
            trial.arcs.push_back(l);
            Action next=action(trial,model);
            ratio=log_add_ratio(state,cached,next,l,model);
            if (accept(ratio)) { state=std::move(trial); cached=next; ++accepted[type]; }
            return;
        }
        if (type==1) {
            if (n==0) return;
            const unsigned j=index(n);
            const Arc l=trial.arcs[j];
            trial.arcs.erase(trial.arcs.begin()+j);
            Action next=action(trial,model);
            ratio=log_remove_ratio(state,cached,next,l,model);
            if (accept(ratio)) { state=std::move(trial); cached=next; ++accepted[type]; }
            return;
        }
        if (type==2) {
            if (n==0) return;
            trial.arcs[index(n)].q=2*pi*uniform()-pi;
        } else if (type==3) {
            if (n==0) return;
            auto& l=trial.arcs[index(n)];
            if (uniform()<.5) l.u=uniform()*l.v;
            else l.v=l.u+uniform()*(state.tau-l.u);
            if (!(l.u<l.v && l.v<state.tau)) return;
        } else {
            if (n==0) {
                const double rate=-2*model.t*std::cos(model.k)-model.mu;
                if (std::abs(rate)<1e-12) trial.tau=tau_max*uniform();
                else if (rate>0) trial.tau=-std::log1p(uniform()*std::expm1(-rate*tau_max))/rate;
                else trial.tau=tau_max+std::log1p(uniform()*std::expm1(rate*tau_max))/(-rate);
                state=std::move(trial); cached=action(state,model); ++accepted[type];
                return;
            }
            const bool logarithmic=uniform()<.5;
            trial.tau=logarithmic ? state.tau*std::exp(1.2*(uniform()-.5)) : tau_max*uniform();
            if (!(trial.tau>0 && trial.tau<tau_max)) return;
            const double r=trial.tau/state.tau;
            for (auto& l : trial.arcs) { l.u*=r; l.v*=r; }
            // The 2n endpoint Jacobian is essential when stretching the diagram.
            ratio=(2*n+(logarithmic ? 1 : 0))*std::log(r)-cached.value*(r-1);
            if (accept(ratio)) { state=std::move(trial); cached=action(state,model); ++accepted[type]; }
            return;
        }
        const Action next=action(trial,model);
        if (accept(cached.value-next.value)) { state=std::move(trial); cached=next; ++accepted[type]; }
    }
    double energy_estimator() const {
        return (cached.kinetic+cached.phonon-2*state.arcs.size())/state.tau;
    }
    unsigned crossing_pairs() const {
        unsigned count=0;
        for (unsigned i=0;i<state.arcs.size();++i)
            for (unsigned j=i+1;j<state.arcs.size();++j) {
                Arc a=state.arcs[i], b=state.arcs[j];
                if (a.u>b.u) std::swap(a,b);
                count += a.u<b.u && b.u<a.v && a.v<b.v;
            }
        return count;
    }
};
}
