#pragma once
#include "conditional.hpp"
#include <limits>

// Continuous-time expansion in hopping and electron-phonon vertices.
// The static potential is exact. All electron paths live on the infinite chain.
namespace defect {
struct Model { double t=1, omega=1, g=std::sqrt(.5), U=1, mu=-2.65; int origin=0; };
struct Hop { double time; int direction; };
struct Arc { double u,v; };
struct State { double tau=1; std::vector<Hop> hops; std::vector<Arc> arcs; };
struct Action {
    double value=0, residence=0, phonon_time=0;
    int midpoint=0, nph_mid=0, extent=0;
    bool valid=true;
};
inline unsigned vertices(const State& s) { return s.hops.size()+2*s.arcs.size(); }
inline void sort_hops(State& s) {
    std::sort(s.hops.begin(),s.hops.end(),[](const Hop& a,const Hop& b){return a.time<b.time;});
}
inline Action action(const State& s,const Model& m) {
    Action a; int position=m.origin; double previous=0;
    std::vector<int> positions{position}; positions.reserve(s.hops.size()+1);
    for (const auto& h:s.hops) {
        if (!(h.time>previous && h.time<s.tau) || std::abs(h.direction)!=1) {a.valid=false;return a;}
        if (position==0) a.residence+=h.time-previous;
        position+=h.direction; positions.push_back(position);
        a.extent=std::max(a.extent,std::abs(position)); previous=h.time;
    }
    if (position!=m.origin) {a.valid=false;return a;}
    if (position==0) a.residence+=s.tau-previous;
    auto at=[&](double time) {
        const auto j=std::upper_bound(s.hops.begin(),s.hops.end(),time,
            [](double t,const Hop& h){return t<h.time;})-s.hops.begin();
        return positions[j];
    };
    a.midpoint=at(s.tau/2);
    for (const auto& l:s.arcs) {
        if (!(0<l.u && l.u<l.v && l.v<s.tau) || at(l.u)!=at(l.v)) {a.valid=false;return a;}
        a.phonon_time+=l.v-l.u;
        a.nph_mid+=(l.u<s.tau/2 && l.v>s.tau/2);
    }
    a.value=m.omega*a.phonon_time-m.U*a.residence-m.mu*s.tau;
    return a;
}
inline double hop_add_ratio(const State& old,const Action& a,const Action& b,const Model& m) {
    // Ordered times from two independent uniform draws, then direction +/-.
    // q_add=1/tau^2; q_remove=1/(number_of_positive_hops)^2.
    return 2*std::log(m.t*old.tau/(old.hops.size()/2.+1)) + a.value-b.value;
}
inline double hop_remove_ratio(const State& old,const Action& a,const Action& b,const Model& m) {
    return -2*std::log(m.t*old.tau/(old.hops.size()/2.)) + a.value-b.value;
}
inline double arc_add_ratio(const State& old,const Action& a,const Action& b,const Arc& l,const Model& m) {
    return 2*std::log(m.g)+a.value-b.value+std::log(old.tau)
        +std::log(-std::expm1(-m.omega*(old.tau-l.u)))-std::log(m.omega)
        +m.omega*(l.v-l.u)-std::log(old.arcs.size()+1.);
}
inline double arc_remove_ratio(const State& old,const Action& a,const Action& b,const Arc& l,const Model& m) {
    return -2*std::log(m.g)+a.value-b.value-std::log(old.tau)
        -std::log(-std::expm1(-m.omega*(old.tau-l.u)))+std::log(m.omega)
        -m.omega*(l.v-l.u)+std::log(static_cast<double>(old.arcs.size()));
}

class Sampler {
public:
    Model model; State state; Action cached;
    double tau_max; unsigned max_arcs,max_hop_pairs;
    std::mt19937_64 rng;
    std::array<std::uint64_t,7> attempted{},accepted{};
    std::uint64_t arc_cap_attempts=0,hop_cap_attempts=0,invalid_proposals=0;
    Sampler(Model m,double T,unsigned cap,unsigned hopcap,std::uint64_t seed)
        :model(m),tau_max(T),max_arcs(cap),max_hop_pairs(hopcap),rng(seed) {
        state.tau=std::min(1.,T/2);cached=action(state,model);
    }
    double uniform() {double u;do {u=std::generate_canonical<double,53>(rng);}while(u==0);return u;}
    unsigned index(unsigned n) {return std::uniform_int_distribution<unsigned>(0,n-1)(rng);}
    void reset() {attempted.fill(0);accepted.fill(0);arc_cap_attempts=hop_cap_attempts=invalid_proposals=0;}
    void step() {
        // Symmetric, state-independent add/remove probabilities (including nulls).
        const double c=uniform();
        const unsigned type=c<.16?0:c<.32?1:c<.48?2:c<.64?3:c<.76?4:c<.88?5:6;
        ++attempted[type]; State trial=state;double ratio=0;
        if (type==0) {
            if (model.t==0) return;
            if (state.hops.size()/2==max_hop_pairs) {++hop_cap_attempts;return;}
            double u=uniform()*state.tau,v=uniform()*state.tau;if(u>v)std::swap(u,v);
            const int direction=uniform()<.5?1:-1;
            trial.hops.push_back({u,direction});trial.hops.push_back({v,-direction});sort_hops(trial);
        } else if (type==1) {
            if(state.hops.empty())return;
            std::vector<unsigned> positive,negative;
            for(unsigned i=0;i<state.hops.size();++i)(state.hops[i].direction>0?positive:negative).push_back(i);
            unsigned i=positive[index(positive.size())],j=negative[index(negative.size())];if(i<j)std::swap(i,j);
            trial.hops.erase(trial.hops.begin()+i);trial.hops.erase(trial.hops.begin()+j);
        } else if(type==2) {
            if(model.g==0)return;
            if(state.arcs.size()==max_arcs){++arc_cap_attempts;return;}
            const double u=uniform()*state.tau;
            const double v=u-std::log1p(-uniform()*(-std::expm1(-model.omega*(state.tau-u))))/model.omega;
            trial.arcs.push_back({u,v});
        } else if(type==3) {
            if(state.arcs.empty())return;
            trial.arcs.erase(trial.arcs.begin()+index(trial.arcs.size()));
        } else if(type==4) {
            if(state.hops.empty())return;
            auto& h=trial.hops[index(trial.hops.size())];
            h.time=uniform()<.5?uniform()*state.tau:h.time+2*(uniform()-.5);
            if(!(h.time>0 && h.time<state.tau))return;
            sort_hops(trial);
        } else if(type==5) {
            if(state.arcs.empty())return;
            auto& l=trial.arcs[index(trial.arcs.size())];
            if(uniform()<.5)l.u=uniform()*l.v;else l.v=l.u+uniform()*(state.tau-l.u);
        } else {
            if(vertices(state)==0) {
                const double rate=-model.mu-(model.origin==0?model.U:0);
                trial.tau=-std::log1p(-uniform()*(-std::expm1(-rate*tau_max)))/rate;
                state=std::move(trial);cached=action(state,model);++accepted[type];return;
            }
            const bool logarithmic=uniform()<.7;
            trial.tau=logarithmic?state.tau*std::exp(.8*(uniform()-.5)):tau_max*uniform();
            if(!(trial.tau>0 && trial.tau<tau_max))return;
            const double r=trial.tau/state.tau;
            for(auto& h:trial.hops)h.time*=r;
            for(auto& l:trial.arcs){l.u*=r;l.v*=r;}
            ratio=(vertices(state)+(logarithmic?1:0))*std::log(r);
        }
        const Action next=action(trial,model);
        if(!next.valid){++invalid_proposals;return;}
        if(type==0)ratio=hop_add_ratio(state,cached,next,model);
        if(type==1)ratio=hop_remove_ratio(state,cached,next,model);
        if(type==2)ratio=arc_add_ratio(state,cached,next,trial.arcs.back(),model);
        if(type==3) {
            // Identify the deleted arc; vector order is retained by erase.
            unsigned j=0;while(j<trial.arcs.size() && trial.arcs[j].u==state.arcs[j].u && trial.arcs[j].v==state.arcs[j].v)++j;
            ratio=arc_remove_ratio(state,cached,next,state.arcs[j],model);
        }
        if(type>=4)ratio+=cached.value-next.value;
        if(ratio>=0 || std::log(uniform())<ratio){state=std::move(trial);cached=next;++accepted[type];}
    }
};

// Exact integration over the external time at fixed scaled vertex positions.
// Density proportional to tau^N exp(-a*tau), N=N_hop+2*N_arc.
inline void conditional(const State& s,const Action& a,const Model& m,double T,
                        std::vector<double>& probability,std::vector<double>& energy) {
    const unsigned N=vertices(s),shape=N+1;
    const double rate=a.value/s.tau,norm=diagmc::gamma_p(shape,rate*T);
    if(!(rate>0 && norm>std::numeric_limits<double>::min()))throw std::runtime_error("Invalid conditional time normalization");
    double previous=0,previous_lower=0;
    for(unsigned i=0;i<probability.size();++i){
        const double edge=T*(i+1)/probability.size();
        const double next=i+1==probability.size()?1:diagmc::gamma_p(shape,rate*edge)/norm;
        probability[i]=std::max(0.,next-previous);
        const double lower=N?diagmc::gamma_p(N,rate*edge)/norm:0;
        energy[i]=(rate+m.mu)*probability[i]-rate*(lower-previous_lower);
        previous=next;previous_lower=lower;
    }
}
}
