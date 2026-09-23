#pragma once
#include "defect.hpp"

// Open electron paths. Both endpoints lie in the observation window, while
// intermediate positions are unrestricted on the infinite integer chain.
namespace defect_momentum {
using defect::State;
using defect::Model;
using defect::Arc;
using defect::Hop;
struct Action : defect::Action { int endpoint=0; };

inline Action action(const State& s,const Model& m,int window,std::vector<int>& positions) {
    Action a;int position=m.origin;double previous=0;
    if(std::abs(position)>window){a.valid=false;return a;}
    positions.clear();positions.reserve(s.hops.size()+1);positions.push_back(position);
    a.extent=std::abs(position);
    for(const auto& h:s.hops) {
        if(!(h.time>previous && h.time<s.tau) || std::abs(h.direction)!=1){a.valid=false;return a;}
        if(position==0)a.residence+=h.time-previous;
        position+=h.direction;positions.push_back(position);
        a.extent=std::max(a.extent,std::abs(position));previous=h.time;
    }
    a.endpoint=position;
    if(std::abs(position)>window){a.valid=false;return a;}
    if(position==0)a.residence+=s.tau-previous;
    auto at=[&](double time){
        auto j=std::upper_bound(s.hops.begin(),s.hops.end(),time,
            [](double t,const Hop& h){return t<h.time;})-s.hops.begin();
        return positions[j];
    };
    for(const auto& l:s.arcs) {
        if(!(0<l.u && l.u<l.v && l.v<s.tau) || at(l.u)!=at(l.v)){a.valid=false;return a;}
        a.phonon_time+=l.v-l.u;
    }
    a.value=m.omega*a.phonon_time-m.U*a.residence-m.mu*s.tau;
    return a;
}

inline double add_single_ratio(const State& s,const Action& a,const Action& b,const Model& m) {
    return std::log(2*m.t*s.tau/(s.hops.size()+1.))+a.value-b.value;
}
inline double remove_single_ratio(const State& s,const Action& a,const Action& b,const Model& m) {
    return -std::log(2*m.t*s.tau/s.hops.size())+a.value-b.value;
}
inline unsigned positives(const State& s) {
    return std::count_if(s.hops.begin(),s.hops.end(),[](const Hop& h){return h.direction>0;});
}
inline double add_pair_ratio(const State& s,const Action& a,const Action& b,const Model& m) {
    const unsigned np=positives(s),nm=s.hops.size()-np;
    return 2*std::log(m.t*s.tau)-std::log((np+1.)*(nm+1.))+a.value-b.value;
}
inline double remove_pair_ratio(const State& s,const Action& a,const Action& b,const Model& m) {
    const unsigned np=positives(s),nm=s.hops.size()-np;
    return -2*std::log(m.t*s.tau)+std::log(static_cast<double>(np*nm))+a.value-b.value;
}

class Sampler {
public:
    Model model;State state,scratch;Action cached;
    std::vector<int> positions;
    std::vector<double> residence,weights;
    double tau_max;int window;unsigned max_arcs,max_hops;
    std::mt19937_64 rng;
    std::array<std::uint64_t,10> attempted{},accepted{};
    std::uint64_t arc_cap_attempts=0,hop_cap_attempts=0,invalid_proposals=0;
    Sampler(Model m,double T,int W,unsigned cap,unsigned hopcap,std::uint64_t seed)
      :model(m),residence(2*W+1),weights(2*W+1),tau_max(T),window(W),max_arcs(cap),max_hops(hopcap),rng(seed){
        state.tau=std::min(1.,T/2);cached=action(state,model,window,positions);
    }
    double uniform(){double u;do{u=std::generate_canonical<double,53>(rng);}while(u==0);return u;}
    unsigned index(unsigned n){return std::uniform_int_distribution<unsigned>(0,n-1)(rng);}
    void reset(){attempted.fill(0);accepted.fill(0);arc_cap_attempts=hop_cap_attempts=invalid_proposals=0;}
    void translate_origin(){
        // Exact heat bath over all origins compatible with the endpoints.
        // At fixed relative path only the defect residence changes.
        const int d=cached.endpoint-model.origin;
        const int lo=std::max(-window,-window-d),hi=std::min(window,window-d);
        std::fill(residence.begin(),residence.end(),0.);
        int relative=0;double previous=0;
        auto add=[&](double duration){
            int origin=-relative;
            if(origin>=lo && origin<=hi)residence[origin+window]+=duration;
        };
        for(const auto& h:state.hops){add(h.time-previous);relative+=h.direction;previous=h.time;}
        add(state.tau-previous);
        double largest=0.;for(int o=lo;o<=hi;++o)largest=std::max(largest,model.U*residence[o+window]);
        double sum=0.;for(int o=lo;o<=hi;++o){weights[o+window]=std::exp(model.U*residence[o+window]-largest);sum+=weights[o+window];}
        double u=uniform()*sum;int chosen=hi;
        for(int o=lo;o<=hi;++o){u-=weights[o+window];if(u<=0){chosen=o;break;}}
        model.origin=chosen;cached=action(state,model,window,positions);
    }
    void step(){
        const double c=uniform();
        const unsigned type=c<.12?0:c<.24?1:c<.32?2:c<.40?3:c<.55?4:c<.70?5:c<.77?6:c<.84?7:c<.94?8:9;
        ++attempted[type];
        if(type==9){translate_origin();++accepted[type];return;}
        State& trial=scratch;trial=state;double ratio=0.;unsigned deleted=0;
        if(type==0){
            if(model.t==0)return;
            if(state.hops.size()>=max_hops){++hop_cap_attempts;return;}
            trial.hops.push_back({uniform()*state.tau,uniform()<.5?1:-1});defect::sort_hops(trial);
        }else if(type==1){
            if(state.hops.empty())return;
            trial.hops.erase(trial.hops.begin()+index(trial.hops.size()));
        }else if(type==2){
            if(model.t==0)return;
            if(state.hops.size()+2>max_hops){++hop_cap_attempts;return;}
            double u=uniform()*state.tau,v=uniform()*state.tau;if(u>v)std::swap(u,v);
            int direction=uniform()<.5?1:-1;
            trial.hops.push_back({u,direction});trial.hops.push_back({v,-direction});defect::sort_hops(trial);
        }else if(type==3){
            const unsigned np=positives(state),nm=state.hops.size()-np;
            if(!np || !nm)return;
            unsigned ip=index(np),im=index(nm),i=0,j=0,p=0,n=0;
            for(unsigned k=0;k<state.hops.size();++k){
                if(state.hops[k].direction>0){if(p++==ip)i=k;}else{if(n++==im)j=k;}
            }
            if(i<j)std::swap(i,j);
            trial.hops.erase(trial.hops.begin()+i);trial.hops.erase(trial.hops.begin()+j);
        }else if(type==4){
            if(model.g==0)return;
            if(state.arcs.size()==max_arcs){++arc_cap_attempts;return;}
            double u=uniform()*state.tau;
            double v=u-std::log1p(-uniform()*(-std::expm1(-model.omega*(state.tau-u))))/model.omega;
            trial.arcs.push_back({u,v});
        }else if(type==5){
            if(state.arcs.empty())return;
            deleted=index(state.arcs.size());trial.arcs.erase(trial.arcs.begin()+deleted);
        }else if(type==6){
            if(state.hops.empty())return;
            auto& h=trial.hops[index(state.hops.size())];
            h.time=uniform()<.5?uniform()*state.tau:h.time+2*(uniform()-.5);
            if(!(h.time>0 && h.time<state.tau))return;
            defect::sort_hops(trial);
        }else if(type==7){
            if(state.arcs.empty())return;
            auto& l=trial.arcs[index(state.arcs.size())];
            if(uniform()<.5)l.u=uniform()*l.v;else l.v=l.u+uniform()*(state.tau-l.u);
        }else{
            if(defect::vertices(state)==0){
                double rate=-model.mu-(model.origin==0?model.U:0);
                trial.tau=-std::log1p(-uniform()*(-std::expm1(-rate*tau_max)))/rate;
                std::swap(state,trial);cached=action(state,model,window,positions);++accepted[type];return;
            }
            bool logarithmic=uniform()<.7;
            trial.tau=logarithmic?state.tau*std::exp(.8*(uniform()-.5)):tau_max*uniform();
            if(!(trial.tau>0 && trial.tau<tau_max))return;
            double r=trial.tau/state.tau;
            for(auto& h:trial.hops)h.time*=r;
            for(auto& l:trial.arcs){l.u*=r;l.v*=r;}
            ratio=(defect::vertices(state)+(logarithmic?1:0))*std::log(r);
        }
        Action next=action(trial,model,window,positions);
        if(!next.valid){++invalid_proposals;return;}
        if(type==0)ratio=add_single_ratio(state,cached,next,model);
        if(type==1)ratio=remove_single_ratio(state,cached,next,model);
        if(type==2)ratio=add_pair_ratio(state,cached,next,model);
        if(type==3)ratio=remove_pair_ratio(state,cached,next,model);
        if(type==4)ratio=defect::arc_add_ratio(state,cached,next,trial.arcs.back(),model);
        if(type==5)ratio=defect::arc_remove_ratio(state,cached,next,state.arcs[deleted],model);
        if(type>=6)ratio+=cached.value-next.value;
        if(ratio>=0 || std::log(uniform())<ratio){std::swap(state,trial);cached=next;++accepted[type];}
    }
};
}
