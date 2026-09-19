#include "diagmc.hpp"
#include "conditional.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>

int main(int argc,char** argv) {
    try {
        std::map<std::string,std::string> opt;
        for (int i=1;i<argc;i+=2) {
            if (i+1>=argc || std::string(argv[i]).rfind("--",0)!=0) throw std::runtime_error("Arguments are --name value pairs");
            if (!opt.emplace(argv[i]+2,argv[i+1]).second) throw std::runtime_error("Duplicate argument");
        }
        const std::vector<std::string> allowed={"t","omega","g","k","mu","tau-max","max-order","bins","blocks","steps-per-block","warmup","thin","rb-thin","seed","out"};
        for (const auto& p:opt) if (std::find(allowed.begin(),allowed.end(),p.first)==allowed.end()) throw std::runtime_error("Unknown argument: "+p.first);
        auto number=[&](std::string key,double fallback) { return opt.count(key) ? std::stod(opt.at(key)) : fallback; };
        auto integer=[&](std::string key,std::uint64_t fallback) { return opt.count(key) ? std::stoull(opt.at(key)) : fallback; };
        diagmc::Model model{number("t",1),number("omega",1),number("g",1),number("k",0),number("mu",-2.6)};
        const double taumax=number("tau-max",24);
        const auto cap=integer("max-order",96), bins=integer("bins",96), blocks=integer("blocks",200);
        const auto steps=integer("steps-per-block",50000), warmup=integer("warmup",200000), thin=integer("thin",5), seed=integer("seed",1729);
        const auto rbthin=integer("rb-thin",0);
        if (model.t<0 || model.omega<=0 || model.g<0 || taumax<=0 || bins<1 || blocks<2 || steps<1 || thin<1 || steps%thin || cap>100000 || bins>100000 || !opt.count("out"))
            throw std::runtime_error("Invalid parameters or missing --out; steps-per-block must be a multiple of thin");
        for (double x : {model.t,model.omega,model.g,model.k,model.mu,taumax}) if (!std::isfinite(x)) throw std::runtime_error("Non-finite parameter");
        if (rbthin && (rbthin%thin || steps%rbthin || model.mu>=-2*model.t))
            throw std::runtime_error("rb-thin must divide steps, be a multiple of thin, and requires mu < -2t");
        const std::filesystem::path out=opt.at("out");
        if (std::filesystem::exists(out)) throw std::runtime_error("Refusing to overwrite output directory: "+out.string());
        std::filesystem::create_directories(out);
        auto metadata=[&](bool complete,double seconds,const diagmc::Sampler* s) {
            std::ofstream f(out/"run.json"); f<<std::setprecision(17);
            f<<"{\n  \"complete\": "<<(complete ? "true" : "false")<<",\n  \"method\": \"bare_G_continuous_tau_momentum\",\n"
              <<"  \"t\": "<<model.t<<", \"omega\": "<<model.omega<<", \"g\": "<<model.g<<", \"k\": "<<model.k<<", \"mu\": "<<model.mu<<",\n"
              <<"  \"tau_max\": "<<taumax<<", \"max_order\": "<<cap<<", \"bins\": "<<bins<<", \"blocks\": "<<blocks<<",\n"
              <<"  \"steps_per_block\": "<<steps<<", \"warmup\": "<<warmup<<", \"thin\": "<<thin<<", \"seed\": "<<seed<<",\n"
              <<"  \"rb_thin\": "<<rbthin<<",\n"
              <<"  \"seconds\": "<<seconds;
            if (s) {
                f<<",\n  \"cap_attempts\": "<<s->cap_attempts<<",\n  \"attempted\": [";
                for(unsigned j=0;j<5;++j) f<<(j ? "," : "")<<s->attempted[j];
                f<<"],\n  \"accepted\": [";
                for(unsigned j=0;j<5;++j) f<<(j ? "," : "")<<s->accepted[j];
                f<<"]";
            }
            f<<"\n}\n";
        };
        metadata(false,0,nullptr);
        const auto start=std::chrono::steady_clock::now();
        diagmc::Sampler sampler(model,taumax,cap,seed);
        for(std::uint64_t i=0;i<warmup;++i) sampler.step();
        sampler.reset_counters();
        std::ofstream data(out/"blocks.csv"); data<<std::setprecision(17)<<"block,reference,samples,crossings";
        for (std::string field : {"count","energy","current","mass"})
            for (unsigned j=0;j<bins;++j) data<<","<<field<<"_"<<j;
        data<<"\n";
        std::ofstream rbdata;
        if (rbthin) {
            rbdata.open(out/"rb_blocks.csv"); rbdata<<std::setprecision(17)<<"block,reference,samples,crossings";
            for (unsigned j=0;j<bins;++j) rbdata<<",count_"<<j;
            rbdata<<"\n";
        }
        std::vector<std::uint64_t> orders(cap+1,0);
        for(unsigned block=0;block<blocks;++block) {
            std::vector<std::uint64_t> counts(bins,0);
            std::vector<double> energy(bins,0),current(bins,0),mass(bins,0);
            std::uint64_t reference=0,measurements=0,crossings=0;
            std::vector<double> rbcounts(bins,0);
            std::uint64_t rbreference=0,rbmeasurements=0;
            for(std::uint64_t i=0;i<steps;++i) {
                sampler.step();
                if ((i+1)%thin) continue;
                const auto n=sampler.state.arcs.size();
                const unsigned j=std::min<unsigned>(bins-1,static_cast<unsigned>(sampler.state.tau/taumax*bins));
                ++counts[j]; ++measurements; ++orders[n]; reference+=(n==0);
                energy[j]+=sampler.energy_estimator();
                current[j]+=sampler.cached.current;
                mass[j]+=sampler.cached.current*sampler.cached.current-sampler.cached.curvature;
                if (n==2) crossings+=sampler.crossing_pairs();
                if (rbthin && (i+1)%rbthin==0) {
                    ++rbmeasurements; rbreference+=(n==0);
                    diagmc::conditional_histogram(n,sampler.cached.value/sampler.state.tau,taumax,rbcounts);
                }
            }
            data<<block<<","<<reference<<","<<measurements<<","<<crossings;
            for(auto x:counts) data<<","<<x;
            for(const auto& values : {energy,current,mass}) for(auto x:values) data<<","<<x;
            data<<"\n"; data.flush();
            if (rbthin) {
                rbdata<<block<<","<<rbreference<<","<<rbmeasurements<<",0";
                for (auto x:rbcounts) rbdata<<","<<x;
                rbdata<<"\n"; rbdata.flush();
            }
        }
        std::ofstream histogram(out/"orders.csv"); histogram<<"order,count\n";
        for(unsigned n=0;n<orders.size();++n) histogram<<n<<","<<orders[n]<<"\n";
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        metadata(true,seconds,&sampler);
        std::cout<<out.string()<<": "<<blocks*steps<<" steps in "<<seconds<<" s\n";
        return 0;
    } catch (const std::exception& e) { std::cerr<<e.what()<<"\n"; return 1; }
}
