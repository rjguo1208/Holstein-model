#include "defect.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>

int main(int argc,char** argv){try{
    std::map<std::string,std::string> opt;
    for(int i=1;i<argc;i+=2){
        if(i+1>=argc || std::string(argv[i]).rfind("--",0)!=0 || !opt.emplace(argv[i]+2,argv[i+1]).second)throw std::runtime_error("Use unique --name value arguments");
    }
    const std::vector<std::string> allowed={"t","omega","g","U","mu","origin","tau-max","max-arcs","max-hop-pairs","bins","blocks","steps-per-block","warmup","thin","radius","seed","out","observables"};
    for(const auto& p:opt)if(std::find(allowed.begin(),allowed.end(),p.first)==allowed.end())throw std::runtime_error("Unknown argument: "+p.first);
    auto number=[&](std::string k,double d){return opt.count(k)?std::stod(opt.at(k)):d;};
    auto integer=[&](std::string k,std::uint64_t d){if(opt.count(k)&&opt.at(k).front()=='-')throw std::runtime_error("Negative integer option");return opt.count(k)?std::stoull(opt.at(k)):d;};
    defect::Model m{number("t",1),number("omega",1),number("g",std::sqrt(.5)),number("U",1),number("mu",-2.65),static_cast<int>(number("origin",0))};
    const double T=number("tau-max",48);
    const std::string observables=opt.count("observables")?opt.at("observables"):"all";
    if(observables!="all" && observables!="green")throw std::runtime_error("observables must be all or green");
    const bool projectors=observables=="all";
    const auto bins=integer("bins",48),blocks=integer("blocks",160),steps=integer("steps-per-block",100000),thin=integer("thin",100),warmup=integer("warmup",500000),radius=integer("radius",24),seed=integer("seed",2901001),cap=integer("max-arcs",128),hopcap=integer("max-hop-pairs",256);
    if(m.t<0||m.omega<=0||m.g<0||m.U<0||m.mu>=-m.U||T<=0||bins<8||bins%4||bins>1024||blocks<2||!steps||!thin||steps%thin||radius>1000||cap>10000||hopcap>10000||!opt.count("out"))throw std::runtime_error("Invalid parameters; require mu < -U, bins divisible by 4 and thin divides steps");
    for(double x:{m.t,m.omega,m.g,m.U,m.mu,T,number("origin",0)})if(!std::isfinite(x))throw std::runtime_error("Non-finite argument");
    if(number("origin",0)!=m.origin)throw std::runtime_error("origin must be an integer");
    const std::filesystem::path out=opt.at("out");
    if(std::filesystem::exists(out))throw std::runtime_error("Refusing to overwrite output");
    std::filesystem::create_directories(out);
    defect::Sampler sampler(m,T,cap,hopcap,seed);
    std::uint64_t maximum_arcs=0,maximum_hops=0,maximum_extent=0;
    auto metadata=[&](bool complete,double seconds){
        std::ofstream f(out/"run.json");f<<std::setprecision(17);
        f<<"{\n\"complete\":"<<(complete?"true":"false")<<",\n\"method\":\"real_space_bare_hopping_and_phonon_diagrams\",\n"
         <<"\"observables\":\""<<observables<<"\",\n"
         <<"\"t\":"<<m.t<<",\"omega\":"<<m.omega<<",\"g\":"<<m.g<<",\"U\":"<<m.U<<",\"mu\":"<<m.mu<<",\"origin\":"<<m.origin<<",\n"
         <<"\"tau_max\":"<<T<<",\"bins\":"<<bins<<",\"blocks\":"<<blocks<<",\"steps_per_block\":"<<steps<<",\"thin\":"<<thin<<",\"warmup\":"<<warmup<<",\"seed\":"<<seed<<",\n"
         <<"\"radius\":"<<radius<<",\"spatial_boundary\":\"none; radius is a reporting window with overflow\",\"max_arcs\":"<<cap<<",\"max_hop_pairs\":"<<hopcap<<",\n"
         <<"\"arc_cap_attempts\":"<<sampler.arc_cap_attempts<<",\"hop_cap_attempts\":"<<sampler.hop_cap_attempts<<",\"invalid_proposals\":"<<sampler.invalid_proposals<<",\n"
         <<"\"maximum_arcs\":"<<maximum_arcs<<",\"maximum_hops\":"<<maximum_hops<<",\"maximum_extent\":"<<maximum_extent<<",\"seconds\":"<<seconds<<",\n\"attempted\":[";
        for(unsigned j=0;j<7;++j)f<<(j?",":"")<<sampler.attempted[j];
        f<<"],\"accepted\":[";
        for(unsigned j=0;j<7;++j)f<<(j?",":"")<<sampler.accepted[j];
        f<<"]\n}\n";
    };
    metadata(false,0);const auto start=std::chrono::steady_clock::now();
    for(std::uint64_t i=0;i<warmup;++i)sampler.step();
    sampler.reset();
    std::ofstream csv(out/"blocks.csv"),prof;if(projectors)prof.open(out/"profiles.csv");
    csv<<std::setprecision(17);prof<<std::setprecision(17);
    csv<<"block,reference,samples";
    for(auto name:{"count","energy","zero_ph","nph","p0","r2"}) {
        if(!projectors && std::string(name)!="count")continue;
        for(unsigned b=0;b<bins;++b)csv<<","<<name<<"_"<<b;
    }
    csv<<"\n";
    // Three projection windows [T/4,T/2], [T/2,3T/4], [3T/4,T].
    const unsigned sites=2*radius+2;
    if(projectors) {
        prof<<"block";
        for(unsigned w=0;w<3;++w)for(auto name:{"p","z"})for(unsigned x=0;x<sites;++x)prof<<","<<name<<"_"<<w<<"_"<<x;
        prof<<"\n";
    }
    for(unsigned block=0;block<blocks;++block){
        std::vector<std::vector<double>> fields(projectors?6:1,std::vector<double>(bins,0));
        std::vector<double> profiles(projectors?3*2*sites:0,0),probability(bins),energy(projectors?bins:0);
        std::uint64_t reference=0,samples=0;
        for(std::uint64_t step=0;step<steps;++step){
            sampler.step();if((step+1)%thin)continue;
            const auto& s=sampler.state;const auto& a=sampler.cached;
            reference+=defect::vertices(s)==0;++samples;
            maximum_arcs=std::max<std::uint64_t>(maximum_arcs,s.arcs.size());maximum_hops=std::max<std::uint64_t>(maximum_hops,s.hops.size());maximum_extent=std::max<std::uint64_t>(maximum_extent,a.extent);
            defect::conditional(s,a,m,T,probability,energy);
            for(unsigned b=0;b<bins;++b){
                const double p=probability[b];fields[0][b]+=p;
                if(projectors){fields[1][b]+=energy[b];fields[2][b]+=p*(a.nph_mid==0);fields[3][b]+=p*a.nph_mid;fields[4][b]+=p*(a.midpoint==0);fields[5][b]+=p*a.midpoint*a.midpoint;}
            }
            const unsigned x=std::abs(a.midpoint)<=static_cast<int>(radius)?a.midpoint+radius:sites-1;
            for(unsigned w=0;projectors && w<3;++w){
                double p=0;for(unsigned b=(w+1)*bins/4;b<(w+2)*bins/4;++b)p+=probability[b];
                profiles[(w*2)*sites+x]+=p;
                if(a.nph_mid==0)profiles[(w*2+1)*sites+x]+=p;
            }
        }
        csv<<block<<","<<reference<<","<<samples;for(const auto& f:fields)for(double x:f)csv<<","<<x;csv<<"\n";csv.flush();
        if(projectors){prof<<block;for(double x:profiles)prof<<","<<x;prof<<"\n";prof.flush();}
    }
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();metadata(true,seconds);
    std::cout<<out.string()<<": "<<blocks*steps<<" production steps in "<<seconds<<" s\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}}
