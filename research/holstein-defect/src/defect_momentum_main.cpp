#include "defect_momentum.hpp"
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
    const std::vector<std::string> allowed={"t","omega","g","U","mu","window","tau-max","max-arcs","max-hops","bins","blocks","steps-per-block","warmup","thin","seed","out"};
    for(const auto& p:opt)if(std::find(allowed.begin(),allowed.end(),p.first)==allowed.end())throw std::runtime_error("Unknown argument: "+p.first);
    auto number=[&](std::string k,double d){return opt.count(k)?std::stod(opt.at(k)):d;};
    auto integer=[&](std::string k,std::uint64_t d){
        if(opt.count(k)&&opt.at(k).front()=='-')throw std::runtime_error("Negative integer option");
        return opt.count(k)?std::stoull(opt.at(k)):d;
    };
    defect::Model m{number("t",1),number("omega",1),number("g",std::sqrt(.5)),number("U",1),number("mu",-3.),0};
    const double T=number("tau-max",12);
    const auto W=integer("window",16),bins=integer("bins",192),blocks=integer("blocks",320),
      steps=integer("steps-per-block",100000),thin=integer("thin",100),warmup=integer("warmup",500000),
      seed=integer("seed",109001),cap=integer("max-arcs",256),hopcap=integer("max-hops",1024);
    if(m.t<0||m.omega<=0||m.g<0||m.U<0||m.mu>=-m.U||T<=0||bins<8||bins%4||bins>1024||blocks<2||!steps||!thin||steps%thin||W>256||cap>10000||hopcap>10000||!opt.count("out"))throw std::runtime_error("Invalid parameters");
    for(double x:{m.t,m.omega,m.g,m.U,m.mu,T})if(!std::isfinite(x))throw std::runtime_error("Nonfinite parameter");
    const std::filesystem::path out=opt.at("out");
    if(std::filesystem::exists(out))throw std::runtime_error("Refusing to overwrite output");
    std::filesystem::create_directories(out);
    defect_momentum::Sampler sampler(m,T,W,cap,hopcap,seed);
    std::uint64_t maximum_arcs=0,maximum_hops=0,maximum_extent=0,nonlocal_samples=0,odd_displacement_samples=0;
    auto metadata=[&](bool complete,double seconds){
        std::ofstream f(out/"run.json");f<<std::setprecision(17);
        f<<"{\n\"complete\":"<<(complete?"true":"false")<<",\n\"method\":\"real-space bare diagrams with open electron endpoints\",\n"
         <<"\"t\":"<<m.t<<",\"omega\":"<<m.omega<<",\"g\":"<<m.g<<",\"U\":"<<m.U<<",\"mu\":"<<m.mu<<",\"window\":"<<W<<",\n"
         <<"\"tau_max\":"<<T<<",\"bins\":"<<bins<<",\"blocks\":"<<blocks<<",\"steps_per_block\":"<<steps<<",\"thin\":"<<thin<<",\"warmup\":"<<warmup<<",\"seed\":"<<seed<<",\n"
         <<"\"max_arcs\":"<<cap<<",\"max_hops\":"<<hopcap<<",\"arc_cap_attempts\":"<<sampler.arc_cap_attempts<<",\"hop_cap_attempts\":"<<sampler.hop_cap_attempts<<",\n"
         <<"\"invalid_proposals\":"<<sampler.invalid_proposals<<",\"maximum_arcs\":"<<maximum_arcs<<",\"maximum_hops\":"<<maximum_hops<<",\"maximum_extent\":"<<maximum_extent<<",\n"
         <<"\"nonlocal_samples\":"<<nonlocal_samples<<",\"odd_displacement_samples\":"<<odd_displacement_samples<<",\"seconds\":"<<seconds<<",\n"
         <<"\"binary_layout\":\"native-endian float64: blocks x [reference, samples, counts(abs(endpoint-origin), tau_bin) in C order]\",\n"
         <<"\"rows\":"<<blocks<<",\"columns\":"<<2+(2*W+1)*bins<<",\"byte_order\":\"little\",\n"
         <<"\"spatial_boundary\":\"none; only initial and final states are windowed\",\n\"attempted\":[";
        for(unsigned j=0;j<10;++j)f<<(j?",":"")<<sampler.attempted[j];
        f<<"],\"accepted\":[";for(unsigned j=0;j<10;++j)f<<(j?",":"")<<sampler.accepted[j];f<<"]\n}\n";
    };
    const std::uint16_t byte_test=1;
    if(*reinterpret_cast<const unsigned char*>(&byte_test)!=1 || sizeof(double)!=8)throw std::runtime_error("Binary format needs little-endian float64");
    metadata(false,0);const auto start=std::chrono::steady_clock::now();
    for(std::uint64_t i=0;i<warmup;++i)sampler.step();
    sampler.reset();
    std::ofstream binary(out/"blocks.f64",std::ios::binary);
    std::vector<double> fields(2+(2*W+1)*bins),probability(bins),unused;
    for(unsigned block=0;block<blocks;++block){
        std::fill(fields.begin(),fields.end(),0.);
        for(std::uint64_t step=0;step<steps;++step){
            sampler.step();if((step+1)%thin)continue;
            const auto& s=sampler.state;const auto& a=sampler.cached;
            fields[0]+=defect::vertices(s)==0;++fields[1];
            unsigned d=std::abs(a.endpoint-sampler.model.origin);
            if(d>2*W)throw std::runtime_error("Endpoint escaped window");
            nonlocal_samples+=d>0;odd_displacement_samples+=d%2;
            maximum_arcs=std::max<std::uint64_t>(maximum_arcs,s.arcs.size());
            maximum_hops=std::max<std::uint64_t>(maximum_hops,s.hops.size());
            maximum_extent=std::max<std::uint64_t>(maximum_extent,a.extent);
            defect::conditional(s,a,sampler.model,T,probability,unused);
            auto dest=fields.begin()+2+d*bins;
            for(unsigned b=0;b<bins;++b)dest[b]+=probability[b];
        }
        binary.write(reinterpret_cast<const char*>(fields.data()),fields.size()*sizeof(double));
        binary.flush();if(!binary)throw std::runtime_error("Failed writing blocks");
    }
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();metadata(true,seconds);
    std::cout<<out.string()<<": "<<blocks*steps<<" production steps in "<<seconds<<" s\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}}
