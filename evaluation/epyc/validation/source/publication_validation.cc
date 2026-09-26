// Small machine-readable validation driver. Energy traversal shares the solver's
// energy code; publication_validation.py supplies the external ViennaRNA check.
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/resource.h>
#include "codon_table.h"
#include "dfa.h"
#include "energy.h"
#include "fold_turner.h"

using namespace ldclean;
static std::string quote(const std::string& s) {
    std::string out="\"";
    for (char c:s) { if(c=='\\'||c=='\"') out+='\\'; if(c=='\n') out+="\\n"; else out+=c; }
    return out+'\"';
}
int main(int argc, char** argv) {
  try {
    if(argc<2) throw std::runtime_error("mode required: design, fixed, eval, special");
    const std::string mode=argv[1];
    std::cout<<std::setprecision(17);
    if(mode=="special") {
      std::cout<<"["; bool first=true;
      for(int n:{3,4,6}) for(const auto& x:energy::special_hairpins(n)) {
        if(!first)std::cout<<","; first=false;
        std::cout<<"{\"sequence\":"<<quote(x.first)<<",\"loop_size\":"<<n<<"}";
      }
      std::cout<<"]\n";return 0;
    }
    if(mode=="eval") {
      if(argc!=4)throw std::runtime_error("eval RNA STRUCTURE");
      std::cout<<"{\"status\":\"ok\",\"energy_centikcal\":"<<turner_eval(argv[2],argv[3])<<"}\n";return 0;
    }
    if(argc<5)throw std::runtime_error("design|fixed TABLE PROTEIN_OR_RNA LAMBDA [MOTIFS...]");
    CodonTable table(argv[2]);const std::string input=argv[3];
    const double lambda=std::stod(argv[4])*100.0;
    std::vector<std::string> motifs;for(int i=5;i<argc;++i)motifs.push_back(argv[i]);
    const auto start=std::chrono::steady_clock::now();
    std::unique_ptr<DFA> lattice;
    if(mode=="fixed")lattice.reset(new DFA(input));
    else if(mode=="design")lattice.reset(new DFA(input,table,motifs));
    else throw std::runtime_error("unknown mode");
    const DFA& dfa=*lattice;
    std::size_t edges=0,width=0;
    for(int i=0;i<dfa.num_nodes();++i)edges+=dfa.out_edges(i).size();
    for(int i=0;i<=dfa.rna_length();++i)width=std::max(width,dfa.nodes_at(i).size());
    TurnerSparseStats stats;
    TurnerDesign d=design_turner(dfa,table,lambda,3,1,&stats);
    const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    const int energy=turner_eval(d.mrna,d.structure);
    double penalty=0;std::string translated;
    if(d.mrna.size()%3==0)for(std::size_t i=0;i<d.mrna.size();i+=3){
      translated+=table.aa_of(d.mrna.substr(i,3));penalty-=std::log(table.w_of(d.mrna.substr(i,3)));
    }
    struct rusage usage;getrusage(RUSAGE_SELF,&usage);
    std::cout<<"{\"status\":\"ok\",\"rna\":"<<quote(d.mrna)
      <<",\"structure\":"<<quote(d.structure)<<",\"translated\":"<<quote(translated)
      <<",\"dp_centikcal\":"<<d.cost<<",\"mfe_kcal\":"<<d.mfe_kcal
      <<",\"energy_centikcal\":"<<energy<<",\"codon_penalty\":"<<penalty
      <<",\"cai\":"<<d.cai<<",\"objective_recomputed_centikcal\":"
      <<energy+(mode=="fixed"?0:lambda*penalty)
      <<",\"nodes\":"<<dfa.num_nodes()<<",\"edges\":"<<edges<<",\"width\":"<<width
      <<",\"candidates\":"<<stats.candidates<<",\"direct_intervals\":"<<stats.direct_intervals
      <<",\"instrumented_seconds\":"<<elapsed<<",\"peak_rss_kib\":"<<usage.ru_maxrss;
    // Small cases use both retained dense decompositions. Large cases use
    // external fixed-sequence refolding, not a cubic dense lattice replay.
    if(dfa.rna_length()<=360&&width<=16){
      std::cout<<",\"dense_reference_centikcal\":"<<turner_min_cost_dense_reference(dfa,lambda)
        <<",\"dense_wavefront_centikcal\":"<<turner_min_cost_dense_wavefront(dfa,lambda,3,1)
        <<",\"dense_right_normal_centikcal\":"<<turner_min_cost_dense_right_normal(dfa,lambda,3,1);
    }
    std::cout<<"}\n";return 0;
  } catch(const std::exception& e) {
    std::cout<<"{\"status\":\"error\",\"message\":"<<quote(e.what())<<"}\n";return 2;
  }
}
