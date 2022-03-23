// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// we only include RcppEigen.h which pulls Rcpp.h in for us

#include "cpp11.hpp"
#include "chain_file.h"






[[cpp11::register]]
cpp11::data_frame lift_over(std::string path,cpp11::strings chroms, cpp11::integers starts){

    const size_t num_pos = chroms.size();
    cpp11::writable::integers from_id_vec;
    from_id_vec.reserve(num_pos);
    cpp11::writable::strings to_chrom_vec;
    cpp11::writable::integers to_pos_vec;
    cpp11::writable::logicals to_strand_vec;


    to_chrom_vec.reserve(num_pos);
    to_pos_vec.reserve(num_pos);
    to_strand_vec.reserve(num_pos);

    std::map<std::string, liftover::Target> targets = liftover::open_chainfile(path);

    // search for a given coordinate
    std::string chrom;
    for(int i=0; i<num_pos; i++){
        chrom = chroms[i];
        int pos = starts[i];
        for (const auto x : targets[chrom][pos]) {
            from_id_vec.push_back(i);
            to_chrom_vec.push_back(x.contig);
            to_pos_vec.push_back(x.pos);
            to_strand_vec.push_back(x.fwd_strand);
        }
    }
    using namespace cpp11::literals;
    return cpp11::writable::data_frame({
            "index"_nm = from_id_vec,
            "to_chrom"_nm = to_chrom_vec,
            "to_pos"_nm = to_pos_vec,
            "to_strand"_nm = to_strand_vec});
}
