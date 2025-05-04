
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
DataFrame hit_to_lead(
  const int n_compounds,
  DataFrame tcs,
  DataFrame model) {

  // requires tcs to be sorted by ref_cid
  IntegerVector tcs_ref_cids = tcs["ref_cid"];
  IntegerVector tcs_query_cids = tcs["query_cid"];

  // each row in model represents a hit at target <tid>
  // for which there are <n_leads> additional leads
  IntegerVector model_tids = model["tid"];
  IntegerVector model_hit_cids = model["hit_cid"];
  IntegerVector model_n_leads = model["n_leads"];


  // return data.frame with columns (tid, hit_cid, lead_cid)
  // for each tid, hit_cids are sampled without replacement from 1:n_compounds
  // for each (tid, hit_cid) let candidate_hits = tcs[tcs$ref_cid == hit_cid,"query_cid"]
  // if n_leads < len(candidate_leads) then sample n_leads from candidate_hits
  // otherwise take all candidate_leads


  // Warning, R uses 1-based indexing but Rcpp uses 0-based indexing


  // // 1) for each target pick hits
  // IntegerVector model_hit_cids(model_tids.size());
  // {
  //   int first_hit_for_tid = 0;
  //   int cur_model_tid = 1;
  //   int n_hits = 0;
  //   for(int ii=0; ii < model_tids.size(); ii++){
  //     std::cout
  //       << "ii: " << ii
  //       << " first_hit_for_tid: " << first_hit_for_tid
  //       << " cur_model_tid: " << cur_model_tid
  //       << " n_hits: " << n_hits;
  //
  //     if(cur_model_tid == model_tids[ii] && ii+1 <= model_tids.size()){
  //       std::cout << " n_hits++" << std::endl;
  //       n_hits++;
  //     } else {
  //       arma::uvec hit_cids(n_hits);
  //       RcppArmadillo::SampleNoReplace(
  //         hit_cids,
  //         n_compounds,
  //         n_hits);
  //       std::cout << std::endl;
  //       for(int jj = 0; jj < n_hits; jj++){
  //         // SampleNoReplace samples in 0:(n_compounds-1) so add 1 for R 1-based indexing
  //         int hit_cid = hit_cids[jj] + 1;
  //         model_hit_cids[first_hit_for_tid + jj] = hit_cid;
  //         std::cout << "    " << "model_hit_cids[" << first_hit_for_tid << " + " << jj << "] = " << hit_cid << std::endl;
  //
  //       }
  //
  //       first_hit_for_tid = ii;
  //       n_hits = 1;
  //       cur_model_tid = model_tids[ii];
  //     }
  //   }
  // }


  std::cout << "Preprocess tcs to able to rapidly pick leads for a given hit ..." << std::endl;
  std::vector<int> lead_start(n_compounds);
  {
    // lead_start indexes into first query_cid (if any) for ref_cid
    int ii = -1;  // index into tcs_lookup
    for(int jj=0; jj < tcs_ref_cids.size(); jj++){
      while(ii+1 < tcs_ref_cids[jj]){
        ii++;
        lead_start[ii] = jj;
	if( ii % 100000 == 0){
	  std::cout << "  ii:" << ii << " jj:" << jj << " tcs_ref_cids[jj]:" << tcs_ref_cids[jj] << std::endl;
	}
      }
    }
  }

  std::cout << "For each hit pick the leads ... " << std::endl;
  IntegerVector m_tids;
  IntegerVector m_hit_cids;
  IntegerVector m_lead_cids;
  for(int ii = 0; ii < model_hit_cids.size(); ii++){
    if( ii % 500 == 0 ) {
      std::cout << "    ii= " << ii << " of " << model_hit_cids.size() << std::endl;
    }
    int const tid = model_tids[ii];
    int const hit_cid = model_hit_cids[ii];
    int const n_leads = model_n_leads[ii];

    int start = lead_start[hit_cid-1];
    int n_candidate_leads;
    if(hit_cid == (n_compounds)){
      n_candidate_leads = tcs_query_cids.size() - start;
    } else {
      n_candidate_leads = lead_start[(hit_cid-1)+1] - start;
    }

    // std::cout
    //   << "ii: " << ii
    //   << " tid: " << tid
    //   << " hit_cid: " << hit_cid
    //   << " n_leads: " << n_leads
    //   << " start: " << start
    //   << " n_candidate_leads: " << n_candidate_leads
    //   << std::endl;

    if(n_candidate_leads < n_leads){
      for(int jj = 0; jj < n_candidate_leads; jj++){
        m_tids.push_back(tid);
        m_hit_cids.push_back(hit_cid);
        m_lead_cids.push_back(tcs_query_cids[start + jj]);
      }
    } else {
      arma::uvec which_candidate_leads(n_leads);
      RcppArmadillo::SampleNoReplace(
        which_candidate_leads,
        n_candidate_leads,
        n_leads);

      for(int jj = 0; jj < n_leads; jj++){
        m_tids.push_back(tid);
        m_hit_cids.push_back(hit_cid);
        m_lead_cids.push_back(tcs_query_cids[start + which_candidate_leads[jj]]);
      }
    }
  }

  // std::cout << "m_tids.size():" << m_tids.size() << std::endl;
  // std::cout << "m_hit_cids.size():" << m_hit_cids.size() << std::endl;
  // std::cout << "m_lead_cids.size():" << m_lead_cids.size() << std::endl;

  return DataFrame::create(
    _["tid"]=m_tids,
    _["hit_cid"]=m_hit_cids,
    _["lead_cid"]=m_lead_cids);
}


