#ifndef K_NEAR_POINTS_H
#define K_NEAR_POINTS_H

#include <zjucad/matrix/matrix.h>
#include <ANN/ANN.h>
#include <vector>

class K_near_points{
public:
  K_near_points(const zjucad::matrix::matrix<double> & node)
    :node_(node){init();}
  void query_k_near_points(
      const zjucad::matrix::matrix<double> & query_pt,
      const size_t k, std::vector<int> & kpts, std::vector<double> & dis,
      const double eps = 1e-5){
    ANNpoint query = const_cast<ANNpoint>(&query_pt[0]);
    kpts.resize(k);
    dis.resize(k);
    ANNidxArray nnidx = const_cast<ANNidxArray>(&kpts[0]);
    ANNdistArray dists = const_cast<ANNdistArray>(&dis[0]);
    kdTree_->annkSearch(query, k, nnidx, dists, eps);
  }
  virtual ~K_near_points(){
    delete []datapts_;
    delete kdTree_;
    annClose();
  }
private:
  void init(){
    datapts_ = new ANNpoint[node_.size(2)];
    for(size_t pi = 0; pi < node_.size(2); ++pi){
        *(datapts_+pi) = const_cast<ANNpoint>(&node_(0,pi));
      }
    kdTree_ = new ANNkd_tree(datapts_, node_.size(2), node_.size(1));
  }
private:
  const zjucad::matrix::matrix<double> &node_;
  ANNpointArray datapts_;
  ANNkd_tree * kdTree_;
};

#endif // K_NEAR_POINTS_H
