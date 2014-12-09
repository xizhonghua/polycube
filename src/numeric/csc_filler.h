#ifndef CSC_FILLER_H
#define CSC_FILLER_H

#include <hjlib/sparse/operation.h>
#include <hjlib/sparse/sparse.h>

template <typename val_type>
class csc_filler
{
public:
  csc_filler(size_t row_num, size_t col_num):row_num_(row_num), col_num_(col_num){}
  void reset(size_t row_num, size_t col_num){
    row_num_ = row_num;
    col_num_ = col_num;
    entry.clear();
  }
  void add(size_t row, size_t col, val_type v) {
    auto it_col = entry.find(col);
    if(it_col == entry.end()) {
        entry[col][row] = v;
      }else{
        std::map<size_t,val_type> & one_col = it_col->second;
        auto it_row = one_col.find(row);
        if(it_row == one_col.end()){
            one_col[row] = v;
          }else{
            it_row->second += v;
          }
      }
  }
  void out(hj::sparse::csc<val_type,int32_t> &csc){
    size_t nnz = 0;
    for(const auto & one_entry : entry) nnz += one_entry.second.size();
    csc.resize(row_num_, col_num_, nnz);

    size_t col_idx = 0;
    for(const auto & one_entry : entry){
        if(one_entry.first != col_idx){
            assert(one_entry.first > col_idx);
            for(; col_idx < one_entry.first; ++col_idx){
                csc.ptr()[col_idx+1] = csc.ptr()[col_idx];
              }
          }
        assert(one_entry.first == col_idx);
        csc.ptr()[col_idx+1] = csc.ptr()[col_idx] + one_entry.second.size();
        size_t one_row_idx = 0;
        for(const auto & one_row_item : one_entry.second){
            csc.idx()[csc.ptr()[col_idx]+one_row_idx] = one_row_item.first;
            csc.val()[csc.ptr()[col_idx]+one_row_idx] = one_row_item.second;
            ++one_row_idx;
          }
        ++col_idx;
      }
  }
private:
  csc_filler(){}
  csc_filler & operator = (const csc_filler &){}
  size_t row_num_;
  size_t col_num_;
  std::map<size_t, std::map<size_t, val_type> > entry; // col, row, v
};

#endif // CSC_FILLER_H
