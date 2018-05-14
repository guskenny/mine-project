#ifndef __MY_MATH_H__
#define __MY_MATH_H__

#include <math.h>
#include <vector>

namespace my_math{
  template <typename T>
  double calc_mean(const std::vector<T> &vals){
    double sum = 0;
    for (int i = 0; i < vals.size(); ++i){
      sum += vals[i];
    }
    return sum/(double)vals.size();
  }

  template <typename T>
  double calc_std_dev(const std::vector<T> &vals){
    double mean = calc_mean(vals);
    std::vector<double> sqr_diff(vals.size());
    for (int i = 0; i < vals.size(); ++i){
      sqr_diff[i] = pow(vals[i]-mean,2);
    }
    double mean_sqr_diff = calc_mean(sqr_diff);
    return sqrt(mean_sqr_diff);
  }
}

#endif
