#ifndef __SetObj_H__
#define __SetObj_H__

#include <set>
#include <vector>
#include <random>

class SetObj{
  public:
    int num_elements;
    std::vector<int> idx_data;
    std::vector<int> set_data;

    SetObj(int num_elements) : num_elements(num_elements){
      idx_data = std::vector<int>(num_elements, -1);
      set_data = std::vector<int>();
    };

    void clear(){
      idx_data = std::vector<int>(num_elements, -1);
      set_data.clear();
    };

    void clear(int num_elements_){
      num_elements = num_elements_;
      idx_data = std::vector<int>(num_elements, -1);
      set_data.clear();
    };

    bool empty(){
      return set_data.empty();
    }

    bool is_element(int idx){
      return (idx_data[idx] > -1);
    }

    int size(){
      return set_data.size();
    }

    int getRandomElement(std::mt19937 &rng){
      if (set_data.empty()){
        return -1;
      }
      std::uniform_int_distribution<int> uni(0,set_data.size()-1);
      int idx = uni(rng);
      return set_data[idx];
    };

    void addElement(const int idx){
      if (idx_data[idx] < 0){
        idx_data[idx] = set_data.size();
        set_data.push_back(idx);
      }
    };

    void removeElement(const int idx){
      if (idx_data[idx] > -1){
        // replace current element with back element
        set_data[idx_data[idx]] = set_data.back();
        // change idx_data for back element
        idx_data[set_data.back()] = idx_data[idx];
        // remove from idx data
        idx_data[idx] = -1;
        // remove back element
        set_data.pop_back();
      }
    };

    void makeUnion(const SetObj &src){
      for (int i = 0; i < src.set_data.size(); ++i){
        addElement(src.set_data[i]);
      }
    };

    void setMinus(const SetObj &src){
      for (int i = 0; i < src.set_data.size(); ++i){
        removeElement(src.set_data[i]);
      }
    };

    SetObj getIntersection(const SetObj &src){
      SetObj temp(num_elements);

      if (src.size() > set_data.size()){
        for (int i = 0; i < src.set_data.size(); ++i){
          if (is_element(src.set_data[i])){
            temp.addElement(src.set_data[i]);
          }
        }
      }
      else{
        for (int i = 0; i < set_data.size(); ++i){
          if (src.is_element(set_data[i])){
            temp.addElement(set_data[i]);
          }
        }
      }
      return temp;
    };

    ~SetObj(){};
};
#endif
