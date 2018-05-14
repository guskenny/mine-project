#ifndef __SetObj_H__
#define __SetObj_H__

#include <set>
#include <vector>
#include <random>

class SetObj{
  public:
    int num_elements;
    std::vector<bool> bool_data;
    std::set<int> set_data;

    SetObj(std::vector<bool> bool_input, std::set<int> set_input) :
        bool_data(bool_input), set_data(set_input),
        num_elements(bool_input.size()) {};

    SetObj(std::vector<int> int_input, int num_elements) :
        num_elements(num_elements){
      bool_data = std::vector<bool>(num_elements, false);
      for (int i = 0; i < int_input.size(); ++i){
        addElement(int_input[i]);
      }
    };

    SetObj(std::set<int> int_input, int num_elements) :
        num_elements(num_elements) {
      bool_data = std::vector<bool>(num_elements, false);
      set_data = int_input;
      for (auto it = int_input.begin(); it != int_input.end(); ++it){
          bool_data[*it] = true;
      }
    };

    SetObj(int num_elements) : num_elements(num_elements){
      bool_data = std::vector<bool>(num_elements, false);
      set_data = std::set<int>();
    };

    void clear(){
      bool_data = std::vector<bool>(num_elements, false);
      set_data.clear();
    };

    void clear(int num_elements_){
      num_elements = num_elements_;
      bool_data = std::vector<bool>(num_elements, false);
      set_data.clear();
    };

    bool is_element(int idx){
      return bool_data[idx];
    }

    int get_set_size(){
      return set_data.size();
    }

    const int getRandomElement(std::mt19937 rng){
      if (set_data.empty()){
        return -1;
      }
      std::uniform_int_distribution<int> uni(0,set_data.size()-1);
      int idx = uni(rng);
      auto it = set_data.begin();
      std::advance(it, idx);
      return *it;
    };

    void addElement(const int idx){
      if (!bool_data[idx]){
        bool_data[idx] = true;
        set_data.insert(idx);
      }
    };

    void removeElement(const int idx){
      if (bool_data[idx]){
        bool_data[idx] = false;
        set_data.erase(idx);
      }
    };

    void getUnion(const SetObj &src){
      for (auto it=src.set_data.begin(); it!=src.set_data.end(); ++it){
        bool_data[*it] = true;
        set_data.insert(*it);
      }
    };

    void getDiff(const SetObj &src){
      for (auto it=src.set_data.begin(); it!=src.set_data.end(); ++it){
        bool_data[*it] = false;
        set_data.erase(*it);
      }
    };


    ~SetObj(){};
};
#endif
