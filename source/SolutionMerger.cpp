#include "SolutionMerger.h"

void SolutionMerger::fullMerge(const std::vector<Sol_Int>&sols,
  const std::vector<int> &include, std::vector<int> &fixed, std::vector<std::vector<int> > &groups, std::vector<int> &group_map){

    std::cout << "merging " << sols.size() << " CPIT solutions using fullMerge" << std::endl;

    qol::CpuTimer merge_timer;

    int nB = sols[0].x.size();
    int nT = sols[0].nT;

    fixed = std::vector<int>(nB*nT, 0);

    std::vector<int> var_cnt(3,0);

    std::vector<std::array<int, 2> > intervals(nB, std::array<int,2>{{nT,0}});

    for (size_t b=0; b < nB; ++b){
      if (!include[b]){
        intervals[b][1] = nT;
        continue;
      }
      for (size_t sol=0; sol < sols.size(); ++sol){
        if (sols[sol].x[b] < intervals[b][0]){ // fixed to 0
          intervals[b][0] = sols[sol].x[b];
        }
        if (sols[sol].x[b] > intervals[b][1]){ // fixed to 1
          intervals[b][1] = sols[sol].x[b];
        }
      }
    }

    for (int b = 0; b<nB; ++b){
      for (int t = 0; t < nT; ++t){
        if (t >= intervals[b][0]){
          fixed[t*nB + b] = -1;
          if (t > intervals[b][1]){
            fixed[t*nB + b] = 1;
          }
        }
      }
    }

    std::cout << "initialising group data structures... " << std::flush;

    groups = std::vector<std::vector<int> >();
    group_map = std::vector<int>(nB*nT,0);

    // establish vector for storing groups
    groups.push_back(std::vector<int>());
    groups[0].reserve(nB*nT);
    for(int j = 0; j < nB*nT; ++j){
      groups[0].push_back(j);
    }

    std::cout << "done!" << std::endl;
    std::cout << "getting groups... " << std::flush;

    // iterate over all solutions
    for (int sol = 0; sol < sols.size(); ++sol){
      // in case groups are added, have fixed iteration stopping point
      int groups_size = groups.size();
      // iterate over all groups
      for (int group = 0; group < groups_size; ++group){
        std::vector<int> temp_group;
        // value for group
        int group_val = (sols[sol].x[groups[group][0]%nB] <= groups[group][0]/nB);
        // iterate over every element in the group to check if same still
        for (int i = 1; i < groups[group].size(); ++i){
          // groups[group][i]] is current variable
          // if sol is less than t, sol = 1
          if (group_val != (sols[sol].x[groups[group][i]%nB] <= groups[group][i]/nB)) { // test if same value
            // change value in map to new group
            group_map[groups[group][i]] = groups.size();
            // push current variable index to new group
            temp_group.push_back(groups[group][i]);
            // remove current variable from current group
            // by overwriting variable with last variable in group
            groups[group][i] = groups[group].back();
            // remove copy of last variable in group
            groups[group].pop_back();
            i--; // decrement i to test newly inserted variable
          }
        }
        // if new group created, add to list of groups, but no need to
        // test it
        if (!temp_group.empty()){
          groups.push_back(temp_group);
        }
      }
    }

    std::cout << "done!" << std::endl;
    //
    // std::cout << "sorting groups... " << std::flush;
    // vector to keep track of what has been added already
    // std::vector<bool> already_done(nB, false);


    // std::cout << "testing groups.. " << std::flush;
    //
    // for (int sol = 0; sol < sols.size(); ++sol){
    //   for (int group = 0; group < groups.size(); ++group){
    //     int group_val = (groups[group][0]/nB >= sols[sol].x[groups[group][0]%nB]);
    //     for (int member_idx = 1; member_idx < groups[group].size(); ++member_idx){
    //       int member = groups[group][member_idx];
    //       if (group_map[member] != group){
    //         std::cout << "\nmember " << member << " in group " << group << " not " << group_map[member] << "!" <<std::flush;
    //       }
    //       if ((groups[group][member_idx]/nB >= sols[sol].x[groups[group][member_idx]%nB]) != group_val){
    //         std::cout << "\nmember " << member_idx << " in group " << group << " different value in sol " << sol << std::endl;
    //       }
    //     }
    //   }
    // }
    //
    // // for (int var = 0; var < nB * nT; ++var){
    // //   int found = -1;
    // //   for (int group = 0; group < groups.size(); ++group){
    // //     for (int member = 0; member < groups[group].size(); ++member){
    // //       if (groups[group][member] == var){
    // //         if (found < 0){
    // //           found = group;
    // //         }
    // //         else{
    // //           std::cout << "ERROR: variable " << var << " found in groups " << found << " and " << group << std::endl;
    // //         }
    // //       }
    // //     }
    // //   }
    // //   if (found < 0){
    // //     std::cout << "ERROR: variable " << var << " not found!" << std::endl;
    // //   }
    // // }
    //
    // std::cout << "done!" << std::endl;

    std::cout << "Finding fixed groups..." << std::flush;
    std::vector<int> fixed_groups;
    for (int group = 0; group < groups.size(); ++group){
      bool fixed_group = true;
      int base_val = (sols[0].x[groups[group][0]%nB] <= groups[group][0]/nB);
      for (int sol = 1; sol < sols.size(); ++sol){
        if((sols[sol].x[groups[group][0]%nB] <= groups[group][0]/nB) != base_val){
          // not fixed
          fixed_group = false;
          break;
        }
      }
      if (fixed_group){
        fixed_groups.push_back(group);
      }
    }
    std::cout << " done!" << std::endl;

    std::cout << "\nFixed groups: "<< std::endl;
    for (int i = 0; i < fixed_groups.size(); ++i){
      std::cout << "group " << fixed_groups[i] << " - value: " << (sols[0].x[groups[fixed_groups[i]][0]%nB] <= groups[fixed_groups[i]][0]/nB) << ", size: " << groups[fixed_groups[i]].size() << std::endl;
    }
    std::cout << std::endl;


    int groups_count = 0;
    // int total_count = 0;

    for (int group = 0; group < groups.size(); ++group){
      if (groups[group].size() > 1){
        groups_count++;
      }
    }

    int zero_count = 0;

    for (int i = 0; i < nB*nT; ++i){
      if (group_map[i] == 0){
        zero_count++;
      }
    }


    std::cout << "groups[0].size(): " << groups[0].size() << ", zero_count: " << zero_count << std::endl;

    // std::vector<int> sizes(groups.size(), 0);
    //
    //
    // int max_group_size = 0;
    // int second_max_size = 0;
    // int third_max_size = 0;
    //
    // for (int group = 0; group < groups.size(); ++group){
    //
    //   max_group_size = std::max(max_group_size, (int)groups[group].size());
    //   if ((int)groups[group].size() < max_group_size){
    //     second_max_size = std::max(second_max_size, (int)groups[group].size());
    //   }
    //   if ((int)groups[group].size() < second_max_size){
    //     third_max_size = std::max(third_max_size, (int)groups[group].size());
    //   }
    // }

    int total_vars = 0;


    std::priority_queue<std::pair<double, int>> q;
    for (int i = 0; i < groups.size(); ++i) {
      total_vars += groups[i].size();
      q.push(std::pair<double, int>(groups[i].size(), i));
    }
    int k = 5; // number of indices we need
    std::cout << "\nTop " << k << " group sizes:\n------------------" << std::endl;
    for (int i = 0; i < k; ++i) {
      int ki = q.top().second;
      double ks = q.top().first;
      std::cout << i+1 << ": group[" << ki << "] = " << ks << std::endl;
      q.pop();
    }

    std::cout << "\ntotal_vars: " << total_vars << ", group_map.size(): " << group_map.size() << std::endl;
    //
    // std::cout << "max_group_size: " << max_group_size << std::endl;
    // std::cout << "second_max_size: " << second_max_size << std::endl;
    // std::cout << "third_max_size: " << third_max_size << std::endl;

    // std::cout << "testing fixed.. " << std::flush;
    //
    // int zeroes = -1;
    // int ones = -1;
    //
    // for (int i = 0; i < fixed.size(); ++i){
    //   if (fixed[i] == 1 && ones < 0){
    //     ones = group_map[i];
    //   }
    //   if (fixed[i] == 0 && zeroes < 0){
    //     zeroes = group_map[i];
    //   }
    //   if (fixed[i] == 1 && group_map[i] != ones){
    //     std::cout << "multiple groups fixed to one: " << group_map[i] << std::endl;
    //     // break;
    //   }
    //   if (fixed[i] == 0 && group_map[i] != zeroes){
    //     std::cout << "multiple groups fixed to zero" << group_map[i] << std::endl;
    //     // break;
    //   }
    // }
    // std::cout << "done! group " << zeroes << " fixed to zero and group " << ones << " fixed to one" << std::endl;

    // std::cout << "Testing groups to see if blocks exist in different periods and groups.." << std::endl;
    //
    // for (int b = 0; b < nB; ++b){
    //   if (groups[group_map[b]].size() == 1){
    //     continue;
    //   }
    //   for (int t = 1; t < nT; ++t){
    //     if (group_map[b] != group_map[b+(nB*t)] && groups[group_map[b]].size() > 1){
    //       std::cout << "block " << b << " in group " << group_map[b] << " and " << group_map[b+(nB*t)] << "!" << std::endl;
    //       break;
    //     }
    //   }
    // }
    //
    // std::cout << "done!" << std::endl;

    std::cout << "\nMerge time: " << merge_timer.elapsedSeconds() << std::endl << std::endl;

    std::cout << "merge finished!" << std::endl << nB*nT << " variables found in " << groups.size() << " groups and " << groups_count << " groups with size greater than 1"<<std::endl;

    std::cout << "reduction in problem size: " << (double)groups.size()/(nB*nT) << std::endl;
    // std::cin.get();

    std::vector<int> var_cnt2(3,0);

    for (int i = 0; i < fixed.size(); ++i){
      var_cnt2[fixed[i]+1]++;
    }

    std::cout << "validation:" << std::endl;
    std::cout << "variables fixed to 0: " << var_cnt2[1] <<std::endl;
    std::cout << "variables fixed to 1: " << var_cnt2[2] <<std::endl;
    std::cout << "variables unfixed: " << var_cnt2[0] <<std::endl << std::endl;
}

void SolutionMerger::simpleMerge(const std::vector<Sol_Int>&sols,
  const std::vector<int> &include, std::vector<int> &fixed){

  std::cout << "merging " << sols.size() << " CPIT solutions using simpleMerger" << std::endl;

  qol::CpuTimer merge_timer;

  int nB = sols[0].x.size();
  int nT = sols[0].nT;

  fixed = std::vector<int>(nB*nT, 0);

  std::vector<int> var_cnt(3,0);

  std::vector<std::array<int, 2> > intervals(nB, std::array<int,2>{{nT,0}});

  for (size_t b=0; b < nB; ++b){
    if (!include[b]){
      intervals[b][1] = nT;
      continue;
    }
    for (size_t sol=0; sol < sols.size(); ++sol){
      if (sols[sol].x[b] < intervals[b][0]){ // fixed to 0
        intervals[b][0] = sols[sol].x[b];
      }
      if (sols[sol].x[b] > intervals[b][1]){ // fixed to 1
        intervals[b][1] = sols[sol].x[b];
      }
    }
  }

  for (int b = 0; b<nB; ++b){
    for (int t = 0; t < nT; ++t){
      if (t >= intervals[b][0]){
        fixed[t*nB + b] = -1;
        if (t > intervals[b][1]){
          fixed[t*nB + b] = 1;
        }
      }
    }
  }

  //   for (size_t t = 0; t < nT; ++t){
  //     for (size_t b = 0; b < nB; ++b){
  //     // std::cout << " " << var_cnt[0] << " " << var_cnt[1] << " " << var_cnt[2] << std::endl;
  //     int first_sol = (sols[0].x[b] <= t);
  //     bool unfixed = false;
  //     // std::cout << (sols[0].x[b] < t) << " ";
  //     for (size_t sol = 1; sol < sols.size(); ++sol){
  //       // std::cout << (sols[sol].x[b] < t) << " ";
  //       if (first_sol == (sols[sol].x[b] <= t)){
  //         continue;
  //       }
  //       else{
  //         var_cnt.back()++;
  //         unfixed = true;
  //         // std::cout << "unfixed";
  //         break;
  //       }
  //     }
  //     if (!unfixed){
  //       // if it makes it through the loop then it is fixed
  //       fixed[t*nB + b] = first_sol;
  //       // std::cout << "fixed to " << first_sol;
  //       //update variable counts
  //       var_cnt[first_sol]++;
  //     }
  //   }
  // }

  std::cout << "\nMerge time: " << merge_timer.elapsedSeconds() << std::endl << std::endl;

  std::cout << "merge finished!" << std::endl;

  std::vector<int> var_cnt2(3,0);

  for (int i = 0; i < fixed.size(); ++i){
    var_cnt2[fixed[i]+1]++;
  }

  std::cout << "validation:" << std::endl;
  std::cout << "variables fixed to 0: " << var_cnt2[1] <<std::endl;
  std::cout << "variables fixed to 1: " << var_cnt2[2] <<std::endl;
  std::cout << "variables unfixed: " << var_cnt2[0] <<std::endl << std::endl;

  for (int i = 0; i < var_cnt.size()-1; ++i){
    std::cout << var_cnt[i] << " variables fixed to " << i << std::endl;
  }
  std::cout << var_cnt.back() << " variables unfixed" << std::endl;
}

void SolutionMerger::mergeCPIT(const std::vector<Sol_Int>&sols,
  std::vector<int> &fixed, std::vector<std::vector<int> > &groups, std::vector<int> &group_map){
  std::cout << "merging " << sols.size() << " CPIT solutions" << std::endl;

  qol::CpuTimer merge_timer;

  int nB = sols[0].x.size();
  int nT = sols[0].nT;

  // std::cout << "getting fixed... " << std::flush;

  fixed = std::vector<int>(nB*nT, 0);

  // find all variables that are always 0 or always 1
  std::vector<std::array<int, 2> > intervals(nB, std::array<int,2>{{nT,0}});

  for (size_t b=0; b < nB; ++b){
    for (size_t sol=0; sol < sols.size(); ++sol){
      if (sols[sol].x[b] < intervals[b][0]){ // fixed to 0
        intervals[b][0] = sols[sol].x[b];
      }
      if (sols[sol].x[b] > intervals[b][1]){ // fixed to 1
        intervals[b][1] = sols[sol].x[b];
      }
    }
  }

  for (int b = 0; b<nB; ++b){
    for (int t = 0; t < nT; ++t){
      if (t >= intervals[b][0]){
        fixed[t*nB + b] = -1;
        if (t > intervals[b][1]){
          fixed[t*nB + b] = 1;
        }
      }
    }
  }

  // std::cout << "done!" << std::endl;

  // std::cout << "initialising group data structures... " << std::flush;

  groups = std::vector<std::vector<int> >();
  group_map = std::vector<int>(nB*nT, 0);

  // establish vector for storing groups
  groups.push_back(std::vector<int>());
  groups[0].reserve(nB*nT);
  for(int j = 0; j < nB*nT; ++j){
    groups[0].push_back(j);
  }

  // std::cout << "done!" << std::endl;
  // std::cout << "getting groups... " << std::flush;

  // iterate over all solutions
  for (int sol = 0; sol < sols.size(); ++sol){
    // in case groups are added, have fixed iteration stopping point
    int groups_size = groups.size();
    // iterate over all groups
    for (int group = 0; group < groups_size; ++group){
      std::vector<int> temp_group;
      // iterate over every element in the group to check if same still
      for (int i = 1; i < groups[group].size(); ++i){
        // if sol is less than t, sol = 1
        if ((sols[sol].x[groups[group][0]%nB] <= i/nB) != (sols[sol].x[groups[group][i]%nB] <= i/nB)) {
          // change value in map
          group_map[groups[group][i]] = groups.size();
          // push current variable index to new group
          temp_group.push_back(groups[group][i]);
          // remove current variable from current group
          groups[group][i] = groups[group].back();
          groups[group].pop_back();
          i--; // decrement i to test newly inserted variable
        }
      }
      if (!temp_group.empty()){
        groups.push_back(temp_group);
      }
    }
  }

  // std::cout << "done!" << std::endl;
  //
  // std::cout << "sorting groups... " << std::flush;
  // vector to keep track of what has been added already
  // std::vector<bool> already_done(nB, false);

  int groups_count = 0;
  // int total_count = 0;

  for (int group = 0; group < groups.size(); ++group){
    if (groups[group].size() > 1){
      groups_count++;
    }
  }

  int zero_count = 0;

  for (int i = 0; i < nB*nT; ++i){
    if (group_map[i] == 0){
      zero_count++;
    }
  }


  std::cout << "groups[0].size(): " << groups[0].size() << ", zero_count: " << zero_count << std::endl;

  int total_vars = 0;
  int max_group_size = 0;

  for (int group = 0; group < groups.size(); ++group){
    total_vars += groups[group].size();
    max_group_size = std::max(max_group_size, (int)groups[group].size());
  }

  std::cout << "total_vars: " << total_vars << ", group_map.size(): " << group_map.size() << std::endl;

  std::cout << "max_group_size: " << max_group_size << std::endl;

  std::cout << "Testing groups to see if blocks exist in different periods and groups.." << std::endl;

  // for (int b = 0; b < nB; ++b){
  //   if (groups[group_map[b]].size() == 1){
  //     continue;
  //   }
  //   for (int t = 1; t < nT; ++t){
  //     if (group_map[b] != group_map[b+(nB*t)] && groups[group_map[b]].size() > 1){
  //       std::cout << "block " << b << " in group " << group_map[b] << " and " << group_map[b+(nB*t)] << "!" << std::endl;
  //       break;
  //     }
  //   }
  // }

  // std::cout << "done!" << std::endl;

  std::cout << "\nMerge time: " << merge_timer.elapsedSeconds() << std::endl << std::endl;

  std::cout << "merge finished!" << std::endl << nB*nT << " variables found in " << groups.size() << " groups and " << groups_count << " groups with size greater than 1"<<std::endl;

  std::cout << "reduction in problem size: " << (double)groups.size()/(nB*nT) << std::endl;
  // std::cin.get();
}

void SolutionMerger::mergePCPSP(const std::vector<Sol_Int>&sols, BranchNode_info &merged){
  std::cout << "merging PCPSP solutions" << std::endl;

  int nB = sols[0].x.size();
  int nT = sols[0].nT;

  std::vector<std::array<int, 2> > intervals(nB, std::array<int,2>{{nT,0}});

  // for (size_t sol = 0; sol < sols.size(); ++sol){
  //   for(size_t i=0;i<nB;++i){
  //     std::cout << sols[sol].x[i] << " ";
  //   }
  //   std::cout << std::endl;
  // }

  for (size_t b=0; b < nB; ++b){
    for (size_t sol=0; sol < sols.size(); ++sol){
      if (sols[sol].x[b]-1 < intervals[b][0]){ // fixed to 0
        intervals[b][0] = sols[sol].x[b]-1;
      }
      if (sols[sol].x[b] > intervals[b][1]){ // fixed to 1
        intervals[b][1] = sols[sol].x[b];
      }
    }
  }

  // std::cout << std::endl << "merged solutions:" << std::endl;
  //
  // for (size_t i = 0; i<intervals[0].size();++i){
  //   for (size_t b=0; b<nB; ++b){
  //     std::cout << intervals[b][i] << " ";
  //   }
  //   std::cout << std::endl;
  // }

  merged.init(nB,nT,sols[0].y[0].size(),sols[0].y[0].size());
  merged.fixed = intervals;
}

// void SolutionMerger::buildNewProblem(std::vector<std::vector<int> > &groups, std::vector<int> &group_map, const SinglePModel *old_model, SinglePModel *new_model){
//
//   std::vector<Block> * blocks=model->getBlock();
//   std::vector<Block> blocks;
//
// }

void SolutionMerger::writePCPSPmerged2D(const std::vector<Sol_Int> &sols){
  int nSols = sols.size();
  int nB = sols[0].x.size();
  int nT = sols[0].nT;

  std::ofstream outfile;
  outfile.open ("merged_sols2D.xsol");
  for (size_t t=0; t < nT; ++t){
    for (size_t b=0; b < nB; ++b){
      double total = 0.0;
      for (size_t sol=0; sol < nSols; ++sol){
        if (sols[sol].x[b] <= t){
          total++;
        }
      }
      outfile << total/nSols << std::endl;
    }
  }
  outfile << "END" << std::endl;
  outfile.close();
}

void SolutionMerger::writePCPSPmerged3D(const std::vector<Sol_Int> &sols){
  int nSols = sols.size();
  int nB = sols[0].x.size();
  int nT = sols[0].nT;

  std::ofstream outfile;
  outfile.open ("merged_sols3D.xsol");
  for (size_t b=0; b < nB; ++b){
    for (size_t t=0; t < nT; ++t){
      double total = 0.0;
      for (size_t sol=0; sol < nSols; ++sol){
        if (sols[sol].x[b] <= t){
          total++;
        }
      }
      outfile << total/nSols << " ";
    }
    outfile << std::endl;
  }
  outfile << "END" << std::endl;
  outfile.close();
}
