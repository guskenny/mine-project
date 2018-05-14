#include "MineProblem.h"

MineProblem::MineProblem(const char* path, int id, int count)
  : problemPath(path), id(id), total(count), updated(0), update_ready(false), current(NULL)
{
   start_time = std::chrono::high_resolution_clock::now();

  //Load the problem file...

  data = new Daten(path,'p');

  //Just create one ProblemData object for each proc for testing
  for (int p=0; p<count; p++)
  {
    ProblemData* pdata = new ProblemData();
    worklist.push_back(pdata);
  }

  //Initial problem split
  split();
}

MineProblem::~MineProblem()
{
  delete data;

  for (ProblemData* problem : worklist)
    delete problem;
}

void MineProblem::split() 
{
  //Split loaded problem, assign processes to worklist items
  // for now just divide the problem count among processes as evenly as possible
  int rounded = worklist.size() / total;
  int remain = worklist.size() % total;
  int workitem = 0;
  //For each proc...
  for (int proc=0; proc<total; proc++)
  {
    int assigned = rounded;
    //Add one of remaining units to each proc until there are none
    if (remain > 0) 
    {
      assigned++;
      remain--;
    }
    //Set assigned proc
    for (int i=0; i<assigned; i++, workitem++)
    {
      worklist[workitem]->assignedProc = proc;
      if (id==0) std::cout << "Assigned item " << workitem << " to proc " << proc << std::endl;
    }
  }
}

void MineProblem::next()
{
  //Get next assigned work item
  current = NULL;
  for (int i=0; i<worklist.size(); i++)
  {
    if (worklist[i]->assignedProc == id && !worklist[i]->solved)
    {
      current = worklist[i];
      break;
    }
  }
}

void MineProblem::monitor()
{
  assert(current);
  std::cout << "Started monitor thread: " << std::this_thread::get_id() << std::endl; 
  //mutex.lock();
  while (!solved() && !timeout())
  {
    auto start = std::chrono::high_resolution_clock::now();
    //Unlock the solver and sleep for a bit to let it work
    //mutex.unlock();
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));

    //Check if updated
    bool status = sync();

    //Lock and show status
    mutex.lock();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end-start;
#ifdef DEBUG
    std::cout << " -- Updated and Synched? " << status << ". Waited " << elapsed.count() << " ms, on proc " << (id+1) << std::endl;
#endif
    //Clear updated flag
    updated = false;
    mutex.unlock();
  }
  mutex.unlock();
  std::cout << "Monitor thread done " << std::this_thread::get_id() << std::endl; 
}

bool MineProblem::sync()
{
  std::lock_guard<std::mutex> guard(mutex);

  //Check if update required
  MPI_Allreduce(MPI_IN_PLACE, &updated, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (updated == 0)
    return false;

  //if (setsolved)
  //  current->solved = true;

  double bestobj = current->sol_int.obj;

  //Update obj via MPI
  MPI_Allreduce(MPI_IN_PLACE, &bestobj, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  //printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
  //std::cout << "Synchronising data, proc " << (id+1) << ", current obj: " << current->sol_int.obj << " best available: " << bestobj << std::endl;
  //Get solution obj from each procs to array on all procs
  double obj[total];
  //MPI_Gather(&current->sol_int.obj, 1, MPI_DOUBLE, obj, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Allgather(&current->sol_int.obj, 1, MPI_DOUBLE, obj, 1, MPI_DOUBLE, MPI_COMM_WORLD);

  //Find first process with best solution, send this solution to other processes
  double orig = current->sol_int.obj;
  for (int p=0; p<total; p++)
  {
    //printf("ON PROC %d : PROC %d/%d OBJ %f BEST %f\n", id+1, p+1, total, obj[p], bestobj);
    if (obj[p] == bestobj)
    {
      //std::cout << "Broadcasting best solution from PROC " << p << std::endl;
      MPI_Bcast(&current->sol_int.obj, 1, MPI_DOUBLE, p, MPI_COMM_WORLD);
      MPI_Bcast(current->sol_int.x.data(), current->sol_int.x.size(), MPI_INT, p, MPI_COMM_WORLD);
      //This could be slow, alternative is to use a contiguous container for Y data (same for dest/w_start below)
      for (auto row : current->sol_int.y)
        MPI_Bcast(row.data(), row.size(), MPI_DOUBLE, p, MPI_COMM_WORLD);

      //Also send the branch node info
      MPI_Bcast(&current->node.ub, 1, MPI_DOUBLE, p, MPI_COMM_WORLD);
      MPI_Bcast(current->node.time.data(), current->node.time.size()*2, MPI_INT, p, MPI_COMM_WORLD);
      for (auto row : current->node.dest)
        MPI_Bcast(row.data(), row.size(), MPI_DOUBLE, p, MPI_COMM_WORLD);
      for (auto row : current->node.w_start)
        MPI_Bcast(row.data(), row.size(), MPI_DOUBLE, p, MPI_COMM_WORLD);

      //Flag update received and waiting to be applied on receiving procs
      if (id != p)
        update_ready = true;
      break;
    }
  }

  //std::cout << "OBJECTIVE orig: " << orig << " updated: " << current->sol_int.obj << " should equal == best: " << bestobj << std::endl;
  //printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

  return true;
}


