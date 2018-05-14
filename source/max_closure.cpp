/* LP Relaxation: Open Pit Mining(PCPSP)
 * author: Davaa Baatar
 * date: 19.06.2016
*/
#include <max_closure.h>
MaxClosure::MaxClosure(Daten & data): _data(data), _tails(0), _heads(0){
	
	const int nBlocks = data.getNBlock();
	const int t_max = data.getNPeriod();
	const int d_max = data.getnDestination();
	_nVertices = nBlocks*t_max*d_max + 1;  // one dummy node(0)
	_nArcs = data.getnArcs();				 // block precedence arcs
	_nArcs += nBlocks*d_max;				 // destination: chain + linking(incl. dummy)
	_nArcs *= t_max;						 // dublicate for each period
	_nArcs += nBlocks;						 // return arcs for the dummy
	
	_tails = new int[_nArcs]();
	_heads = new int[_nArcs]();
	_obj = new double[_nArcs]();
	_sol = std::vector<int> (_nVertices-1, 0);
	
	int current_index = 0;
	for(int t=0; t<t_max; t++)
		for(int d=0; d<d_max; d++){
			std::vector<Block> * blocks = data.getBlock();
			for(int b=0; b<nBlocks; b++){
				int index = convert_triplet_to_index(data, b, t, d) + 1; // index in nNodes
				// add block precedence arcs
				if(d==d_max-1){
					int n = (*blocks)[b].getNumPred();
					const std::vector<int> * pred = (*blocks)[b].getPreds();
					for(int i=0; i<n; i++){
						int a = (*pred)[i];
						int a_index = convert_triplet_to_index(data, a, t, d) + 1;
						_heads[current_index] = index;
						_tails[current_index] = a_index;
						current_index++;
					}
				
					// add linking arcs
					if(t<t_max-1){
						int tail_index = convert_triplet_to_index(data, b, t+1, 0) + 1;
						_heads[current_index] = index;
						_tails[current_index] = tail_index;
						current_index++;
					}else {// arcs from dummy
						_heads[current_index] = index;
						_tails[current_index] = 0;
						_obj[current_index] = 1.0;
						current_index++;
					}
				}
				// add chain precedence
				if(d>0){
					int head_index = convert_triplet_to_index(data, b, t, d-1) + 1;
					_heads[current_index] = head_index;
					_tails[current_index] = index;
					current_index++;
				}else{
					if(t==0){ // arcs to dummy
						_heads[current_index] = 0; 
						_tails[current_index] = index;
						current_index++;
					}
				}
			}
		}
		// initialize cplex
		_env = NULL;
		_net = NULL;	
		int        	  status = 0;
		
		_env = CPXopenCPLEX (&status);
		if ( _env == NULL ) {
			  char  errmsg[CPXMESSAGEBUFSIZE];
			  fprintf (stderr, "Could not open CPLEX environment.\n");
			  CPXgeterrorstring (_env, status, errmsg);
			  fprintf (stderr, "%s", errmsg);
			  throw(-1);
		}
		
		status = CPXsetintparam (_env, CPX_PARAM_SCRIND, CPX_OFF);
	    if ( status ) {
			  fprintf (stderr, 
					   "Failure to turn on screen indicator, error %d.\n", status);
			  throw(-1);
		}
		status = CPXsetintparam (_env, CPX_PARAM_CLOCKTYPE, 1);
	    if ( status ) {
			  fprintf (stderr, 
					   "Failure to set cpu timer, error %d.\n", status);
			  throw(-1);
		}
		_net = CPXNETcreateprob (_env, &status, "maxcut");
		if ( _net == NULL ) {
			  fprintf (stderr, "Failed to create network object.\n");
			  throw(-1);
		}
		if ( CPXNETgetnumnodes (_env, _net) > 0 ) {
			  status = CPXNETdelnodes (_env, _net, 0, CPXNETgetnumnodes (_env, _net)-1);
			  if ( status ) throw(-1);
		}
		
	    status = CPXNETchgobjsen (_env, _net, CPX_MIN);
	    if ( status ) throw(-1);
		
		status = CPXNETaddnodes (_env, _net, _nVertices, NULL, NULL);
		if ( status ) throw(-1);
		
		status = CPXNETaddarcs (_env, _net, _nArcs, _tails, _heads, NULL, NULL, _obj, NULL);
		if ( status ) throw(-1);
}

MaxClosure::~MaxClosure(){
	// free memory
	if(_tails != NULL) delete [] _tails;
	if(_heads != NULL) delete [] _heads;
	if(_obj != NULL)   delete [] _obj;
	int  status = 0;
	if ( _net != NULL ) {
		status = CPXNETfreeprob (_env, &_net);
		if ( status ) {
			fprintf (stderr, "CPXNETfreeprob failed, error code %d.\n", status);
		}
    }
    if (_env != NULL ) {
		  status = CPXcloseCPLEX (&_env);
		  if ( status ) {
		  char  errmsg[CPXMESSAGEBUFSIZE];
			 fprintf (stderr, "Could not close CPLEX environment.\n");
			 CPXgeterrorstring (_env, status, errmsg);
			 fprintf (stderr, "%s", errmsg);
		  }
   }
}

void MaxClosure::solve(double * supply)
{
	clock_t start, end;
	start = clock();
	
	int  status = 0;
		
	int * indices = NULL;
	indices = new int[_nVertices];
	for(int i =0; i< _nVertices; i++) indices[i] = i;
/* scalling the supplies*/
	
for(int i=0; i< _nVertices; i++) // SCALLING
		supply[i] *= 1e-8;		 // SCALLING
	status =  CPXNETchgsupply( _env, _net, _nVertices, indices, supply );
	if ( status ) {
	  fprintf (stderr, 
		   "Failure to update the supplies, error %d.\n", status);
	  throw(-1);
	}
	delete[] indices; indices = NULL;
   
   /* Optimize the problem and obtain solution. */
   status = CPXNETprimopt (_env, _net);
   if ( status ) {
      fprintf (stderr, "Failed to optimize network.\n");
      throw(-1);
   }

   /* get network dimensions */

   int nnodes = CPXNETgetnumnodes (_env, _net);
   if(nnodes!=_nVertices){
	  fprintf (stderr, 
		   "Number of Nodes does not match, error %d, %d.\n", _nVertices, nnodes);
	  throw(-1);
   }
   double   * pi    = NULL;
   pi = new double[nnodes];
   status = CPXNETsolution (_env, _net, &_solstat, &_objVal, NULL, pi, NULL, NULL);
   _objVal *= 1e8; // RESCALLING
   if ( status ) {
      fprintf (stderr, "Failed to obtain dual values.\n");
      throw(-1);
   }
   
   /* get solution(max cut) */
   bool reduce = 0;
   double val;
   if(pi[0] > 1e6 || pi[0] < -1e-6){ 
		reduce = 1;
		val = pi[0];
   }	
   for(int i=1; i<_nVertices; i++){
		if(reduce){
			_sol[i-1] = abs( (int) (pi[i] - val)); 
		}else _sol[i-1] = abs(pi[i]);
   }	
   end = clock();
   _t_cpu = ((double)(end-start)) / (double) CLOCKS_PER_SEC;
   
   /* Free memory */
   delete[] pi; 
   pi = NULL;
}
double MaxClosure::getObjVal(){return _objVal;}
double MaxClosure::getTime(){return _t_cpu;}
bool MaxClosure::updatePartitioning(std::vector<std::vector<int> > & Cpart){
	bool stop = 1;
	int n = (int) Cpart.size();
	std::vector<std::vector<int> > temp;
	if(n==0){
		stop = 0;
		temp.resize(2);
		temp[0].reserve(_nVertices-1);
		temp[1].reserve(_nVertices-1);
		for(int j=0; j< _nVertices-1; j++){
			if(_sol[j] > 0.5){
				temp[1].push_back(j);
			}else temp[0].push_back(j);
		}
	}else{
		temp.resize(2*n);
		
		int index = 0;
		for(int i=0; i< n; i++){
			int m = (int) Cpart[i].size();
			std::vector<int> I; I.reserve(m);
			std::vector<int> O; O.reserve(m);
			for(int j=0; j<m; j++){
				if(_sol[Cpart[i][j]]>0.5) I.push_back(Cpart[i][j]);
				else O.push_back(Cpart[i][j]);
			}
			int overlaps=0;
			if(I.size()>0){ temp[index]=I; index++; overlaps++;}
			if(O.size()>0){ temp[index]=O; index++; overlaps++;}
			/// check stopping criteria
			if(stop && overlaps==2) stop = 0;
		}
		temp.resize(index);
	}
	Cpart = temp;
	
	return stop;
}
