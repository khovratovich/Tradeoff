// Lyra-test.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <vector>
#include <algorithm>
#include <stdint.h>
using namespace std;


 #include "MersenneTwister.h"
#define WIDTH_TEST 4096

#define START_LATENCY 0
#define RAND32 (mtr.randInt())

uint64_t Q1  =  0xaaaaaaaaaaaaaaaa;
uint64_t Q2  =  0x5555555555555555;
uint64_t Q3 = 0xcccccccccccccccc;
uint64_t Q4 = 0x3333333333333333;

FILE* fp;


inline unsigned wt2(uint64_t v)//Hamming weight
{
	uint64_t a;
	uint64_t v1 = v&Q1;
	uint64_t v2 = v&Q2;
	a = (v1>>1) + v2;
	v1 = a&Q3;
	v2 = a&Q4;
	a = (v1>>2) + v2;
	a = (a&0xf0f0f0f0f0f0f0f0)>>4 + (a&0x0f0f0f0f0f0f0f0f);
	a = (a&0xff00ff00ff00ff00)>>8 + (a&0x00ff00ff00ff00ff);
	a = (a&0xffff0000ffff0000)>>16+ (a&0x0000ffff0000ffff);
	a = (a>>32)+ (a&0xffffffff);
	return a;
}

inline unsigned wt(uint64_t v)//Hamming weight
{
	v = ((v&0xaaaaaaaaaaaaaaaa)>>1) + (v&0x5555555555555555);
	v = ((v&0xcccccccccccccccc)>>2) + (v&0x3333333333333333);
	v = ((v&0xf0f0f0f0f0f0f0f0)>>4) + (v&0x0f0f0f0f0f0f0f0f);
	v = ((v&0xff00ff00ff00ff00)>>8) + (v&0x00ff00ff00ff00ff);
	v = ((v&0xffff0000ffff0000)>>16)+ (v&0x0000ffff0000ffff);
	v = (v>>32)+ (v&0xffffffff);
	return v;
}

struct intervals_st
{
	vector<uint64_t> masks;
	intervals_st(){};
	intervals_st(vector<uint64_t> &m){masks=m;};
	intervals_st operator|(const intervals_st &r) 
	{
		vector<uint64_t> tmp;
		for(unsigned i=0; i<min(masks.size(),r.masks.size()); ++i)
			tmp.push_back(masks[i] | r.masks[i]);
		return intervals_st(tmp);

	};

	unsigned weight()
	{
		unsigned a = 0;
		for(unsigned i=0; i<masks.size(); ++i)
			a += wt(masks[i]);
		return a;
	}
};

struct row{
	
	uint64_t latency;
	bool stored;
	unsigned extra_row;
	uint64_t access;
	uint64_t read;
	row(bool b=false, unsigned e=0, uint64_t a=0)
	{stored = b; extra_row = e; access=a;latency=START_LATENCY; read=1;}
};

struct ext_row: public row{
	intervals_st intervals;
	ext_row(bool b=false, unsigned e=0)
	{stored = b; extra_row = e; access=0;latency=START_LATENCY;}
};

struct tradeoff
{
	float save;
	float penalty;
	float latency;
	float intervals;
	float read;
	tradeoff(float s=1, float p=1, float l =1 , float i=1, float r=1){save = s; penalty=p;latency = l;intervals = i; read =r;};
	bool operator<(const tradeoff& r)
	{
		if(save < r.save)
			return true;
		else if(save == r.save)
		{
			if(penalty < r.penalty)
				return true;
			else /*if(penalty == r.penalty)
				return true; 
			else */return false;
		}
		else return false;
	}
};



struct RandomGraph
{	
	unsigned width;
	vector<row> rows;
	RandomGraph(unsigned w){rows.resize(w); width = w; visited.resize(w,false);};
	vector<bool> visited;

	void clean(){visited.assign(width,false);}
	row& operator[](unsigned i){return rows[i%width];}
	void Traverse(unsigned position, unsigned time);
	unsigned count(){unsigned a=0;
					for(unsigned j=0; j<width; ++j)
						if(visited[j]) 
							a++; 
					return a;};
};

struct RandomGraphInt
{
	unsigned width;
	vector<ext_row> rows;
	RandomGraphInt(unsigned w){rows.resize(w); width = w; };
	
	ext_row& operator[](unsigned i){return rows[i%width];}
};


void RandomGraph::Traverse(unsigned  position, unsigned time)
{
	if(time>width)
		time = width;
	if(time==0)
		return;
	else
	{
		if(time >position)//Normal update
		{
			for(int j=time-1; j>position; j--)//Check if the row has been updated as extra_row
			{
				if(rows[j].extra_row == position)
					Traverse(j-1,j);
			}
		
			visited[position] = true;
			if(!rows[position].stored)
			{
				Traverse((position==0)?(width-1):(position-1),position); //Previous
				Traverse(rows[position].extra_row, position);
			}
		}
		if(!rows[position].stored)
		{
			for(int j=time-1; j>=0; j--) //Earlier updates as extra row
			{
				if(rows[j].extra_row == position)
						Traverse((j==0)?(width-1):(j-1),j);
			}
		}
	}
}

void LyraTest(unsigned width, unsigned q)
{
	//Fill
	MTRand mtr;
	RandomGraph Lyra(width);
	vector<vector<unsigned> > columns; 
	columns.resize(width);
	for(unsigned i=0; i< width; ++i)
	{
		unsigned extra_row = RAND32%width;
		if(i%q==0)
			Lyra.rows[i] = row(true, extra_row);
		else
			Lyra.rows[i] = row(false, extra_row);
	}

	//Test
	unsigned tests = 32;
	float average=0;
	for(unsigned t=0; t<tests; ++t)
	{
		unsigned pos = RAND32%width;
		Lyra.Traverse(pos, pos+1);
		unsigned depends = Lyra.count();
		average += (float)depends/(tests*width);
		printf("Position %d (of %d), depends on %d\n", pos, width, depends);
		Lyra.clean();
	}
	printf("Average proportion %2.2f Step %d \n", average,q);

}

tradeoff SimplePenalty(unsigned q)
{
	MTRand mtr;
	vector<float> online;
	online.resize(q);
	vector<float> access;
	access.resize(q);
	
	online[0] = 1;
	access[0] = 0;
	float sum=online[0];
	for(unsigned i=1; i< q; ++i)
	{
		float tmp=1 + 1/(float)i;
		for(unsigned j=1; j<i; j++)
			tmp+=access[j]/j;
		online[i] = tmp;
		access[i] = (i>q/2)?(q*online[i]/2):(i*online[i]);
		sum+=online[i];
	}
	tradeoff a;
	a.penalty = sum/q;
	a.save = 2*q/(3+log((double)q-1)/log((double)2));
	return a;
}



tradeoff LyraCompute(unsigned width, unsigned q)  //Wandering phase
{
	//Fill
	MTRand mtr;
	unsigned input_cost = 40;
	RandomGraph Lyra(width);
	vector<vector<unsigned> > columns;  //Alternative storage
	columns.resize(width);
	unsigned extra_row =1;
	for(unsigned i=0; i< width; ++i)
	{
		extra_row = RAND32%width;
		columns[extra_row].push_back(i);
		if(i%q==0)
		{
			Lyra[i] = row(true, extra_row);
		}
		else
		{
			Lyra[i] = row(false, extra_row);

		}		
		Lyra[i].access=input_cost;
	}

	uint64_t layer_access=0;
	uint64_t layer_compute=0;
	uint64_t layer_stored = 0;
	uint64_t total_stored=0;
	uint64_t total_compute=0;
	uint64_t total_extra_latency=0;
	
	//Test
	float prev_access=1;  
	float prev_compute=1;
	unsigned layer_length = width/q;
	unsigned top_size = layer_length;
	vector<uint64_t> top;  //The largest access cost values, sorted, top[0] is the highest
	top.resize(top_size,input_cost); 
	for(unsigned i=1; i< width; ++i)
	{
		uint64_t extra_access = Lyra[Lyra[i].extra_row].access;  //Cost of extra row
		uint64_t extra_latency = Lyra[Lyra[i].extra_row].latency;  //Cost of extra row
		layer_compute += extra_access+1;
		total_compute += extra_access+1;
		total_extra_latency += extra_latency;
		
		
		if(extra_access> top[top_size-1]) //Among the most expensive rows => need to store output of F
		{
			Lyra[i].stored = true;
			layer_stored++;
			total_stored++;
			//Access cost is not changed
			//Latency is not changed
		}
		else if(Lyra[i].stored)  
		{
			layer_stored++;
			total_stored++;
			//Access cost is not changed
			//Latency is not changed
		}
		else
		{
			uint64_t new_access = Lyra[i-1].access + Lyra[Lyra[i].extra_row].access +1;  //New access cost
			Lyra[Lyra[i].extra_row].access  = Lyra[i].access = new_access;
			uint64_t new_latency = max(Lyra[i-1].latency,Lyra[Lyra[i].extra_row].latency) +1; //New latency
			Lyra[Lyra[i].extra_row].latency  = Lyra[i].latency = new_latency;
			
			unsigned rank=top_size;  //Check if the new access complexity is in the top
			while(rank>0)
			{
				if(top[rank-1]<new_access)
					rank--;
				else
					break;
			}
			
			top.insert(top.begin()+rank,new_access);
			top.insert(top.begin()+rank,new_access);
			top.pop_back();
			top.pop_back();
		}

		if((i+1)%layer_length==0)//End of layer
		{
			/*printf("Layer %d Average compute: %2.2f (Factor %2.2f) Stored %2.2f\n", i/(width/q),
				(float)layer_compute/layer_length,(float)layer_compute/(layer_length*prev_compute), (float)layer_stored/layer_length);*/
			
			prev_compute = (float)layer_compute/layer_length;
			layer_compute = layer_stored =  0;
		}
	}
	float save = (float) width/total_stored;
	float penalty = (float)total_compute/(width-1);
	float latency = (float)total_extra_latency/(width-1);
	//printf("Memory saving %2.2f ", (float) width/total_stored);
	//double log_penalty = log((double)total_compute/width)/log((double)2);
//	printf("Log-Total compute: %2.2f Log-Penalty: %2.2f Memory factor:%d\n", log((double)total_compute)/log((double)2),log((double)total_compute/width)/log((double)2),q/2);
	return tradeoff(save,penalty,latency);

}

tradeoff SimpleRandom(unsigned width, unsigned q)  //Attacking random graph with simple technique
{
	//Fill
	MTRand mtr;
	unsigned input_cost = 0;
	RandomGraph Lyra(width);
	vector<vector<unsigned> > columns;  //Alternative storage
	columns.resize(width);
	unsigned extra_row =1;
	for(unsigned i=0; i< width; ++i)
	{
		extra_row = (i==0)?0: RAND32%i;
		columns[extra_row].push_back(i);
		if(i%q==0)
		{
			Lyra[i] = row(true, extra_row);
		}
		else
		{
			Lyra[i] = row(false, extra_row);

		}		
		Lyra[i].access=input_cost;
	}

	uint64_t layer_access=0;
	uint64_t layer_compute=0;
	uint64_t layer_stored = 0;
	uint64_t total_stored=0;
	uint64_t total_compute=0;
	uint64_t total_extra_latency=0;
	
	//Test
	float prev_access=1;  
	float prev_compute=1;
	unsigned layer_length = width/q;
	unsigned top_size = layer_length;
	//vector<uint64_t> top;  //The largest access cost values, sorted, top[0] is the highest
	//top.resize(top_size,input_cost); 
	for(unsigned i=1; i< width; ++i)
	{
		extra_row = Lyra[i].extra_row;
		uint64_t extra_access = Lyra[extra_row].access;  //Cost of extra row
		uint64_t extra_latency = Lyra[extra_row].latency;  //Cost of extra row
		total_extra_latency += extra_latency;
		
		unsigned curr_layer = i/layer_length;
		unsigned extra_layer = extra_row/layer_length;

		if(extra_layer==curr_layer)
		{
			Lyra[i].access=0;
			total_compute += 1;
			total_stored++;
			Lyra[i].latency =1;
		}
		else
		{
			if(Lyra[i].stored)  
			{
				Lyra[i].access=0;
				total_stored++;
				Lyra[i].latency =1;
			}
			else
			{
				Lyra[i].latency =Lyra[extra_row].latency+1;
				Lyra[i].access = Lyra[i-1].access + extra_access+1;
			}
			total_compute += extra_access+1;

		}

		
	}
	float save = (float) width/total_stored + 1/(float)q;
	float penalty = (double)total_compute/(width-1);
	float latency = (double)total_extra_latency/(width-1);
	//printf("Memory saving %2.2f ", (float) width/total_stored);
	//double log_penalty = log((double)total_compute/width)/log((double)2);
//	printf("Log-Total compute: %2.2f Log-Penalty: %2.2f Memory factor:%d\n", log((double)total_compute)/log((double)2),log((double)total_compute/width)/log((double)2),q/2);
	return tradeoff(save,penalty,latency);

}

tradeoff RankingLyra(unsigned q, unsigned top_size)  //Ranking method with parameter <q> for Lyra graph, where F is stored
{
		//Fill
	MTRand mtr;
	unsigned input_cost = 0;
	uint32_t width = WIDTH_TEST;
	RandomGraph Graph(width);
	unsigned extra_row =1;
	for(unsigned i=0; i< width; ++i)
	{
		extra_row = (i==0)?0: RAND32%i;
		Graph[i] = row(false, extra_row,input_cost);
	}

	uint64_t layer_access=0;
	uint64_t layer_compute=0;
	uint64_t layer_stored = 0;
	uint64_t total_stored=0;
	uint64_t total_compute=0;
	uint64_t total_latency=0;
	uint64_t total_read=0;
	
	//Test
	float prev_access=1;  
	float prev_compute=1;
	unsigned layer_length = width/q;
	vector<uint64_t> top;  //The largest access cost values, sorted, top[0] is the highest
	top.resize(top_size,input_cost); 
	for(unsigned i=1; i< width; ++i)
	{
		extra_row = Graph[i].extra_row;
		total_compute += Graph[i].access+Graph[extra_row].access+1;
		total_latency += Graph[i].latency+ Graph[extra_row].latency+1;
		total_read +=  Graph[i].read+Graph[extra_row].read;

		if(Graph[extra_row].access> top[top_size-1]) //Among the most expensive rows => need to store output of F
		{
			Graph[i].stored = true;
			//Access cost remains unchanged for both nodes
			total_stored++;
			//Latency  remains unchanged for both nodes
			Graph[i].read++;
			Graph[extra_row].read++;
		}
		else
		{
			if(i%q==0)//Segment store rule
			{
				Graph[i].stored = true;
				//Access cost remains unchanged for both nodes
				total_stored++;
				//Latency  remains unchanged for both nodes
				Graph[i].read++;
				Graph[extra_row].read++;
			}
			else
			{
				Graph[extra_row].latency = Graph[i].latency = max(Graph[extra_row].latency, max(Graph[i-1].latency, Graph[i].latency))+1; //both nodes have new latency
				Graph[extra_row].access = Graph[i].access = Graph[i].access + Graph[i-1].access + Graph[extra_row].access+1; //both nodes have new access cost
				Graph[extra_row].read = Graph[i].read = Graph[i-1].read + Graph[i].read+ Graph[extra_row].read; ////both nodes have new red values
				unsigned rank=top_size;  //Check if the new access complexity is in the top
				while(rank>0)//adding new value to the top
				{
					if(top[rank-1]<Graph[i].access)
						rank--;
					else break;
				}			
				top.insert(top.begin()+rank,Graph[i].access);
				top.pop_back(); 
			}
		}		
	}
	float save = (float) width/total_stored + 1/(float)96;
	float penalty = (double)total_compute/width;
	float latency = (double)total_latency/width;
	float read = (double)total_read/width;
	return tradeoff(save,penalty,latency,0,read);

}

tradeoff RankingRandom(unsigned q, unsigned top_size)  //Ranking method with parameter <q> for random graph
{
		//Fill
	MTRand mtr;
	unsigned input_cost = 0;
	uint32_t width = WIDTH_TEST;
	RandomGraph Graph(width);
	vector<vector<unsigned> > columns;  //Alternative storage
	columns.resize(width);
	unsigned extra_row =1;
	for(unsigned i=0; i< width; ++i)
	{
		extra_row = (i==0)?0: RAND32%i;
		columns[extra_row].push_back(i);
		if(i%q==0)
		{
			Graph[i] = row(true, extra_row);
		}
		else
		{
			Graph[i] = row(false, extra_row);

		}		
		Graph[i].access=input_cost;
	}

	uint64_t layer_access=0;
	uint64_t layer_compute=0;
	uint64_t layer_stored = 0;
	uint64_t total_stored=0;
	uint64_t total_compute=0;
	uint64_t total_latency=0;
	uint64_t total_read=0;
	
	//Test
	float prev_access=1;  
	float prev_compute=1;
	unsigned layer_length = width/q;
	vector<uint64_t> top;  //The largest access cost values, sorted, top[0] is the highest
	top.resize(top_size,input_cost); 
	for(unsigned i=1; i< width; ++i)
	{
		extra_row = Graph[i].extra_row;
		total_compute += Graph[extra_row].access+1;
		total_latency += Graph[extra_row].latency+1;
		total_read +=  Graph[extra_row].read;

		if(Graph[extra_row].access> top[top_size-1]) //Among the most expensive rows => need to store output of F
		{
			Graph[i].stored = true;
			Graph[i].access=0;
			total_stored++;
			Graph[i].latency =0;
			Graph[i].read =1;
		}
		else
		{
			if(Graph[i].stored)  
			{
				Graph[i].access=0;
				total_stored++;
				Graph[i].latency =0;
				Graph[i].read =1;
			}
			else
			{
				Graph[i].latency = max(Graph[extra_row].latency, Graph[i-1].latency)+1;
				Graph[i].access = Graph[i-1].access + Graph[extra_row].access+1;
				Graph[i].read = Graph[i-1].read +  Graph[extra_row].read;
				unsigned rank=top_size;  //Check if the new access complexity is in the top
				while(rank>0)//adding new value to the top
				{
					if(top[rank-1]<Graph[i].access)
						rank--;
					else break;
				}			
				top.insert(top.begin()+rank,Graph[i].access);
				top.pop_back(); 
			}
		}		
	}
	float save = (float) width/total_stored;
	float penalty = (double)total_compute/(width-1);
	float latency = (double)total_latency/(width-1);
	float read = (double)total_read/(width-1);
	return tradeoff(save,penalty,latency,0,read);

}

tradeoff AdvancedMultiRandom(unsigned width, unsigned q, unsigned k)  //Attacking random k-pass graph with advanced technique
{
	//Fill
	MTRand mtr;
	unsigned input_cost = 0;
	RandomGraph Graph(width);
	vector<vector<unsigned> > columns;  //Alternative storage
	columns.resize(width);
	unsigned extra_row =1;
	for(unsigned i=0; i< width; ++i)
	{
		if(i==0)
			extra_row=0;
		else if(i< width/k)
			extra_row = RAND32%i;
		else extra_row = i- RAND32%(width/k);
		columns[extra_row].push_back(i);
		if(i%q==0)
		{
			Graph[i] = row(true, extra_row);
		}
		else
		{
			Graph[i] = row(false, extra_row);

		}		
		Graph[i].access=input_cost;
	}

	uint64_t layer_access=0;
	uint64_t layer_compute=0;
	uint64_t layer_stored = 0;
	uint64_t total_stored=0;
	uint64_t total_compute=0;
	uint64_t total_latency=0;
	uint64_t total_read=0;
	
	//Test
	float prev_access=1;  
	float prev_compute=1;
	unsigned layer_length = width/q;
	unsigned top_size = layer_length;
	vector<uint64_t> top;  //The largest access cost values, sorted, top[0] is the highest
	top.resize(top_size,input_cost); 
	for(unsigned i=1; i< width; ++i)
	{
		extra_row = Graph[i].extra_row;
		total_compute += Graph[extra_row].access+1;
		total_latency += Graph[extra_row].latency+1;
		total_read +=  Graph[extra_row].read;

		if(Graph[extra_row].access> top[top_size-1]) //Among the most expensive rows => need to store output of F
		{
			if((i< width/k) || (!Graph[i].stored))//first pass or the same element was not stored at the first pass
			{
				Graph[i].stored = true;
				total_stored++;
			}

			Graph[i].access=0;
			Graph[i].latency =0;
			Graph[i].read =1;
		}
		else
		{
			if(i%q==0) //Store always, do not set the .stored field to true
			{
				Graph[i].access=0;
				if(i< width/k)//first pass
					total_stored++;
				Graph[i].latency =0;
				Graph[i].read =1;
			}
			else
			{
				if((i> width/k) && (Graph[i].stored))//do not store it anymore
				{
					Graph[i].stored = false;
					total_stored--;
				}
				Graph[i].latency = max(Graph[extra_row].latency, Graph[i-1].latency)+1;
				Graph[i].access = Graph[i-1].access + Graph[extra_row].access+1;
				Graph[i].read = Graph[i-1].read +  Graph[extra_row].read;
				unsigned rank=top_size;  //Check if the new access complexity is in the top
				while(rank>0)//adding new value to the top
				{
					if(top[rank-1]<Graph[i].access)
						rank--;
					else break;
				}			
				top.insert(top.begin()+rank,Graph[i].access);
				top.pop_back(); 
			}
		}		
	}
	float save = (float) width/total_stored + 1/(float)32;
	float penalty = (double)total_compute/(width-1);
	float latency = (double)total_latency/(width-1);
	float read = (double)total_read/(width-1);
	//printf("Memory saving %2.2f ", (float) width/total_stored);
	//double log_penalty = log((double)total_compute/width)/log((double)2);
//	printf("Log-Total compute: %2.2f Log-Penalty: %2.2f Memory factor:%d\n", log((double)total_compute)/log((double)2),log((double)total_compute/width)/log((double)2),q/2);
	return tradeoff(save,penalty,latency,0,read);

}

tradeoff SetupTradeoff(unsigned	q, unsigned l)
{
	float save =1/( (1-pow((float)2,-(float)l))/q +  pow((float)2,-(float)l));
	float penalty = ((float)q-1)*(l)/2;
	return tradeoff(save, penalty,q);
}


tradeoff LyraWandInt(unsigned width, unsigned q, unsigned q_setup)  //Wandering phase with intervals 
{
	//Fill
	MTRand mtr;
	unsigned input_cost = 0;
	RandomGraphInt Lyra(width);
	vector<vector<unsigned> > columns;  //Alternative storage
	columns.resize(width);
	unsigned extra_row =1;
	for(unsigned i=0; i< width; ++i)
	{
		extra_row = RAND32%width;
		columns[extra_row].push_back(i);
		if(i%q==0)
		{
			Lyra[i] = ext_row(true, extra_row);
		}
		else
		{
			Lyra[i] = ext_row(false, extra_row);

		}		
		Lyra[i].access=input_cost;
		Lyra[i].intervals.masks.resize((width/q_setup)/64+1,0);
		Lyra[i].intervals.masks[(i/q_setup)/64] =((uint64_t)1)<< ((i/q_setup)%64);
	}

	uint64_t layer_access=0;
	uint64_t layer_compute=0;
	uint64_t layer_stored = 0;
	uint64_t total_stored=0;
	uint64_t total_compute=0;
	uint64_t total_extra_latency=0;
	uint64_t total_intervals = 0;
	
	//Test
	float prev_access=1;  
	float prev_compute=1;
	unsigned layer_length = width/q;
	unsigned top_size = layer_length;
	vector<uint64_t> top;  //The largest access cost values, sorted, top[0] is the highest
	top.resize(top_size,input_cost); 
	for(unsigned i=1; i< width; ++i)
	{
		uint64_t extra_access = Lyra[Lyra[i].extra_row].access;  //Cost of extra row
		uint64_t extra_latency = Lyra[Lyra[i].extra_row].latency;  //Cost of extra row
		layer_compute += extra_access+1;
		total_compute += extra_access+1;
		total_extra_latency += extra_latency;
		intervals_st extra_intervals  = Lyra[Lyra[i].extra_row].intervals;
		total_intervals += extra_intervals.weight();
		
		if(extra_access> top[top_size-1]) //Among the most expensive rows => need to store output of F
		{
			Lyra[i].stored = true;
			layer_stored++;
			total_stored++;
			//Access cost is not changed
			//Latency is not changed
			//Intervals are not changed
		}
		else if(Lyra[i].stored)  
		{
			layer_stored++;
			total_stored++;
			//Access cost is not changed
			//Latency is not changed
			//Intervals are not changed
		}
		else
		{
			uint64_t new_access = Lyra[i-1].access + Lyra[Lyra[i].extra_row].access +1;  //New access cost
			Lyra[Lyra[i].extra_row].access  = Lyra[i].access = new_access;
			uint64_t new_latency = max(Lyra[i-1].latency,Lyra[Lyra[i].extra_row].latency) +1; //New latency
			Lyra[Lyra[i].extra_row].latency  = Lyra[i].latency = new_latency;
			intervals_st new_intervals = Lyra[i-1].intervals | Lyra[Lyra[i].extra_row].intervals;
			Lyra[i].intervals = Lyra[Lyra[i].extra_row].intervals = new_intervals;
			
			unsigned rank=top_size;  //Check if the new access complexity is in the top
			while(rank>0)
			{
				if(top[rank-1]<new_access)
					rank--;
				else
					break;
			}
			
			top.insert(top.begin()+rank,new_access);
			top.insert(top.begin()+rank,new_access);
			top.pop_back();
			top.pop_back();
		}

		if((i+1)%layer_length==0)//End of layer
		{
			/*printf("Layer %d Average compute: %2.2f (Factor %2.2f) Stored %2.2f\n", i/(width/q),
				(float)layer_compute/layer_length,(float)layer_compute/(layer_length*prev_compute), (float)layer_stored/layer_length);*/
			
			prev_compute = (float)layer_compute/layer_length;
			layer_compute = layer_stored =  0;
		}
	}
	float save = (float) width/total_stored;
	float penalty = (double)total_compute/(width-1);
	float latency = (double)total_extra_latency/(width-1);
	float intervals = (double)total_intervals/(width-1);
	//printf("Memory saving %2.2f ", (float) width/total_stored);
	//double log_penalty = log((double)total_compute/width)/log((double)2);
//	printf("Log-Total compute: %2.2f Log-Penalty: %2.2f Memory factor:%d\n", log((double)total_compute)/log((double)2),log((double)total_compute/width)/log((double)2),q/2);
	return tradeoff(save,penalty,latency,intervals);

}

tradeoff LyraSetup(unsigned log_width, unsigned q, unsigned depth)
{
	//Fill
	MTRand mtr;
	unsigned input_cost = 0;
	unsigned width = 1<<log_width;
	RandomGraph Lyra(width);
	vector<vector<unsigned> > columns;  //Alternative storage
	columns.resize(width);
	Lyra[0] = row(true,0);
	Lyra[1] = row(true,0);

	unsigned extra_row =1;
	unsigned level=1;
	unsigned curr_q = q<<(log_width-1);
	for(unsigned i=2; i< width; ++i)
	{
		
		extra_row = (extra_row==0)?(i-2):(extra_row-1);//RAND32%width;
		columns[extra_row].push_back(i);
		if((i%q==0)||(i<(1<<(log_width-depth))))
		{
			Lyra[i] = row(true, extra_row);
		}
		else
		{
			Lyra[i] = row(false, extra_row);

		}		
		Lyra[i].access=input_cost;
	}

	uint64_t layer_access =0;
	uint64_t layer_compute=0;
	uint64_t layer_stored =0;
	uint64_t total_stored =0;
	uint64_t total_compute=0;
	uint64_t total_latency=0;
	
	//Test
	float prev_access=1;  
	float prev_compute=1;
	unsigned layer_length = width/q;
	unsigned top_size = layer_length;
	vector<uint64_t> top;  //The largest access cost values, sorted, top[0] is the highest
	top.resize(top_size,input_cost); 
	for(unsigned i=1; i< width; ++i)
	{
		uint64_t extra_access = Lyra[Lyra[i].extra_row].access;  //Cost of extra row
		uint64_t extra_latency = Lyra[Lyra[i].extra_row].latency;  //Cost of extra row
		layer_compute += extra_access+1;
		total_compute += extra_access+1;
		//total_latency += extra_latency+1;
		
		
		if(extra_access> top[top_size-1]) //Among the most expensive rows => need to store output of F
		{
			Lyra[i].stored = true;
			layer_stored++;
			total_stored++;
			//Access cost is not changed
			//Latency is not changed
		}
		else if(Lyra[i].stored)  
		{
			layer_stored++;
			total_stored++;
			//Access cost is not changed
			//Latency is not changed
		}
		else
		{
			uint64_t new_access = Lyra[i-1].access + Lyra[Lyra[i].extra_row].access +1;  //New access cost
			Lyra[Lyra[i].extra_row].access  = Lyra[i].access = new_access;
			uint64_t new_latency = max(Lyra[i-1].latency,Lyra[Lyra[i].extra_row].latency) +1; //New latency
			Lyra[Lyra[i].extra_row].latency  = Lyra[i].latency = new_latency;
			
			unsigned rank=top_size;  //Check if the new access complexity is in the top
			while(rank>0)
			{
				if(top[rank-1]<new_access)
					rank--;
				else
					break;
			}
			
			top.insert(top.begin()+rank,new_access);
			top.insert(top.begin()+rank,new_access);
			top.pop_back();
			top.pop_back();
		}

		if((i+1)%layer_length==0)//End of layer
		{
			/*printf("Layer %d Average compute: %2.2f (Factor %2.2f) Stored %2.2f\n", i/(width/q),
				(float)layer_compute/layer_length,(float)layer_compute/(layer_length*prev_compute), (float)layer_stored/layer_length);*/
			
			prev_compute = (float)layer_compute/layer_length;
			layer_compute = layer_stored =  0;
		}
	}
	float save = (float) width/total_stored;
	float penalty = (double)total_compute/(width);
	for(unsigned i=0; i<width; ++i)
		total_latency += Lyra[i].latency;
	float latency = (double)total_latency/(width);
	//printf("Memory saving %2.2f ", (float) width/total_stored);
	//double log_penalty = log((double)total_compute/width)/log((double)2);
//	printf("Log-Total compute: %2.2f Log-Penalty: %2.2f Memory factor:%d\n", log((double)total_compute)/log((double)2),log((double)total_compute/width)/log((double)2),q/2);
	return tradeoff(save,penalty,latency);

}


double LyraOld(unsigned width, unsigned q)
{
	//Fill
	MTRand mtr;
	RandomGraph Lyra(width);
	vector<vector<unsigned> > columns; 
	columns.resize(width);
	for(unsigned i=0; i< width; ++i)
	{
		unsigned extra_row = RAND32%width;
		if(i%q==0)
		{
			Lyra[i] = row(true, extra_row);
		}
		else
		{
			Lyra[i] = row(false, extra_row);

		}		
		Lyra[i].access=0;
	}

	uint64_t layer_access=0;
	uint64_t layer_compute=0;
	uint64_t layer_stored = 0;
	uint64_t total_stored=0;
	uint64_t total_compute=0;
	//Test
	//printf("Width %d Step %d\n", width, q);
	float prev_access=1;
	float prev_compute=1;
	unsigned top_size = (width/q);
	vector<uint64_t> top;
	top.resize(top_size,0);
	for(unsigned i=1; i< width; ++i)
	{
		unsigned curr_layer = i/(width/q);
		unsigned extra_layer = Lyra[i].extra_row/(width/q);
		uint64_t extra_access = Lyra[Lyra[i].extra_row].access;
		layer_compute += extra_access+1;
		total_compute += extra_access+1;
		
		if(extra_access> top[top_size-1]) //need to store
		{
			Lyra[i].stored = true;
			layer_stored++;
			total_stored++;
		}
		/*if(curr_layer == extra_layer)
		{
			Lyra[i].stored = true;
			//Lyra[i].access = 0;
			layer_stored++;
		}*/
		else if(Lyra[i].stored){layer_stored++;}
			//Lyra[i].access = 0;
		else
		{
			uint64_t new_access = Lyra[i-1].access + Lyra[Lyra[i].extra_row].access +1;
			Lyra[Lyra[i].extra_row].access  = Lyra[i].access = new_access;
			unsigned rank=top_size;
			while(rank>0)
			{
				if(top[rank-1]<new_access)
					rank--;
				else
					break;
			}
			
			top.insert(top.begin()+rank,new_access);
			top.insert(top.begin()+rank,new_access);
			top.pop_back();
			top.pop_back();
		}
		layer_access += Lyra[i].access;

		if((i+1)%(width/q)==0)//Start of layer
		{
			/*printf("Layer %d Average access: %2.2f Average compute: %2.2f (Factor %2.2f) Stored %2.2f\n", i/(width/q),(float)layer_access/(width/q), 
				(float)layer_compute/(width/q),(float)layer_compute/(width*prev_compute/q), (float)layer_stored/(width/q));*/
			
			prev_compute = (float)layer_compute/(width/q);
			layer_access=layer_compute = layer_stored =  0;
		}
	}
	printf("Memory save %2.2f ", (float) width/total_stored);
	/*for(unsigned i=width; i>0; --i)
	{
		unsigned curr_layer = i/(width/q);
		unsigned extra_layer = Lyra[i].extra_row/(width/q);
		layer_compute += Lyra[Lyra[i].extra_row].access +1;
		total_compute += Lyra[Lyra[i].extra_row].access +1;
		if(curr_layer == extra_layer)
		{
			Lyra[i].stored = true;
			//Lyra[i].access = 0;
			layer_stored++;
		}
		else if(Lyra[i].stored){}
			//Lyra[i].access = 0;
		else
		{
			Lyra[Lyra[i].extra_row].access  = Lyra[i].access = Lyra[i-1].access + Lyra[Lyra[i].extra_row].access +1;
		}
		layer_access += Lyra[i].access;

		if((i+1)%(width/q)==0)//Start of layer
		{
			prev_compute = (float)layer_compute/(width/q);
			layer_access=layer_compute = layer_stored =  0;
		}
	}*/
	//double log_penalty = log((double)total_compute/width)/log((double)2);
//	printf("Log-Total compute: %2.2f Log-Penalty: %2.2f Memory factor:%d\n", log((double)total_compute)/log((double)2),log((double)total_compute/width)/log((double)2),q/2);
	return (double)total_compute/width;

}
void LyraStat(unsigned width, unsigned q)
{
	//Fill
	MTRand mtr;
	RandomGraph Lyra(width);
	vector<vector<unsigned> > columns; 
	columns.resize(width);
	for(unsigned i=0; i< width; ++i)
	{
		unsigned extra_row = RAND32%width;
		if(i%q==0)
			Lyra.rows[i] = row(true, extra_row);
		else
			Lyra.rows[i] = row(false, extra_row);
	}

	//Test
	printf("Width %d Step %d\n", width, q);
	//By layer elements
		vector<float> average_compute;
		
		average_compute.resize(q,0);
	for(unsigned l=0; l<width; l+=width/q) //By layers
	{
		unsigned current_layer = l/(width/q);
		vector<unsigned> ref_main;
		vector<unsigned> ref_extra;
		ref_main.resize(current_layer+1,0);
		ref_extra.resize(current_layer+1,0);

		for(unsigned j=0; j< width/q; ++j)
		{
			unsigned index = l+j;
			unsigned extra_index = Lyra.rows[index].extra_row;

			int index_layer=-1; //Layer when index was updated
			int extra_layer=-1; //Layer when extra_row was updated

			for(unsigned k=0; k<index; k++)
			{
				if(Lyra.rows[k].extra_row == index)
					index_layer = k/(width/q);
				if(Lyra.rows[k].extra_row == extra_index)
					extra_layer = k/(width/q);
				if(k==extra_index)
					extra_layer = k/(width/q);
			}
			if(index_layer >=0)
				ref_main[index_layer]++;
			if(extra_layer >=0)
				ref_extra[extra_layer]++;
		}
		printf("\nLayer %d: \n",current_layer);
		for(unsigned s=0; s<=current_layer; ++s)
		{
			float portion1 = (float)ref_main[s]/(width/q);
			float portion2 = (float)ref_extra[s]/(width/q);
			printf("Index updated at layer %d: %2.3f\n", s,portion1);
			printf("Extra row updated at layer %d: %2.3f\n", s,portion2);
			if(s!=current_layer)
				average_compute[current_layer] += (portion1+portion2)*average_compute[s]*q/4;
			else
				average_compute[current_layer] += 1;
		}
		printf("Layer %d complexity: %2.2f\n", current_layer,average_compute[current_layer]);

		
		
	}

}

vector<tradeoff> TestSetup()
{

	vector<tradeoff> tradeoffs;
	unsigned depth=4;
	for(unsigned step = 2; step<=24; step++)
	{
		for(unsigned d=2; d<=10; d++)
		{
			tradeoffs.push_back(SetupTradeoff(step,d));
		}
	}
	std::sort(tradeoffs.begin(), tradeoffs.end());

	/*for(int i=0; i<tradeoffs.size(); ++i)
	{
		printf(" Penalty for fraction %2.2f: %2.2f  Latency: %2.2f\n", tradeoffs[i].save, tradeoffs[i].penalty, tradeoffs[i].latency);
	}*/
	
	
	vector<tradeoff> filtered;
	filtered.push_back(tradeoffs[tradeoffs.size()-1]);
	unsigned filter_index=0;
	for(unsigned i=tradeoffs.size()-1; i>0; i--)
	{
		if(filtered[filter_index].penalty>tradeoffs[i-1].penalty)
		{
			filtered.push_back(tradeoffs[i-1]);
			filter_index++;
			printf("Save %2.2f Penalty %2.2f Latency: %2.2f\n", tradeoffs[i-1].save, tradeoffs[i-1].penalty, tradeoffs[i-1].latency);
		}
	}
	return filtered;
}

vector<tradeoff> TestWanderInt(unsigned q_setup)
{
	unsigned tests=100;
	vector<tradeoff> tradeoffs;
	for(unsigned t=0; t<tests; ++t)
	{
		for(unsigned q=2; q<=32; q++)
			tradeoffs.push_back(LyraWandInt(1<<12,q,q_setup));
	}

	std::sort(tradeoffs.begin(), tradeoffs.end());

	vector<float> averages(100);
	vector<float> av_latency(100);
	vector<float> av_intervals(100);
	vector<unsigned> counts(100,0);

	for(int i=0; i<tradeoffs.size(); ++i)
	{
		int close = floor(tradeoffs[i].save);
		if(tradeoffs[i].save-close <0.1)
		{
			averages[close]+= tradeoffs[i].penalty;
			av_latency[close]+=tradeoffs[i].latency;
			av_intervals[close]+=tradeoffs[i].intervals;
			counts[close]++;
		}
		else if (close+1 - tradeoffs[i].save <0.1)
		{
			averages[close+1]+= tradeoffs[i].penalty;
			av_latency[close+1]+=tradeoffs[i].latency;
			av_intervals[close+1]+=tradeoffs[i].intervals;
			counts[close+1]++;
		}
	}
	vector<tradeoff> results;
	for(unsigned i=1; i<100; ++i)
	{
		if(counts[i]>0)
		{
			printf("Average Penalty for fraction %d (%d close values): %2.2f  Latency: %2.2f %d-intervals: %2.2f\n", i, counts[i], averages[i]/counts[i], av_latency[i]/counts[i],
				q_setup, av_intervals[i]/counts[i]);
			fprintf(fp,"%2.2f, ",av_intervals[i]/counts[i]);
			//results.push_back(tradeoff(i,averages[i]/counts[i], av_latency[i]/counts[i]));
		}
	}
	fprintf(fp,"\n");
	
	
	

	return results;
}

vector<tradeoff> TestEmulate(tradeoff TradeoffFunc(unsigned,unsigned))//Computing average tradeoff for <TradeoffFunc> method with one parameter
{
	unsigned tests=100;
	vector<tradeoff> tradeoffs;
	for(unsigned t=0; t<tests; ++t)
	{
		for(unsigned q=2; q<=80; q++)
			tradeoffs.push_back(TradeoffFunc(q,WIDTH_TEST/(2*q)));
	}

	std::sort(tradeoffs.begin(), tradeoffs.end());

	vector<float> averages(100,0);
	vector<float> av_latency(100,0);
	vector<float> av_read(100,0);
	vector<unsigned> counts(100,0);

	for(int i=0; i<tradeoffs.size(); ++i)
	{
		int close = floor(tradeoffs[i].save);
		if(tradeoffs[i].save-close <0.1)
		{
			averages[close]+= tradeoffs[i].penalty;
			av_latency[close]+=tradeoffs[i].latency;
			av_read[close]+=tradeoffs[i].read;
			counts[close]++;
		}
		else if (close+1 - tradeoffs[i].save <0.1)
		{
			averages[close+1]+= tradeoffs[i].penalty;
			av_latency[close+1]+=tradeoffs[i].latency;
			av_read[close+1]+=tradeoffs[i].read;
			counts[close+1]++;
		}
	}
	vector<tradeoff> results;
	FILE* fp = fopen("ranking-lyra.log","w+");
	for(unsigned i=1; i<100; ++i)
	{
		if(counts[i]>0)
		{
			fprintf(fp,"Average Penalty for fraction %d (%d close values): %2.2f  Latency: %2.2f Read: %2.2f\n", i, counts[i], averages[i]/counts[i], av_latency[i]/counts[i],
				 av_read[i]/counts[i]);
			results.push_back(tradeoff(i,averages[i]/counts[i], av_latency[i]/counts[i],0,av_read[i]/counts[i]));
		}
	}
	fclose(fp);
	
	

	return results;
}

vector<tradeoff> TestEmulatePass(tradeoff TradeoffFunc(unsigned, unsigned, unsigned), unsigned pass)
{
	unsigned tests=100;
	vector<tradeoff> tradeoffs;
	for(unsigned t=0; t<tests; ++t)
	{
		for(unsigned q=2; q<=80; q++)
			tradeoffs.push_back(TradeoffFunc(1<<12,q,pass));
	}

	std::sort(tradeoffs.begin(), tradeoffs.end());

	vector<float> averages(100,0);
	vector<float> av_latency(100,0);
	vector<float> av_read(100,0);
	vector<unsigned> counts(100,0);

	for(int i=0; i<tradeoffs.size(); ++i)
	{
		int close = floor(tradeoffs[i].save/pass);
		if(tradeoffs[i].save-close <0.1)
		{
			averages[close]+= tradeoffs[i].penalty;
			av_latency[close]+=tradeoffs[i].latency;
			av_read[close]+=tradeoffs[i].read;
			counts[close]++;
		}
		else if (close+1 - tradeoffs[i].save <0.1)
		{
			averages[close+1]+= tradeoffs[i].penalty;
			av_latency[close+1]+=tradeoffs[i].latency;
			av_read[close+1]+=tradeoffs[i].read;
			counts[close+1]++;
		}
	}
	vector<tradeoff> results;
	for(unsigned i=1; i<100; ++i)
	{
		if(counts[i]>0)
		{
			printf("Average Penalty for fraction %d (%d close values): %2.2f  Latency: %2.2f Read: %2.2f\n", i, counts[i], averages[i]/counts[i], av_latency[i]/counts[i],
				 av_read[i]/counts[i]);
			results.push_back(tradeoff(i,averages[i]/counts[i], av_latency[i]/counts[i],0,av_read[i]/counts[i]));
		}
	}
	
	
	

	return results;
}

vector<tradeoff> TestWander()
{
	unsigned tests=100;
	vector<tradeoff> tradeoffs;
	for(unsigned t=0; t<tests; ++t)
	{
		for(unsigned q=4; q<=16; q++)
			tradeoffs.push_back(LyraCompute(1<<3,q));
	}

	std::sort(tradeoffs.begin(), tradeoffs.end());

	vector<float> averages(100);
	vector<float> av_latency(100);
	vector<unsigned> counts(100,0);

	for(int i=0; i<tradeoffs.size(); ++i)
	{
		int close = floor(tradeoffs[i].save);
		if(tradeoffs[i].save-close <0.1)
		{
			averages[close]+= tradeoffs[i].penalty;
			av_latency[close]+=tradeoffs[i].latency;
			counts[close]++;
		}
		else if (close+1 - tradeoffs[i].save <0.1)
		{
			averages[close+1]+= tradeoffs[i].penalty;
			av_latency[close+1]+=tradeoffs[i].latency;
			counts[close+1]++;
		}
	}
	vector<tradeoff> results;
	for(unsigned i=1; i<100; ++i)
	{
		if(counts[i]>0)
		{
			printf("Average Penalty for fraction %d (%d close values): %2.2f  Latency: %2.2f\n", i, counts[i], averages[i]/counts[i], av_latency[i]/counts[i]);
			results.push_back(tradeoff(i,averages[i]/counts[i], av_latency[i]/counts[i]));
		}
	}
	
	
	

	return results;
}

void GetPenaltyOld()
{
	vector<tradeoff> trade_setup = 	TestSetup();
	printf("Setup finished: %d tradeoffs\n", trade_setup.size());
	vector<tradeoff> trade_wander = 	TestWander();
	printf("Wandering finished: %d tradeoffs\n", trade_wander.size());

	vector<tradeoff> trade_tmp;

	for(unsigned i=0; i<trade_setup.size(); ++i)
		for(unsigned j=0; j<trade_wander.size(); ++j)
			trade_tmp.push_back(tradeoff(trade_setup[i].save*trade_wander[j].save/(trade_setup[i].save+trade_wander[j].save), trade_setup[i].penalty*trade_wander[j].penalty/2,
			trade_setup[i].latency+trade_wander[j].latency));

	vector<tradeoff> filtered;
	filtered.push_back(trade_tmp[trade_tmp.size()-1]);
	unsigned filter_index=0;
	for(unsigned i=trade_tmp.size()-1; i>0; i--)
	{
		if(filtered[filter_index].penalty>trade_tmp[i-1].penalty)
		{
			filtered.push_back(trade_tmp[i-1]);
			filter_index++;
			printf("FINAL: Save %2.2f Penalty %2.2f Latency: %2.2f\n", trade_tmp[i-1].save, trade_tmp[i-1].penalty, trade_tmp[i-1].latency);
		}
	}
}

void list_simple_penalty()
{
	for(unsigned q=2; q<40; q++)
	{
		tradeoff t = SimplePenalty(q);
		printf("Q %d Reduction %2.2f, Penalty %2.2f\n",q,t.save,t.penalty);
	}
}

void list_simple_emul_penalty()
{
	for(unsigned q=4; q<40; q++)
	{
		tradeoff t = SimpleRandom(1<<12,q);
		printf("Q %d Reduction %2.2f, Penalty %2.2f\n",q,t.save,t.penalty);
	}
}


void test_wt()
{
	MTRand mtr;
	for(unsigned i=0; i<100; ++i)
	{
		uint64_t a = (((uint64_t)RAND32)<<32)+ RAND32;
		unsigned w1 = wt(a);
		unsigned w2 = 0;
		while(a!=0)
		{
			if(a&((uint64_t)1))
				w2++;
			a>>=1;
		}
		if(w1!=w2)
			printf("ERROR\n");

	}
	printf("WT tested\n");
}

int _tmain(int argc, _TCHAR* argv[])
{
	//`test_wt();
	/*fp = fopen("total.log","w+");
	TestWanderInt(4);
	TestWanderInt(5);
	TestWanderInt(6);
	TestWanderInt(7);
	TestWanderInt(10);
	TestWanderInt(9);
	TestWanderInt(11);
	TestWanderInt(13);
	TestWanderInt(12);
	TestWanderInt(14);
	TestWanderInt(15);
	TestWanderInt(17);
	TestWanderInt(18);
	TestWanderInt(20);
	TestWanderInt(22);
	TestWanderInt(20);
	TestWanderInt(22);
	TestWanderInt(23);
	TestWanderInt(25);
	TestWanderInt(26);
	TestWanderInt(28);
	TestWanderInt(29);
	TestWanderInt(31);
	TestWanderInt(32);
	TestWanderInt(30);
	TestWanderInt(32);
	TestWanderInt(31);
	TestWanderInt(32);
	TestWanderInt(32);
	fclose(fp);*/

	//list_simple_penalty();
	//list_simple_emul_penalty();
	//TestEmulatePass(AdvancedMultiRandom,2);
	//TestEmulate(RankingRandom);
	TestEmulate(RankingLyra);
	//TestWander();
	//TestSetup();
	return 0;

}

