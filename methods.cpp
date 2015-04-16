
#include <vector>
#include <algorithm>
#include <stdint.h>
using namespace std;


#include "./MersenneTwister.h"
#define WIDTH_TEST 4000

#define START_LATENCY 0
#define RAND32 (mtr.randInt())

//Fill
MTRand mtr;

/**
* p2floor(x):
* Largest power of 2 not greater than argument.
*/
static uint32_t
p2floor(uint32_t x)
{
	uint32_t y;
	while ((y = x & (x - 1)))
		x = y;
	return x;
}

/**
* wrap(x, i):
* Wrap x to the range 0 to i-1.
*/
static uint32_t
wrap(uint32_t x, uint32_t i)
{
	uint32_t n = p2floor(i);
	return (x & (n - 1)) + (i - n);
}


struct tradeoff
{
	float save;    //Memory reduction
	float penalty; //Computational penalty
	float latency; //Depth penalty
	float intervals;
	float read;
	tradeoff(float s = 1, float p = 1, float l = 1, float i = 1, float r = 1){ save = s; penalty = p; latency = l; intervals = i; read = r; };
	bool operator<(const tradeoff& r)  //Tradeoff comparison operator
	{
		if (save < r.save)
			return true;
		else if (save == r.save)
		{
			if (penalty < r.penalty)
				return true;
			else /*if(penalty == r.penalty)
				 return true;
				 else */return false;
		}
		else return false;
	}
};

struct node{

	uint64_t latency;
	bool stored;
	unsigned extra_node;
	uint64_t access;
	uint64_t read;
	node(bool b = false, unsigned e = 0, uint64_t a = 0)
	{
		stored = b; extra_node = e; access = a; latency = START_LATENCY; read = 1;
	}
};

struct RandomGraph
{
	unsigned width;
	vector<node> nodes;
	RandomGraph(unsigned w){ nodes.resize(w); width = w; visited.resize(w, false); };
	vector<bool> visited;

	void clean(){ visited.assign(width, false); }
	node& operator[](unsigned i){ return nodes[i%width]; }
	void Traverse(unsigned position, unsigned time);
	unsigned count(){
		unsigned a = 0;
		for (unsigned j = 0; j<width; ++j)
			if (visited[j])
				a++;
		return a;
	};
};


tradeoff RankingLyra(unsigned q, unsigned top_size)  //Ranking method for Lyra/Wandering graph, where F is stored.  <q> is the segment length, and <top_size> is number of highest complexities stored
{
	//Fill
	unsigned input_cost = 0;
	uint32_t width = WIDTH_TEST;
	RandomGraph Graph(width);
	unsigned extra_node = 1;
	for (unsigned i = 0; i< width; ++i)
	{
		extra_node = (i == 0) ? 0 : RAND32%i;
		Graph[i] = node(false, extra_node, input_cost);
	}

	uint64_t layer_access = 0;
	uint64_t layer_compute = 0;
	uint64_t layer_stored = 0;
	uint64_t total_stored = 0;
	uint64_t total_compute = 0;
	uint64_t total_latency = 0;
	uint64_t total_read = 0;

	//Test
	float prev_access = 1;
	float prev_compute = 1;
	unsigned layer_length = width / q;
	vector<uint64_t> top;  //The largest access cost values, sorted, top[0] is the highest
	top.resize(top_size, input_cost);
	for (unsigned i = 1; i< width; ++i)
	{
		extra_node = Graph[i].extra_node;
		total_compute += Graph[i].access + Graph[extra_node].access + 1;
		total_latency += Graph[i].latency + Graph[extra_node].latency + 1;
		total_read += Graph[i].read + Graph[extra_node].read;

		if (Graph[extra_node].access> top[top_size - 1]) //Among the most expensive nodes => need to store output of F
		{
			Graph[i].stored = true;
			//Access cost remains unchanged for both nodes
			total_stored++;
			//Latency  remains unchanged for both nodes
			Graph[i].read++;
			Graph[extra_node].read++;
		}
		else
		{
			if (i%q == 0)//Segment store rule
			{
				Graph[i].stored = true;
				//Access cost remains unchanged for both nodes
				total_stored++;
				//Latency  remains unchanged for both nodes
				Graph[i].read++;
				Graph[extra_node].read++;
			}
			else
			{
				Graph[extra_node].latency = Graph[i].latency = max(Graph[extra_node].latency, max(Graph[i - 1].latency, Graph[i].latency)) + 1; //both nodes have new latency
				Graph[extra_node].access = Graph[i].access = Graph[i].access + Graph[i - 1].access + Graph[extra_node].access + 1; //both nodes have new access cost
				Graph[extra_node].read = Graph[i].read = Graph[i - 1].read + Graph[i].read + Graph[extra_node].read; ////both nodes have new red values
				unsigned rank = top_size;  //Check if the new access complexity is in the top
				while (rank>0)//adding new value to the top
				{
					if (top[rank - 1]<Graph[i].access)
						rank--;
					else break;
				}
				top.insert(top.begin() + rank, Graph[i].access);
				top.pop_back();
			}
		}
	}
	float save = (float)width / total_stored + 1 / (float)96;
	float penalty = (double)total_compute / width;
	float latency = (double)total_latency / width;
	float read = (double)total_read / width;
	return tradeoff(save, penalty, latency, 0, read);

}

tradeoff RankingYescrypt(unsigned q, unsigned top_size)  //Ranking method with parameter <q> for yescrypt graph, initial pass
{
	//Fill
	unsigned input_cost = 0;
	uint32_t width = WIDTH_TEST;
	RandomGraph Graph(width);
	vector<vector<unsigned> > columns;  //Alternative storage
	columns.resize(width);
	unsigned extra_node = 1;
	for (unsigned i = 0; i< width; ++i)
	{
		extra_node = (i == 0) ? 0 : RAND32;
		if(i!=0)
			extra_node = wrap(extra_node, i);
		columns[extra_node].push_back(i);
		if (i%q == 0)
		{
			Graph[i] = node(true, extra_node);//supposed to store
		}
		else
		{
			Graph[i] = node(false, extra_node);

		}
		Graph[i].access = input_cost;
	}

	uint64_t total_stored = 0;
	uint64_t total_compute = 1; //First element skipped
	uint64_t total_latency = 1;
	uint64_t total_read = 0;

	//Test
	float prev_access = 1;
	float prev_compute = 1;
	vector<uint64_t> top;  //The largest access cost values, sorted, top[0] is the highest
	top.resize(top_size, input_cost);
	for (unsigned i = 1; i< width; ++i)
	{
		extra_node = Graph[i].extra_node;
		total_compute += Graph[extra_node].access + 1;
		total_latency += Graph[extra_node].latency + 1;
		total_read += Graph[extra_node].read;

		if ((Graph[extra_node].access> top[top_size - 1]) ||( Graph[i].stored)) //Among the most expensive nodes or every q-th => need to store 
		{
			Graph[i].stored = true;
			Graph[i].access = 0;
			total_stored++;
			Graph[i].latency = 0;
			Graph[i].read = 1;
		}
		else
		{
			Graph[i].latency = max(Graph[extra_node].latency, Graph[i - 1].latency) + 1;
			Graph[i].access = Graph[i - 1].access + Graph[extra_node].access + 1;
			Graph[i].read = Graph[i - 1].read + Graph[extra_node].read;
			unsigned rank = top_size;  //Check if the new access complexity is in the top
			while (rank>0)//adding new value to the top
			{
				if (top[rank - 1]<Graph[i].access)
					rank--;
				else break;
			}
			top.insert(top.begin() + rank, Graph[i].access);
			top.pop_back();
		}
	}
	float save = (float)width / total_stored;
	float penalty = (float)total_compute / width;
	float latency = (float)total_latency / width;
	float read = (float)total_read / width;
	return tradeoff(save, penalty, latency, 0, read);

}

tradeoff RankingYescryptLatency(unsigned q, unsigned top_size)  
//Ranking method with parameter <q> for yescrypt graph, initial pass. Top latency values are stored, not top access values
{
	//Fill
	unsigned input_cost = 0;
	uint32_t width = WIDTH_TEST;
	RandomGraph Graph(width);
	vector<vector<unsigned> > columns;  //Alternative storage
	columns.resize(width);
	unsigned extra_node = 1;
	for (unsigned i = 0; i< width; ++i)
	{
		extra_node = (i == 0) ? 0 : RAND32;
		if (i != 0)
			extra_node = wrap(extra_node, i);
		columns[extra_node].push_back(i);
		if (i%q == 0)
		{
			Graph[i] = node(true, extra_node);
		}
		else
		{
			Graph[i] = node(false, extra_node);

		}
		Graph[i].access = input_cost;
	}

	uint64_t total_stored = 0;
	uint64_t total_compute = 1; //First element skipped
	uint64_t total_latency = 1;
	uint64_t total_read = 0;

	//Test
	float prev_access = 1;
	float prev_compute = 1;
	unsigned layer_length = width / q;
	vector<uint64_t> top;  //The largest access cost values, sorted, top[0] is the highest
	top.resize(top_size, input_cost);
	for (unsigned i = 1; i< width; ++i)
	{
		extra_node = Graph[i].extra_node;
		total_compute += Graph[extra_node].access + 1;
		total_latency += Graph[extra_node].latency + 1;
		total_read += Graph[extra_node].read;

		if (Graph[extra_node].latency> top[top_size - 1]) //Among the most expensive nodes => need to store output of F
		{
			Graph[i].stored = true;
			Graph[i].access = 0;
			total_stored++;
			Graph[i].latency = 0;
			Graph[i].read = 1;
		}
		else
		{
			if (Graph[i].stored)
			{
				Graph[i].access = 0;
				total_stored++;
				Graph[i].latency = 0;
				Graph[i].read = 1;
			}
			else
			{
				Graph[i].latency = max(Graph[extra_node].latency, Graph[i - 1].latency) + 1;
				Graph[i].access = Graph[i - 1].access + Graph[extra_node].access + 1;
				Graph[i].read = Graph[i - 1].read + Graph[extra_node].read;
				unsigned rank = top_size;  //Check if the new access complexity is in the top
				while (rank>0)//adding new value to the top
				{
					if (top[rank - 1]<Graph[i].latency)
						rank--;
					else break;
				}
				top.insert(top.begin() + rank, Graph[i].latency);
				top.pop_back();
			}
		}
	}
	float save = (float)width / total_stored;
	float penalty = (double)total_compute / width;
	float latency = (double)total_latency / width;
	float read = (double)total_read / width;
	return tradeoff(save, penalty, latency, 0, read);

}


tradeoff RankingRandom(unsigned q, unsigned top_size)  //Ranking method with parameter <q> for random graph
{
	//Fill
	unsigned input_cost = 0;
	uint32_t width = WIDTH_TEST;
	RandomGraph Graph(width);
	vector<vector<unsigned> > columns;  //Alternative storage
	columns.resize(width);
	unsigned extra_node = 1;
	for (unsigned i = 0; i< width; ++i)
	{
		extra_node = (i == 0) ? 0 : RAND32%i;
		columns[extra_node].push_back(i);
		if (i%q == 0)
		{
			Graph[i] = node(true, extra_node);
		}
		else
		{
			Graph[i] = node(false, extra_node);

		}
		Graph[i].access = input_cost;
	}

	uint64_t total_stored = 0;
	uint64_t total_compute = 1; //First element skipped
	uint64_t total_latency = 1;
	uint64_t total_read = 0;

	//Test
	float prev_access = 1;
	float prev_compute = 1;
	unsigned layer_length = width / q;
	vector<uint64_t> top;  //The largest access cost values, sorted, top[0] is the highest
	top.resize(top_size, input_cost);
	for (unsigned i = 1; i< width; ++i)
	{
		extra_node = Graph[i].extra_node;
		total_compute += Graph[extra_node].access + 1;
		total_latency += Graph[extra_node].latency + 1;
		total_read += Graph[extra_node].read;

		if (Graph[extra_node].access> top[top_size - 1]) //Among the most expensive nodes => need to store output of F
		{
			Graph[i].stored = true;
			Graph[i].access = 0;
			total_stored++;
			Graph[i].latency = 0;
			Graph[i].read = 1;
		}
		else
		{
			if (Graph[i].stored)
			{
				Graph[i].access = 0;
				total_stored++;
				Graph[i].latency = 0;
				Graph[i].read = 1;
			}
			else
			{
				Graph[i].latency = max(Graph[extra_node].latency, Graph[i - 1].latency) + 1;
				Graph[i].access = Graph[i - 1].access + Graph[extra_node].access + 1;
				Graph[i].read = Graph[i - 1].read + Graph[extra_node].read;
				unsigned rank = top_size;  //Check if the new access complexity is in the top
				while (rank>0)//adding new value to the top
				{
					if (top[rank - 1]<Graph[i].access)
						rank--;
					else break;
				}
				top.insert(top.begin() + rank, Graph[i].access);
				top.pop_back();
			}
		}
	}
	float save = (float)width / total_stored;
	float penalty = (double)total_compute / (width - 1);
	float latency = (double)total_latency / (width - 1);
	float read = (double)total_read / (width - 1);
	return tradeoff(save, penalty, latency, 0, read);

}

tradeoff AdvancedMultiRandom(unsigned width, unsigned q, unsigned k)  //Attacking random k-pass graph with advanced technique
{
	//Fill
	unsigned input_cost = 0;
	RandomGraph Graph(width);
	vector<vector<unsigned> > columns;  //Alternative storage
	columns.resize(width);
	unsigned extra_node = 1;
	for (unsigned i = 0; i< width; ++i)
	{
		if (i == 0)
			extra_node = 0;
		else if (i< width / k)
			extra_node = RAND32%i;
		else extra_node = i - RAND32 % (width / k);
		columns[extra_node].push_back(i);
		if (i%q == 0)
		{
			Graph[i] = node(true, extra_node);
		}
		else
		{
			Graph[i] = node(false, extra_node);

		}
		Graph[i].access = input_cost;
	}

	uint64_t layer_access = 0;
	uint64_t layer_compute = 0;
	uint64_t layer_stored = 0;
	uint64_t total_stored = 0;
	uint64_t total_compute = 0;
	uint64_t total_latency = 0;
	uint64_t total_read = 0;

	//Test
	float prev_access = 1;
	float prev_compute = 1;
	unsigned layer_length = width / q;
	unsigned top_size = layer_length;
	vector<uint64_t> top;  //The largest access cost values, sorted, top[0] is the highest
	top.resize(top_size, input_cost);
	for (unsigned i = 1; i< width; ++i)
	{
		extra_node = Graph[i].extra_node;
		total_compute += Graph[extra_node].access + 1;
		total_latency += Graph[extra_node].latency + 1;
		total_read += Graph[extra_node].read;

		if (Graph[extra_node].access> top[top_size - 1]) //Among the most expensive nodes => need to store output of F
		{
			if ((i< width / k) || (!Graph[i].stored))//first pass or the same element was not stored at the first pass
			{
				Graph[i].stored = true;
				total_stored++;
			}

			Graph[i].access = 0;
			Graph[i].latency = 0;
			Graph[i].read = 1;
		}
		else
		{
			if (i%q == 0) //Store always, do not set the .stored field to true
			{
				Graph[i].access = 0;
				if (i< width / k)//first pass
					total_stored++;
				Graph[i].latency = 0;
				Graph[i].read = 1;
			}
			else
			{
				if ((i> width / k) && (Graph[i].stored))//do not store it anymore
				{
					Graph[i].stored = false;
					total_stored--;
				}
				Graph[i].latency = max(Graph[extra_node].latency, Graph[i - 1].latency) + 1;
				Graph[i].access = Graph[i - 1].access + Graph[extra_node].access + 1;
				Graph[i].read = Graph[i - 1].read + Graph[extra_node].read;
				unsigned rank = top_size;  //Check if the new access complexity is in the top
				while (rank>0)//adding new value to the top
				{
					if (top[rank - 1]<Graph[i].access)
						rank--;
					else break;
				}
				top.insert(top.begin() + rank, Graph[i].access);
				top.pop_back();
			}
		}
	}
	float save = (float)width / total_stored + 1 / (float)32;
	float penalty = (double)total_compute / (width - 1);
	float latency = (double)total_latency / (width - 1);
	float read = (double)total_read / (width - 1);
	//printf("Memory saving %2.2f ", (float) width/total_stored);
	//double log_penalty = log((double)total_compute/width)/log((double)2);
	//	printf("Log-Total compute: %2.2f Log-Penalty: %2.2f Memory factor:%d\n", log((double)total_compute)/log((double)2),log((double)total_compute/width)/log((double)2),q/2);
	return tradeoff(save, penalty, latency, 0, read);

}


tradeoff RankingYescryptMulti(unsigned q, unsigned top_size, float k)  //Attacking random k-pass graph with advanced technique
{
	unsigned input_cost = 0;
	uint32_t width = WIDTH_TEST;
	RandomGraph Graph(width);
	vector<vector<unsigned> > columns;  //Alternative storage
	columns.resize(width);
	unsigned extra_node = 1;
	for (unsigned i = 0; i< width; ++i)
	{
		extra_node = (i == 0) ? 0 : RAND32;
		if( (i != 0) && (i< width / k))
			extra_node = wrap(extra_node, i);
		else extra_node %= (int)(width / k);
		columns[extra_node].push_back(i);
		if (i%q == 0)
		{
			Graph[i] = node(true, extra_node); //supposed to store
		}
		else
		{
			Graph[i] = node(false, extra_node);

		}
		Graph[i].access = input_cost;
	}

	uint64_t total_stored = 0;
	uint64_t total_compute = 1;//First element skipped
	uint64_t total_latency = 1;
	uint64_t total_read = 0;

	//Test
	float prev_access = 1;
	float prev_compute = 1;
	vector<uint64_t> top;  //The largest access cost values, sorted, top[0] is the highest
	top.resize(top_size, input_cost);
	for (unsigned i = 1; i< width; ++i)
		//Strategy at the second pass: overwritten blocks are considered as new ones. However, for calculation simplicity, we overwrite the access complexity and 
	{
		extra_node = Graph[i].extra_node;
		total_compute += Graph[extra_node].access + 1;
		total_latency += Graph[extra_node].latency + 1;
		total_read += Graph[extra_node].read;

		if ((Graph[extra_node].access> top[top_size - 1]) || Graph[i].stored) //Among the most expensive nodes or every q-th block => need to store it
		{
			if ((i < width / k))//first pass 
			{
				total_stored++;
				Graph[i].access = 0;
				Graph[i].latency = 0;
				Graph[i].read = 1;
			}
			else  //second pass. Overwrite the access complexity of extra node
			{
				total_stored++;
				Graph[extra_node].access = Graph[i].access = 0;
				Graph[extra_node].latency = Graph[i].access = 0;
				Graph[extra_node].read = Graph[i].access = 1;
			}
		}
		else
		{
			Graph[i].latency = max(Graph[extra_node].latency, Graph[i - 1].latency) + 1;
			Graph[i].access = Graph[i - 1].access + Graph[extra_node].access + 1;
			Graph[i].read = Graph[i - 1].read + Graph[extra_node].read;
			unsigned rank = top_size;  //Check if the new access complexity is in the top
			while (rank>0)//adding new value to the top
			{
				if (top[rank - 1]<Graph[i].access)
					rank--;
				else break;
			}
			top.insert(top.begin() + rank, Graph[i].access);
			top.pop_back();
		}
		
	}
	total_stored += (width - total_stored) / 16; //storing output of Salsa20
	float save = ((float)width/k) / total_stored;
	float penalty = (double)total_compute / (width);
	float latency = (double)total_latency / (width);
	float read = (double)total_read / (width);
	return tradeoff(save, penalty, latency, 0, read);

}

vector<tradeoff> OutputPenalties(FILE* fp,vector<tradeoff> tradeoffs)
{
	std::sort(tradeoffs.begin(), tradeoffs.end());

	vector<float> averages(100, 0);
	vector<float> av_latency(100, 0);
	vector<float> av_read(100, 0);
	vector<unsigned> counts(100, 0);

	for (int i = 0; i<tradeoffs.size(); ++i)
	{
		int close = floor(tradeoffs[i].save);
		if (tradeoffs[i].save - close <0.1)
		{
			averages[close] += tradeoffs[i].penalty;
			av_latency[close] += tradeoffs[i].latency;
			av_read[close] += tradeoffs[i].read;
			counts[close]++;
		}
		else if (close + 1 - tradeoffs[i].save <0.1)
		{
			averages[close + 1] += tradeoffs[i].penalty;
			av_latency[close + 1] += tradeoffs[i].latency;
			av_read[close + 1] += tradeoffs[i].read;
			counts[close + 1]++;
		}
	}
	vector<tradeoff> results;
	
	for (unsigned i = 1; i<100; ++i)
	{
		if (counts[i]>0)
		{
			fprintf(fp, "Average Penalty for fraction %d (%d close values): %2.2f  Latency: %2.2f Read: %2.2f\n", i, counts[i], averages[i] / counts[i], av_latency[i] / counts[i],
				av_read[i] / counts[i]);
			results.push_back(tradeoff(i, averages[i] / counts[i], av_latency[i] / counts[i], 0, av_read[i] / counts[i]));
		}
	}



	return results;
}

vector<tradeoff> ComputePenalties(tradeoff TradeoffFunc(unsigned, unsigned))//Computing average tradeoff for <TradeoffFunc> method called with (memory reduction,
{
	unsigned tests = 100;
	vector<tradeoff> tradeoffs;
	for (unsigned t = 0; t<tests; ++t)
	{
		for (unsigned q = 2; q <= 80; q++)
			tradeoffs.push_back(TradeoffFunc(q, WIDTH_TEST / (3*q)));
	}
	FILE* fp = fopen("penalties-yes.log", "w+");
	tradeoffs =  OutputPenalties(fp,tradeoffs);
	fclose(fp);
	return tradeoffs;
}


vector<tradeoff> ComputePenaltiesMulti(tradeoff TradeoffFunc(unsigned, unsigned,float), float passes)//Computing average tradeoff for <TradeoffFunc> method called with (memory reduction,
{
	unsigned tests = 100;
	vector<tradeoff> tradeoffs;
	for (unsigned t = 0; t<tests; ++t)
	{
		for (unsigned q = 2; q <= 80; q++)
			tradeoffs.push_back(TradeoffFunc(q, WIDTH_TEST / (3 * q),passes));
	}

	FILE* fp = fopen("penalties-yes-multi.log", "w+");
	tradeoffs = OutputPenalties(fp, tradeoffs);
	fclose(fp);
	return tradeoffs;
}


int main(unsigned argc, void** argv)
{
	mtr.seed(time(0));
	//ComputePenalties(RankingLyra);
	//ComputePenalties(RankingRandom);
	//ComputePenalties(RankingYescryptLatency);
	vector<tradeoff> tradeoff1 = ComputePenaltiesMulti(RankingYescryptMulti, 1.33);
	vector<tradeoff> tradeoff2 = ComputePenalties(RankingYescrypt);
	return 0;
}