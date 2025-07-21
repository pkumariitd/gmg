#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <iterator>
#include <iomanip>
#include <ios>
#include <sstream>
#include <fstream>
#include <list>
#include <map>
#include <set>
#include <math.h>
#include <random>
#include <cmath>
#include <ctime>
#include <chrono>
using namespace std;
typedef pair<long, long> edge;
typedef set<long> group;
void rand_seed()
{
    int seed = static_cast<int>(time(0));
    srand(seed);
}
long rand(long, long);
float random_number(float , float );
void set_parameters(long&, long&, long&, bool&, float&, int, char*[], int& );
void set_fill(set<long>&, long, long); //fills the set with all the longegers in the given range
void print_set(set<long>&);
void print_map(map<int,map<int,float>>&);
void print_edge_vector(vector<edge>&);
void print_vector(vector<long>& );
void print_vector(vector<float>&);
map<int, map<int, float> > weight_in,weight_cross,weight;
void write_graph(string&, map<int, map<int, float>>&, bool, vector<long>& );
void write_modules(string&, vector<group>&);
void usage();
vector<long> generate_community_sizes(int, int, int, double);
vector<long> generate_ov(long, long, double);
void calculate_weights(map<int, map<int, float>>&, const vector<long>&, long);
void make_distribution(vector<long>&, vector<long>&);
int preferential_module_selection(vector<long>&);



int main(int argc, char* argv[]){
    std::chrono::time_point<std::chrono::system_clock> start_time, end_time;
    start_time = std::chrono::system_clock::now();  //time starts
	rand_seed();
	int k=0;                  //overlapping nodes
	long n_min = 3, l;
	long n_max = n_min;
	float mu=0, alpha = 2.5;
	bool weighted = false;
	if(argc < 3 || argc > 12)
		usage();
	set_parameters(l, n_min, n_max, weighted, mu, argc, argv, k);
	cout<<"Setting parameters as :"<<endl;
	cout<<"l = "<<l << ", nmin = "<<n_min<<", nmax = "<<n_max<<", mu = "<<mu<<", ov = "<<k<<", weighted = "<<weighted<<endl;
	vector<group> module(l);

	vector<long> module_size = generate_community_sizes(l, n_min, n_max, alpha);
    vector<long> ov_nodes = generate_ov(l, k, alpha);

	set_fill(module.at(0), 1, module_size.at(0));
	long last = module_size.at(0);
    long n = module_size.at(0);
    int curr_overlap = 0;
    for(int i = 1; i < l; i++)
	{
		    set_fill(module[i], last-ov_nodes[i-1]+1, last+module_size.at(i)-ov_nodes[i-1]);
            n +=module_size.at(i)-ov_nodes[i-1];
            last = last + module[i].size()-ov_nodes[i-1];
	}
	long e_in_max = 0;
	int j;
	for(int i = 0; i < l; ++i)
		e_in_max += module.at(i).size()*(module.at(i).size()-1)/2;
	long required_within_edges = long(e_in_max*(1-mu));
	long required_crossing_edges = long(e_in_max*(mu));
	set<long>::iterator si, sj;

	vector<long> deg(n+1);
    int i = 0;
    vector<long> dist;

    make_distribution(dist, module_size);

	while(i < required_within_edges)
	{
		int j = preferential_module_selection(dist);
		int u = rand(*module.at(j).begin(), *module.at(j).begin()+module.at(j).size()-1);
        int v;
        do
            v = rand(*module.at(j).begin(), *module.at(j).begin()+module.at(j).size()-1);
        while(v == u);

        if(weight[u].find(v) != weight[u].end() ||weight[v].find(u) != weight[v].end())
            continue;
        if(weighted == false)
            weight[u][v] = 1;
        else
         weight[u][v] = (1-mu)*random_number(5,10);
        deg.at(u)++;
		deg.at(v)++;
        i++;

	}

    int s=0;
	while(s < required_crossing_edges)
	{
	    int j = rand(0, l-1);
	    int k;
	    do
            k = rand(0, l-1);
        while(k == j);

        int u = rand(*module.at(j).begin(), *module.at(j).begin()+module.at(j).size()-1);
        int v = rand(*module.at(k).begin(), *module.at(k).begin()+module.at(k).size()-1);
        if(weight[u].find(v) != weight[u].end() || weight[v].find(u) != weight[v].end())
            continue;
         weight[u][v] = mu*random_number(0,5);
         deg.at(u)++;
		 deg.at(v)++;
         s++;
    }


	ostringstream str_mu;
    str_mu.setf(ios::fixed, ios::floatfield);
    str_mu.precision(2);
    str_mu<<mu;
	ostringstream graph_file;
	graph_file<<"./gmg_graph"<<"_mu"<<str_mu.str()<<".txt";
    string c_graph_file = graph_file.str();
    write_graph(c_graph_file, weight,weighted, deg);
	ostringstream module_file;
	module_file << "./gmg_modules"<<"_mu"<<str_mu.str()<<".txt";
	string c_module_file = module_file.str();
	write_modules(c_module_file, module);
	 end_time = std::chrono::system_clock::now();  //time ends
    std::chrono::duration<double> elapsed_seconds = end_time-start_time;
	cout<<"-------------------------------------------------------------"<<endl;
	cout<<setw(30)<<"Total number of vertices = "<<n<<endl;
    cout<<setw(30)<<"total number of edges = "<<e_in_max<<endl;
	cout<<setw(30)<<"intra-module edges = "<<required_within_edges<<endl;
	cout<<setw(30)<<"inter-module edges = "<<required_crossing_edges<<endl;
	cout<<setw(30)<<"Time elapsed(seconds) = "<<elapsed_seconds.count()<<endl;
	cout<<"-------------------------------------------------------------"<<endl;
	return 0;
}

void set_parameters(long& l, long& n_min, long& n_max, bool& weighted, float& mu, int argc, char* argv[], int &k)
{
	int i=1;
	while(i <= argc-1)
	{
      string arg = string(argv[i]);
		if(arg == "-l")
		{
			istringstream is(argv[i+1]);
			is>>l;
			if( l < 2)
				usage();
			i += 2;
		}
        else
			if(arg == "-nmin")
			{
				istringstream is(argv[i+1]);
				is>>n_min;
				if( n_min < 1)
					usage();
				i += 2;
			}
			else
				if(arg == "-nmax")
				{
					istringstream is(argv[i+1]);
					is>>n_max;
					if( n_max < n_min)
						usage();
					i += 2;
				}
				else
					if(arg == "-mu")
					{
						istringstream is(argv[i+1]);
						is>>mu;
						if( mu < 0 or mu > 1)
							usage();
						i += 2;
					}
                    else
                        if(arg == "-ov")
                        {
                            istringstream is(argv[i+1]);
                            is>>k;
                            if( k < 0 or k > n_max)
                                usage();
                            i += 2;
                        }
                        else
                            if(arg == "-w")
                            {
                                weighted = true;
                                i++;
                            }
                            else
                                usage();
	}
}


void write_graph(string& file, map<int, map<int, float>>& w, bool weighted, vector<long>& deg)
{
    ofstream fout(file);
    if (!fout.is_open())
    {
        cout << "Destination file for communities could not be opened." << endl;
        exit(1);
    }

    // Iterate over the nested map `w`
    for (const auto& it : w)
        for (const auto& jt : it.second)
        {
            fout << it.first << " " << jt.first; // Write the source and target
            if (weighted)
                fout << " " << jt.second; // Write the weight if `weighted` is true
            fout << endl;
        }
    // Handle nodes with degree 0
    for (int v = 1; v < deg.size(); ++v)
        if (deg.at(v) == 0)
        {
            fout << v << " " << v; // Write self-loop for isolated nodes
            if (weighted)
                fout << " " << 0; // Write weight 0 if `weighted` is true
            fout << endl;
        }
}


void usage()
{
    cout<<"Please, follow the syntax as given below:"<<endl;
    cout<<"-----------------------------------------------------------"<<endl;
    cout<<"./gmg -l [option] -nmin [option] -nmax [option] -mu [option] -ov [option] -w"<<endl;
    cout<<"-----------------------------------------------------------"<<endl;
    cout<<"Here, l is the number of groups in the network."<<endl;
    cout<<"nmin and nmax are the minimum and maximum sizes of group, respectively."<<endl;
    cout<<"mu (varying 0 to 1), mu = 0 means complete anti-community structure, mu = 1 means complete community structure"<<endl;
    cout<<"Use flag w if your network is weighted."<<endl;
    exit(1);
}
long rand(long a, long b)
{
    long r = a + rand()%(b-a+1);
    return r;
}

float random_number(float a, float b)
{
    float r = a + static_cast<float>(rand()) / (static_cast<float>(RAND_MAX) / (b - a));
    return r;
}

vector<long> generate_community_sizes(int l, int nmin, int nmax, double alpha)
 {
    vector<long> sizes;
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dist(0.0, 1.0);

    double a = 1.0 - alpha;

    for (int i = 0; i < l; ++i)
    {
        double u = dist(gen);
        double x = pow((pow(nmax, a) - pow(nmin, a)) * u + pow(nmin, a), 1/a);
        int size = max(nmin, min(nmax, static_cast<int>(round(x))));
        sizes.push_back(size);
    }
    return sizes;
 }

vector<long> generate_ov(long subsets, long on, double exponent)
{
    vector<long> values(subsets - 1);
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dist(0.0, 1.0);

    double totalWeight = 0;
    for (int i = 0; i < subsets - 1; ++i)
    {
        values[i] = pow(dist(gen), -1.0 / (exponent - 1));
        totalWeight += values[i];
    }

    for (long &value : values)
        value = static_cast<long>((value / totalWeight) * on);

    long currentSum = accumulate(values.begin(), values.end(), 0L);
    while (currentSum != on)
    {
        for (int i = 0; i < subsets - 1; ++i)
        {
            if (currentSum < on)
            {
                values[i]++;
                currentSum++;
            }
            else
                if (currentSum > on && values[i] > 0)
                {
                    values[i]--;
                    currentSum--;
                }
            if (currentSum == on)
                break;
        }
    }

    return values;
}

void set_fill(set<long>& s, long a, long b)
{
	for(long i = a; i <= b; ++i)
		s.insert(i);
}
void write_modules(string& module_file, vector<group>& module)
{
	ofstream fout(module_file);
   if(!fout.is_open())
    {
        cout<<"Destination file for communities could not be opened."<<endl;
        exit(1);
    }
	for(long i = 0; i < module.size(); ++i)
	{
		copy(module.at(i).begin(), module.at(i).end(), ostream_iterator<long>(fout, " "));
		fout<<endl;
	}
}
void calculate_weights(map<int, map<int, float>>& weight, const vector<long>& deg, long max_deg)
{
    for (auto mi = weight.begin(); mi != weight.end(); ++mi)
    {
          int u = mi->first;
          for (auto jt = mi->second.begin(); jt != mi->second.end(); ++jt)
           {
             int v = jt->first;
             float w = sqrt(deg.at(u) * deg.at(v)) / float(max_deg);
             jt->second *= w;
           }
    }
}

void make_distribution(vector<long>& dist, vector<long>& module_size)
{
    int l=module_size.size();
    for(int i=0;i<l;++i)
        for(int j=0;j<module_size[i];j++)
            dist.push_back(i);
    random_shuffle(dist.begin(),dist.end());

}

int preferential_module_selection(vector<long>& dist)
{
    int r=rand(0,dist.size()-1);
    return dist.at(r);
}
