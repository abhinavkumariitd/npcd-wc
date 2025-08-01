#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <sstream>
#include <string>
#include <list>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <iterator>
#include <set>
#include <fstream>
#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <ios>
#include <string>
#include <set>
#include <map>
#include <list>
#include <numeric>
#include<unordered_map>
using namespace std;
class Graph
    {
    public:
        set<int> V; //only vertices with positive degree are stored
        map<int, float> strength;
         map<int, set<int> > neighbours;
        map<int, map<int, float> > weight;
        map<int,float>Centnode,Centnode1;
        void read_edgelist(string&, bool, bool);
        inline int order(){ return V.size(); }
        int ecount();
        float get_weight(int, int);
        void print_graph();
        float proximity(int, int);
        float get_sn_nbrs_proximity(int, map<int, float>&);
        friend float proximity(Graph&, int, int);
        friend float overlap(Graph&, set<int>&, set<int>&);
        friend float overlapping_weighted_modularity(Graph&, map<int, set<int> >&, int, map<int, set<int> >& );
    };
//other function declarations
void merge_communities(Graph&, map<int, set<int> >&, map<int, set<int> >&, float, int);
void write_seeds(set<int>&, string&);
float overlap(Graph& g, set<int>&, set<int>&);
void usage(void);
void read_communities(string& comfile, map<int, set<int> >&, map<int, set<int> >&);
void print_communities(map<int, set<int> >&);
float overlapping_weighted_modularity(Graph&, map<int, set<int> >&, map<int, set<int> >& );
void print_set(set<int>&);
void print_vector(vector<int>& );
map<int, set<int> > coms;
map<int, set<int> > memberships;

void Graph::read_edgelist(string& edgefile, bool weighted, bool directed)
{
    ifstream fin;
    fin.open(edgefile);
    if(!fin.is_open())
    {
        cout<<"The file containing the edgelist could not be opened."<<endl;
        exit(1);
    }
    string line;
    while ( getline(fin, line ) )
    {
        istringstream is( line );
        int u, v;
        is>>u;
        is>>v;
        if(u != v and neighbours[u].find(v) == neighbours[u].end())
        {
            neighbours[u].insert(v);
            neighbours[v].insert(u);
        }
        if(u == v || weight[u].find(v) != weight[u].end())
            continue;
        float w;
        if(weighted == true)
        {
            if(is.str() == " " || is.eof())
            {   cout<<endl<<"The edge list has missing weights."<<endl;
                exit(1);
            }
            is>>w;
        }
        else
            w = 1;
        weight[u][v] = w;
        if(strength.find(u) == strength.end())
            strength[u] = w;
        else
            strength[u] = strength[u]+w;
        if(directed == false)
        {   weight[v][u] = w;
            if(strength.find(v) == strength.end())
                strength[v] = w;
            else
                strength[v] = strength[v]+w;
        }
        V.insert(u);
        V.insert(v);
    }

     set<int>::iterator si,sj,st;
    for(st=V.begin();st!=V.end();++st)
    {
        if(strength[*st]>0)
            {
              int i=*st;
              float s=0;
              for(si=neighbours[i].begin();si!=neighbours[i].end();++si)
                 {
                   for(sj=neighbours[i].begin();sj!=neighbours[i].end();++sj)
                    {
                       if(si!=sj&&neighbours[*si].find(*sj)!=neighbours[*si].end())
                       s=s+weight[*si][*sj];
                    }
                  }
              Centnode.insert({i,(strength[i]+(s/2))});
           }
     }
  /*    cout<<endl<<"centnode"<<endl;
   map<int,float>::iterator ih;
   for(ih=Centnode.begin();ih!=Centnode.end();++ih){
    cout<<ih->first<<"    "<<ih->second<<endl;
    }
    cout<<"centnode size"<<Centnode.size();*/
}

int Graph::ecount()
{
    map<int, map<int, float> >::iterator mi;
    int degree_sum = 0;
    for(mi = weight.begin(); mi != weight.end(); ++mi)
        degree_sum += mi->second.size();
    return degree_sum/2;
}

void print_set(set<int>& s)
{
    set<int>::iterator sitr1;
    for(sitr1 = s.begin(); sitr1 != s.end(); ++sitr1)
            cout<<*sitr1<<" ";
    cout<<endl;
}

float Graph::get_weight(int u, int v)
{
    if(weight[u].find(v) != weight[u].end())
        return weight[u][v];
    else
        return 0;
}

void Graph::print_graph()
{
    cout<<"No. of vertices: "<<order()<<endl;
    cout<<"No. of edges: "<<ecount()<<endl;
    cout<<"\n"<<"vertex"<<setw(8)<<"degree"<<endl;
    map<int, float>::iterator i;
    for(i = strength.begin(); i != strength.end(); ++i)
        cout<<i->first<<setw(8)<<i->second<<endl;
    cout<<"\n"<<"node"<<setw(10)<<"nbrs"<<setw(5)<<"weight"<<endl;
    map<int, map<int, float> >::iterator it;
    for(it=weight.begin(); it != weight.end(); ++it)
    {
        cout<<it->first<<setw(5)<<"--- "<<endl;
        map<int, float>::iterator jt;
        for(jt = it->second.begin(); jt != it->second.end(); ++jt)
            cout<<setw(10)<<jt->first<<setw(5)<<jt->second<<endl;
    }
    map<int,set<int>>::iterator ci;
    for(ci=neighbours.begin();ci!=neighbours.end();++ci)
    {
         cout<<ci->first<<"  ";
        print_set(neighbours[ci->first]);

    }

}

float Graph::proximity(int u, int v)
{
    float min_weight_sum = 0, common_nbrs_strength = 0, weight_sum = 0, W=0;
    map<int, float>::iterator i, j;
   if(weight[u].size() > weight[v].size()) //to ensure that next loop runs minimum times.
    {
        int temp = u;
        u = v;
        v = temp;
    }
    int summ=0;
    set <float> nbrs_weight;
    set<int> common_nbrs;
    set<int>::iterator si;
    for(i = weight[u].begin(); i != weight[u].end(); ++i)
    {
        nbrs_weight.insert(i->second);
        j = weight[v].find(i->first);
        if(j != weight[v].end())
        {
            min_weight_sum += min(i->second, j->second);
            weight_sum += i->second + j->second;
            common_nbrs_strength += strength[i->first];
            for(si = common_nbrs.begin(); si != common_nbrs.end(); ++si)
                weight_sum += 2*get_weight(i->first, *si);
            common_nbrs.insert(i->first);
           // cout<<"curr common vertex "<<i->first<<" its strength ="<<g.strength[i->first]<<endl;
        }
    }

  /*  for(i = weight[v].begin(); i != weight[v].end(); ++i)
    {
        nbrs_weight.insert(i->second);

    }
    W = *max_element(nbrs_weight.begin(),nbrs_weight.end());*/
   // cout<<endl<<"min weight sum = "<<min_weight_sum<<" weight sum = "<<weight_sum<<" common nbrs strength = "<<common_nbrs_strength<<endl;

    float prox;

    float A = weight_sum/(common_nbrs_strength);
    float B = float (common_nbrs.size()+1)/min(neighbours[u].size(), neighbours[v].size());
    if(common_nbrs_strength == 0)
        prox = get_weight(u,v)/min(strength[u],strength[v]);
    else
        prox =pow((A*B),(1-(get_weight(u,v)/min(strength[u],strength[v]))));


    return prox;
}

void merge_communities(Graph& g, map<int, set<int> >& coms, map<int, set<int> >& members, float given_max_ov, int minc)
{
    //merges communities with overlap >= given_max_ov
    map<int, set<int> >::iterator ci, cj;
    set<int>::iterator si, sj, sk;
    map<int, float>::iterator mi;
    set<int> labels;
    for(ci = coms.begin(); ci != coms.end(); ++ci)
        labels.insert(ci->first);
    do
    {
        //cout<<"all coms"<<endl; print_community(coms);
        si = labels.begin();
        int r = rand()%labels.size();
        for(int i = 1; i<=r; ++i)
            ++si;
        //cout<<"curr com"<<endl<<*si<<"->"; print_set(coms[*si]);
        set<int> neighboring_com_labels;
        for(sj = coms[*si].begin(); sj != coms[*si].end(); ++sj)
            for(mi = g.weight[*sj].begin(); mi != g.weight[*sj].end(); ++mi)
                if(coms[*si].find(mi->first) == coms[*si].end())
                {
                    //cout<<"members["<<*sk<<"]->"; print_set(members[*sk]);
                    for(sk = members[mi->first].begin(); sk != members[mi->first].end(); ++sk)
                        if(coms[*sk].size() >= coms[*si].size())  //neighboring coms that are bigger than coms[*si]
                            neighboring_com_labels.insert(*sk);
                }
        //cout<<"neighboring coms "<<endl; print_set(neighboring_com_labels);
        float max_ov = 0;
        set<int> high_overlapping_coms;
        for(sj = neighboring_com_labels.begin(); sj != neighboring_com_labels.end(); ++sj)
        {
            float curr_ov = overlap(g, coms[*si], coms[*sj]);
            //cout<<"curr ov = "<<curr_ov<<endl;
            if(coms[*si].size() < minc)
            {
                //cout<<"This comm has size smaller than min"<<endl;
                if(curr_ov > max_ov)
                {
                    max_ov = curr_ov;
                    high_overlapping_coms.clear();
                    high_overlapping_coms.insert(*sj);
                }
                else
                    if(curr_ov == max_ov)
                        high_overlapping_coms.insert(*sj);
            }
            else
            {
                if( curr_ov >= given_max_ov)
                   high_overlapping_coms.insert(*sj);
            }
        }
        /*cout<<"high overlapping coms"<<endl;
        for(sj = high_overlapping_coms.begin(); sj != high_overlapping_coms.end(); ++sj)
        {
            cout<<*sj<<"->";
            print_set(coms[*sj]);
        } */
        for(sj = high_overlapping_coms.begin(); sj != high_overlapping_coms.end(); ++sj)
        {   //cout<<"curr ov = "<<curr_ov<<endl;
            for(sk = coms[*si].begin(); sk != coms[*si].end(); ++sk)
            {
                coms[*sj].insert(*sk);
                members[*sk].erase(*si);
                members[*sk].insert(*sj);
            }
        }
        /* cout<<"modified coms"<<endl;
        for(sj = high_overlapping_coms.begin(); sj != high_overlapping_coms.end(); ++sj)
        {
            cout<<*sj<<"->";
            print_set(coms[*sj]);
        } */
        if(!high_overlapping_coms.empty())
            coms.erase(*si);
        labels.erase(*si);
    }while(!labels.empty());
}

float overlap(Graph& g, set<int>& C1, set<int>& C2)
{
    set<int>::iterator si, sj;
    float intercom_weight = 0, deg1 = 0, deg2 = 0;
    for(si = C1.begin(); si != C1.end(); ++si)
    {
        for(sj = C2.begin(); sj != C2.end(); ++sj)
            intercom_weight += g.get_weight(*si, *sj);
        deg1 += g.strength[*si];
    }
    for(sj = C2.begin(); sj != C2.end(); ++sj)
        deg2 += g.strength[*sj];
    return intercom_weight/min(deg1, deg2);
}




void read_communities(string& comfile, map<int, set<int> >& coms, map<int, set<int> >& memberships)
{
    std::ifstream fin(comfile);
    string line;
    if(!fin.is_open())
    {   cout<<"Community file could not be opened."<<endl;
        exit(1);
    }
    int i = 1;
    while ( std::getline(fin, line ) )
    {
        istringstream is( line );
        int v;
        while(is>>v)
        {
            coms[i].insert(v);
            memberships[v].insert(i);
        }
        i++;
   }
}

pair<int, int>LargestValue(map<int, float> sampleMap)
{
	pair<int, float> MaxValue = *sampleMap.begin();
	map<int, float>::iterator currentEntry;
	for (currentEntry = sampleMap.begin();currentEntry != sampleMap.end();++currentEntry) {
		if (currentEntry->second> MaxValue.second) {
             MaxValue= make_pair(currentEntry->first,currentEntry->second);
		   }
	    }


    return MaxValue;
}
void print_communities(map<int, set<int> >& C)
{
    cout<<"community"<<endl;
    map<int, set<int> >::iterator sitr;
    for(sitr = C.begin(); sitr != C.end(); ++sitr)
    {
        copy((*sitr).second.begin(), (*sitr).second.end(), ostream_iterator<int>(cout, " "));
        cout<<endl;
    }
}


float Graph::get_sn_nbrs_proximity(int u, map<int, float>& prox)
{
    map<int, float>::iterator mi, mj;
    float max_prox = 0;
    for(mi = weight[u].begin(); mi != weight[u].end(); ++mi)
    {
        prox[mi->first] = proximity(u, mi->first);
        if(prox[mi->first] > max_prox)
            max_prox = prox[mi->first];
        for(mj = weight[mi->first].begin(); mj != weight[mi->first].end(); ++mj)
            if(mj->first != u && prox.find(mj->first) == prox.end())
            {
                prox[mj->first] = proximity(u, mj->first);
                if(prox[mj->first] > max_prox)
                    max_prox = prox[mj->first];
            }
    }
    //cout<<"max prox = "<<max_prox<<endl;
    return max_prox;
}


int main(int argc, char* argv[])
{
    std::chrono::time_point<std::chrono::system_clock> start_time, end_time;
    start_time = std::chrono::system_clock::now();  //time starts
     float max_overlap = 0.4, rho_0 = 0.4;
    bool weighted = false;
    Graph g;
    map<int, set<int> > coms;
    string network_file ;
 
     if(argc < 2 || argc > 7)
     {
	cout<<"Please see README file."<<endl;
	exit(1);
     }
    if(argv[1][0] == '-')
    {
	cout<<"Please see README file."<<endl;
	exit(1);
     }
    else
    {
        network_file = string(argv[1]);
    }
    int i=2;
    while(i < argc)
    {
        string arg = string(argv[i]);
        if(arg == "-w")
        {
            weighted = true;
            i++;
        }
        else
            if(arg == "-ov")
            {
                istringstream is(argv[i+1]);
                is>>max_overlap;
                if( max_overlap < 0 || max_overlap > 0.5)
		{
		    cout<<"Please see README file."<<endl;
		    exit(1);
		}
                i += 2;
            }
            else
                if(arg == "-rh")
                {
                    istringstream is(argv[i+1]);
                    is>>rho_0;
                    if( rho_0 < 0 || rho_0 >= 1)
		    {    
			cout<<"Please see README file."<<endl;
			exit(1);
		    }
                    i += 2;
                }
                else
                {
                    cout<<"Please see README file."<<endl;
                    exit(1);
                  }
}
    g.read_edgelist(network_file ,weighted, false);
 
    cout<<"-------------------------------------------------------------------"<<endl;
    cout<<"Neighbourhood Proximity based Community Detection using Weighted Centrality (NPCD-WC) "<<endl;
    cout<<"Authors: Abhinav Kumar, Pawan Kumar and Ravins Dohare"<<endl;
    cout<<"Email: abhinavkumar080395@gmail.com, pkumariitd@gmail.com, ravinsdohare@gmail.com"<<endl;
    cout<<"-------------------------------------------------------------------"<<endl;
  
   
    map<int, set<int> >::iterator ci;
    set<int> uncovered = g.V;
    //cout<<"initializing communities..."<<endl;
    set<int> seeds;
    map<int, set<int> > membership;
    set<int>::iterator si, sj;
    //cout<<"expanding communities..."<<endl;
    int com_count = 1;
      while(!g.Centnode.empty())
    {
      //  cout<<"g.centnode size = "<<g.Centnode.size()<<endl;
        pair<int,int> yup=LargestValue(g.Centnode);
        int r=yup.first;
       // cout<<"curr uncovered = "<<r<<endl;
        map<int, float> prox;
        float max_prox;
        max_prox = g.get_sn_nbrs_proximity(r, prox);
        map<int, float>::iterator mi;
        //cout<<"nbrs proximities"<<endl;
       //for(mi = prox.begin(); mi != prox.end(); ++mi)
        //cout<<mi->first<<"--"<<mi->second<<endl;
        //cout<<"max prox = "<<max_prox<<endl;
        for(mi = prox.begin(); mi != prox.end(); ++mi)
        {
            if(mi->second > rho_0 || mi->second >= max_prox)
                membership[r].insert(membership[mi->first].begin(), membership[mi->first].end());
        }
        if(membership[r].empty())
        {
            membership[r].insert(com_count);
            com_count++;
        }
        for(mi = prox.begin(); mi != prox.end(); ++mi)
             if(mi->second > rho_0 && membership[mi->first].empty())
             {
                 membership[mi->first].insert(membership[r].begin(), membership[r].end());
                 g.Centnode.erase(mi->first);
               // cout<<"membership "<<mi->first<<" = ";
                //print_set(membership[mi->first]);
             }
          g.Centnode.erase(r);
     }
      for(si = g.V.begin(); si != g.V.end(); ++si)
    {
        for(sj = membership[*si].begin(); sj != membership[*si].end(); ++sj)
            coms[*sj].insert(*si);
    }
    //cout<<"communities before  merging = "<<coms.size()<<endl;
    int avg_csize = 0;
    for(ci = coms.begin(); ci != coms.end(); ++ci)
        {   //cout<<ci->second.size()<<", ";
            avg_csize += ci->second.size();
        }
    avg_csize /= coms.size();
    //cout<<"average size of initial communities = "<<avg_csize<<endl;
    //cout<<"merging communities..."<<endl;
    merge_communities(g, coms, membership, max_overlap, avg_csize);
    //cout<<"writing final communities..."<<endl;
    cout<<endl;
    ostringstream str_o;
    str_o.setf(ios::fixed, ios::floatfield);
    str_o.precision(2);
    str_o<<max_overlap;

    ostringstream comfile;
    comfile<<"./npcd_wc-coms.txt";
    ofstream fout(comfile.str());
    if(!fout.is_open())
    {
        cout<<"Destination file for communities could not be opened.";
        exit(1);
    }
    for(ci = coms.begin(); ci != coms.end(); ++ci)
    {
        copy(ci->second.begin(), ci->second.end(), ostream_iterator<int>(fout, " "));
        fout<<endl;
    }
    fout.close();
    
   end_time = std::chrono::system_clock::now();  //time ends
    std::chrono::duration<double> elapsed_seconds = end_time-start_time;
   
    cout.setf(ios::left, ios::adjustfield);  
    cout<<"-------------------------------------------------------------------"<<endl;
    cout<<setw(33)<<"Network file"<<"= "<<setw(20)<<network_file<<endl;
    cout<<setw(33)<<"Network order"<<"= "<<setw(20)<<g.order()<<endl;
    cout<<setw(33)<<"No. of edges"<<"= "<<setw(20)<<g.ecount()<<endl;
    cout<<setw(33)<<"Total non-singleton communities"<<"= "<<coms.size()<<endl;
    //cout<<setw(33)<<"Seed file"<<"= "<<setw(20)<<seedfile<<endl;
    cout<<setw(33)<<"Community file"<<"= "<<setw(20)<<comfile.str()<<endl;
    //cout<<setw(33)<<"No. of uncovered vertices"<<"= "<<setw(20)<<uncovered_vertices.size()<<endl;
  cout<<setw(33)<<"Time elapsed"<<"= "<<elapsed_seconds.count()<<"s\n";
    cout<<"-------------------------------------------------------------------"<<endl;
    return 0;
}
