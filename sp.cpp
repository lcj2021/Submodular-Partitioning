#include <ctime>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <unordered_map>
#include <set>
#include <cstring>
#include <string>
#include <algorithm>
using namespace std;
clock_t st, ed;
double endtime;

vector<string> 
split(string textline, string tag)
{
    // textline = textline.substr(0, textline.find_last_of(" "));
	vector<string> res;
	if (tag == "\t")
	{
		std::size_t pre_pos = 0;
		std::size_t pos = textline.find(tag);

		//依次处理line中的每个tag数据
		while (pos != std::string::npos)
		{
			string curStr = textline.substr(pre_pos, pos - pre_pos);
			curStr.erase(0, curStr.find_first_not_of("\r\t\n "));
			curStr.erase(curStr.find_last_not_of("\r\t\n ") + 1);

			if (strcmp(curStr.c_str(), "") != 0)
				res.push_back(curStr);
			pre_pos = pos + tag.size();
			pos = textline.find(tag, pre_pos);
		}

		string curStr = textline.substr(pre_pos, pos - pre_pos);
		curStr.erase(0, curStr.find_first_not_of("\r\t\n "));
		curStr.erase(curStr.find_last_not_of("\r\t\n ") + 1);
		if (strcmp(curStr.c_str(), "") != 0)
			res.push_back(curStr);
	}
	else
	{
		string curStr = textline;
		string subjectStr = curStr.substr(0, curStr.find(" "));
		curStr = curStr.substr(curStr.find(" ") + 1);
		string predicateStr = curStr.substr(0, curStr.find(" "));
		curStr = curStr.substr(curStr.find(" ") + 1);
		string objectStr = curStr;
		// cout << objectStr << 1 << endl;

		res.push_back(subjectStr);
		res.push_back(predicateStr);
		res.push_back(objectStr);
	}

	return res;
}

void 
trim(std::string &textline)
{
	if (textline.empty())	return ;
	textline.erase(0, textline.find_first_not_of(" "));
	textline.erase(textline.find_last_not_of(" ") + 1);
}

set<string>
s_inter(set<string> a, set<string> b)
{
	set<string> res;
	set_intersection(a.begin(), a.end(), b.begin(), b.end(), inserter(res, res.begin()));
	return res;
}

set<string>
s_union(vector<set<string> > _v)
{
	set<string> res = _v.front();
	for (auto it : _v)
	set_union(res.begin(), res.end(), it.begin(), it.end(), inserter(res, res.begin()));
	return res;
}

class graph
{
public:
	long long triples;

	// entityCnt : entity count
	long long entity_cnt;
	long long strEntity_cnt;

	// predicate count
	int pre_cnt;

    set<string> predicate;
	unordered_map<string, int> entity2ID;
	vector<string> ID2Entity;
	unordered_map<string, int> strEntity2ID;
	vector<string> ID2StrEntity;
	unordered_map<string, int> predicate2ID;
	vector<string> ID2Predicate;
	vector<vector<pair<int, int>>> edge;

    void init();
    void load(int part_cnt, string db_name, string tag, int db_id);
};

double 
f(int i, set<string> miu, graph *G)
{
    auto F = set<string>(G[i].ID2Entity.begin(), G[i].ID2Entity.end());
    // for (auto it : F)   cout << it << endl;
	set<string> intersect;
	intersect = s_inter(F, miu);	
    // for (auto it : intersect)   cout << it << endl;	puts("");
	return (double)intersect.size();
}

double 
F(int part_cnt, double c, vector<set<string> > V, vector<set<int> > pi, graph *G)
{
	double ret = 0;
	for (int i = 1; i <= part_cnt; ++ i)
	{
		set<string > miu_i;
		for (auto it : pi[i])	
			miu_i.insert(V[it].begin(), V[it].end());
		ret += min((double)f(i, miu_i, G), c);
	}
	ret /= part_cnt;
	return ret;
}

vector<set<int> >
greed_swp(int part_cnt, double c, vector<set<string> > V, graph * G)
{
	vector<set<int> > ret(part_cnt + 1);		// return val;
	vector<set<string> > A(part_cnt + 1);
	vector<set<string> > R = V;
	set<int> R_map;
	for (int i = 1; i <= V.size() - 1; ++ i)	R_map.insert(i);

	while (R_map.size())
	{
		int j_V = *R_map.begin(), j_A = 1;
		double delta = -1;

		cout << "round " << V.size() - R_map.size() << endl;

		for (int i = 1; i <= part_cnt; ++ i)
		{
			double item2 = min(f(i, A[i], G), c), item1 = -1;
			for (auto r : R_map)
			{
				cout << "i=" << i << " r=" << r << "\t" << "max(min{fi(Ai U r), c})=" << min(f(i, s_union({A[i], V[r]}), G), c) << "\t" << "min{fi(Ai), c}=" << item2; 
				item1 = max(min(f(i, s_union({A[i], V[r]}), G), c), item1);
				double val = item1 - item2;
				if (val > delta)
					delta = val, j_V = r, j_A = i;
				cout << "\t" << "delta=" << val << endl;
			}
		}
		ret[j_A].insert(j_V);
		A[j_A].insert(V[j_V].begin(), V[j_V].end());
		R_map.erase(j_V);

		
		cout << "=============== J_A:J_V " << j_A << ":" << j_V << endl;
	}
	return ret;
}

vector<set<int> >
greed_sat(int part_cnt, double epsilon, double alpha, vector<set<string> > V, graph * G)
{
	vector<set<int> > ret(part_cnt + 1);		// return val;
	set<string> V_set = s_union(V);
	
	double c_min = 0, c_max = 1e10;
	for (int i = 1; i <= part_cnt; ++ i)
		c_max = min(c_max, f(i, V_set, G));
	// cout << c_min << " " << c_max << endl;

	while (c_max - c_min >= epsilon)
	{
		double c = (c_max + c_min) / 2;
		vector<set<int> > pi_c = greed_swp(part_cnt, c, V, G);
		double F_c = F(part_cnt, c, V, pi_c, G);

		if (F_c < alpha * c)	c_max = c, puts("left");
		else					c_min = c, ret = pi_c, puts("right");
		
		for (int i = 1; i <= pi_c.size() - 1; ++ i)
		{
			for (auto it : pi_c[i])
				cout << it << "\t";	
			puts("");
		}
		cout << "c=" << c << "\t" << "c_max - c_min=" << c_max - c_min << endl;
	}
	return ret;
}

//[0]./bin/sp   [1]db_name   [2]tag    [3]k   [4]pattern_template  [5]pattern_result
int main(int argc, char *argv[])
{
	string db_name = argv[1];
	string tag = (string(argv[2]) == "1") ? " " : "\t";
	int part_cnt = atoi(argv[3]);
	// string pattern_path = argv[4];

	vector<set<string> > id2miu;
	// unordered_map<set<string>, int> miu2id;

	st = clock();
		graph *test = new graph[part_cnt + 1];
        for (int i = 1; i <= part_cnt; ++ i)
        {
            test[i].init();
            test[i].load(part_cnt, db_name, tag, i);
        }
		// test->part_cnt = part_cnt;


	ed = clock();
	endtime = (double)(ed - st) / CLOCKS_PER_SEC;
	cout << "partition : " << endtime << " s" << endl;
    
    // cout << test[3].triples << endl;
    // for (auto it : test[1].predicate)           cout << it << endl;
    // for (auto it : test->edge[test->predicate2ID["<b>"]])     cout << it.first << "\t" << it.second << endl;
    // for (auto it : test->ID2StrEntity)          cout << it << endl;
    // set<string> s;
    // set_intersection(test[1].predicate.begin(), test[1].predicate.end(), test[2].predicate.begin(), test[2].predicate.end(), inserter(s, s.begin()));   puts("");
    // for (auto it : s)   cout << it << endl;

    set<string> miu1{"<C>", "<A>", "<B>", "<E>"};
    set<string> miu2{"<J>", "<C>", "<D>", "<I>"};
    set<string> miu3{"<H>", "<G>", "<D>", "<I>"};
    set<string> miu4{"<H>", "<G>", "<F>", "<E>"};
    set<string> miu5{"<H>", "<A>", "<C>", "<E>"};
    set<string> miu6{"<G>", "<A>", "<D>", "<B>"};
    set<string> miu7{"<H>", "<G>", "<D>", "<G>"};
    set<string> miu8{"<A>", "<B>", "<E>", "<F>"};
    // cout << f(1, miu3, test) << endl;
	// cout << f(1, s_union({miu1, miu3}), test) << endl;
	// cout << f(1, s_union({miu2, miu3}), test) << endl;
	// cout << f(1, s_union({miu2, miu3, miu1}), test) << endl;

	// vector<set<string> > miu = {set<string>(), miu1};
	// vector<set<string> > miu = {set<string>(), miu1, miu2};
	// vector<set<string> > miu = {set<string>(), miu1, miu2, miu3};
	// vector<set<string> > miu = {set<string>(), miu1, miu2, miu3, miu4};
	// vector<set<string> > miu = {set<string>(), miu1, miu2, miu3, miu4, miu5};
	// vector<set<string> > miu = {set<string>(), miu1, miu2, miu3, miu4, miu5, miu6};
	vector<set<string> > miu = {set<string>(), miu1, miu2, miu3, miu4, miu5, miu6, miu7};
	// vector<set<string> > miu = {set<string>(), miu1, miu2, miu3, miu4, miu5, miu6, miu7, miu8};

	// vector<set<int> > swp1 = greed_swp(4, 2, miu, test);
	// for (int i = 1; i <= part_cnt; ++ i)
	// {
	// 	for (auto it : swp1[i])
	// 		cout << it << "\t";	
	// 	puts("");
	// }

	// cout << F(4, 2, miu, swp1, test) << endl;

	vector<set<int> > pi = greed_sat(part_cnt, 1, 0.6321, miu, test);
	for (int i = 1; i <= part_cnt; ++ i)
	{
		for (auto it : pi[i])
			cout << it << "\t";	
		puts("");
	}


	// delete test;
	return 0;

}

void
graph::load(int part_cnt, string db_name, string tag, int db_id)
{
    string line, data_name = db_name + to_string(db_id);
    ifstream in(data_name);
    vector<pair<int, int>> tmp;
    cout << data_name << "========" << endl;
    while (getline(in, line))
    {
        if (triples % 10000 == 0)
            cout << "loading triples : " << triples << endl;
        triples ++;
        // line.resize(line.length() - 2);

        vector<string> s;
        if (line.find_last_of(".") != string::npos)
            line = line.substr(0, line.find_last_of("."));
        trim(line);
        s = split(line, tag);

        predicate.insert(s[1]);

        for (int i = 0; i < 3; i += 2)
        {
            if ((s[i][0] == '<' || s[i][0] == '_') && entity2ID.count(s[i]) == 0)
            {
                entity2ID[s[i]] = ++entity_cnt;
                ID2Entity.push_back(s[i]);
                // entityTriples.push_back(0);
            }
            else if ((s[i][0] == '"') && strEntity2ID.count(s[i]) == 0)
            {
                strEntity2ID[s[i]] = ++strEntity_cnt;
                ID2StrEntity.push_back(s[i]);
                // entityTriples.push_back(0);
            }
        }

        // edge_cnt[s[1]]++;
        int a = entity2ID[s[0]];

        // entityTriples[a]++;

        if ((s[0][0] == '<' || s[0][0] == '_') && (s[2][0] == '<' || s[2][0] == '_'))
        {
            if (predicate2ID.count(s[1]) == 0)
            {
                predicate2ID[s[1]] = ++pre_cnt;
                ID2Predicate.push_back(s[1]);
                edge.push_back(tmp);
            }
            int b = entity2ID[s[2]];
            edge[predicate2ID[s[1]]].push_back({a, b});

            // entityTriples[b]++;
        }
    }
    in.close();
}

void 
graph::init()
{
	vector<pair<int, int>> tmp;
	edge.push_back(tmp);

	// entityTriples.push_back(0);

	ID2Entity.push_back("");
	ID2Predicate.push_back("");
	pre_cnt = entity_cnt = triples = 0;
}

// g++ sp.cpp -o sp -std=c++11
// ./sp test_data 2 4 