// uagd_nsga.cpp
// Urban Air–Ground Cooperative Delivery — NSGA-II (single-seed C++)
// Primary objective:  minimize makespan (longest per-truck completion time)
// Secondary objective: minimize total movement time (truck + drone flight)
//
// 与原 Python 单种子版本等价：
// - 支持“同一车边多无人机并发”的解码（partition -> chains）
// - mask 仅作起点提示（不强制整段）
// - 2-opt 平滑；增益驱动的链接受；轻量强化 intensify_once
// - 约束：仅在顾客点汇合；顾客各一次；无人机容量/航程；每条车边并发≤车载无人机数；车 L1/机 L2 距离；
//        车一次出发终回仓；等待计入完工时间
//
// Usage:
//   uagd_nsga.exe --input 6.txt --pop 80 --gens 240 --seed 42 --out Instance6.txt
//
// 依赖：C++17，无第三方库。建议以 -O3 -march=native 编译。

#include <bits/stdc++.h>
using namespace std;

struct Instance {
    string file_path, name;
    int drones_per_truck = 0;
    int trucks = 0;
    double truck_capacity = 0.0;
    double drone_capacity = 0.0;
    double truck_max_dist = 0.0;
    double drone_max_dist = 0.0;
    double truck_speed = 1.0;
    double drone_speed = 1.0;
    int num_nodes_total = 0; // depot + customers
    int num_customers = 0;

    vector<double> xs; // [0..n] 0==depot
    vector<double> ys;
    vector<double> demands; // [0..n]
    // 预计算
    vector<double> D1; // L1 距离矩阵 (N*N)
    vector<double> DE; // 欧氏距离矩阵 (N*N)
    inline int N() const { return num_customers + 1; }
};

// --------- 工具 ----------
static inline void logline(const string &s){ cout<<s<<endl; }
static inline double l1(double xa,double ya,double xb,double yb){
    return fabs(xa-xb)+fabs(ya-yb);
}
static inline double euc(double xa,double ya,double xb,double yb){
    double dx=xa-xb, dy=ya-yb; return sqrt(dx*dx+dy*dy);
}
static inline double rnd01(std::mt19937_64 &rng){
    static uniform_real_distribution<double> U(0.0,1.0); return U(rng);
}
template<class T> static inline T rnd_int(std::mt19937_64 &rng, T lo, T hi){
    uniform_int_distribution<long long> U((long long)lo,(long long)hi);
    return (T)U(rng);
}

// --------- 解析输入 ----------
static int first_number(const string &s, double &out, bool &is_int){
    // 找到第一个数字（int/float），返回找到位置（>=0），并设置 out/is_int
    regex re("[-+]?\\d*\\.?\\d+(?:[eE][-+]?\\d+)?");
    smatch m;
    if(regex_search(s, m, re)){
        string v = m.str();
        out = stod(v);
        is_int = (v.find('.')==string::npos && v.find('e')==string::npos && v.find('E')==string::npos);
        return (int)m.position();
    }
    return -1;
}
static Instance parse_instance(const string &path){
    ifstream fin(path);
    if(!fin) throw runtime_error("Cannot open file: "+path);
    string raw((istreambuf_iterator<char>(fin)), istreambuf_iterator<char>());
    // 去掉 C 风格块注释
    {
        string out; out.reserve(raw.size());
        for(size_t i=0;i<raw.size();){
            if(i+1<raw.size() && raw[i]=='/' && raw[i+1]=='*'){
                i+=2;
                while(i+1<raw.size() && !(raw[i]=='*' && raw[i+1]=='/')) ++i;
                i+= (i+1<raw.size()?2:1);
            }else out.push_back(raw[i++]);
        }
        raw.swap(out);
    }
    // 分行+trim+过滤空行
    vector<string> lines; lines.reserve(4096);
    {
        string cur;
        for(char c: raw){
            if(c=='\r') continue;
            if(c=='\n'){ 
                auto s=cur; 
                auto l=s.find_first_not_of(" \t"); if(l==string::npos) s.clear();
                else { auto r=s.find_last_not_of(" \t"); s=s.substr(l,r-l+1); }
                if(!s.empty()) lines.push_back(s);
                cur.clear();
            }else cur.push_back(c);
        }
        if(!cur.empty()){
            auto s=cur; auto l=s.find_first_not_of(" \t");
            if(l!=string::npos){ auto r=s.find_last_not_of(" \t"); s=s.substr(l,r-l+1); }
            if(!s.empty()) lines.push_back(s);
        }
    }
    auto atnum=[&](int idx)->double{
        double v=0; bool is_int=false; int pos=first_number(lines[idx], v, is_int);
        if(pos<0) throw runtime_error("parse number failed at line "+to_string(idx));
        return v;
    };

    int idx=0;
    Instance ins;
    ins.file_path=path;
    ins.drones_per_truck=(int)atnum(idx++); 
    ins.trucks=(int)atnum(idx++);
    ins.truck_capacity=atnum(idx++);
    ins.drone_capacity=atnum(idx++);
    ins.truck_max_dist=atnum(idx++);
    ins.drone_max_dist=atnum(idx++);
    ins.truck_speed=atnum(idx++);
    ins.drone_speed=atnum(idx++);
    ins.num_nodes_total=(int)atnum(idx++);
    ins.num_customers=ins.num_nodes_total-1;

    // depot
    {
        regex re("[-+]?\\d*\\.?\\d+(?:[eE][-+]?\\d+)?");
        smatch m; string s=lines[idx++];
        vector<double> vals;
        auto it=s.cbegin();
        while(regex_search(it, s.cend(), m, re)){ vals.push_back(stod(m.str())); it=m.suffix().first; }
        if(vals.size()<2) throw runtime_error("depot line need xy");
        double depot_x=vals[0], depot_y=vals[1];
        ins.xs.resize(ins.num_nodes_total);
        ins.ys.resize(ins.num_nodes_total);
        ins.demands.resize(ins.num_nodes_total);
        ins.xs[0]=depot_x; ins.ys[0]=depot_y; ins.demands[0]=0.0;
    }
    for(int k=0;k<ins.num_customers;k++){
        string s=lines[idx+k];
        regex re("[-+]?\\d*\\.?\\d+(?:[eE][-+]?\\d+)?|loc\\d+|depot|^0$");
        smatch m; vector<string> toks;
        auto it=s.cbegin();
        while(regex_search(it, s.cend(), m, re)){ toks.push_back(m.str()); it=m.suffix().first; }
        if(toks.size()<3) throw runtime_error("customer line need x y dem");
        double x=stod(toks[0]), y=stod(toks[1]); double dem=stod(toks.back());
        ins.xs[1+k]=x; ins.ys[1+k]=y; ins.demands[1+k]=dem;
    }

    ins.name = std::filesystem::path(path).filename().string();

    // 预计算距离矩阵
    int N=ins.N();
    ins.D1.assign((size_t)N*(size_t)N, 0.0);
    ins.DE.assign((size_t)N*(size_t)N, 0.0);
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            ins.D1[(size_t)i*N+j] = l1(ins.xs[i],ins.ys[i],ins.xs[j],ins.ys[j]);
            ins.DE[(size_t)i*N+j] = euc(ins.xs[i],ins.ys[i],ins.xs[j],ins.ys[j]);
        }
    }
    return ins;
}

// --------- 几何/时间 ----------
static inline double truck_leg_time(const Instance& ins, int a, int b){
    int N=ins.N(); return ins.D1[(size_t)a*N+b] / ins.truck_speed;
}
static double truck_path_time_segment(const Instance& ins, const vector<int>& nodes){
    if(nodes.size()<2) return 0.0;
    double sum=0.0; int N=ins.N();
    for(size_t i=0;i+1<nodes.size();++i){
        int a=nodes[i], b=nodes[i+1];
        sum += ins.D1[(size_t)a*N+b] / ins.truck_speed;
    }
    return sum;
}
static pair<double,double> chain_time_and_length(const Instance& ins, int start_node, const vector<int>& chain, int end_node){
    int N=ins.N();
    double e=0.0;
    if(chain.empty()){
        e = ins.DE[(size_t)start_node*N+end_node];
    }else{
        e += ins.DE[(size_t)start_node*N+chain[0]];
        for(size_t i=0;i+1<chain.size();++i) e += ins.DE[(size_t)chain[i]*N+chain[i+1]];
        e += ins.DE[(size_t)chain.back()*N+end_node];
    }
    return {e / ins.drone_speed, e};
}

static double simulate_truck_time_with_chains(const Instance& ins, const vector<int>& route,
                                              const map<pair<int,int>, vector<vector<int>>>& chains){
    double tot=0.0;
    for(size_t k=0;k+1<route.size();++k){
        int a=route[k], b=route[k+1];
        double t_leg = truck_leg_time(ins,a,b);
        double max_drone = 0.0;
        auto it=chains.find({a,b});
        if(it!=chains.end()){
            for(const auto& ch: it->second){
                auto dt = chain_time_and_length(ins,a,ch,b).first;
                if(dt>max_drone) max_drone=dt;
            }
        }
        double wait = max(0.0, max_drone - t_leg);
        tot += t_leg + wait;
    }
    return tot;
}
static double total_movement_time(const Instance& ins, const vector<int>& route,
                                  const map<pair<int,int>, vector<vector<int>>>& chains){
    double tot=0.0;
    for(size_t k=0;k+1<route.size();++k){
        int a=route[k], b=route[k+1];
        tot += truck_leg_time(ins,a,b);
    }
    for(const auto& kv: chains){
        auto a=kv.first.first, b=kv.first.second;
        for(const auto& ch: kv.second){
            tot += chain_time_and_length(ins,a,ch,b).first;
        }
    }
    return tot;
}

// 最近邻 + 2-opt（L1）
static vector<int> nearest_neighbor_tsp_l1(const Instance& ins, const vector<int>& nodes){
    if(nodes.empty()) return {0,0};
    auto l1d=[&](int A,int B){ return l1(ins.xs[A],ins.ys[A],ins.xs[B],ins.ys[B]); };
    unordered_set<int> unvisited(nodes.begin(),nodes.end());
    int cur = *min_element(unvisited.begin(),unvisited.end(),
                           [&](int j1,int j2){return l1d(j1,0)+0 < l1d(j2,0)+0;});
    vector<int> route={0,cur}; unvisited.erase(cur);
    while(!unvisited.empty()){
        int best=-1; double bestd=1e100;
        for(int j: unvisited){
            double d=l1d(j,cur);
            if(d<bestd){bestd=d;best=j;}
        }
        route.push_back(best); unvisited.erase(best); cur=best;
    }
    route.push_back(0);
    return route;
}
static vector<int> two_opt_improve_l1(const Instance& ins, vector<int> route, int max_iter=60){
    auto l1d=[&](int A,int B){ return l1(ins.xs[A],ins.ys[A],ins.xs[B],ins.ys[B]); };
    bool improved=true; int it=0;
    while(improved && it<max_iter){
        improved=false; ++it;
        for(size_t i=1;i+2<route.size();++i){
            for(size_t k=i+1;k+1<route.size();++k){
                int a=route[i-1], b=route[i], c=route[k], d=route[k+1];
                double curd=l1d(a,b)+l1d(c,d);
                double newd=l1d(a,c)+l1d(b,d);
                if(curd>newd+1e-9){
                    reverse(route.begin()+i, route.begin()+k+1);
                    improved=true;
                }
            }
        }
    }
    return route;
}

// ---------- KMeans for seeding ----------
static vector<int> kmeans_xy(const Instance& ins, int K, int max_iter, std::mt19937_64 &rng){
    int n=ins.num_customers;
    vector<array<double,2>> coords(n);
    for(int i=0;i<n;i++) coords[i]={ins.xs[i+1],ins.ys[i+1]};
    vector<array<double,2>> centers;
    centers.reserve(K);
    if(n==0){ return vector<int>(0); }
    // choose K unique random indexes
    vector<int> idx(n); iota(idx.begin(),idx.end(),0);
    shuffle(idx.begin(),idx.end(),rng);
    for(int i=0;i<min(K,n);++i) centers.push_back(coords[idx[i]]);
    while((int)centers.size()<K) centers.push_back(centers[0]);

    vector<int> lab(n,0);
    for(int it=0; it<max_iter; ++it){
        // assign
        for(int i=0;i<n;i++){
            double best=1e100; int bestj=0;
            for(int j=0;j<K;j++){
                double dx=coords[i][0]-centers[j][0];
                double dy=coords[i][1]-centers[j][1];
                double d=dx*dx+dy*dy;
                if(d<best){best=d;bestj=j;}
            }
            lab[i]=bestj;
        }
        // update
        vector<double> sx(K,0), sy(K,0); vector<int> cnt(K,0);
        for(int i=0;i<n;i++){ sx[lab[i]]+=coords[i][0]; sy[lab[i]]+=coords[i][1]; cnt[lab[i]]++; }
        for(int j=0;j<K;j++){
            if(cnt[j]>0){ centers[j][0]=sx[j]/cnt[j]; centers[j][1]=sy[j]/cnt[j]; }
        }
    }
    return lab;
}

static vector<char> heuristic_mask_guess(const Instance& ins, const vector<int>& assign){
    int nCust=ins.num_customers, m=ins.trucks;
    vector<vector<int>> per_truck(m);
    for(int cid=1; cid<=nCust; ++cid){
        int t=assign[cid-1];
        per_truck[t].push_back(cid);
    }
    vector<char> mask(nCust, 0);
    for(int t=0;t<m;++t){
        const auto &base = per_truck[t];
        if(base.empty()) continue;
        vector<int> route; route.reserve(base.size()+2);
        route.push_back(0);
        route.insert(route.end(), base.begin(), base.end());
        route.push_back(0);
        for(size_t k=1;k+1<route.size();++k){
            int a=route[k-1], b=route[k], c=route[k+1];
            double with_i = truck_path_time_segment(ins, {a,b,c});
            double without_i = truck_path_time_segment(ins, {a,c});
            if(with_i - without_i > 1e-6) mask[b-1]=1;
        }
    }
    return mask;
}

static void init_population_smart(const Instance& ins, int pop_size, int seed,
                                  vector<vector<int>>& assign,
                                  vector<vector<int>>& perms,
                                  vector<vector<char>>& masks){
    std::mt19937_64 rng(seed);
    int nCust=ins.num_customers, m=ins.trucks;
    auto labels = kmeans_xy(ins, m, 25, rng);
    int half = pop_size/2;

    assign.assign(pop_size, vector<int>(nCust));
    for(int i=0;i<half;i++){
        auto a = labels;
        for(int j=0;j<nCust;j++){
            if(rnd01(rng)<0.10) a[j]=rnd_int<int>(rng,0,m-1);
        }
        assign[i]=move(a);
    }
    for(int i=half;i<pop_size;i++){
        for(int j=0;j<nCust;j++) assign[i][j]=rnd_int<int>(rng,0,m-1);
    }

    perms.assign(pop_size, vector<int>(nCust));
    vector<int> base(nCust); iota(base.begin(),base.end(),1);
    for(int i=0;i<pop_size;i++){
        perms[i]=base;
        shuffle(perms[i].begin(), perms[i].end(), rng);
    }

    masks.assign(pop_size, vector<char>(nCust,0));
    for(int i=0;i<half;i++){
        masks[i]=heuristic_mask_guess(ins, assign[i]);
        for(int j=0;j<nCust;j++){
            if(rnd01(rng)<0.10) masks[i][j] = !masks[i][j];
        }
    }
    for(int i=half;i<pop_size;i++){
        for(int j=0;j<nCust;j++){
            masks[i][j] = (rnd01(rng)<0.30)?1:0;
        }
    }
}

// --------- 多机链分割 ----------
static bool partition_segment_into_chains(const Instance& ins, int start,
                                          const vector<int>& segment, int end,
                                          int max_chains, vector<vector<int>>& out_chains){
    out_chains.clear();
    const double dcap = ins.drone_capacity;
    const double dmax = ins.drone_max_dist;
    double cur_dem=0.0;
    vector<int> cur;

    auto feasible = [&](const vector<int>& nodes)->bool{
        auto p = chain_time_and_length(ins,start,nodes,end);
        return p.second <= dmax + 1e-9;
    };

    for(int cid: segment){
        bool need_new=false;
        if(!cur.empty() && (cur_dem + ins.demands[cid] > dcap + 1e-9)){
            need_new=true;
        }else{
            if(ins.demands[cid] > dcap + 1e-9) return false;
            vector<int> test=cur; test.push_back(cid);
            if(!feasible(test)) need_new=true;
        }
        if(need_new){
            if(cur.empty()) return false; // 单点不可行
            out_chains.push_back(cur);
            if((int)out_chains.size()>=max_chains) return false;
            cur.clear();
            cur.push_back(cid);
            cur_dem = ins.demands[cid];
            if(!feasible(cur)) return false;
        }else{
            cur.push_back(cid);
            cur_dem += ins.demands[cid];
        }
    }
    if(!cur.empty()){
        out_chains.push_back(cur);
        if((int)out_chains.size()>max_chains) return false;
    }
    return true;
}

// --------- 解码（支持同一车边多无人机并发） ----------
static void decode_individual(const Instance& ins,
                              const vector<int>& assign,
                              const vector<int>& perm,
                              const vector<char>& mask,
                              vector<vector<int>>& truck_routes,
                              vector<map<pair<int,int>, vector<vector<int>>>>& chains_all){
    int nCust=ins.num_customers, m=ins.trucks;

    vector<vector<int>> per_truck_customers(m);
    for(int cid: perm){
        int t=assign[cid-1];
        per_truck_customers[t].push_back(cid);
    }
    truck_routes.clear(); truck_routes.resize(m);
    for(int t=0;t<m;t++){
        const auto &base = per_truck_customers[t];
        if(base.empty()){ truck_routes[t]={0,0}; continue; }
        vector<int> route; route.reserve(base.size()+2);
        route.push_back(0);
        route.insert(route.end(), base.begin(), base.end());
        route.push_back(0);
        route = two_opt_improve_l1(ins, route, 40);
        truck_routes[t]=move(route);
    }

    chains_all.clear(); chains_all.resize(m);
    const int MAX_WINDOW=20;

    for(int t=0;t<m;t++){
        auto &route = truck_routes[t];
        auto &chains_per_leg = chains_all[t];
        if(route.size()<=2) continue;

        size_t i=1;
        while(i+0<route.size()-1){
            int start = route[i-1];
            int cur = route[i];
            if(start==0 || !mask[cur-1]){ i++; continue; }

            double best_gain_ms=-1e100, best_tie_mov=-1e100;
            int best_r=-1; vector<vector<int>> best_chains;

            int max_r = (int)min(route.size()-2, (size_t)i + MAX_WINDOW - 1);
            for(int r=max_r; r>= (int)i; --r){
                int end = route[r+1];
                if(end==0) continue;
                vector<int> segment(route.begin()+i, route.begin()+r+1);
                vector<vector<int>> chains;
                if(!partition_segment_into_chains(ins, start, segment, end, ins.drones_per_truck, chains)) continue;

                double T_seg = truck_path_time_segment(ins, [&]{
                    vector<int> tmp; tmp.reserve(segment.size()+2);
                    tmp.push_back(start); tmp.insert(tmp.end(), segment.begin(), segment.end()); tmp.push_back(end); return tmp;}());
                double T_dir = truck_leg_time(ins,start,end);
                vector<double> drone_times; drone_times.reserve(chains.size());
                for(auto &ch: chains) drone_times.push_back(chain_time_and_length(ins,start,ch,end).first);
                double wait = max(0.0, *max_element(drone_times.begin(),drone_times.end()) - T_dir);
                double gain_ms = T_seg - (T_dir + wait);
                double tie_gain_mov = T_seg - (T_dir + accumulate(drone_times.begin(),drone_times.end(),0.0));

                if( (gain_ms > 1e-12) ||
                    (fabs(gain_ms)<=1e-12 && tie_gain_mov > 1e-12) ){
                    if( (best_r==-1) ||
                        (gain_ms > best_gain_ms + 1e-12) ||
                        (fabs(gain_ms-best_gain_ms)<=1e-12 && tie_gain_mov > best_tie_mov + 1e-12) ){
                        best_r=r; best_gain_ms=gain_ms; best_tie_mov=tie_gain_mov; best_chains=move(chains);
                    }
                }
            }
            if(best_r!=-1){
                int end=route[best_r+1];
                for(auto &ch: best_chains) chains_per_leg[{start,end}].push_back(ch);
                // remove i..best_r, connect start->end
                vector<int> newroute; newroute.reserve(route.size() - (best_r - (int)i + 1));
                newroute.insert(newroute.end(), route.begin(), route.begin()+i);
                newroute.push_back(end);
                newroute.insert(newroute.end(), route.begin()+best_r+2, route.end());
                route.swap(newroute);
                // continue after end
                auto it=find(route.begin(),route.end(), end);
                i = (size_t)( (it==route.end()) ? route.size()-1 : (it - route.begin()) ) + 1;
            }else{
                i++;
            }
        }
    }
}

// --------- 目标 ----------
static tuple<double,double, vector<vector<int>>, vector<map<pair<int,int>, vector<vector<int>>>>>
evaluate_objectives(const Instance& ins,
                    const vector<int>& assign,
                    const vector<int>& perm,
                    const vector<char>& mask){
    vector<vector<int>> truck_routes;
    vector<map<pair<int,int>, vector<vector<int>>>> chains_all;
    decode_individual(ins, assign, perm, mask, truck_routes, chains_all);
    vector<double> per_truck(ins.trucks, 0.0);
    for(int t=0;t<ins.trucks;t++){
        per_truck[t] = simulate_truck_time_with_chains(ins, truck_routes[t], chains_all[t]);
    }
    double f1 = per_truck.empty()?0.0:*max_element(per_truck.begin(),per_truck.end());
    double f2 = 0.0;
    for(int t=0;t<ins.trucks;t++) f2 += total_movement_time(ins, truck_routes[t], chains_all[t]);
    return {f1, f2, truck_routes, chains_all};
}

// --------- NSGA-II 组件 ----------
static vector<int> tournament_select(const vector<int>& ranks, const vector<double>& crowd,
                                     int size, std::mt19937_64 &rng){
    int N=(int)ranks.size();
    vector<int> idx; idx.reserve(size);
    for(int s=0;s<size;s++){
        int a=rnd_int<int>(rng,0,N-1), b=rnd_int<int>(rng,0,N-1);
        if( (ranks[a]<ranks[b]) || (ranks[a]==ranks[b] && crowd[a]>crowd[b]) ) idx.push_back(a);
        else idx.push_back(b);
    }
    return idx;
}

static vector<int> ox_crossover(const vector<int>& p1, const vector<int>& p2, std::mt19937_64 &rng){
    int n=(int)p1.size();
    if(n<=1) return p1;
    int a=rnd_int<int>(rng,0,n-1), b=rnd_int<int>(rng,0,n-1);
    if(a>b) swap(a,b);
    vector<int> child(n,-1);
    vector<char> used(n+1,0);
    for(int i=a;i<=b;i++){ child[i]=p1[i]; used[p1[i]]=1; }
    int j=0;
    for(int i=0;i<n;i++){
        if(child[i]==-1){
            while(j<n && used[p2[j]]) j++;
            child[i]=p2[j++];
        }
    }
    return child;
}
static void swap_mutation(vector<int>& arr, double p, std::mt19937_64 &rng){
    int n=(int)arr.size();
    if(n>=2 && rnd01(rng)<p){
        int i=rnd_int<int>(rng,0,n-1), j=rnd_int<int>(rng,0,n-1);
        if(i!=j) swap(arr[i],arr[j]);
    }
}
static pair<vector<int>, vector<int>> uniform_crossover_assign(const vector<int>& a, const vector<int>& b, std::mt19937_64 &rng){
    int n=(int)a.size(); vector<int> c1(n), c2(n);
    for(int i=0;i<n;i++){
        if(rnd01(rng)<0.5){ c1[i]=a[i]; c2[i]=b[i]; } else { c1[i]=b[i]; c2[i]=a[i]; }
    }
    return {c1,c2};
}
static pair<vector<char>, vector<char>> uniform_crossover_mask(const vector<char>& a, const vector<char>& b, std::mt19937_64 &rng){
    int n=(int)a.size(); vector<char> c1(n), c2(n);
    for(int i=0;i<n;i++){
        if(rnd01(rng)<0.5){ c1[i]=a[i]; c2[i]=b[i]; } else { c1[i]=b[i]; c2[i]=a[i]; }
    }
    return {c1,c2};
}
static void bitflip(vector<char>& x, double p, std::mt19937_64 &rng){
    for(auto &v: x) if(rnd01(rng)<p) v = !v;
}
static void mutate_assign(vector<int>& a, int m, double p, std::mt19937_64 &rng){
    for(auto &v: a) if(rnd01(rng)<p) v = rnd_int<int>(rng,0,m-1);
}

static pair<vector<int>, vector<vector<int>>> fast_non_dominated_sort(const vector<double>& F1, const vector<double>& F2){
    int N=(int)F1.size();
    vector<vector<int>> S(N);
    vector<int> n(N,0), rank(N,0);
    vector<vector<int>> fronts(1);
    for(int p=0;p<N;p++){
        vector<int> Sp;
        int np_count=0;
        for(int q=0;q<N;q++){
            if(p==q) continue;
            bool pdomq = (F1[p]<=F1[q] && F2[p]<=F2[q]) && (F1[p]<F1[q] || F2[p]<F2[q]);
            bool qdomp = (F1[q]<=F1[p] && F2[q]<=F2[p]) && (F1[q]<F1[p] || F2[q]<F2[p]);
            if(pdomq) Sp.push_back(q);
            else if(qdomp) np_count++;
        }
        S[p].swap(Sp); n[p]=np_count;
        if(n[p]==0){ rank[p]=0; fronts[0].push_back(p); }
    }
    int i=0;
    while(!fronts[i].empty()){
        vector<int> Q;
        for(int p: fronts[i]){
            for(int q: S[p]){
                n[q]--;
                if(n[q]==0){ rank[q]=i+1; Q.push_back(q); }
            }
        }
        i++; fronts.push_back(Q);
    }
    fronts.pop_back();
    return {rank, fronts};
}
static vector<double> crowding_distance(const vector<double>& F1, const vector<double>& F2, const vector<int>& front){
    if(front.empty()) return {};
    int m=2;
    vector<double> dist(front.size(), 0.0);
    vector<int> order = front;

    // for F1
    sort(order.begin(), order.end(), [&](int a,int b){return F1[a]<F1[b];});
    dist[ distance(front.begin(), find(front.begin(), front.end(), order.front())) ] = numeric_limits<double>::infinity();
    dist[ distance(front.begin(), find(front.begin(), front.end(), order.back())) ] = numeric_limits<double>::infinity();
    double fmin=F1[order.front()], fmax=F1[order.back()];
    if(fmax - fmin > 1e-12){
        for(size_t j=1;j+1<order.size();++j){
            int prev=order[j-1], next=order[j+1], cur=order[j];
            dist[ distance(front.begin(), find(front.begin(), front.end(), cur)) ] += (F1[next]-F1[prev])/(fmax-fmin);
        }
    }
    // for F2
    order = front;
    sort(order.begin(), order.end(), [&](int a,int b){return F2[a]<F2[b];});
    dist[ distance(front.begin(), find(front.begin(), front.end(), order.front())) ] = numeric_limits<double>::infinity();
    dist[ distance(front.begin(), find(front.begin(), front.end(), order.back())) ] = numeric_limits<double>::infinity();
    fmin=F2[order.front()], fmax=F2[order.back()];
    if(fmax - fmin > 1e-12){
        for(size_t j=1;j+1<order.size();++j){
            int prev=order[j-1], next=order[j+1], cur=order[j];
            dist[ distance(front.begin(), find(front.begin(), front.end(), cur)) ] += (F2[next]-F2[prev])/(fmax-fmin);
        }
    }
    return dist;
}

// 轻量强化
static void intensify_once(const Instance& ins, vector<int>& assign, const vector<int>& perm, const vector<char>& mask, std::mt19937_64 &rng){
    int nCust=ins.num_customers, m=ins.trucks;
    if(m<=1) return;
    auto [f1_cur, f2_cur, tr_cur, ch_cur] = evaluate_objectives(ins, assign, perm, mask);

    int longest_t=-1; double longest_val=-1.0;
    for(int t=0;t<m;t++){
        double val = simulate_truck_time_with_chains(ins, tr_cur[t], ch_cur[t]);
        if(val>longest_val){longest_val=val;longest_t=t;}
    }
    vector<int> cand;
    for(int c=1;c<=nCust;c++) if(assign[c-1]==longest_t) cand.push_back(c);
    if(cand.empty()) return;
    int c = cand[rnd_int<int>(rng,0,(int)cand.size()-1)];
    vector<int> others;
    for(int t=0;t<m;t++) if(t!=longest_t) others.push_back(t);
    if(others.empty()) return;
    int to_t = others[rnd_int<int>(rng,0,(int)others.size()-1)];

    auto new_assign = assign; new_assign[c-1]=to_t;
    auto [nf1, nf2, _tr, _ch] = evaluate_objectives(ins, new_assign, perm, mask);
    if( (nf1 < f1_cur - 1e-9) || (fabs(nf1-f1_cur)<=1e-9 && nf2 < f2_cur - 1e-9) ){
        assign.swap(new_assign);
    }
}

// --------- NSGA-II 主过程（单种子） ----------
struct BestResult {
    vector<vector<int>> truck_routes;
    vector<map<pair<int,int>, vector<vector<int>>>> chains_all;
    double f1=0, f2=0;
};
static BestResult nsga2_optimize(const Instance& ins, int pop_size, int generations, int seed,
                                 bool verbose=true, int log_every=10){
    std::mt19937_64 rng(seed);
    int nCust=ins.num_customers, m=ins.trucks;

    if(verbose){
        stringstream ss;
        ss<<"[NSGA-II:C++] instance="<<ins.name<<" customers="<<nCust<<" trucks="<<m
          <<" pop="<<pop_size<<" gens="<<generations<<" seed="<<seed;
        logline(ss.str());
    }

    vector<vector<int>> assign, perms;
    vector<vector<char>> masks;
    init_population_smart(ins, pop_size, seed, assign, perms, masks);

    vector<double> F1(pop_size,0.0), F2(pop_size,0.0);

    auto eval_pop=[&](){
        for(int i=0;i<pop_size;i++){
            auto [f1,f2,_,__] = evaluate_objectives(ins, assign[i], perms[i], masks[i]);
            F1[i]=f1; F2[i]=f2;
        }
    };

    auto t0 = chrono::high_resolution_clock::now();
    eval_pop();
    auto [ranks, fronts] = fast_non_dominated_sort(F1,F2);
    if(verbose){
        int best0 = int(min_element(F1.begin(),F1.end()) - F1.begin());
        stringstream ss;
        ss<<"[Gen "<<0<<"] best_f1="<<F1[best0]<<" best_f2="<<F2[best0]<<" front1="<<fronts[0].size();
        auto t1=chrono::high_resolution_clock::now();
        ss<<" elapsed="<<chrono::duration<double>(t1-t0).count()<<"s";
        logline(ss.str());
    }

    for(int gen=1; gen<=generations; ++gen){
        vector<double> crowd(pop_size,0.0);
        {
            auto F = fronts;
            for(const auto &fr: F){
                if(fr.empty()) continue;
                auto dist = crowding_distance(F1,F2,fr);
                for(size_t i=0;i<fr.size();++i) crowd[fr[i]]=dist[i];
            }
        }
        auto sel_idx = tournament_select(ranks,crowd,pop_size,rng);
        vector<int> a_idx(sel_idx.begin(), sel_idx.begin()+pop_size/2);
        vector<int> b_idx(sel_idx.begin()+pop_size/2, sel_idx.end());

        vector<vector<int>> child_assign; child_assign.reserve(pop_size);
        vector<vector<int>> child_perm; child_perm.reserve(pop_size);
        vector<vector<char>> child_mask; child_mask.reserve(pop_size);

        for(size_t z=0; z<a_idx.size(); ++z){
            int ai=a_idx[z], bi=b_idx[z];
            auto AB = uniform_crossover_assign(assign[ai], assign[bi], rng);
            auto MB = uniform_crossover_mask(masks[ai], masks[bi], rng);
            auto CP  = ox_crossover(perms[ai], perms[bi], rng);
            auto CP2 = ox_crossover(perms[bi], perms[ai], rng);
            swap_mutation(CP, 0.25, rng);
            swap_mutation(CP2, 0.25, rng);
            mutate_assign(AB.first,  m, 0.06, rng);
            mutate_assign(AB.second, m, 0.06, rng);
            bitflip(MB.first,  0.06, rng);
            bitflip(MB.second, 0.06, rng);

            child_assign.push_back(move(AB.first));
            child_assign.push_back(move(AB.second));
            child_perm.push_back(move(CP));
            child_perm.push_back(move(CP2));
            child_mask.push_back(move(MB.first));
            child_mask.push_back(move(MB.second));
        }
        child_assign.resize(pop_size);
        child_perm.resize(pop_size);
        child_mask.resize(pop_size);

        // 轻量强化 top 20%
        int topK = max(1, pop_size/5);
        vector<int> order(pop_size); iota(order.begin(),order.end(),0);
        // 按 (F1,F2) 词典序
        sort(order.begin(),order.end(),[&](int i,int j){
            if(F1[i]!=F1[j]) return F1[i]<F1[j];
            return F2[i]<F2[j];
        });
        for(int ii=0; ii<topK; ++ii){
            int ti=order[ii];
            intensify_once(ins, child_assign[ti], child_perm[ti], child_mask[ti], rng);
        }

        // 评估子代
        vector<double> off_F1(pop_size,0.0), off_F2(pop_size,0.0);
        for(int i=0;i<pop_size;i++){
            auto [f1,f2,_,__] = evaluate_objectives(ins, child_assign[i], child_perm[i], child_mask[i]);
            off_F1[i]=f1; off_F2[i]=f2;
        }

        // 环境选择
        vector<vector<int>> all_assign = assign; all_assign.insert(all_assign.end(), child_assign.begin(), child_assign.end());
        vector<vector<int>> all_perm   = perms;  all_perm.insert(all_perm.end(),   child_perm.begin(),   child_perm.end());
        vector<vector<char>> all_mask  = masks;  all_mask.insert(all_mask.end(),   child_mask.begin(),   child_mask.end());
        vector<double> all_F1=F1; all_F1.insert(all_F1.end(), off_F1.begin(), off_F1.end());
        vector<double> all_F2=F2; all_F2.insert(all_F2.end(), off_F2.begin(), off_F2.end());

        auto [r2, fronts2] = fast_non_dominated_sort(all_F1, all_F2);

        vector<int> new_idx; new_idx.reserve(pop_size);
        for(auto &fr: fronts2){
            if((int)new_idx.size() + (int)fr.size() <= pop_size){
                new_idx.insert(new_idx.end(), fr.begin(), fr.end());
            }else{
                auto dist = crowding_distance(all_F1, all_F2, fr);
                vector<int> ord(fr.size()); iota(ord.begin(),ord.end(),0);
                sort(ord.begin(),ord.end(),[&](int i,int j){return dist[i]>dist[j];});
                int need = pop_size - (int)new_idx.size();
                for(int i=0;i<need;i++) new_idx.push_back(fr[ord[i]]);
                break;
            }
        }

        vector<vector<int>> new_assign(pop_size), new_perms(pop_size);
        vector<vector<char>> new_masks(pop_size);
        vector<double> new_F1(pop_size), new_F2(pop_size);
        for(int i=0;i<pop_size;i++){
            int k=new_idx[i];
            new_assign[i].swap(all_assign[k]);
            new_perms[i].swap(all_perm[k]);
            new_masks[i].swap(all_mask[k]);
            new_F1[i]=all_F1[k]; new_F2[i]=all_F2[k];
        }
        assign.swap(new_assign); perms.swap(new_perms); masks.swap(new_masks); F1.swap(new_F1); F2.swap(new_F2);

        if(verbose && ((gen%log_every==0) || (gen==generations))){
            int best = int(min_element(F1.begin(),F1.end()) - F1.begin());
            double best_f1=F1[best], best_f2=F2[best];
            double mean_f1=accumulate(F1.begin(),F1.end(),0.0)/pop_size;
            double mean_f2=accumulate(F2.begin(),F2.end(),0.0)/pop_size;
            auto [_,frs] = fast_non_dominated_sort(F1,F2);
            auto t1=chrono::high_resolution_clock::now();
            stringstream ss;
            ss<<"[Gen "<<gen<<"] best_f1="<<best_f1<<" best_f2="<<best_f2
              <<" mean=("<<mean_f1<<","<<mean_f2<<") front1="<<(frs.empty()?0:frs[0].size())
              <<" elapsed="<<chrono::duration<double>(t1-t0).count()<<"s";
            logline(ss.str());
        }
    }

    int best = 0;
    for(int i=1;i<pop_size;i++){
        if( (F1[i] < F1[best]-1e-12) ||
            (fabs(F1[i]-F1[best])<=1e-12 && F2[i] < F2[best]-1e-12) ) best=i;
    }
    auto [bf1,bf2, tr, ch] = evaluate_objectives(ins, assign[best], perms[best], masks[best]);
    BestResult br; br.f1=bf1; br.f2=bf2; br.truck_routes=move(tr); br.chains_all=move(ch);
    return br;
}

// --------- 构建输出 ----------
struct Submission {
    vector<vector<int>> truck_route;
    vector<vector<vector<int>>> drone_routes; // size = m*ndr
};
static Submission build_submission(const Instance& ins,
                                   const vector<vector<int>>& best_tr,
                                   const vector<map<pair<int,int>, vector<vector<int>>>>& best_ch){
    int m=ins.trucks, ndr=ins.drones_per_truck;
    Submission S; S.truck_route=best_tr;
    S.drone_routes.assign(m*ndr, {});
    auto gdi=[&](int t,int d){ return t*ndr+d; };

    for(int t=0;t<m;t++){
        const auto &reduced = best_tr[t];
        unordered_map<int,int> idx_map; idx_map.reserve(reduced.size());
        for(int i=0;i<(int)reduced.size();++i) idx_map[reduced[i]]=i;
        int cur=0;
        for(const auto& kv: best_ch[t]){
            int a=kv.first.first, b=kv.first.second;
            auto ia=idx_map.find(a), ib=idx_map.find(b);
            if(ia==idx_map.end() || ib==idx_map.end()) continue;
            int st_idx=ia->second, ed_idx=ib->second;
            for(const auto& chain: kv.second){
                vector<int> sortie; sortie.reserve(chain.size()+2);
                sortie.push_back(st_idx);
                for(int x: chain) sortie.push_back(x);
                sortie.push_back(ed_idx);
                S.drone_routes[gdi(t,cur)].push_back(sortie);
                cur = (cur+1)%ndr;
            }
        }
    }
    return S;
}
static void write_output_file(const string& out_fp, const Submission& subm){
    ofstream fout(out_fp);
    if(!fout) throw runtime_error("Cannot write file: "+out_fp);
    fout<<"truck_route = [";
    for(size_t i=0;i<subm.truck_route.size();++i){
        if(i) fout<<", ";
        fout<<"[";
        for(size_t j=0;j<subm.truck_route[i].size();++j){
            if(j) fout<<", ";
            fout<<subm.truck_route[i][j];
        }
        fout<<"]";
    }
    fout<<"]\n";
    for(size_t k=0;k<subm.drone_routes.size();++k){
        fout<<"drone_routes"<<(k+1)<<" = [";
        for(size_t i=0;i<subm.drone_routes[k].size();++i){
            if(i) fout<<", ";
            fout<<"[";
            for(size_t j=0;j<subm.drone_routes[k][i].size();++j){
                if(j) fout<<", ";
                fout<<subm.drone_routes[k][i][j];
            }
            fout<<"]";
        }
        fout<<"]\n";
    }
}

// --------- CLI ----------
struct Args {
    string input_file;
    string out_file; // 可选
    int pop=80;
    int gens=240;
    int seed=42;
    bool verbose=true;
};
static Args parse_args(int argc, char** argv){
    Args a;
    for(int i=1;i<argc;i++){
        string s=argv[i];
        auto need=[&](const string& k){ return s==k && i+1<argc; };
        if(s=="--input" && i+1<argc){ a.input_file=argv[++i]; }
        else if(s=="--out" && i+1<argc){ a.out_file=argv[++i]; }
        else if(s=="--pop" && i+1<argc){ a.pop=stoi(argv[++i]); }
        else if(s=="--gens" && i+1<argc){ a.gens=stoi(argv[++i]); }
        else if(s=="--seed" && i+1<argc){ a.seed=stoi(argv[++i]); }
        else if(s=="--quiet"){ a.verbose=false; }
        else if(s=="-h"||s=="--help"){
            cout<<"Usage: uagd_nsga --input <instance.txt> [--out InstanceX.txt] [--pop 80] [--gens 240] [--seed 42] [--quiet]\n";
            exit(0);
        }
    }
    if(a.input_file.empty()){ cerr<<"ERROR: --input <file> is required\n"; exit(1); }
    if(a.out_file.empty()){
        // 缺省命名：尝试从输入名提取数字，否则用 base 名
        string base = std::filesystem::path(a.input_file).filename().string();
        string digits; for(char c: base) if(isdigit((unsigned char)c)) digits.push_back(c);
        if(!digits.empty()) a.out_file = "Instance"+digits+".txt";
        else a.out_file = "Instance_output.txt";
    }
    return a;
}

int main(int argc, char** argv){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    Args args = parse_args(argc, argv);
    Instance ins = parse_instance(args.input_file);

    auto best = nsga2_optimize(ins, args.pop, args.gens, args.seed, args.verbose, 10);
    auto subm = build_submission(ins, best.truck_routes, best.chains_all);

    write_output_file(args.out_file, subm);

    // 控制台输出摘要（JSON风格）
    cout<<std::fixed<<setprecision(6);
    cout<<"{\n";
    cout<<"  \"input\": \""<<ins.name<<"\",\n";
    cout<<"  \"output\": \""<<args.out_file<<"\",\n";
    cout<<"  \"makespan\": "<<best.f1<<",\n";
    cout<<"  \"total_movement_time\": "<<best.f2<<",\n";
    cout<<"  \"pop\": "<<args.pop<<", \"gens\": "<<args.gens<<", \"seed\": "<<args.seed<<"\n";
    cout<<"}\n";
    return 0;
}
