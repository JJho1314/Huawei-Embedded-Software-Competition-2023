#include <algorithm>
#include <cstring>
#include <iostream>
#include <iterator>
#include <map>
#include <queue>
#include <string>
#include <vector>
// #include "tools_mem_used.h"
using namespace ::std;
#define max_N 1350
#define max_P 80
#define INF 0x3f3f3f3f // 无穷大

// 是节点数量N、连边数量M、业务数量T、单边通道数量P、最大衰减距离D。
//  (2 ≤ N, M ≤ 5000; 2 ≤ T ≤ 10,000; 2 ≤ P ≤ 80; 2 ≤ D ≤ 1000)
int M; // 连边数量
int P; // 单边通道数量
int D; // 最大衰减距离D
int N; // 节点数量

//-------------------------------------------------------------------------
// 存储边
class edge
{
public:
    int src;            // 边的起点
    int dst;            // 边的终点
    int weight;         // 边的权重(距离)
    int idx;            // 边的idx，按输入顺序排，从0开始
    int ch_occ_num = 0; // 通道占用数量
    vector<int> stat;
    edge(int s, int d, int w, int i, int p_num)
        : src(s), dst(d), weight(w), idx(i)
    {
        for (int i = 0; i < p_num; i++)
        {
            stat.push_back(0);
        }
    }
    edge(int s, int d, int w, int i) : src(s), dst(d), weight(w), idx(i) {}
    bool operator<(const edge a) const { return weight < a.weight; }
};

// 无向图
vector<edge> Graph[max_N];

// 根据idx，存入对于的边，方便调用
vector<edge> edge_idx;

//-----------------------------------业务类--------------------------------------
class transaction
{
public:
    int begin;             // 业务起点
    int end;               // 业务终点
    int channel_idx = -1;  // 业务占的通道,初始时刻为-1
    int magnifier_num = 0; // 放大器数量
    vector<int> verticle;  // 存储业务经过的顶点(包括源点)    (end,...,begin)
    vector<int> path;      // 存储业务经过的边的index  (边也是逆序)
    vector<int> magnifier; // 用来存储放置放大器的结点的index

    int idx; // 业务index
    transaction(int b, int e, int i) : begin(b), end(e), idx(i) {}
};

//---------------------------------用于存储业务的-----------------------------------
vector<transaction> tranx;

class Add_edge
{
public:
    int num = 0;
    vector<pair<int, int>> begin_end; // 存储新边的起始结点
};

Add_edge add_edge;

//-------------------------------------------------------------------------
void init(); // 输入初始化
void dijkstra(transaction &tran, int begin, int end);
// vector<int> channel_is_spare(int edge_index, vector<edge> edge_tmp); // 判断某条边的通道是否占满（若已满返回1，否则返回0）
// vector<int> vectors_intersection(vector<int> v1, vector<int> v2);    // 计算通道交集

int isPathGo(vector<int> path, vector<edge> edge_tmp)
{
    for (int i = 0; i < P; i++)
    {
        int sum = path.size();
        for (int j = 0; j < path.size(); j++)
        {
            if (edge_tmp[path[j]].stat[i] == -1)
            {
                sum--;
            }
        }

        if (sum == path.size())
            return i; // channel_idx
    }
    return -1;
}

//----------------------------------------------------------按照通道数量进行排序-------------------------------------------------------
vector<int> sortedByChannelNumber(vector<int> path)
{
    for (int i = 0; i < path.size() - 1; i++)
    {
        for (int j = i + 1; j < path.size(); j++)
        {
            if (edge_idx[path[i]].ch_occ_num < edge_idx[path[j]].ch_occ_num)
            {
                int tmp = 0;
                tmp = path[j];
                path[j] = path[i];
                path[i] = tmp;
            }
        }
    }
    return path;
}

//-------------------------------------------------------------判断best_path中index是否都在edge_idx中----------------------------------------------
int is_add_edge(vector<int> best_path)
{
    int max_index = -1;
    for (auto index : best_path)
    {
        max_index = max(max_index, index);
    }
    return max_index - (M - 1);
}

void init()
{
    int T; // 业务数量

    cin >> N >> M >> T >> P >> D;

    for (int i = 0; i < M; i++)
    {
        int src, dst, weight;
        cin >> src >> dst >> weight;
        Graph[src].push_back(edge(src, dst, weight, i, P)); // edges列表的index为图顶点
        Graph[dst].push_back(edge(dst, src, weight, i, P));
        edge_idx.push_back(edge(src, dst, weight, i, P));
    }

    for (int i = 0; i < T; i++)
    {
        int begin, end;
        cin >> begin >> end;
        tranx.push_back(transaction(begin, end, i));
    }
}

//-----------------------------------通道是否空闲-----------------------------------
// vector<int> channel_is_spare(int edge_index, vector<edge> edge_tmp)
// {
//     vector<int> ret;
//     for (int i = 0; i < P; i++)
//     {
//         if (edge_tmp[edge_index].stat[i] == 0)
//         // if (edge_channel[edge_idx].stat[i] == 0)
//         {
//             ret.push_back(i);
//         }
//     }
//     return ret;
// }

//-----------------------------------路径上是否有相交的通道出现--------------------------------------
// vector<int> vectors_intersection(vector<int> v1, vector<int> v2)
// {
//     vector<int> v;
//     set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(),
//                      back_inserter(v)); // 求交集
//     return v;
// }

//-------------------------------------------------------------bfs-----------------------------------------------------------------------

vector<vector<int>> all_path_node_v;
vector<vector<int>> all_path_v;

vector<int> best_verticle;
vector<int> best_path;
int best_channel;
int add_edge_num;
int add_magnifier_num;
vector<pair<int, int>> add_edge_verticle;

void new_bfs(int begin, int end)
{
    queue<vector<int>> q;
    queue<vector<bool>> q_st;
    queue<vector<int>> q_path;
    q.push({begin});
    q_st.push(vector<bool>(max_N, false));

    while (q.size())
    {
        // 获取队列中的第一条路径
        vector<int> path = q.front();
        q.pop();
        vector<bool> st = q_st.front();
        q_st.pop();
        vector<int> new_ret;
        if (q_path.size() > 0)
        {
            new_ret = q_path.front();
            q_path.pop();
        }
        // 路径的最后一个元素
        int last = path.back();
        st[last] = true;

        // 若路径的最后一个元素与 end 相等，说明已经找到一条路径，加入 all
        if (last == end)
        {

            // 把src和end对应的边加入到result_path中
            all_path_v.push_back(new_ret);

            all_path_node_v.push_back(path);

            if (isPathGo(new_ret, edge_idx) != -1)
            {
                best_path = new_ret;
                break;
            }
        }

        for (auto &e : Graph[last])
        {
            if (!st[e.dst] && e.weight != INF)
            {
                vector<int> next_path(path);
                next_path.push_back(e.dst);
                vector<bool> next_st(st);
                next_st[e.dst] = true;
                // 添加该条路径，等待后续 bfs
                q.push(next_path);
                q_st.push(next_st);
                if (q_path.size() > 0)
                {
                    vector<int> nest_edge(new_ret);
                    nest_edge.push_back(e.idx);
                    q_path.push(nest_edge);
                }
                else
                {
                    vector<int> nest_edge;
                    nest_edge.push_back(e.idx);
                    q_path.push(nest_edge);
                }
            }
        }
    }
}

//-------------------------------------------------------输入边和输出加边后的边------------------------------------------------

void addedge(vector<int> current_path, vector<edge> &after_egde_idx)
{
    // 按照通道数进行排序
    vector<int> new_path = sortedByChannelNumber(current_path);
    vector<edge> tem_edge_idx = after_egde_idx;
    int current_M = M;
    // 优先从占用通道数最多的开始加边
    for (int i = 0; i < new_path.size(); i++)
    {
        // 需要加边
        // add_edge_idx.push_back(new_path[i]);
        // 边数+1
        add_edge_num++;
        auto it = find(current_path.begin(), current_path.end(), new_path[i]) - current_path.begin();

        current_path[it] = current_M;
        current_M = current_M + 1;
        tem_edge_idx.push_back(edge(tem_edge_idx[new_path[i]].src,
                                    tem_edge_idx[new_path[i]].dst,
                                    tem_edge_idx[new_path[i]].weight,
                                    current_path[it], P));

        if (isPathGo(current_path, tem_edge_idx) != -1)
        {
            best_path = current_path;
            after_egde_idx = tem_edge_idx;
            break;
        }
    }
}

//------------------------------------------------------------花费计算------------------------------------------------------------------

void best_path_verticle(transaction &tran, int src, int dst,
                        vector<int> &best_path)
{
    tran.verticle.clear();
    for (int j = 0; j < best_path.size(); j++)
    {
        if (src == edge_idx[best_path[j]].src)
        {
            tran.verticle.push_back(edge_idx[best_path[j]].src);
            src = edge_idx[best_path[j]].dst;
        }
        else if (src == edge_idx[best_path[j]].dst)
        {
            tran.verticle.push_back(edge_idx[best_path[j]].dst);
            src = edge_idx[best_path[j]].src;
        }
    }
    tran.verticle.push_back(dst);
    // 得到的tran.verticle是顺序存入
}

// u是源点，求的是源点到其他顶点最短路径
void dijkstra(transaction &tran, int begin, int end)
{
    // 初始化距离数组
    vector<int> dist(N + 1, INF);
    // 初始化路径数组
    vector<int> path(N + 1, -1);
    // 初始化优先队列
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<>> q; // 队首元素距离dist一定是最小的
    // 初始化距离数组
    dist[begin] = 0;
    // 初始化路径数组
    path[begin] = begin;
    // 初始化优先队列
    q.push({0, begin});
    // Dijkstra算法
    while (!q.empty())
    {
        // 取出队首元素
        // const int v = -1;//顶点
        // const int d = -1;//距离
        auto ret = q.top();
        int d, v; // d为距离，v为顶点
        d = ret.first;
        v = ret.second;
        if (v == end)
        {
            break;
        }
        q.pop();
        // 如果当前距离大于已知距离，则跳过
        if (d > dist[v])
        {
            // cout<<v<<endl;
            continue;
        }
        // 遍历当前节点的所有邻接点
        for (auto &e : Graph[v])
        {
            // 如果当前距离加上邻接点的权重大于已知距离，则跳过
            if (dist[v] + e.weight > dist[e.dst])
                continue;

            if (dist[v] + e.weight <= dist[e.dst])
            {
                // 更新距离数组
                dist[e.dst] = dist[v] + e.weight;
                // 更新路径数组
                path[e.dst] = v;
                // 将邻接点加入优先队列
                q.push({dist[e.dst], e.dst});
            }
        }
    }

    //--------------------------将路径经过的边的index和顶点加入Class
    // transaction---------------------------------------
    int dst = end;
    int src = path[dst];
    int edge_weight = 0;
    // 遍历边，添加路径经过的边的index和顶点

    while (dst != begin)
    {
        edge_weight = dist[dst] - dist[src];
        for (int i = Graph[src].size() - 1; i >= 0; i--)
        {
            // 有多条边edge_weight相同
            if (Graph[src][i].dst == dst && Graph[src][i].weight == edge_weight)
            {
                tran.path.push_back(Graph[src][i].idx);
                break;
            }
        }
        tran.verticle.push_back(dst);
        // 节点前移
        dst = src;
        src = path[src];
    }

    tran.verticle.push_back(begin);

    //-----------------------------------------------选择通道-----------------------------------------------
    // 得到path以后 开始处理通道
    // 遍历所用边的空余通道，最后取交集
    vector<edge> tem_edge_idx;
    tem_edge_idx = edge_idx;

    // vector<int> first = channel_is_spare(tran.path[tran.path.size() - 1], tem_edge_idx);
    // vector<int> second;
    // vector<int> intersection;

    // for (int k = tran.path.size() - 2; k >= 0; k--)
    // {
    //     second = channel_is_spare(tran.path[k], tem_edge_idx);

    //     intersection = vectors_intersection(first, second);

    //     first = intersection;
    // }

    // first.clear();
    // first.shrink_to_fit();
    // second.clear();
    // second.shrink_to_fit();

    // 有空余通道
    int channel_idx = isPathGo(tran.path, tem_edge_idx);
    if (channel_idx != -1)
    {
        tran.channel_idx = channel_idx;

        for (auto &k : tran.path)
        {
            edge_idx[k].stat[channel_idx] = -1;
            edge_idx[k].ch_occ_num += 1;
        }

        //------------------------------------------------------放大器--------------------------------------------
        // 计算加放大器的数量
        // 考虑路径 src----(weight1)-----end
        int weight = 0;
        int D_current = D;
        for (int k = tran.verticle.size() - 1; k > 0; k--)
        {
            weight = dist[tran.verticle[k - 1]] - dist[tran.verticle[k]];
            D_current -= weight;
            if (D_current < 0)
            {
                tran.magnifier_num++;
                tran.magnifier.push_back(tran.verticle[k]);
                D_current = D - weight;
            }
        }
    }
    // 无同一类型的空余通道
    else
    {
        //------------------------加边-------------------------------------
        // 加边的情况：
        // 不能在一条路径上的所有边中找到同一个类型的通道

        // 利用DFS遍历begin到end的所有路径并存储，查找是该路径是否有同一类型的空余通道
        // 若有，找一条边数最短的路径进行切换，不需要加边:
        //      tran.path.swap(new_path);
        //      tran.verticle.swap(new_verticle);
        // int visit[N] = {0};
        best_path.clear();
        all_path_node_v.clear();
        // all_node_v.clear();
        // result_verticle.clear();

        // path verticle的存储顺序与之前不同 顺序src——>......——>end

        // 获得all_path,all_verticle
        // new_bfs(begin, end);

        // 获得best_path
        // cost(all_path_v, tem_edge_idx);
        reverse(tran.path.begin(), tran.path.end());
        addedge(tran.path, tem_edge_idx);

        // 获得best_path以后 更新edge_idx
        edge_idx = tem_edge_idx;

        // vector<int> first = channel_is_spare(best_path[0], tem_edge_idx);
        // vector<int> second;
        // vector<int> intersection;

        // for (int k = 1; k < best_path.size(); k++)
        // {
        //     second = channel_is_spare(best_path[k], tem_edge_idx);

        //     intersection = vectors_intersection(first, second);

        //     first = intersection;
        // }

        int channel_idx = isPathGo(tran.path, tem_edge_idx);

        // 有空余通道,只需修改路径 不需要加边
        // 判断best_path对应的index是否在原始的path里面，返回0表示没有增加新的边
        if (channel_idx != -1 && is_add_edge(best_path) == 0)
        {
            // int channel_idx = intersection[0];
            tran.channel_idx = channel_idx;
            for (auto &k : best_path)
            {
                edge_idx[k].stat[channel_idx] = -1;
                // edge_channel[k].stat[channel_idx] = 1;
                edge_idx[k].ch_occ_num += 1;
            }

            // 将best_path路径上所经过的顶点push进入tran.verticle
            best_path_verticle(tran, begin, end, best_path);

            // 计算放大器数量
            add_magnifier_num = 0;
            int weight = 0;
            int D_current = D;

            for (int i = 0; i < best_path.size(); i++)
            {
                weight = edge_idx[best_path[i]].weight;
                D_current -= weight;
                if (D_current < 0)
                {
                    // add_magnifier_num++;
                    tran.magnifier_num++;
                    tran.magnifier.push_back(tran.verticle[i]);
                    D_current = D - weight;
                }
            }

            tran.path.clear();
            // reverse(best_path.begin(), best_path.end());
            tran.path = best_path;
        }
        else
        { // 表示增加边
            // 边上的通道数量
            int edge_spare_channel_cnt[P] = {0};

            // 遍历路径上的所有边
            // 记录每个通道的空闲数
            for (int i = 0; i < best_path.size(); i++)
            {
                for (int j = 0; j < P; j++)
                {
                    if (edge_idx[best_path[i]].stat[j] == 0)
                    {
                        edge_spare_channel_cnt[j]++;
                    }
                }
            }

            // 初始化
            int max_cnt = edge_spare_channel_cnt[0];
            int channel_idx = 0;
            int cnt;
            // 选择共通空闲最多的通道
            for (int i = 1; i < P; i++)
            {
                cnt = edge_spare_channel_cnt[i];
                if (cnt > max_cnt)
                {
                    max_cnt = cnt;
                    channel_idx = i;
                }
            }

            tran.channel_idx = channel_idx;

            // 循环最好的路径，找到对应的空闲通道
            int max_index = -1;

            for (int i = 0; i < best_path.size(); i++)
            {
                // 得到最好路径中最大的index
                if (best_path[i] >= M)
                {
                    // 得到这条新加边的两个端点
                    int ver1 = edge_idx[best_path[i]].src;
                    int ver2 = edge_idx[best_path[i]].dst;
                    // 加边
                    // 同时图graph和edge_idx也要改变,还要修改best_path中边的index
                    add_edge.num++;
                    add_edge.begin_end.push_back(make_pair(ver1, ver2));
                    // 更新全局graph
                    Graph[ver1].push_back(edge(ver1, ver2,
                                               edge_idx[best_path[i]].weight,
                                               M - 1 + add_edge.num, P));
                    Graph[ver2].push_back(edge(ver2, ver1,
                                               edge_idx[best_path[i]].weight,
                                               M - 1 + add_edge.num, P));

                    // edge_idx.push_back(edge(ver1, ver2,
                    // edge_idx[best_path[i]].weight, M - 1 + add_edge.num, P));
                    // edge_idx.push_back(
                    //     edge(ver2, ver1, edge_idx[best_path[i]].weight, M - 1
                    //     + add_edge.num));
                    // best_path[i] = M - 1 + add_edge.num;
                }
                // 更新channel
                edge_idx[best_path[i]].stat[channel_idx] = -1;
            }

            best_path_verticle(tran, begin, end, best_path);

            // 计算放大器数量
            add_magnifier_num = 0;
            int weight = 0;
            int D_current = D;
            for (int i = 0; i < best_path.size(); i++)
            {
                weight = edge_idx[best_path[i]].weight;
                D_current -= weight;
                if (D_current < 0)
                {
                    // add_magnifier_num++;
                    tran.magnifier_num++;
                    tran.magnifier.push_back(tran.verticle[i]);
                    D_current = D - weight;
                }
            }
            // 更新M长度
            M = edge_idx.size();
        }

        // 清空之前存储的path
        tran.path.clear();
        // 将best_path 逆置
        reverse(best_path.begin(), best_path.end());
        tran.path = best_path;
    }
}

int main()
{
    init();

    // 对每个顶点的边按距离进行排序
    for (int i = 0; i < tranx.size(); i++)
    {
        int begin = tranx[i].begin;
        int end = tranx[i].end;
        dijkstra(tranx[i], begin, end);
    }

    cout << add_edge.num << endl;

    if (add_edge.begin_end.size() > 0)
    {
        for (int i = 0; i < add_edge.begin_end.size(); i++)
        {
            // 输出加边的源点和终点
            cout << add_edge.begin_end[i].first << " "
                 << add_edge.begin_end[i].second << endl;
        }
    }

    // 输出交易的通道编号+边数+放大器数+经过边的index+依次经过的放大器所在节点的编号
    for (int i = 0; i < tranx.size(); i++)
    {
        cout << tranx[i].channel_idx << " " << tranx[i].path.size() << " "
             << tranx[i].magnifier_num << " ";

        for (int j = tranx[i].path.size() - 1; j >= 0; j--)
        {
            cout << tranx[i].path[j] << " ";
        }

        for (auto &k : tranx[i].magnifier)
        {
            cout << k << " ";
        }
        cout << endl;
    }
    // float mem_used_mb = (Common_tools::get_RSS_Mb());
    // printf("Memory used (Mb): %f Mb \n", mem_used_mb);
    return 0;
}