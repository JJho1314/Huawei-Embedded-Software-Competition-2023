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
int M;     // 连边数量
int P;     // 单边通道数量
int D;     // 最大衰减距离D
int N;     // 节点数量
int raw_M; // 初始边数量

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

void dijkstra(int begin, int end, int i);
void addedge(vector<int> current_path, vector<edge> &after_egde_idx, int num);

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

//-------------------------------------------------------------bfs-----------------------------------------------------------------------

// vector<int> best_verticle;
vector<int> best_path;
vector<int> best_node;

//-------------------------------------------------------输入边和输出加边后的边------------------------------------------------

void addedge(vector<int> current_path, vector<edge> &after_egde_idx, vector<int> &pre_result)
{
    // 按照通道数进行排序
    vector<int> new_path = sortedByChannelNumber(current_path);
    vector<edge> tem_edge_idx = after_egde_idx;
    int add_edge_num;
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

        int channel_idx = isPathGo(current_path, tem_edge_idx);
        if (channel_idx != -1)
        {
            best_path = current_path;
            // cout << channel_idx << ' ';
            pre_result.push_back(channel_idx);
            after_egde_idx = tem_edge_idx;
            break;
        }
    }
}

void new_bfs(int begin, int end, vector<edge> &edge_tmp, vector<int> &pre_result)
{
    queue<vector<int>> q;
    queue<vector<bool>> q_st;
    queue<vector<int>> q_path;
    vector<vector<int>> all_path_node_v;
    vector<vector<int>> all_path_v;
    all_path_v.clear();
    all_path_node_v.clear();

    q.push({begin});
    q_st.push(vector<bool>(max_N, false));
    int path_num = 0;
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

            int channel_idx = isPathGo(new_ret, edge_tmp);

            if (channel_idx != -1)
            {
                best_path = new_ret;
                best_node = path;
                // std::cout << channel_idx << ' ';
                pre_result.push_back(channel_idx);
                for (auto &k : best_path)
                {
                    edge_tmp[k].stat[channel_idx] = -1;
                    edge_tmp[k].ch_occ_num += 1;
                }
                break;
            }

            if (path_num > 0)
            {
                best_path = all_path_v[1];
                best_node = all_path_node_v[1];
                addedge(best_path, edge_tmp, pre_result);
                channel_idx = isPathGo(best_path, edge_tmp);
                for (auto &k : best_path)
                {
                    edge_tmp[k].stat[channel_idx] = -1;
                    edge_tmp[k].ch_occ_num += 1;
                }
                break;
            }

            path_num = path_num + 1;
        }

        int number = 0;
        for (auto &e : Graph[last])
        {
            if (number >= 4)
            {
                break;
            }
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
            number++;
        }
    }
}

// u是源点，求的是源点到其他顶点最短路径
void dijkstra(int begin, int end, vector<int> &pre_result)
{
    //-----------------------------------------------选择通道-----------------------------------------------
    // 得到path以后 开始处理通道
    // 遍历所用边的空余通道，最后取交集

    vector<int> tran_magnifier;
    tran_magnifier.clear();
    tran_magnifier.shrink_to_fit();

    best_path.clear();
    best_node.clear();

    // path verticle的存储顺序与之前不同 顺序src——>......——>end

    // 获得all_path,all_verticle
    new_bfs(begin, end, edge_idx, pre_result);

    // 计算放大器数量
    int weight = 0;
    int D_current = D;
    for (int i = 0; i < best_path.size(); i++)
    {
        weight = edge_idx[best_path[i]].weight;
        D_current -= weight;
        if (D_current < 0)
        {
            // add_magnifier_num++;
            tran_magnifier.push_back(best_node[i]);
            D_current = D - weight;
        }
    }
    M = edge_idx.size();
    pre_result.push_back(best_path.size());
    pre_result.push_back(tran_magnifier.size());
    // std::cout << best_path.size() << ' ' << tran_magnifier.size() << ' ';

    for (int i = 0; i < best_path.size(); i++)
    {
        pre_result.push_back(best_path[i]);
        // std::cout << best_path[i] << ' ';
    }

    for (int i = 0; i < tran_magnifier.size(); i++)
    {
        pre_result.push_back(tran_magnifier[i]);
        // std::cout << tran_magnifier[i] << ' ';
    }
}

int main()
{
    int T; // 业务数量

    cin >> N >> M >> T >> P >> D;
    raw_M = M;
    vector<pair<int, int>> yewu;
    vector<vector<int>> result;
    vector<int> pre_result;
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
        yewu.emplace_back(make_pair(begin, end));
    }

    // 对每个顶点的边按距离进行排序
    for (int i = 0; i < yewu.size(); i++)
    {
        vector<int>().swap(pre_result);
        int begin = yewu[i].first;
        int end = yewu[i].second;
        dijkstra(begin, end, pre_result);
        result.push_back(pre_result);
    }

    cout << M - raw_M << endl;

    for (int i = raw_M; i < edge_idx.size(); i++)
    {
        cout << edge_idx[i].src << ' ' << edge_idx[i].dst << endl;
    }

    for (int i = 0; i < result.size(); i++)
    {
        for (int j = 0; j < result[i].size(); j++)
        {
            cout << result[i][j] << " ";
        }
        cout << endl;
    }

    // float mem_used_mb = (Common_tools::get_RSS_Mb());
    // printf("Memory used (Mb): %f Mb \n", mem_used_mb);
    return 0;
}