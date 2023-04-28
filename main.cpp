#include <algorithm>
#include <cstring>
#include <iostream>
#include <iterator>
#include <map>
#include <queue>
#include <string>
#include <vector>
using namespace ::std;
#define max_N 1350
#define max_P 80
#define INF 0x3f3f3f3f // 无穷大

short int M; // 连边数量
short int P; // 单边通道数量
short int D; // 最大衰减距离D
short int N; // 节点数量

//-------------------------------------------------------------------------
// 存储边
class edge
{
public:
    short int src;            // 边的起点
    short int dst;            // 边的终点
    short int weight;         // 边的权重(距离)
    short int idx;            // 边的idx，按输入顺序排，从0开始
    short int ch_occ_num = 0; // 通道占用数量
    vector<short int> stat;
    edge(short int s, short int d, short int w, short int i, short int p_num)
        : src(s), dst(d), weight(w), idx(i)
    {
        for (short int i = 0; i < p_num; i++)
        {
            stat.emplace_back(0);
        }
    }
    edge(short int s, short int d, short int w, short int i) : src(s), dst(d), weight(w), idx(i) {}
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
    short int begin;             // 业务起点
    short int end;               // 业务终点
    short int channel_idx = -1;  // 业务占的通道,初始时刻为-1
    short int magnifier_num = 0; // 放大器数量
    vector<short int> verticle;  // 存储业务经过的顶点(包括源点)    (end,...,begin)
    vector<short int> path;      // 存储业务经过的边的index  (边也是逆序)
    vector<short int> magnifier; // 用来存储放置放大器的结点的index

    short int idx; // 业务index
    transaction(short int b, short int e, short int i) : begin(b), end(e), idx(i) {}
};

//---------------------------------用于存储业务的-----------------------------------
vector<transaction> tranx;

class Add_edge
{
public:
    short int num = 0;
    vector<pair<short int, short int>> begin_end; // 存储新边的起始结点
};

Add_edge add_edge;

void init();
void dijkstra(transaction &tran, short int begin, short int end);

short int isPathGo(vector<short int> path, vector<edge> &edge_tmp)
{
    for (short int i = 0; i < P; i++)
    {
        short int sum = 0;
        for (short int j = 0; j < path.size(); j++)
        {
            if (edge_tmp[path[j]].stat[i] == -1)
            {
                break;
            }
            sum = sum + 1;
        }

        if (sum == path.size())
            return i; // channel_idx
    }
    return -1;
}

int partition(vector<short int> &path, int left, int right)
{
    int pivot = path[right];
    int i = left - 1;
    for (int j = left; j < right; j++)
    {
        if (edge_idx[path[j]].ch_occ_num >= edge_idx[pivot].ch_occ_num)
        {
            i++;
            swap(path[i], path[j]);
        }
    }
    swap(path[i + 1], path[right]);
    return i + 1;
}

void quickSort(vector<short int> &path, int left, int right)
{
    if (left < right)
    {
        int pivot = partition(path, left, right);
        quickSort(path, left, pivot - 1);
        quickSort(path, pivot + 1, right);
    }
}

vector<short int> sortedByChannelNumber(vector<short int> path)
{
    quickSort(path, 0, path.size() - 1);
    return path;
}

short int is_add_edge(vector<short int> &best_path)
{
    short int max_index = -1;
    for (auto index : best_path)
    {
        max_index = max(max_index, index);
    }
    return max_index - (M - 1);
}

void init()
{
    short int T; // 业务数量

    cin >> N >> M >> T >> P >> D;

    for (short int i = 0; i < M; i++)
    {
        short int src, dst, weight;
        cin >> src >> dst >> weight;
        Graph[src].emplace_back(edge(src, dst, weight, i, P)); // edges列表的index为图顶点
        Graph[dst].emplace_back(edge(dst, src, weight, i, P));
        edge_idx.emplace_back(edge(src, dst, weight, i, P));
    }

    for (short int i = 0; i < T; i++)
    {
        short int begin, end;
        cin >> begin >> end;
        tranx.emplace_back(transaction(begin, end, i));
    }
}
vector<vector<short int>> all_path_node_v;
vector<vector<short int>> all_path_v;

vector<short int> best_verticle;
vector<short int> best_path;
short int best_channel;
short int add_edge_num;
short int add_magnifier_num;
vector<pair<short int, short int>> add_edge_verticle;

void new_bfs(short int begin, short int end, transaction &tran)
{
    queue<vector<short int>> q;
    queue<vector<bool>> q_st;
    queue<vector<short int>> q_path;
    q.push({begin});
    q_st.push(vector<bool>(max_N, false));

    while (q.size())
    {
        // 获取队列中的第一条路径
        vector<short int> path = q.front();
        q.pop();
        vector<bool> st = q_st.front();
        q_st.pop();
        vector<short int> new_ret;
        if (q_path.size() > 0)
        {
            new_ret = q_path.front();
            q_path.pop();
        }
        // 路径的最后一个元素
        short int last = path.back();
        st[last] = true;

        // 若路径的最后一个元素与 end 相等，说明已经找到一条路径，加入 all
        if (last == end)
        {
            // 把src和end对应的边加入到result_path中
            all_path_v.emplace_back(new_ret);

            all_path_node_v.emplace_back(path);

            if (isPathGo(new_ret, edge_idx) != -1)
            {
                best_path = new_ret;
                tran.verticle = path;
                tran.path = new_ret;
                break;
            }
        }

        short int number = 0;
        for (auto &e : Graph[last])
        {
            if (number > 4)
            {
                break;
            }
            if (!st[e.dst] && e.weight != INF)
            {
                vector<short int> next_path(path);
                next_path.emplace_back(e.dst);
                vector<bool> next_st(st);
                next_st[e.dst] = true;
                // 添加该条路径，等待后续 bfs
                q.push(next_path);
                q_st.push(next_st);
                if (q_path.size() > 0)
                {
                    vector<short int> nest_edge(new_ret);
                    nest_edge.emplace_back(e.idx);
                    q_path.push(nest_edge);
                }
                else
                {
                    vector<short int> nest_edge;
                    nest_edge.emplace_back(e.idx);
                    q_path.push(nest_edge);
                }
            }
            number++;
        }
    }
}

void addedge(vector<short int> current_path, vector<edge> &after_egde_idx)
{
    // 按照通道数进行排序
    vector<short int> new_path = current_path;
    short int current_M = M;
    // 优先从占用通道数最多的开始加边
    for (short int i = 0; i < new_path.size(); i++)
    {
        add_edge_num++;
        auto it = find(current_path.begin(), current_path.end(), new_path[i]) - current_path.begin();

        current_path[it] = current_M;
        current_M = current_M + 1;
        after_egde_idx.emplace_back(edge(after_egde_idx[new_path[i]].src,
                                         after_egde_idx[new_path[i]].dst,
                                         after_egde_idx[new_path[i]].weight,
                                         current_path[it], P));

        if (isPathGo(current_path, after_egde_idx) != -1)
        {
            best_path = current_path;
            break;
        }
    }
    vector<short int>().swap(new_path);
}

//------------------------------------------------------------花费计算------------------------------------------------------------------

void best_path_verticle(transaction &tran, short int src, short int dst,
                        vector<short int> &best_path)
{
    // tran.verticle.clear();
    vector<short int>().swap(tran.verticle);
    for (short int j = 0; j < best_path.size(); j++)
    {
        if (src == edge_idx[best_path[j]].src)
        {
            tran.verticle.emplace_back(edge_idx[best_path[j]].src);
            src = edge_idx[best_path[j]].dst;
        }
        else if (src == edge_idx[best_path[j]].dst)
        {
            tran.verticle.emplace_back(edge_idx[best_path[j]].dst);
            src = edge_idx[best_path[j]].src;
        }
    }
    tran.verticle.emplace_back(dst);
}

void dijkstra(transaction &tran, short int begin, short int end)
{
    vector<short int> dist(N + 1, INF);
    vector<short int> path(N + 1, -1);
    priority_queue<pair<short int, short int>, vector<pair<short int, short int>>, greater<>> q; // 队首元素距离dist一定是最小的
    dist[begin] = 0;
    path[begin] = begin;
    q.push({0, begin});
    while (!q.empty())
    {
        // 取出队首元素
        auto ret = q.top();
        short int d, v; // d为距离，v为顶点
        d = ret.first;
        v = ret.second;
        if (v == end)
        {
            break;
        }
        q.pop();
        if (d > dist[v])
        {
            continue;
        }
        for (auto &e : Graph[v])
        {
            if (dist[v] + e.weight > dist[e.dst])
                continue;

            if (dist[v] + e.weight <= dist[e.dst])
            {
                dist[e.dst] = dist[v] + e.weight;
                path[e.dst] = v;
                q.push({dist[e.dst], e.dst});
            }
        }
    }
    short int dst = end;
    short int src = path[dst];
    short int edge_weight = 0;

    while (dst != begin)
    {
        edge_weight = dist[dst] - dist[src];
        for (short int i = Graph[src].size() - 1; i >= 0; i--)
        {
            // 有多条边edge_weight相同
            if (Graph[src][i].dst == dst)
            {
                tran.path.emplace_back(Graph[src][i].idx);
                Graph[src][i].ch_occ_num++;
                break;
            }
        }
        tran.verticle.emplace_back(dst);
        // 节点前移
        dst = src;
        src = path[src];
    }

    tran.verticle.emplace_back(begin);

    // 有空余通道
    short int channel_idx = isPathGo(tran.path, edge_idx);
    if (channel_idx != -1)
    {
        tran.channel_idx = channel_idx;

        for (auto &k : tran.path)
        {
            edge_idx[k].stat[channel_idx] = -1;
            edge_idx[k].ch_occ_num += 1;
        }

        short int weight = 0;
        short int D_current = D;
        for (short int k = tran.verticle.size() - 1; k > 0; k--)
        {
            weight = dist[tran.verticle[k - 1]] - dist[tran.verticle[k]];
            D_current -= weight;
            if (D_current < 0)
            {
                tran.magnifier_num++;
                tran.magnifier.emplace_back(tran.verticle[k]);
                D_current = D - weight;
            }
        }
    }

    else
    {
        best_path.clear();
        all_path_node_v.clear();

        reverse(tran.path.begin(), tran.path.end());
        addedge(tran.path, edge_idx);

        channel_idx = isPathGo(best_path, edge_idx);
        tran.channel_idx = channel_idx;

        for (short int i = 0; i < best_path.size(); i++)
        {
            if (best_path[i] >= M)
            {
                // 得到这条新加边的两个端点
                short int ver1 = edge_idx[best_path[i]].src;
                short int ver2 = edge_idx[best_path[i]].dst;

                add_edge.num++;
                add_edge.begin_end.emplace_back(make_pair(ver1, ver2));
                // 更新全局graph
                Graph[ver1].emplace_back(edge(ver1, ver2,
                                              edge_idx[best_path[i]].weight,
                                              M - 1 + add_edge.num, P));
                Graph[ver2].emplace_back(edge(ver2, ver1,
                                              edge_idx[best_path[i]].weight,
                                              M - 1 + add_edge.num, P));
            }
            edge_idx[best_path[i]].stat[channel_idx] = -1;
        }
        best_path_verticle(tran, begin, end, best_path);
        add_magnifier_num = 0;
        short int weight = 0;
        short int D_current = D;
        for (short int i = 0; i < best_path.size(); i++)
        {
            weight = edge_idx[best_path[i]].weight;
            D_current -= weight;
            if (D_current < 0)
            {
                // add_magnifier_num++;
                tran.magnifier_num++;
                tran.magnifier.emplace_back(tran.verticle[i]);
                D_current = D - weight;
            }
        }

        M = edge_idx.size();
        tran.path.clear();
        reverse(best_path.begin(), best_path.end());
        tran.path = best_path;
    }
}

int main()
{
    init();
    for (short int i = 0; i < tranx.size(); i++)
    {
        short int begin = tranx[i].begin;
        short int end = tranx[i].end;
        dijkstra(tranx[i], begin, end);
    }
    cout << add_edge.num << endl;
    if (add_edge.begin_end.size() > 0)
    {
        for (short int i = 0; i < add_edge.begin_end.size(); i++)
        {

            cout << add_edge.begin_end[i].first << " "
                 << add_edge.begin_end[i].second << endl;
        }
    }
    for (short int i = 0; i < tranx.size(); i++)
    {
        cout << tranx[i].channel_idx << " " << tranx[i].path.size() << " "
             << tranx[i].magnifier_num << " ";

        for (short int j = tranx[i].path.size() - 1; j >= 0; j--)
        {
            cout << tranx[i].path[j] << " ";
        }

        for (auto &k : tranx[i].magnifier)
        {
            cout << k << " ";
        }
        cout << endl;
    }

    return 0;
}