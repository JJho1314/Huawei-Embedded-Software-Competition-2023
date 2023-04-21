#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <queue>
#include <iterator>
#include <cstring>
#include <map>
#include <unordered_map>
using namespace ::std;
#define max_N 10001
#define max_P 81
#define max_M 5001
#define INF 0x3f3f3f3f  //无穷大
#define inf 0x3f3f3f

//是节点数量N、连边数量M、业务数量T、单边通道数量P、最大衰减距离D。
// (2 ≤ N, M ≤ 5000; 2 ≤ T ≤ 10,000; 2 ≤ P ≤ 80; 2 ≤ D ≤ 1000)
int N;  //节点数量
int M;  //连边数量
int T;  //业务数量
int P;  //单边通道数量
int D;  //最大衰减距离D

//-------------------------------------------------------------------------
//记录边中的通道数状态
class channel
{
public:
  int stat[max_P] = { 0 };  //初始状态全为0  0：未被占用  1：被占用
};

channel edge_channel[max_M];  //用来保存每条边的通道状态

//-------------------------------------------------------------------------
//存储边
class edge
{
public:
  int src;             //边的起点
  int dst;             //边的终点
  int weight;          //边的权重(距离)
  int idx;             //边的idx，按输入顺序排，从0开始
  int ch_occ_num = 0;  //通道占用数量
  // vector<int> magnifier;  // 放大器

  edge(int s, int d, int w, int i) : src(s), dst(d), weight(w), idx(i)
  {
  }

  bool operator<(const edge a) const
  {
    return weight < a.weight;
  }
};
//无向图
// edge Graph[max_N][max_N];
vector<edge> Graph[max_N];

//根据idx，存入对于的边，方便调用
vector<edge> edge_idx[max_M];

//-------------------------------------------------------------------------
class transaction
{
public:
  int begin;              //业务起点
  int end;                //业务终点
  int channel_idx = -1;   //业务占的通道,初始时刻为-1
  int magnifier_num = 0;  //放大器数量
  vector<int> verticle;   //存储业务经过的顶点(包括源点)    (end,...,begin)
  vector<int> path;       //存储业务经过的边的index  (边也是逆序)
  vector<int> magnifier;  //用来存储放置放大器的结点的index

  // // path_1 vertice_1是经过单条边距离不大于D过滤的
  // vector<int> verticle_1;//存储业务经过的顶点(包括源点)    (end,...,begin)
  // vector<int> path_1; //存储业务经过的边的index  (边也是逆序)

  int idx;  //业务index
  transaction(int b, int e, int i) : begin(b), end(e), idx(i)
  {
  }
};
//存储业务
vector<transaction> tranx;

class Add_edge
{
public:
  int num = 0;
  vector<pair<int, int>> begin_end;  //存储新边的起始结点
};

Add_edge add_edge;

//-------------------------------------------------------------------------

void init();  //输入初始化
void dijkstra(transaction& tran, int begin, int end);
// void Graph_init();//图初始化
// int find_dist(int u,int v);
vector<int> channel_is_spare(int edge_idx);  //判断某条边的通道是否占满（若已满 返回1，否则返回0）
vector<int> vectors_intersection(vector<int> v1, vector<int> v2);  //计算通道交集
void dfs(int begin, int end, int visit[]);
void cost(vector<vector<int>> all_path);
void best_path_verticle(transaction& tran, int src, int dst, vector<int> best_path);
void all_permutations(vector<vector<int>> vec, vector<vector<int>>& choose_path, uint32_t level = 0);

//------------------------------------------------判断输入的path是否畅通----------------------------------------------------
//------------------------------------------------判断输入的path是否畅通----------------------------------------------------
//------------------------------------------------判断输入的path是否畅通----------------------------------------------------
//------------------------------------------------判断输入的path是否畅通----------------------------------------------------

bool isPathGo(vector<int> path)
{
  vector<int> first = channel_is_spare(path[0]);
  vector<int> second;
  vector<int> intersection;
  for (int i = 1; i < path.size(); i++)
  {
    second = channel_is_spare(path[i]);

    intersection = vectors_intersection(first, second);

    first = intersection;
  }
  return intersection.size() != 0;
}

//------------------------------------------------按照每条边占用通道数最多进行降序排序----------------------------------------------------
//------------------------------------------------按照每条边占用通道数最多进行降序排序----------------------------------------------------
bool cmp(int i, int j)
{
  return edge_idx[i][0].ch_occ_num > edge_idx[j][0].ch_occ_num;
}

vector<int> sortedByChannelNumber(vector<int> path)
{
  // sort(path.begin(), path.end(), cmp);
  for (int i = 0; i < path.size() - 1; i++)
  {
    for (int j = i + 1; j < path.size(); j++)
    {
      if (edge_idx[path[i]][0].ch_occ_num < edge_idx[path[j]][0].ch_occ_num)
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
  //节点数量N、连边数量M、业务数量T、单边通道数量P、最大衰减距离D
  cout << "节点数量N、连边数量M、业务数量T、单边通道数量P、最大衰减距离D" << endl;
  cin >> N >> M >> T >> P >> D;

  for (int i = 0; i < M; i++)
  {
    int src, dst, weight;
    cin >> src >> dst >> weight;
    // 无向图
    // Graph[src][dst] = edge(src,dst,weight,i); //edges列表的index为图顶点
    // Graph[dst][src] = edge(dst,dst,weight,i);
    Graph[src].push_back(edge(src, dst, weight, i));  // edges列表的index为图顶点
    Graph[dst].push_back(edge(dst, src, weight, i));
    edge_idx[i].push_back(edge(src, dst, weight, i));
    edge_idx[i].push_back(edge(dst, src, weight, i));
    edge_channel[i] = channel();  //按顺序，存入每条边的通道状态
  }

  for (int i = 0; i < T; i++)
  {
    int begin, end;
    cin >> begin >> end;
    tranx.push_back(transaction(begin, end, i));
  }
}

vector<int> channel_is_spare(int edge_idx)
{
  vector<int> ret;
  for (int i = 0; i < P; i++)
  {
    if (edge_channel[edge_idx].stat[i] == 0)
    {
      ret.push_back(i);
    }
  }
  return ret;
}

vector<int> vectors_intersection(vector<int> v1, vector<int> v2)
{
  vector<int> v;
  sort(v1.begin(), v1.end());
  sort(v2.begin(), v2.end());
  set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(v));  //求交集
  return v;
}

//---------------------------------------------------------递归求路径-------------------------------------------------------------------
vector<int> indexes;  //全局变量用于记录从每个数组中所取元素的索引
void all_permutations(vector<vector<int>> vec, vector<vector<int>>& choose_path, uint32_t level)
{
  std::vector<int> vci;
  if (level < vec.size())
  {
    for (size_t i = 0; i < vec[level].size(); i++)
    {
      indexes[level] = i;
      all_permutations(vec, choose_path, level + 1);
    }
  }
  else
  {
    for (size_t i = 0; i < vec.size(); i++)
    {
      vci.push_back(vec[i][indexes[i]]);
    }
    choose_path.push_back(vci);
  }
}

//------------------------------------------------------------dfs----------------------------------------------------------------------
// vector<vector<int>> choose_path;
vector<vector<int>> all_path;
vector<vector<int>> all_verticle;
vector<int> result_verticle;
// vector<vector<int>> result_path;
// vector<vector<int>> result_path_all;//用来解决 两个顶点间有多条边

void dfs(int begin, int end, int visit[])
{
  result_verticle.push_back(begin);
  visit[begin] = 1;
  if (begin == end)
  {
    //存入顶点

    int tag = 0;  //删除重复出现的verticle
    for (int i = 0; i < all_verticle.size(); i++)
    {
      if (all_verticle[i] == result_verticle)
      {
        tag = 1;
      }
    }
    if (tag == 0)
    {
      all_verticle.push_back(result_verticle);
      int src, dst;

      //根据顶点存入路径
      vector<int> ret;
      vector<vector<int>> result_path;
      for (int i = 0; i < result_verticle.size() - 1; i++)
      {
        src = result_verticle[i];
        dst = result_verticle[i + 1];
        ret.clear();
        // vector<int> ret;
        // ret.clear();
        for (auto& e : Graph[src])
        {
          if (e.dst == dst)
          {
            ret.push_back(e.idx);
          }
        }
        result_path.push_back(ret);
      }
      // result_path.push_back(ret);
      indexes.resize(result_path.size());
      vector<vector<int>> choose_path;

      all_permutations(result_path, choose_path);

      for (int i = 0; i < choose_path.size(); i++)
      {
        all_path.push_back(choose_path[i]);
      }
    }
    // all_verticle.push_back(result_verticle);
    // int src,dst;

    // //根据顶点存入路径
    // vector<int> ret;
    // vector<vector<int>> result_path;
    // for(int i=0;i<result_verticle.size()-1;i++){

    //     src=result_verticle[i];
    //     dst=result_verticle[i+1];
    //     ret.clear();
    //     // vector<int> ret;
    //     // ret.clear();
    //     for(auto &e:Graph[src]){
    //         if(e.dst == dst){
    //             ret.push_back(e.idx);
    //         }
    //     }
    //     result_path.push_back(ret);
    // }
    // // result_path.push_back(ret);
    // indexes.resize(result_path.size());
    // vector<vector<int>> choose_path;

    // all_permutations(result_path,choose_path);

    // for(int i =0; i<choose_path.size() ;i++){
    //     all_path.push_back(choose_path[i]);
    // }

    // for(int k=0;k<choose_path.size() ;k++){
    //     for(int i=0;i<choose_path[k].size();i++){
    //         choose_path[k].pop_back();
    //     }
    // }

    // // result_path清零
    // for(int j=0;j<result_path.size();j++){
    //     for(int i=0;i< result_path[j].size();i++){
    //         result_path[j].pop_back();
    //     }

    // }

    // int flag=0;
    // for(auto &e:Graph[src]){

    //     if(e.dst == dst){
    //         if(flag==0){
    //             result_path.push_back(e.idx);

    //             flag++;
    //         }

    //         //怎么解决两个顶点间有多条边，弹出之前的边的idx，再重新push_back
    //         else{

    //             temp[flag] = result_path;
    //             temp[flag].pop_back();
    //             temp[flag].push_back(e.idx);
    //             flag++;
    //             // result_path_all.push_back(temp);
    //         }

    //     }
    // }

    // all_path.push_back(result_path);
    // for(int j = 0; j< M;j++){
    //     all_path.push_back(result_path_all[j]);
    // }
  }
  for (auto& e : Graph[begin])
  {
    if (visit[e.dst] == 0 && e.weight != INF)
    {
      dfs(e.dst, end, visit);
    }
  }
  //顶点回溯
  result_verticle.pop_back();
  visit[begin] = 0;
}
//------------------------------------------------------------花费计算------------------------------------------------------------------
vector<int> best_verticle;
vector<int> best_path;
int best_channel;
int add_edge_num;
int add_magnifier_num;
vector<pair<int, int>> add_edge_verticle;

void cost(vector<vector<int>> all_path)
{
  //定义当前临时长度
  int current_M = M;
  int min_money = INF;
  int money = 0;
  for (auto& path : all_path)
  {
    // for(int p=0; p< all_path.size(); p++){
    //     vector<int> path = all_path[p];
    vector<int> first = channel_is_spare(path[0]);
    vector<int> second;
    vector<int> intersection;
    for (int i = 1; i < path.size(); i++)
    {
      second = channel_is_spare(path[i]);

      intersection = vectors_intersection(first, second);

      first = intersection;
    }
    //有空余通道,计算花费
    if (intersection.size() > 0)
    {
      //计算放大器数量
      add_magnifier_num = 0;
      int weight = 0;
      int d_sum = 0;

      for (int i = 0; i < path.size(); i++)
      {
        weight = edge_idx[path[i]][0].weight;
        d_sum += weight;
        if (d_sum > D)
        {
          add_magnifier_num++;
          d_sum = 0;
        }
      }

      // 花费
      money = path.size() + add_magnifier_num * 10;
      if (money < min_money)
      {
        min_money = money;
        best_path = path;
      }
    }

    // 加边
    else
    {
      int edge_spare_channel_cnt[P] = { 0 };

      //遍历路径上的所有边
      //记录每个通道的空闲数
      for (int i = 0; i < path.size(); i++)
      {
        for (int j = 0; j < P; j++)
        {
          if (edge_channel[path[i]].stat[j] == 0)
          {
            edge_spare_channel_cnt[j]++;
          }
        }
      }

      //初始化
      int max_cnt = edge_spare_channel_cnt[0];
      int channel_idx = 0;
      int cnt;
      //选择共通空闲最多的通道
      for (int i = 1; i < P; i++)
      {
        cnt = edge_spare_channel_cnt[i];
        if (cnt > max_cnt)
        {
          max_cnt = cnt;
          channel_idx = i;
        }
      }

      //获得选好通道 channel_idx(加边最少)
      //计算加边花费
      vector<int> add_edge_idx;  //用来存储新加边

      // 按照通道数进行排序
      vector<int> new_path = sortedByChannelNumber(path);

      //优先从占用通道数最多的开始加边
      for (int i = 0; i < new_path.size(); i++)
      {
        //需要加边
        add_edge_idx.push_back(new_path[i]);
        //边数+1
        add_edge_num++;
        auto it = find(path.begin(), path.end(), new_path[i]) - path.begin();

        path[it] = current_M;
        current_M = current_M + 1;
        edge_idx[path[it]].push_back(edge(edge_idx[new_path[i]][0].src, edge_idx[new_path[i]][0].dst,
                                          edge_idx[new_path[i]][0].weight, path[it]));

        if (isPathGo(path))
        {
          break;
        }
      }

      //定义新的需要加的放大器数量
      add_magnifier_num = 0;
      int weight = 0;

      //设置当前的D
      int D_current = 0;

      for (int i = 0; i < path.size(); i++)
      {
        //如果该边是重新加入的，则需要重新考虑新加边的权值（取存在边的最小权重）
        // if (find(add_edge_idx.begin(), add_edge_idx.end(), path[i]) != add_edge_idx.end())
        // {
        //   int ver1 = edge_idx[path[i]][0].src;
        //   int ver2 = edge_idx[path[i]][0].dst;
        //   int min_weight = INF;
        //   for (auto& e : Graph[ver1])
        //   {
        //     if (e.dst == ver2)
        //     {
        //       min_weight = (e.weight < min_weight) ? e.weight : min_weight;
        //     }
        //   }
        //   weight = min_weight;
        // }
        // else
        // {
        //   weight = edge_idx[path[i]][0].weight;
        // }

        //当前D减去要走的路径权重
        D_current -= edge_idx[path[i]][0].weight;

        if (D_current < 0)
        {
          //加上一个放大器
          add_magnifier_num++;
          D_current = D - edge_idx[path[i]][0].weight;
        }
      }

      //计算加放大器的开销
      money = add_edge_num * 1000000 + path.size() + add_magnifier_num * 10;
      if (money < min_money)
      {
        min_money = money;
        best_path = path;
      }
    }
  }
  //无空余通道，加边，计算花费
}

void best_path_verticle(transaction& tran, int src, int dst, vector<int> best_path)
{
  tran.verticle.clear();
  for (int j = 0; j < best_path.size(); j++)
  {
    if (src == edge_idx[best_path[j]][0].src)
    {
      tran.verticle.push_back(edge_idx[best_path[j]][0].src);
      src = edge_idx[best_path[j]][0].dst;
    }
    else if (src == edge_idx[best_path[j]][0].dst)
    {
      tran.verticle.push_back(edge_idx[best_path[j]][0].dst);
      src = edge_idx[best_path[j]][0].src;
    }
  }
  tran.verticle.push_back(dst);
  //得到的tran.verticle是顺序存入
}

void dijkstra(transaction& tran, int begin, int end)
{  // u是源点，求的是源点到其他顶点最短路径

  // 初始化距离数组
  vector<int> dist(N + 1, INF);
  // 初始化路径数组
  vector<int> path(N + 1, -1);
  // 初始化优先队列
  priority_queue<pair<int, int>, vector<pair<int, int>>, greater<>> q;  // 队首元素距离dist一定是最小的
  // 初始化距离数组
  dist[begin] = 0;
  // 初始化路径数组
  path[begin] = begin;
  // 初始化优先队列
  q.push({ 0, begin });  //
  // Dijkstra算法
  while (!q.empty())
  {
    // 取出队首元素
    // const int v = -1;//顶点
    // const int d = -1;//距离
    auto ret = q.top();
    int d, v;  // d为距离，v为顶点
    d = ret.first;
    v = ret.second;
    q.pop();
    // 如果当前距离大于已知距离，则跳过
    if (d > dist[v])
    {
      // cout<<v<<endl;
      continue;
    }
    // 遍历当前节点的所有邻接点
    for (auto& e : Graph[v])
    {
      // 如果当前距离加上邻接点的权重大于已知距离，则跳过
      if (dist[v] + e.weight > dist[e.dst])
        continue;
      // vector<int> e_channel = channel_is_spare(e.idx);
      // if (dist[v] + e.weight <= dist[e.dst] && e_channel.size()>0){
      if (dist[v] + e.weight <= dist[e.dst])
      {
        // 更新距离数组
        dist[e.dst] = dist[v] + e.weight;
        // 更新路径数组
        path[e.dst] = v;
        // 将邻接点加入优先队列
        q.push({ dist[e.dst], e.dst });
      }
    }
  }

  //--------------------------将路径经过的边的index和顶点加入Class transaction---------------------------------------
  int dst = end;
  int src = path[dst];
  int edge_weight = 0;
  //遍历边，添加路径经过的边的index和顶点

  while (dst != begin)
  {
    edge_weight = dist[dst] - dist[src];
    for (auto& e : Graph[src])
    {
      if (e.dst == dst && e.weight == edge_weight)
      {
        tran.path.push_back(e.idx);
      }
    }
    tran.verticle.push_back(dst);
    //节点前移
    dst = src;
    src = path[src];
  }

  tran.verticle.push_back(begin);

  //-----------------------------------------------选择通道-----------------------------------------------
  //得到path以后 开始处理通道
  //遍历所用边的空余通道，最后取交集
  vector<int> first = channel_is_spare(tran.path[tran.path.size() - 1]);
  vector<int> second;
  vector<int> intersection;

  for (int k = tran.path.size() - 2; k >= 0; k--)
  {
    second = channel_is_spare(tran.path[k]);

    intersection = vectors_intersection(first, second);

    first = intersection;
  }

  //有空余通道
  if (intersection.size() > 0)
  {
    int channel_idx = intersection[0];
    tran.channel_idx = channel_idx;

    for (auto& k : tran.path)
    {
      edge_channel[k].stat[channel_idx] = 1;
      edge_idx[k][0].ch_occ_num += 1;
    }

    //------------------------------------------------------放大器--------------------------------------------
    //计算加放大器的数量
    //考虑路径 src----(weight1)-----end
    int weight = 0;
    int D_current = D;
    for (int k = tran.verticle.size() - 1; k >= 0; k--)
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
  //无同一类型的空余通道
  else
  {
    //------------------------加边-------------------------------------
    //加边的情况：
    //不能在一条路径上的所有边中找到同一个类型的通道

    // 利用DFS遍历begin到end的所有路径并存储，查找是该路径是否有同一类型的空余通道
    // 若有，找一条边数最短的路径进行切换，不需要加边:
    //      tran.path.swap(new_path);
    //      tran.verticle.swap(new_verticle);
    int visit[N] = { 0 };
    best_path.clear();
    all_path.clear();
    all_verticle.clear();
    result_verticle.clear();
    // result_path.clear();
    // result_path_all.clear();//用来解决 两个顶点间有多条边

    // path verticle的存储顺序与之前不同 顺序src——>......——>end

    //获得all_path,all_verticle
    dfs(begin, end, visit);

    //获得best_path
    cost(all_path);

    vector<int> first = channel_is_spare(best_path[0]);
    vector<int> second;
    vector<int> intersection;
    for (int k = 1; k >= 0; k--)
    {
      second = channel_is_spare(best_path[k]);

      intersection = vectors_intersection(first, second);

      first = intersection;
    }
    // 有空余通道,只需修改路径 不需要加边
    // 判断best_path对应的index是否在原始的path里面
    if (intersection.size() > 0 && is_add_edge(best_path) == 0)
    {
      int channel_idx = intersection[0];
      tran.channel_idx = channel_idx;
      for (auto& k : best_path)
      {
        edge_channel[k].stat[channel_idx] = 1;
        edge_idx[k][0].ch_occ_num += 1;
      }

      //将best_path路径上所经过的顶点push进入tran.verticle
      best_path_verticle(tran, begin, end, best_path);

      //计算放大器数量
      add_magnifier_num = 0;
      int weight = 0;
      int D_current = D;

      for (int i = 0; i < best_path.size(); i++)
      {
        weight = edge_idx[best_path[i]][0].weight;
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
      reverse(best_path.begin(), best_path.end());
      tran.path = best_path;
    }
    else
    {
      int edge_spare_channel_cnt[P] = { 0 };
      //遍历路径上的所有边
      //记录每个通道的空闲数
      for (int i = 0; i < best_path.size(); i++)
      {
        for (int j = 0; j < P; j++)
        {
          if (edge_channel[best_path[i]].stat[j] == 0)
          {
            edge_spare_channel_cnt[j]++;
          }
        }
      }

      //初始化
      int max_cnt = edge_spare_channel_cnt[0];
      int channel_idx;
      int cnt;
      //选择共通空闲最多的通道
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
      for (int i = 0; i < best_path.size(); i++)
      {
        if (edge_channel[best_path[i]].stat[channel_idx] == 1)
        {
          int ver1 = edge_idx[best_path[i]][0].src;
          int ver2 = edge_idx[best_path[i]][0].dst;
          int min_weight = INF;
          for (auto& e : Graph[ver1])
          {
            if (e.dst == ver2)
            {
              min_weight = (e.weight < min_weight) ? e.weight : min_weight;
            }
          }
          //加边 同时图graph和edge_idx也要改变,还要修改best_path中边的index
          add_edge.num++;
          add_edge.begin_end.push_back(make_pair(ver1, ver2));
          Graph[ver1].push_back(edge(ver1, ver2, min_weight, M - 1 + add_edge.num));
          Graph[ver2].push_back(edge(ver2, ver1, min_weight, M - 1 + add_edge.num));
          edge_idx[M - 1 + add_edge.num].push_back(edge(ver1, ver2, min_weight, M - 1 + add_edge.num));
          edge_idx[M - 1 + add_edge.num].push_back(edge(ver2, ver1, min_weight, M - 1 + add_edge.num));
          best_path[i] = M - 1 + add_edge.num;
        }
      }

      best_path_verticle(tran, begin, end, best_path);
      //计算放大器数量
      add_magnifier_num = 0;
      int weight = 0;
      int D_current = D;
      for (int i = 0; i < best_path.size(); i++)
      {
        weight = edge_idx[best_path[i]][0].weight;
        D_current -= weight;
        if (D_current < 0)
        {
          // add_magnifier_num++;
          tran.magnifier_num++;
          tran.magnifier.push_back(tran.verticle[i]);
          D_current = D - weight;
        }
      }

      //清空之前存储的path
      tran.path.clear();
      //将best_path 逆置
      reverse(best_path.begin(), best_path.end());
      tran.path = best_path;
    }
  }
}

int main()
{
  init();

  // 对每个顶点的边按距离进行
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
      //输出加边的源点和终点
      cout << add_edge.begin_end[i].first << " " << add_edge.begin_end[i].second << endl;
    }
  }

  //输出交易的通道编号+边数+放大器数+经过边的index+依次经过的放大器所在节点的编号
  for (int i = 0; i < tranx.size(); i++)
  {
    cout << tranx[i].channel_idx << " " << tranx[i].path.size() << " " << tranx[i].magnifier_num << " ";
    for (int j = tranx[i].path.size() - 1; j >= 0; j--)
    {
      cout << tranx[i].path[j] << " ";
    }
    for (auto& k : tranx[i].magnifier)
    {
      cout << k << " ";
    }
    cout << endl;
  }
  return 0;
}
