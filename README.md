# Huawei-Embedded-Software-Competition-2023
2023年华为嵌入式软件大赛本人代码开源，供大家学习和指正。主要方法是：

1.用一个vector<vector<int>> Graph 保存图，内存占有比较大但是要保存weight和idx和占用通道数量等信息，没有使用邻近矩阵存储，之后为了学习尝试有空会再想想 

2.用dijkstra找最短路径，然后用tranx存下当前业务通过的路径，和放大器数量和放大器的节点 

3.判断不能走通时就将通道占用最多的先加边，Graph对应的边替换边idx，通道数占用数刷新，而不是往Graph增加边减少内存占用 

```
git clone https://github.com/JJho1314/Huawei-Embedded-Software-Competition-2023.git
cd ..
mkdir build
cd build
cmake ..
make
```