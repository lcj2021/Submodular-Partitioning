# Submodular-Partitioning

Question:

1. epsilon的选取?

2. 原algorithm2中 `ret = pi_c;` 何时将计算结果赋给返回值? 

3. db较小时会出现miu_i加入分区时delta全部相等的情况, 这时默认会选取第1个分区加入(顺序遍历)
> ~~若delta全部相等, 随机选取分区加入行不行~~
> 其实是中间结果(c, F函数的值)未用double