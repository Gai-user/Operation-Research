# 一个简单例子，用来理解分支定界算法
"""
max 100 * x1 + 150 * x2
s.t:
    2*x1 + x2 <= 10
    3*x1 + 6x2 <= 40
    x1, x2 >= 0, and integer
"""

"""
分支定界算法伪代码：

1. 初始化模型，创建根节点 S，其模型为对应IP模型的LPr
2. 设置全局上下界：UB - 极大值， LB - 极小值
3. 设置节点集合 Q <-- {S}，当前最优解x* <--Null
4. while Q 非空 且 UB - LB > tol do
5.     选择当前节点 S_i <-- Q.pop()，对应模型为 P_i
6.     stauts <-- 求解模型 P_i的线性松弛模型
7.     更新对偶界限 z_i <-- 线性松弛模型的目标函数
8.     x_i(LPr) <-- 线性松弛模型的解
9.     if status 不是最优 then
10.        根据不可行性剪枝
11.    else if z_i <= LB then
12.        根据界限剪枝
13.    else if x_i(LPr) 是整数解 then
14.        更新全局下界 LB <-- z_i
15.        更新当前最优解 x* <-- x_i(LPr)
16.        根据最优性剪枝
17.    else if status 是最优且 x_i(LPr) 为非整数解 then
18.        将当前解 x_i(LPr) 圆整为整数解 int( x_i(LPr) ) --- 
19：       if 圆整后整数解的目标函数 int( x_i(LPr) ) > LB then
20：           更新LB <-- int( x_i(LPr) ) 对应的目标函数
21：           更新当前最优解 x* <-- int( x_i(LPr) )
22：       end if
23:        选择 x_i(LPr) 中第一个小数变量为分支变量
24:        创建两个子节点 S_i1 和 S_i2 其对应的模型分别为 P_i1 和 P_i2
25：       更新节点的集合 Q <-- Q U {S_i1, S_i2}
26：       tempUB 一剩余叶子节点集合Q 中节点的线性松弛模型的目标函数的最大值
27：       更新 UB - tempUB
28：    end if
29:end While
30:return 最优解 x*
"""

from gurobipy import *
# import matplotlib.pyplot as plt
# 添加项目的根目录到 sys.path
import sys
project_root = "E:/OR/Solver"  # 根据实际情况修改路径
sys.path.append(project_root)
from utils import Node


RLP = Model('relaxed MIP')
x = RLP.addVars(2, lb=0, ub= GRB.INFINITY, vtype=GRB.CONTINUOUS, name = ['x1', 'x2'])
RLP.update()
RLP.setObjective(100*x[0] + 150*x[1], GRB.MAXIMIZE)

RLP.addConstr(2*x[0]+x[1]<=10, name='c1')
RLP.addConstr(3*x[0]+6*x[1]<=10, name='c2')
RLP.optimize()

def Branch_And_Bound(model):
    model.optimize()
    global_UB = model.ObjVal
    global_LB = 0