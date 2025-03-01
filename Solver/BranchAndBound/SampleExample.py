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
import numpy as np


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
    eps = 1e-3
    incumbent_node = None
    Gap = np.inf

    
    """
    Branch and Bound starts
    """
    Queue = []
    node = Node()
    node.local_LB = 0
    node.local_UB = global_UB
    node.model = model.copy()
    node.model.setParam("OutputFlag", 0)
    node.cnt = 0
    Queue.append(node)

    cnt = 0
    global_UB_change = []
    global_LB_change = []

    while (len(Queue)>0 and global_UB - global_LB > eps):
        current_node = Queue.pop()
        cnt += 1
        current_node.model.optimize()
        Solution_status = current_node.model.Status

        is_integer = True
        is_pruned = False
        if Solution_status == 2:  # 有最优解
            for var in current_node.model.getVars():
                current_node.x_sol[var.VarName] = var.x
                print(var.VarName, " = ", var.x)
                current_node.x_int_sol[var.VarName] = (int)(var.x)
                if abs((int)(var.x) - var.x) >= eps:
                    is_integer = False
                    current_node.branch_var_list.append(var.VarName)
            if is_integer: # 该节点的值为整数，更新节点上下界
                is_pruned = True
                current_node.is_integer = True
                current_node.local_LB = current_node.model.ObjVal
                current_node.local_UB = current_node.model.ObjVal
                # 更新全局下界
                if current_node.local_LB > global_LB:
                    global_LB = current_node.local_LB
                    incumbent_node = Node.deepcopy_node(current_node)
            if not is_integer: # 该节点的解为非整数解
                current_node.is_integer = False
                current_node.local_UB = current_node.model.ObjVal
                current_node.local_LB = 0
                for var_name in current_node.x_int_sol.keys():
                    var = current_node.model.getVarByName(var_name)
                    # NOTE 更新节点下界 - 目标函数系数 * 圆整后的值
                    current_node.local_LB += current_node.x_int_sol[var_name] * var.Obj
                if (current_node.local_LB > global_LB) or (current_node.local_LB == global_LB and current_node.is_integer == True):
                    global_LB = current_node.local_LB
                    incumbent_node = Node.deepcopy_node(current_node)
                    incumbent_node.local_LB = current_node.local_LB
                    incumbent_node.local_UB = current_node.local_UB
                    pass
            
            if not is_integer and current_node.local_UB < global_LB:
                is_pruned = True
            
            Gap = round(100*(global_UB - global_LB) / global_LB, 2)
            print('\n ---------------------------------------- \n', cnt, '\t Gap = ', Gap, ' %')
        elif Solution_status != 2:
            is_integer = False
            is_pruned = True
        
        if not is_pruned:
            branch_var_name = current_node.branch_var_list[0]
            left_var_bound = (int)(current_node.x_sol[branch_var_name])
            right_var_bound = (int)(current_node.x_sol[branch_var_name]) + 1

            left_node = Node.deepcopy_node(current_node)
            right_node = Node.deepcopy_node(current_node)

            temp_var = left_node.model.getVarByName(branch_var_name)
            left_node.model.addConstr(temp_var <= left_var_bound, name=f"branch_left_{cnt}")
            left_node.model.setParam("OutputFlag", 0)
            left_node.model.update()
            cnt += 1
            left_node.cnt = cnt

            temp_var = right_node.model.getVarByName(branch_var_name)
            right_node.model.addConstr(temp_var >= right_var_bound, name=f"branch_right_{cnt}")
            right_node.model.setParam("OutputFlag", 0)
            right_node.model.update()
            cnt += 1
            right_node.cnt = cnt

            Queue.append(left_node)
            Queue.append(right_node)

            temp_global_UB = 0
            for node in Queue:
                node.model.optimize()
                if node.model.status == 2:
                    if node.model.ObjVal >= temp_global_UB:
                        temp_global_UB = node.model.ObjVal
            
            global_UB = temp_global_UB
            global_UB_change.append(global_UB)
            global_LB_change.append(global_LB)
    
    Gap = round(100*(global_UB - global_LB) / global_LB, 2)
    print('\n\n\n\n')
    print(' ------------------------------------------- ')
    print (' Branch, and Bound terminates ')
    print(' Optimal solution found  ')
    print ('-------------------------------------------')
    print('\nFinal Gap = ', Gap, ' % ' )
    print('Optimal Solution: ', incumbent_node.x_int_sol)
    print("optimal Obj: ", global_LB)
    return incumbent_node, Gap, global_UB_change, global_LB_change

incumbent_node, Gap, Global_UB_change, Global_LB_change = Branch_And_Bound(RLP)
