"""
列生成求解VRPTW
"""
from DATA import readData
from gurobipy import *


nodeNum = 6
vehicle_num = 5
customer_distance, customer_demand, customer_time_windown, customer_service_time = readData(r'E:\OR\Solomn '
                                                                                            r'Instance\solomon-100\c101.txt', nodeNum)
# 用于存储最优路径
path_set = {}
'''构造主问题'''
# 创建Model类的对象
RMP = Model('RMP')

# 创建变量
y = RMP.addVars(nodeNum-1, lb=0, ub=1, vtype=GRB.CONTINUOUS, name='y')
# 创建目标函数
RMP.setObjective(quicksum((customer_distance[0][i+1] + customer_distance[i+1][0])*y[i] for i in range(nodeNum-1)), GRB.MINIMIZE)
RMP.addConstr(y[0] == 1, name='c1')
RMP.addConstr(y[1] == 1, name='c2')
RMP.addConstr(y[2] == 1, name='c3')
RMP.addConstr(y[3] == 1, name='c4')
RMP.addConstr(y[4] == 1, name='c5')

RMP.update()
RMP.optimize()

RMP.write('RMP_initial.lp')
print('RMP Obj = {}'.format(RMP.ObjVal))
for i in range(1, nodeNum):
    varname = 'y[{}]'.format(i-1)
    path_set[varname] = [[0, i, nodeNum], customer_distance[0][i] + customer_distance[i][0]]

for var in RMP.getVars():
    if var.x > 0.5:
        print(var.VarName, ' = ', var.x, '\t path : ', path_set[var.VarName])
# 获取 各个约束的 对偶变量
RMP_pi = RMP.getAttr('Pi', RMP.getConstrs())
print(RMP_pi)

''' 构造子问题， 找到新的列，加入主函数中 '''
# 创建Model类
SP = Model('subproblem')

# 创建变量
s = {}
for i in range(nodeNum):
    s[i] = SP.addVar(lb=customer_time_windown[i][0], ub=customer_time_windown[i][1], vtype=GRB.CONTINUOUS, name='s[{}]'.format(i))
x = SP.addVars(nodeNum, nodeNum, lb=0, ub=1, vtype=GRB.BINARY, name='x')

# 创建目标函数
SP.setObjective(quicksum((customer_distance[i][j] - RMP_pi[i])*x[i, j] for i in range(nodeNum) for j in range(nodeNum)), GRB.MINIMIZE)

print(customer_time_windown)







