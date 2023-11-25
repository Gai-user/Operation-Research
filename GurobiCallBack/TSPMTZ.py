import copy
from gurobipy import *
import re
import math
import itertools

nodeNum = 11


def ReadData(file, nodeNum):
    posX, posY = [], []
    distance = [[0 for i in range(nodeNum)] for j in range(nodeNum)]
    with open(file) as file:
        readlines = file.readlines()[9:]
        for line in readlines[:nodeNum]:
            location = re.split(r' +', line)
            posX.append(int(location[2]))
            posY.append(int(location[3]))
        for i in range(nodeNum):
            for j in range(nodeNum):
                distance[i][j] = math.sqrt((posX[i] - posX[j]) ** 2 + (posY[i] - posY[j]) ** 2)
    return distance


file = r'E:\PYCHARM\Solomn Instance\solomon-100\c101.txt'
distance = ReadData(file, nodeNum)
# 定义模型
m = Model('TSP')

# 创建变量
x = m.addVars(nodeNum, nodeNum, vtype=GRB.BINARY, name='x')
u = m.addVars(nodeNum, vtype=GRB.INTEGER, name='u')
# 创建目标函数
m.setObjective(quicksum(x[i, j] * distance[i][j] for i in range(nodeNum) for j in range(nodeNum)), GRB.MINIMIZE)

# 添加约束
m.addConstrs((quicksum(x[i, j] for i in range(nodeNum) if i != j) == 1 for j in range(nodeNum)), name='c1')
m.addConstrs((quicksum(x[i, j] for j in range(nodeNum) if j != i) == 1 for i in range(nodeNum)), name='c2')
m.addConstr(quicksum(x[0, j] for j in range(1, nodeNum)) == 1, name='c3')
m.addConstr(quicksum(x[j, 0] for j in range(1, nodeNum)) == 1, name='c4')
m.addConstrs((x[i, i] == 0 for i in range(nodeNum)), name='c5')
m.addConstrs((u[j] - u[i] >= (x[i, j] - 1)*nodeNum + 1 for i in range(1, nodeNum) for j in range(1, nodeNum)), name='c6')

m.optimize()

print('Objective Vlaue:{}'.format(m.ObjVal))
X = m.getAttr('x', x)
for i in range(nodeNum):
    for j in range(nodeNum):
        if X[i, j] > 0.5:
            print('{},{}: {}'.format(i, j, X[i, j]))