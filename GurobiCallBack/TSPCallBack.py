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


def MyCallBack(model, where):
    if where == GRB.Callback.MIPSOL:
        # x是一个字典，下标 - 值
        x = model.cbGetSolution(model._vars)
        tour = subtour(x)
        if len(tour) < nodeNum:
            print('----add sub tour elimination constraints')
            model.cbLazy(quicksum(model._vars[i, j] for i, j in itertools.permutations(tour, 2)) <= len(tour)-1)



def subtour(x: dict):
    unvisited_node = [i for i in range(nodeNum)]
    cycle = [i for i in range(nodeNum)]

    edge = [(i, j) for i in range(nodeNum) for j in range(nodeNum) if x[i, j] > 0.5]
    while unvisited_node:
        thiscycle = []
        neighbors = copy.deepcopy(unvisited_node)
        while neighbors:
            current_node = neighbors[0]
            thiscycle.append(current_node)
            unvisited_node.remove(current_node)
            neighbors1 = [i for i, j in edge if j == current_node and i in unvisited_node]
            neighbors2 = [j for i, j in edge if i == current_node and j in unvisited_node]
            neighbors = [i for i in neighbors2 + neighbors1]
        if len(thiscycle) < len(cycle):
            cycle = [i for i in thiscycle]
    return cycle


file = r'E:\PYCHARM\Solomn Instance\solomon-100\c101.txt'
distance = ReadData(file, nodeNum)
# 定义模型
m = Model('TSP')

# 创建变量
x = m.addVars(nodeNum, nodeNum, vtype=GRB.BINARY, name='x')

# 创建目标函数
m.setObjective(quicksum(x[i, j] * distance[i][j] for i in range(nodeNum) for j in range(nodeNum)), GRB.MINIMIZE)

# 添加约束
m.addConstrs((quicksum(x[i, j] for i in range(nodeNum) if i != j) == 1 for j in range(nodeNum)), name='c1')
m.addConstrs((quicksum(x[i, j] for j in range(nodeNum) if j != i) == 1 for i in range(nodeNum)), name='c2')
m.addConstr(quicksum(x[0, j] for j in range(1, nodeNum)) == 1, name='c3')
m.addConstr(quicksum(x[j, 0] for j in range(1, nodeNum)) == 1, name='c4')
m.addConstrs((x[i, i] == 0 for i in range(nodeNum)), name='c5')
# 获取外部变量
m._vars = x
m.Params.lazyConstraints = 1
m.optimize(MyCallBack)

print('Objective Vlaue:{}'.format(m.ObjVal))
X = m.getAttr('x', x)

for i in range(nodeNum):
    for j in range(nodeNum):
        if X[i, j] > 0.5:
            print('{},{}: {}'.format(i, j, X[i, j]))
