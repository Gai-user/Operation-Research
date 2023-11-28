from DATA import readData
from gurobipy import *

file_path = r'E:\OR\Solomn Instance\solomon-100\c101.txt'
nodeNum = 10
vehicle_number = 4
capacity = 200

node_distance, node_demand, node_time_windown, service_time = readData(file_path, nodeNum)

# 创建模型，进行求解

# creat Model
m = Model('VRPTW')
# creat variables
x = m.addVars(nodeNum + 1, nodeNum + 1, vehicle_number, vtype=GRB.BINARY, name='x')
w = m.addVars(nodeNum + 1, vehicle_number, vtype=GRB.CONTINUOUS, lb=0, ub=1236, name='w')

# creat objective function
m.setObjective(quicksum(
    node_distance[i][j] * x[i, j, k] for i in range(nodeNum + 1) for j in range(nodeNum + 1) for k in
    range(vehicle_number)))
# add constraints
m.addConstrs(
    (quicksum(x[i, j, k] for k in range(vehicle_number) for j in range(nodeNum + 1)) == 1 for i in range(1, nodeNum)),
    name='c1')
m.addConstrs((quicksum(x[0, j, k] for j in range(1, nodeNum + 1)) == 1 for k in range(vehicle_number)), name='c2')
m.addConstrs((quicksum(x[i, h, k] - x[h, j, k] for i in range(nodeNum + 1) for j in range(nodeNum + 1)) == 0 for h in
              range(1, nodeNum) for k in range(vehicle_number)), name='c3')
m.addConstrs((quicksum(x[i, nodeNum, k] for i in range(nodeNum)) == 1 for k in range(vehicle_number)), name='c4')
m.addConstrs(
    (quicksum(node_demand[i] * x[i, j, k] for i in range(1, nodeNum) for j in range(nodeNum + 1)) <= capacity for k in
     range(vehicle_number)), name='c5')
m.addConstrs((w[i, k] + node_distance[i][j] - 1236 * (1 - x[i, j, k]) <= w[j, k] for i in range(nodeNum + 1) for j in
              range(nodeNum + 1) for k in range(vehicle_number)), name='c6')
m.addConstrs((w[i, k] <= node_time_windown[i][1] for i in range(nodeNum + 1) for k in range(vehicle_number)), name='c7')
m.addConstrs((w[i, k] >= node_time_windown[i][0] for i in range(nodeNum + 1) for k in range(vehicle_number)), name='c8')
m.addConstrs((x[i, i, k] == 0 for i in range(nodeNum+1) for j in range(nodeNum+1) for k in range(vehicle_number)), name='c9')
m.update()
m.optimize()
