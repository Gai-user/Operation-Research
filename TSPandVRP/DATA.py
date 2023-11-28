"""
读取数据、处理模型结果
"""
import pandas as pd
import numpy as np

# 显示所有列
pd.set_option('display.max_columns', None)

# 显示所有行
pd.set_option('display.max_rows', None)


def readData(file_name, nodeNum):
    """
    :param file_name: 读取文件的绝对路径
    :param nodeNum: 点个数 - 客户数量
    :return:
    """
    # 使用numpy读取txt文件，将数据转换为ndarray
    ndarray = np.loadtxt(r"E:\OR\Solomn Instance\solomon-100\c101.txt", skiprows=8)
    # 将ndarray转换为DataFrame
    customer_info = pd.DataFrame(ndarray, columns=['CUST NO.', 'XCOORD.', 'YCOORD.', 'DEMAND', 'READY TIME', 'DUE DATE',
                                                   'SERVICE'])
    # 各个顾客之间的距离
    customer_distance = np.zeros((nodeNum + 1, nodeNum + 1))
    customer_demand = np.zeros(nodeNum + 1)
    customer_time_windown = np.zeros((nodeNum+1, 2))
    customer_service_time = np.zeros(nodeNum + 1)
    for i in range(nodeNum):
        customer_demand[i] = customer_info.iloc[i, 3]
        customer_time_windown[i][0], customer_time_windown[i][1] = customer_info.iloc[i, 4], customer_info.iloc[i, 5]
        customer_service_time[i] = customer_info.iloc[i, 6]
        for j in range(nodeNum):
            customer_distance[i][j] = np.sqrt((customer_info.iloc[i, 1] - customer_info.iloc[j, 1]) ** 2 + (
                        customer_info.iloc[i, 2] - customer_info.iloc[j, 2]) ** 2)
    for i in range(nodeNum):
        customer_distance[nodeNum][i] = customer_distance[0][i]
        customer_distance[i][nodeNum] = customer_distance[i][0]
    customer_demand[nodeNum] = customer_demand[0]
    customer_time_windown[nodeNum][0], customer_time_windown[nodeNum][1] = customer_time_windown[0][0], customer_time_windown[0][1]
    customer_service_time[nodeNum] = customer_service_time[0]
    return customer_distance, customer_demand, customer_time_windown, customer_service_time


