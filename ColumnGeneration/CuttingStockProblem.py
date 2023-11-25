from gurobipy import *

# 各钢材的需求长度
TypesDemand = [3, 7, 9, 16]
# 各钢材的需求数量
QuantityDemand = [25, 30, 14, 8]
# 现有钢材的长度
LengthUsable = 20

try:
    # 限制主问题（线性松弛模型）
    MainProbRelax = Model()
    # 子问题
    SubProb = Model()

    # 构造主问题模型，初始切割方案为：每根钢材只切割一种长度
    # 添加变量
    Zp = MainProbRelax.addVars(len(TypesDemand), obj=1.0, vtype=GRB.CONTINUOUS, name='z')
    # 添加约束
    ColumnIndex = MainProbRelax.addConstrs((quicksum(
        Zp[p] * (LengthUsable // TypesDemand[i]) for p in range(len(TypesDemand)) if p == i) >= QuantityDemand[i] for i in
                                            range(len(TypesDemand))), name='c1')
    # 求解模型
    MainProbRelax.optimize()
    # 构造子问题模型
    # 获得主问题约束的对偶值, 这里没有弄清楚，怎么获得对偶值
    # TODO 学习如何获取对偶值
    Dualsolution = MainProbRelax.getAttr(GRB.Attr.Pi, MainProbRelax.getConstrs())
    # NOTE 获取对偶值
    a = MainProbRelax.getAttr('Pi', MainProbRelax.getConstrs())
    print(a)
    Ci = SubProb.addVars(len(TypesDemand), obj=Dualsolution, vtype=GRB.INTEGER, name='c')

    # 添加约束
    SubProb.addConstr(quicksum(Ci[i]*TypesDemand[i] for i in range(len(TypesDemand))) <= LengthUsable, name='s1')
    # 求解子问题
    SubProb.setAttr(GRB.Attr.ModelSense, -1)
    SubProb.optimize()
    print(SubProb.objVal)
    # 判断 Reduced Cost是否小于0
    while SubProb.objVal > 1:
        # 获取变量取值
        columnCoeff = SubProb.getAttr('x', SubProb.getVars())
        # TODO 按列建模？
        column = Column(columnCoeff, MainProbRelax.getConstrs())
        # 添加变量
        MainProbRelax.addVar(obj=1.0, vtype=GRB.CONTINUOUS, name='CG', column=column)
        MainProbRelax.optimize()
        # 修改子问题目标函数系数
        for i in range(len(TypesDemand)):
            Ci[i].obj = ColumnIndex[i].pi
        SubProb.optimize()
    # 将 CG 后的模型转为整数并求解
    for v in MainProbRelax.getVars():
        v.setAttr('Vtype', GRB.INTEGER)
    MainProbRelax.optimize()
    for v in MainProbRelax.getVars():
        if v.X != 0.0:
            print('{} is {}'.format(v.VarName, v.X))
except GurobiError as e:
    print('Error code' + str(e.errno) + ':' + str(e))
except AttributeError:
    print('Encountered an attribute error')
