from gurobipy import *

# 各钢材的需求长度
TypesDemand = [3, 7, 9, 16]
# 各钢材的需求数量
QuantityDemand = [25, 30, 14, 8]
# 现有钢材的长度
LengthUsable = 20


try:
    # 创建模型
    # 主问题模型
    MainProb = Model()
    # 子问题模型
    SubPro = Model()

    # 构建主问题模型
    # 变量
    Zp = MainProb.addVars(len(TypesDemand), vtype=GRB.CONTINUOUS, obj=1.0, name='z')
    # 添加约束
    MainConstrs = MainProb.addConstrs((quicksum(Zp[p]*(LengthUsable // TypesDemand[i]) for p in range(len(TypesDemand)) if p == i) >= QuantityDemand[i] for i in range(len(QuantityDemand))), name='m1')
    # 求解模型
    MainProb.optimize()

    # 获得对偶值
    Dualsolution = MainProb.getAttr('Pi', MainProb.getConstrs())

    # 构造子问题模型
    # 添加变量
    C = SubPro.addVars(len(QuantityDemand), obj=Dualsolution, vtype=GRB.INTEGER, name='c')
    # 添加约束
    SubPro.addConstr(quicksum(C[i]*TypesDemand[i] for i in range(len(QuantityDemand))) <= LengthUsable, name='s1')
    # 求解子问题
    SubPro.setAttr(GRB.Attr.ModelSense, -1)
    SubPro.optimize()

    while SubPro.objVal > 1:
        C_value = SubPro.getAttr('x', SubPro.getVars())
        column = Column(C_value, MainProb.getConstrs())
        # 主问题添加新变量
        MainProb.addVar(obj=1.0, vtype=GRB.CONTINUOUS, column=column)
        # 求解主问题
        MainProb.optimize()
        Dualsolution = MainProb.getAttr('Pi', MainProb.getConstrs())
        for i in range(len(TypesDemand)):
            C[i].obj = Dualsolution[i]
        SubPro.optimize()
    for v in MainProb.getVars():
        v.setAttr('Vtype', GRB.INTEGER)
    MainProb.optimize()
    for v in MainProb.getVars():
        if v.X > 0:
            print('{} is {}'.format(v.VarName, v.X))
    print(MainProb.getVars())
    print(MainProb.getAttr('x', SubPro.getVars()))
except GurobiError as e:
    print(print('Error code' + str(e.errno) + ':' + str(e)))
except AttributeError:
    print('Encountered an attribute error')




