# 分支定价算法流程图

"""
                                初始化
                    初始问题S及其模型P并将其加入NodeList
                                z_lb = -INF （一个非常小的数）
                                incumbent x* = Null
                                |
                                |
        |---------------判断NodeList是否为空 ------是------ 结束，当前最优解x*
        |                       |
        |                       | 否
        |               选择问题S_i，其模型为P_i
        |                       |
        |                       | 否
        |               求解P_i的线性松弛LPr
        |               得到LPr的解，z_i
        |                       |
        |                       | 否
        ----------------如果P_i不可行，根据不可行剪枝
        |                       |
        |                       | 否
        ----------------如果z_i <= z_lb, 根据界限剪枝
        |                       |
        |                       | 否
        |               若是z_i为整数解，且好于当前最好解，
        ----------------更新z_lb = z_i， incumbent x* = z_i
        |               根据最优性剪枝
        |                       |
        |                       | 否
        |----------------返回两个子问题S_i1,S_i2
        
"""