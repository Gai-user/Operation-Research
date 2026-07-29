/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   ProbDataTSP.h
 * @brief  TSP 问题的 C++ 问题数据类
 * @author Timo Berthold
 *
 * @details
 * 在 SCIP 中，问题数据（Problem Data）是存储问题实例具体数据的对象。
 * 对于 TSP，问题数据就是包含节点坐标和边信息的完全图（Graph）。
 *
 * ## SCIP 问题数据的生命周期
 *
 * SCIP 将求解过程分为几个阶段，问题数据在不同的阶段需要不同的处理：
 *
 *   1. **原始阶段 (Original)**：从文件读取后创建，包含所有原始数据
 *   2. **变换阶段 (Transformed)**：经过 presolve 后，数据可能被简化
 *   3. **拷贝**：在并行求解等场景中，需要拷贝问题数据到子 SCIP 实例
 *
 * 对应四个关键回调：
 *   - **scip_trans**：将原始数据变换为 presolved 数据
 *   - **scip_copy**：将数据深拷贝到另一个 SCIP 实例
 *   - **scip_delorig**：释放原始数据
 *   - **scip_deltrans**：释放变换后的数据
 *
 * ## 图的所有权
 *
 * ProbDataTSP 通过 `std::unique_ptr<Graph>` 独有图的所有权。
 * 约束处理器（ConshdlrSubtour）通过 SCIP_ConsData 持有 Graph* 非拥有引用。
 * 这确保了单一所有权，避免内存泄漏和双重释放。
 */

#pragma once

#include <memory>
#include "objscip/objscip.h"
#include "GomoryHuTree.h"

namespace tsp {

/**
 * @brief TSP 问题数据类
 *
 * 封装 TSP 实例的完全图数据，作为 SCIP 问题数据对象。
 */
class ProbDataTSP : public scip::ObjProbData
{
public:
   /**
    * @brief 构造问题数据
    * @param g 图的 unique_ptr，所有权转移给本对象
    */
   explicit ProbDataTSP(std::unique_ptr<Graph> g)
      : graph_(std::move(g))
   {}

   virtual ~ProbDataTSP() = default;

   /** @return 图的可变指针（非拥有，调用方不负责释放） */
   Graph* getGraph() noexcept { return graph_.get(); }
   /** @return 图的只读指针 */
   const Graph* getGraph() const noexcept { return graph_.get(); }

   /**
    * @brief 拷贝问题数据到另一个 SCIP 实例
    *
    * 深拷贝整个图结构（节点、边、邻接表）。
    * 用于并行求解中的子 SCIP 实例，或分支定界中的拷贝。
    */
   virtual SCIP_RETCODE scip_copy(
      SCIP*              scip,
      SCIP*              sourcescip,
      SCIP_HASHMAP*      varmap,
      SCIP_HASHMAP*      consmap,
      scip::ObjProbData** objprobdata,
      SCIP_Bool          global,
      SCIP_RESULT*       result
      ) override;

   /**
    * @brief 释放原始问题数据
    *
    * SCIP 在释放原始问题时调用。释放图中所有边的 SCIP 变量引用，
    * 然后释放图本身。
    */
   virtual SCIP_RETCODE scip_delorig(SCIP* scip) override;

   /**
    * @brief 释放变换后的问题数据
    *
    * 与 scip_delorig 类似，但在求解结束时对变换数据调用。
    */
   virtual SCIP_RETCODE scip_deltrans(SCIP* scip) override;

   /**
    * @brief 将原始问题数据变换为 presolved 数据
    *
    * 创建图的新副本，并将变量映射到变换后的变量空间。
    */
   virtual SCIP_RETCODE scip_trans(
      SCIP*              scip,
      scip::ObjProbData** objprobdata,
      SCIP_Bool*         deleteobject
      ) override;

private:
   std::unique_ptr<Graph> graph_;  ///< TSP 图数据的所有权
};

} // namespace tsp
