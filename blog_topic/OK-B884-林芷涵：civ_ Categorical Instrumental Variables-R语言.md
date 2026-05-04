# OK-B884-林芷涵：civ: Categorical Instrumental Variables-R语言

> **作者：** 林芷涵（首都经济贸易大学）
> 
> **邮箱：** <zhihan_lin0211@163.com>

&emsp;

> **编者按**：本文主要参考了如下资料，特此致谢！
> - Wiemann T (2023). Optimal Categorical Instrumental Variables. *arXiv:2311.17021*. [Link](https://arxiv.org/abs/2311.17021)
> - Belloni A, Chen D, Chernozhukov V, Hansen C (2012). "Sparse Models and Methods for Optimal Instruments With an Application to Eminent Domain." *Econometrica*, 80(6), 2369-2429. [Link](https://doi.org/10.3982/ECTA9493)
> - GitHub: [https://github.com/thomaswiemann/civ](https://github.com/thomaswiemann/civ)


&emsp;

- **Title**: civ-最优分类工具变量估计-R语言
- **Keywords**: 工具变量、类别型工具、弱工具变量、K条件均值、因果推断、R语言

&emsp;

## 1. 简介

### 1.1 什么是工具变量？为什么要"分类"？

在实证研究中，我们经常面临**内生性问题**——即核心解释变量（如"是否接受培训"）与不可观测的因素（如"能力"）相关联，导致普通回归无法准确识别因果效应。

**工具变量（IV）**是解决这一问题的经典方法。其核心思想是找到一个"替身"变量 Z：
- Z 必须与内生变量 D 相关（相关性）
- Z 只能通过 D 影响结果 Y，不能有其他渠道（排他性）

通俗地说：Z 像一根"导管"，通过它我们才能"看到" D 对 Y 的真实影响。

### 1.2 "类别型"工具变量的特殊困难

当工具变量 Z 是**类别型变量**（如法官ID、学校编号、地区代码），且每个类别下的样本量较少时，传统方法会遇到**弱工具变量**问题：

| 场景举例 | 工具变量 Z | 问题 |
|---------|-----------|------|
| 法官固定效应 | 法官ID（可能有上百个） | 每个法官只审少量案件 |
| 学校效应 | 学校编号 | 每所学校样本有限 |
| 地理分组 | 地区代码 | 区域样本不均匀 |

**为什么是问题？** 想象你有一百个班级，每个班只有5个学生，现在要估计"班级"对成绩的影响——样本太分散，估计自然不准确。

这类场景在**司法经济学、教育经济学**等领域非常常见：
- **法官固定效应设计**：法官指派作为工具，分析审前拘留对定罪的影响
- **教育政策评估**：出生季度与州、年份的交互项作为工具，估计教育回报
- **区域分组研究**：地理区域或班级分组作为工具变量

### 1.3 CIV 的解决思路

**civ** 包实现了 Wiemann (2023) 提出的**最优分类工具变量（CIV）估计量**。其核心思想分三步：

1. **假设**：大量类别背后存在少数几个"潜在类型"，这些潜在类型已经包含了预测处理变量的全部信息
2. **聚类**：通过 **KCMeans 算法** 自动将几十甚至上百个类别合并为几个潜在类别
3. **估计**：用合并后的低维工具进行IV估计，获得无偏且有效的处理效应

**一句话概括**：CIV 的智慧在于"化繁为简"——与其费力估计每个类别单独的效果，不如找到类别背后的"共同模式"。

**关键 R 语言命令**：`civ`


## 2. 分类工具变量的理论框架

### 2.1 问题设定：弱工具困境

标准的IV模型包含两个方程：

**第二阶段（结果方程）**：
$$Y_i = D_i \tau_0 + X_i^\top \beta_0 + U_i$$

- $Y_i$：我们关心的结果（如工资、是否定罪）
- $D_i$：内生处理变量（如是否接受培训、审前拘留时长）
- $X_i$：可观测的协变量（年龄、性别等）
- $U_i$：误差项（包含不可观测因素）

**第一阶段（处理方程）**：
$$D_i = m_0(Z_i) + X_i^\top \gamma_0 + V_i$$

- $Z_i$：**类别型工具变量**（如法官ID、学校代码）
- $m_0(Z_i)$：工具变量对处理变量的影响——这是IV估计的关键

**工具变量的核心假设**：$\mathbb{E}[U_i \mid Z_i, X_i] = 0$ 和 $\mathbb{E}[V_i \mid Z_i, X_i] = 0$，即工具与误差项不相关。

**"弱工具"问题**：当 $Z_i$ 的类别很多（如100个法官）但每个类别的样本很少（每个法官只审5个案件）时，第一阶段方程中 $m_0(Z_i)$ 的估计会变得极不精确。这会导致：
- 处理效应的估计值严重偏离真实值
- 置信区间变宽，可靠性下降

### 2.2 CIV 的关键假设：潜在低维工具

Wiemann (2023) 的核心洞察是：虽然观测到的类别很多（如100个法官），但它们**本质上只代表少数几种"风格"或"类型"**。

**假设**：存在一个**潜在的低维工具变量** $Z_i^{(0)}$，满足：
$$|\mathcal{Z}^{(0)}| = K_0 < |\mathcal{Z}|, \quad \text{且} \quad \mathbb{E}[D_i \mid Z_i] = \mathbb{E}[D_i \mid Z_i^{(0)}]$$

用人话说就是：
- 100个法官实际上可以归类为 $K_0$（比如5-10）种"审判风格"
- 这5-10种风格已经包含了法官对案件结果的全部影响
- 多余的类别只是同一风格的重复观测

**这个假设的合理性**：在法官IV设计中，不同法官可能只是"严厉程度"不同，而严厉程度只有几种典型水平——这比记住每个法官的ID更有经济意义。

### 2.3 K条件均值（KCMeans）估计

**核心思想**：将工具变量的"归类"问题转化为一个**优化问题**——找到一种最优方式，把原始的众多类别合并为 K 个潜在类别。

**通俗类比**：这就像把100首歌曲分成5个播放列表，目标是让每个列表内的歌曲"风格一致"——在这里，"风格一致"意味着同一潜在类别的观测点有相似的处理变量 D 值。

**目标函数**：给定允许的类别数 K，KCMeans 通过最小化预测误差来寻找最优合并方案：
$$\hat{m}_K = \arg\min_{m: \mathcal{Z} \to \mathbb{R},\ |m(\mathcal{Z})| \le K} \mathbb{E}_n\left[ (D_i - X_i^\top \hat{\gamma} - m(Z_i))^2 \right]$$

翻译成人话：找到一个合并规则 m，使得合并后的"伪工具" $m(Z_i)$ 对处理变量 D 的预测最准确。

**为什么用动态规划？** 不同于传统K-Means的随机迭代、可能陷入局部最优，KCMeans 利用"类别是有序的"这一特点（法官ID只是编号，本身有顺序），可以用动态规划**精确求解**，保证找到全局最优解。

### 2.4 CIV 估计量

得到 $\hat{m}_K$ 后，就相当于我们有了"降维"后的工具变量 $\hat{m}_K(Z_i)$，用它替换原始工具变量进行IV估计：

$$
\hat{\tau}_{\text{CIV}} = \frac{\mathbb{E}_n\left[ (Y_i - \bar{Y})\,(\hat{m}_K(Z_i) - \bar{\hat{m}}_K) \right]}{\mathbb{E}_n\left[ (D_i - \bar{D})\,(\hat{m}_K(Z_i) - \bar{\hat{m}}_K) \right]}
$$

这个公式的本质是：计算 Y 和"降维工具"的协方差，除以 D 和"降维工具"的协方差——这正是IV估计的标准形式。

### 2.5 渐近性质（理论保证）

**定理**（Wiemann 2023, Theorem 1）：在适当的正则条件下，CIV 估计量具有良好的理论性质：

$$\sqrt{n}\,(\hat{\tau}_{\text{CIV}} - \tau_0) \ \xrightarrow{d}\  N\!\left(0,\ \frac{\mathrm{Var}(U_i\,m_0(Z_i))}{[\mathrm{Cov}(D_i, m_0(Z_i))]^2}\right)$$

翻译：随着样本量增加，CIV 估计量会收敛到真实值 $\tau_0$，且收敛速度为 $\sqrt{n}$（标准速度），服从正态分布。

**关键结论（ practical 意义）**：

| 情况 | 理论性质 | 实际意义 |
|------|---------|---------|
| $K = K_0$ 且同方差 | 达到半参数效率界 | 估计量是最优的，无法做得更好 |
| $K > K_0$（多分了类别） | 仍 $\sqrt{n}$-相合 | 估计正确，但浪费了一些样本 |
| $K < K_0$（少分了类别） | 仍 $\sqrt{n}$-相合，但方差增大 | 估计正确但精度下降 |

**一句话总结**：只要 K 不比真实类别数少，CIV 就能给出正确的因果推断！

## 3. R 语言实现：`civ` 包

### 3.1 安装与加载

```R
# 从 CRAN 安装稳定版
install.packages("civ")

# 从 GitHub 安装最新开发版（推荐）
if (!require("devtools")) install.packages("devtools")
devtools::install_github("thomaswiemann/civ", dependencies = TRUE)

# 加载包
library(civ)
library(AER)  # 用于异方差稳健标准误
```

### 3.2 核心函数：`civ()`

```R
civ(y, D, Z, X = NULL, K, ...)
```

- `y`：结果变量（数值型向量）
- `D`：内生处理变量（数值型向量）
- `Z`：类别型工具变量（因子型向量）
- `X`：外生协变量矩阵（可选）
- `K`：潜在类别数（需要事先设定或通过数据驱动选择）

### 3.3 模拟数据生成

以下模拟遵循 Wiemann (2023) 的数据生成过程，用于验证 CIV 的表现。

**数据生成逻辑**：
- 真实世界存在一个二元潜在工具 Z0（如两种法官风格）
- 观测到的工具 Z 是 Z0 的"嘈杂版本"（40个类别，但本质上只有2类）
- 这模拟了现实中的弱工具场景：类别很多，但背后只有少数几种"模式"

```R
# 设置随机种子（保证结果可复现）
set.seed(51944)

# ========== 参数设定 ==========
nobs <- 800                    # 样本量：800个观测
C <- 0.858                     # 第一阶段：潜在工具对 D 的影响系数
sgm_V <- 0.9                    # 误差项 V 的标准差

# 处理效应随 X 变化：X=1 时效应为 0.5，X=2 时效应为 1.5
tau_X <- c(-0.5, 0.5) + 1      # E[tau(X)] = 1（平均处理效应）

# ========== 生成数据 ==========

# 1. 协变量 X（二元：可理解为两组不同特征的群体）
X <- sample(1:2, nobs, replace = TRUE)

# 2. 观测到的工具变量 Z（40个类别，但本质由 X 和 Z0 决定）
#    model.matrix 创建 X 和潜在类别的交互项，再压缩成单一类别ID
Z <- model.matrix(~ 0 + as.factor(sample(1:20, nobs, replace = TRUE)):as.factor(X))
Z <- Z %*% c(1:ncol(Z))

# 3. 潜在工具 Z0（二元：只反映两种基本类型）
Z0 <- Z %% 2                    # 取模运算：偶数为0，奇数为1

# 4. 相关误差项（U 和 V 相关系数为 0.6）
U_V <- matrix(rnorm(2 * nobs, 0, 1), nobs, 2) %*%
  chol(matrix(c(1, 0.6, 0.6, sgm_V), 2, 2))

# 5. 生成处理变量 D 和结果变量 Y
D <- Z0 * C + U_V[, 2]          # D 受潜在工具 Z0 和误差 V 影响
y <- D * tau_X[X] + U_V[, 1]   # Y 受处理效应 tau_X[X] 和误差 U 影响

# 记录真实处理效应（用于后续比较）
true_effect <- mean(tau_X[X])
true_effect  # 约为 1.0325
```

### 3.4 运行 CIV 估计

```R
# 将 X 转换为因子型（用于回归控制）
X_factor <- as.factor(X)

# 运行 CIV 估计（假设 K=2）
civ_fit <- civ(y = y, D = D, Z = Z, X = X_factor, K = 2)

# 计算稳健标准误
civ_res <- summary(civ_fit, vcov = vcovHC(civ_fit$iv_fit, type = "HC1"))
```

### 3.5 输出结果解读

```R
# CIV 估计结果
c(Estimate = civ_res$coef[2, 1],       # 点估计
  "Std. Error" = civ_res$coef[2, 2],  # 稳健标准误
  "t-val." = abs(civ_res$coef[2, 1] - true_effect) / civ_res$coef[2, 2])
```

输出：
```
  Estimate Std. Error     t-val.
 1.0063143 0.1086868 0.2409285
```

**结果解读**：
- 估计值为 **1.006**，接近真实处理效应 **1.032**
- t 值 = 0.24 < 1.96，说明置信区间覆盖真实值（95% CI 包含真值）
- TSLS 和机器学习方法在该场景下会产生严重偏倚（见下文比较）

### 3.6 聚类结果验证

这一步验证 CIV 的聚类算法是否成功恢复了真实的潜在类别结构。

```R
# 查看每个工具类别的观测数（验证"弱工具"场景）
obs_per_Z <- table(Z)
hist(obs_per_Z, breaks = 20, main = "每个工具类别的观测数",
     xlab = "观测数", ylab = "类别数")

# 查看聚类效果：CIV 估计的潜在类别 vs 真实潜在类别
Z0_hat <- predict(civ_fit$kcmeans_fit, Z, clusters = TRUE) - 1

# 错分样本数（misclassification rate）
missclassified <- sum((Z0 - Z0_hat) != 0)
missclassified  # 仅为 16/800 = 2%
```

**解读**：misclassification rate 只有 2%，说明 KCMeans 成功识别了真实的两类潜在工具！

### 3.7 与其他方法比较

下面比较四种方法在相同数据下的表现，直观展示弱工具变量场景下各方法的优劣。

#### 3.7.1 TSLS（两阶段最小二乘）—— 传统基线

```R
# TSLS 估计（直接用原始类别型工具，不做任何降维处理）
tsls_fit <- ivreg(y ~ D + X_factor | as.factor(Z) + X_factor)
tsls_res <- summary(tsls_fit, vcov = vcovHC(tsls_fit, type = "HC1"))$coef

c(Estimate = tsls_res[2, 1],
  "Std. Error" = tsls_res[2, 2],
  "t-val." = abs(tsls_res[2, 1] - true_effect) / tsls_res[2, 2])
```

输出：
```
   Estimate  Std. Error     t-val.
  1.11620022 0.08664093 0.96605856
```

**分析**：TSLS 偏倚约 8%，尚可接受，但标准误较大。

#### 3.7.2 Post-Lasso IV —— 稀疏学习方法

```R
# Post-Lasso IV（Belloni et al., 2012）：用 Lasso 选择显著的类别
library(hdm)
hdm_fit <- rlassoIV(y ~ D + X_factor | as.factor(Z) + X_factor, select.X = FALSE)
hdm_res <- c(hdm_fit$coefficients, hdm_fit$se)

c(Estimate = hdm_res[1],
  "Std. Error" = hdm_res[2],
  "t-val." = abs(hdm_res[1] - true_effect) / hdm_res[2])
```

输出：
```
  Estimate.D Std. Error.D   t-val..D
   0.7717395    0.1059556    2.4610350
```

**分析**：偏倚高达 **25%**，且 t 值 = 2.46 > 1.96，统计检验会错误地拒绝真值！

#### 3.7.3 随机森林 IV —— 非线性机器学习方法

```R
# 随机森林 IV：用随机森林拟合第一阶段，再用预测值作为工具
library(ranger)
df <- data.frame(D = D, Z = as.factor(Z), X = X_factor)
mhat <- predict(ranger(D ~ Z + X, data = df), df)$predictions
ranger_fit <- ivreg(y ~ D + X_factor | mhat + X_factor)
ranger_res <- summary(ranger_fit, vcov = vcovHC(ranger_fit, type = "HC1"))$coef

c(Estimate = ranger_res[2, 1],
  "Std. Error" = ranger_res[2, 2],
  "t-val." = abs(ranger_res[2, 1] - true_effect) / ranger_res[2, 2])
```

输出：
```
  Estimate Std. Error     t-val.
 1.25887072 0.09053253 2.50043521
```

**分析**：偏倚约 **22%**，且 t 值 > 1.96，同样会错误拒绝真值。

#### 3.7.4 方法比较总结

| 方法 | 估计值 | 偏离真值 | 95% CI 覆盖 | 评判 |
| :--- | :---: | :---: | :---: | :---: |
| **CIV** | **1.006** | **2.5%** | **是** | **最优** |
| TSLS | 1.116 | 8.1% | 是 | 偏倚可接受 |
| Post-Lasso IV | 0.772 | 25.2% | **否** | 失效 |
| 随机森林 IV | 1.259 | 21.9% | **否** | 失效 |

**结论**：在弱工具变量场景下，**只有 CIV** 能够同时保证低偏倚和正确的统计推断！

### 3.8 Oracle 估计（理论基准）

```R
# Oracle 估计（假设已知潜在工具 Z0）
oracle_fit <- ivreg(y ~ D + X_factor | as.factor(Z0) + X_factor)
oracle_res <- summary(oracle_fit, vcov = vcovHC(oracle_fit, type = "HC1"))

c(Estimate = oracle_res$coef[2, 1],
  "Std. Error" = oracle_res$coef[2, 2],
  "t-val." = abs(oracle_res$coef[2, 1] - true_effect) / oracle_res$coef[2, 2])
```

输出：
```
  Estimate Std. Error     t-val.
 1.0184144 0.1065979 0.1321376
```

**结论**：CIV 估计量（1.006）几乎与 Oracle 估计量（1.018）相同，证明 CIV 能够有效恢复潜在最优工具。

## 4. 实际应用：法官固定效应设计

以 Dobbie et al. (2018) 的研究为例，说明如何在实际应用中使用 `civ` 包。

### 4.1 研究背景

该研究分析了审前拘留对定罪的影响。关键挑战是：法官的指派本身可能与被告特征相关，直接使用法官ID作为工具变量会面临弱工具问题（每个法官审理的案件数有限）。

### 4.2 数据结构

```R
# 假设数据结构
# y: 结果变量（定罪与否）
# D: 处理变量（审前拘留时长）
# Z: 工具变量（法官ID，类别数很多）
# X: 协变量（年龄、犯罪历史等）

# 示例数据生成
set.seed(123)
n <- 5000
judge_ids <- sample(1:200, n, replace = TRUE)  # 200个法官
X <- cbind(age = rnorm(n, 35, 10),
           prior = rbinom(n, 1, 0.3))
Z <- as.factor(judge_ids)

# 内生处理变量 D（受法官效应和协变量影响）
D <- 0.3 * (as.numeric(Z) %% 5) + 0.5 * X[, 1] / 10 + rnorm(n)

# 结果变量 Y
Y <- 2 * D - 0.5 * X[, 1] / 10 + rnorm(n)
```

### 4.3 CIV 估计

```R
# 运行 CIV 估计
# K 的选择：可以基于理论（如法官风格数量）或数据驱动（如 BIC/AIC）
civ_result <- civ(y = Y, D = D, Z = Z, X = X, K = 10)

# 输出结果
summary(civ_result, vcov = vcovHC(civ_result$iv_fit, type = "HC1"))
```

### 4.4 K 值选择建议

```R
# 尝试不同的 K 值
K_values <- c(5, 10, 15, 20)
results <- sapply(K_values, function(K) {
  fit <- civ(y = Y, D = D, Z = Z, X = X, K = K)
  return(coef(summary(fit$iv_fit))[2, ])
})

rownames(results) <- c("Estimate", "Std. Error", "t-value")
colnames(results) <- paste("K=", K_values, sep = "")
results
```

**注意**：实际应用中，应结合理论背景和经济含义选择 K，避免过度拟合。

## 5. 模型贡献

- **解决弱工具变量问题**：通过假设潜在低维工具存在，CIV 将复杂的弱工具IV问题转化为可解的低维回归问题，在类别多、每类样本少的场景下表现优异。
- **理论保证的渐近最优性**：在适当正则条件下，CIV 估计量达到半参数效率界，是渐近最优的估计量。
- **计算高效**：KCMeans 算法使用动态规划精确求解，时间复杂度为 $O(|\mathcal{Z}|^2 K)$，避免了迭代 heuristic 方法的收敛问题。
- **实用可靠的推断**：基于稳健标准误的置信区间具有正确的覆盖率，适用于实际政策评估。

## 6. 模型局限

- **K 值选择依赖假设**：需要事先设定潜在类别数 K，K 的选择影响估计精度；虽然可以通过数据驱动方法（如 BIC）选择，但缺乏统一的最优准则。
- **潜在工具假设不可直接检验**：假设存在有限支撑的潜在工具变量，无法由观测数据直接验证，需依赖领域知识。
- **同方差效率假设**：达到半参数效率界需要同方差假设，在异方差情况下效率可能下降。
- **高维协变量场景**：当外生协变量 X 的维度较高时，第一阶段估计的收敛速度可能受限。

## 7. 结论

CIV 估计量为处理弱工具变量问题提供了一个兼具理论严谨性和实践可操作性的解决方案。通过引入潜在低维工具假设和 KCMeans 聚类算法，CIV 能够在传统方法失效的场景下提供可靠的因果推断。

该方法特别适用于：
- 法官/审查员固定效应设计
- 学校/班级作为工具变量的教育研究
- 地理区域分组作为工具的实证研究

`civ` R 包实现了该方法，提供了简洁的接口和完整的估计流程，是实证研究者处理类别型工具变量的有力工具。

## 参考文献

1. Angrist, J. D., & Krueger, A. B. (1991). Does Compulsory Schooling Affect Earnings? *The Quarterly Journal of Economics*, 106(2), 449-479. [Link](https://doi.org/10.2307/2937944)

2. Belloni, A., Chen, D., Chernozhukov, V., & Hansen, C. (2012). Sparse Models and Methods for Optimal Instruments With an Application to Eminent Domain. *Econometrica*, 80(6), 2369-2429. [Link](https://doi.org/10.3982/ECTA9493)

3. Dobbie, W., Goldin, J., & Yang, C. S. (2018). The Effects of Pretrial Detention on Conviction and Outcomes. *American Economic Review*, 108(2), 351-387. [Link](https://doi.org/10.1257/aer.20161559)

4. Wiemann, T. (2023). Optimal Categorical Instrumental Variables. *arXiv:2311.17021*. [Link](https://arxiv.org/abs/2311.17021)

5. Imbens, G. W., & Angrist, J. D. (1994). Identification and Estimation of Local Average Treatment Effects. *Econometrica*, 62(2), 467-475. [Link](https://doi.org/10.2307/2951620)
