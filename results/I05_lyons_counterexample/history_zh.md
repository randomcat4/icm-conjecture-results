# Lyons 的 DPP 熵凹性猜想：提出、发展与本项目的结论

本文说明本仓库所处理问题的来历、几个容易混淆的版本、公开文献中可核实的推进，以及本项目最终证明了什么。文中的“I05 猜想”指 Russell Lyons 关于**有限标号集合上的行列式概率测度之完整子集分布 Shannon 熵**的猜想。它不是 ICM 官方发布的问题；“ICM 猜想”只是因为 Lyons 在 2014 年国际数学家大会邀请报告中再次陈述了它。

## 1. 问题的准确形式

设有限标号集合为 $E$，$Q$ 是 $\ell^2(E;\mathbb C)$ 或 $\ell^2(E;\mathbb R)$ 上的正压缩，即 $Q=Q^*$ 且 $0\preceq Q\preceq I$。相应的行列式概率测度 $\mathbb P^Q$ 由

$$
\mathbb P^Q(A\subseteq X)=\det Q[A]\qquad(A\subseteq E)
$$

确定。这里 $X$ 是取值于 $2^E$ 的随机子集，$Q[A]$ 是由标号坐标 $A$ 取出的主子矩阵。完整分布的 Shannon 熵为

$$
\mathrm{Ent}(Q)
=-\sum_{A\subseteq E}\mathbb P^Q(X=A)\log \mathbb P^Q(X=A).
$$

Lyons 猜想：对任意两个正压缩 $Q_1,Q_2$，

$$
\mathrm{Ent}\!\left(\frac{Q_1+Q_2}{2}\right)
\geq \frac{\mathrm{Ent}(Q_1)+\mathrm{Ent}(Q_2)}2.
$$

中点形式等价于通常的线段凹性。这里研究的是 $2^{|E|}$ 个标号子集的完整概率分布之熵，不是点数 $|X|$ 的熵、核的特征值熵、von Neumann 熵，也不是无限平稳过程的每格熵率。

## 2. 2002–2003：有限猜想的提出

Lyons 的论文 *Determinantal Probability Measures* 于 2002 年 4 月提交 arXiv，2003 年发表于 *Publications Mathématiques de l'IHÉS* 第 98 卷，第 167–212 页。论文主要建立离散 DPP 与拟阵、随机支配、负相关、尾事件和正压缩扩张之间的基础理论；熵凹性出现在第 9 节“Open questions: General theory”中，而不是论文主定理。

在印刷页 201，Lyons 先定义有限 $E$ 上的完整分布熵，随后明确写道数值计算支持该命题，并将它列为 **Conjecture 9.2**。因此最初证据是数值观察；原文没有给出一般证明，也没有报告一个证明计划。

原始资料：

- [期刊页面与完整书目信息](https://pmihes.centre-mersenne.org/articles/10.1007/s10240-003-0016-0/)
- [期刊原始 PDF](https://www.numdam.org/item/10.1007/s10240-003-0016-0.pdf)
- [arXiv:math/0204325 页面](https://arxiv.org/abs/math/0204325)及[原始 PDF](https://arxiv.org/pdf/math/0204325)

同一篇论文在印刷页 209 又讨论了阿贝尔群上的平稳 DPP。若 $f$ 是对偶群上的 $[0,1]$ 值函数，傅里叶变换给出 Toeplitz 型核和一个平稳过程 $\mathbb P^f$。Lyons 指出：相应 Kolmogorov–Sinai 熵率的凹性可由有限 Conjecture 9.2 推出，甚至只需有限猜想在 Toeplitz 核上的限制版。

## 3. 2002–2003：Lyons–Steif 的平稳过程版本

Russell Lyons 与 Jeffrey E. Steif 的 *Stationary Determinantal Processes: Phase Multiplicity, Bernoullicity, Entropy, and Domination* 同样于 2002 年 4 月提交 arXiv，2003 年发表于 *Duke Mathematical Journal* 120(3), 515–575。文章研究 $\mathbb Z^d$ 上的平稳 DPP，证明这些过程是 Bernoulli shifts，并研究熵估计和随机支配。

该文第 9 节把“计算 $H(\mathbb P^f)$”列为 Question 9.1，并把

$$
H\!\left(\mathbb P^{(f+g)/2}\right)
\geq \frac{H(\mathbb P^f)+H(\mathbb P^g)}2
$$

列为 **Conjecture 9.2**。这与 Lyons 单人论文中的有限 Conjecture 9.2 **编号相同但命题不同**：

- 单人论文的 Conjecture 9.2 是任意有限正压缩核的完整子集熵凹性；
- Lyons–Steif 的 Conjecture 9.2 是特殊平稳 Toeplitz 族的熵率凹性。

前者蕴含后者；反方向没有在这些来源中建立。本仓库的五维反例不是 Toeplitz 平稳族的反例，因此不裁决 Lyons–Steif 的平稳熵率猜想。

原始资料：

- [arXiv:math/0204324 页面](https://arxiv.org/abs/math/0204324)
- [原始 PDF](https://arxiv.org/pdf/math/0204324)

## 4. 2012：机器学习文献记录了猜想，但没有推进证明

Alex Kulesza 的宾夕法尼亚大学博士论文 *Learning with Determinantal Point Processes* 把 DPP 用于多样化子集选择和机器学习。论文第 2.5.1 节将 Lyons 的猜想列为 Conjecture 2.1，写道数值模拟有力支持该猜想，但作者不知道证明；第 7.1 节又把熵凹性列为未来工作。

这条记录说明该猜想进入了 DPP 机器学习研究者的视野；论文没有证明新的特殊情形。

原始资料：[论文 PDF](https://www.alexkulesza.com/pubs/thesis.pdf)；[作者的论文目录与书目信息](https://www.alexkulesza.com/)。

## 5. 2014：在 ICM 邀请报告中重新提出

Lyons 的邀请报告 *Determinantal Probability: Basic Properties and Conjectures* 发表于 2014 年国际数学家大会论文集第四卷，第 137–161 页。报告概述行列式概率测度和点过程的基础性质，并集中列出若干问题与猜想。

报告在印刷页 141 重新定义有限 $E$ 上的完整子集熵，把原来的有限 Conjecture 9.2 重编号为 **Conjecture 2.6**，仍以“数值计算支持”为引导语。报告在第 5 节又把平稳熵率版本列为 **Conjecture 5.3**，并明确说它可由 Conjecture 2.6 推出。

因此，“Lyons Conjecture 2.6”不是 2014 年新出现的命题，而是 2003 年有限熵猜想的重述。进入 ICM 论文集提高了它的可见度，却不表示它是 ICM 官方问题或当时已有一批研究者围绕它工作。

原始资料：

- [arXiv:1406.2707 页面](https://arxiv.org/abs/1406.2707)
- [Lyons 作者主页上的已发表 PDF](https://rdlyons.pages.iu.edu/pdf/icm-pub.pdf)
- [arXiv PDF](https://arxiv.org/pdf/1406.2707)

## 6. 2020：Yuzhou Gu 给出三个实质性部分结果

Yuzhou Gu 的 MIT 6.881 课程项目报告 *Entropy of Determinantal Point Processes* 直接以 Lyons 猜想为目标。报告证明：

1. 点数随机变量 $|X(K)|$ 的熵关于 $K$ 凹；
2. 一条线段的一个端点为零核时，完整 DPP 熵满足所需凹性；
3. 当 $\mathrm{rank}(K_1-K_2)=1$ 时，完整 DPP 熵沿该线段凹。

第一条处理的是点数熵而非完整子集熵；后两条是真正的完整熵特殊情形。它们与本仓库的反例兼容，因为这里的反例方向是秩四的纯虚 Hermitian 方向。该报告是课程项目稿，不是经同行评议的期刊论文。

原始资料：[报告 PDF](https://sevenkplus.com/data/dpp.pdf)。

## 7. 2026：本项目构造并严格认证一个五维反例

本项目否定了原猜想量词明确包含的**复 Hilbert 空间分支**。具体地，我们给出：

- 一个有理实对称 $5\times5$ 内部核 $K$；
- 一个纯虚 Hermitian、秩为四的有理方向 $A=iB$；
- 一个固定有理步长 $h=10^{-5}$；
- 两个严格正压缩端点 $K-hA$ 与 $K+hA$，其中点为 $K$；
- 严格不等式

$$
H(K)<\frac{H(K-hA)+H(K+hA)}2.
$$

最终证明不依赖浮点数或继续穷举：

1. Sylvester 判据以正有理数证明中心和端点严格位于 $0$ 与 $I$ 之间；
2. 逐项精确计算五个标号坐标的全部 32 个事件概率，证明每项为正、总和为一，且两个端点逐事件相等；
3. 用带显式几何余项的 80 项 $\mathrm{atanh}$ 展开给出有向有理对数区间；
4. 由此严格证明

   $$
   \frac1{200}<D^2H(K)[A,A]<\frac1{190}
   $$

   以及

   $$
   H(K+hA)-H(K)>10^{-13}.
   $$

后一个不等式和端点逐项相等直接给出有限弦上的中点凹性违反，而正 Hessian 方向给出同一现象的局部二阶证书。公开包包含矩阵、32 项证书、纯标准库复核程序和可重建论文。两个隔离的数学核验分别重构了核心计算，之后另有一次新上下文的论文级复核。

这个结果当前应准确称为“构造并严格认证了一个反例”，还不是“刻画了全部反例”。它没有完成下列事项：

- 没有给出纯实对称方向的反例，因而没有裁决另行限制到实对称核的版本；
- 没有证明维数五最小；
- 没有分类哪些 $K$ 或哪些方向必然保持熵凹性；
- 没有反驳 Lyons–Steif 的平稳 Toeplitz 熵率猜想；

完整公开证明见 [main.pdf](main.pdf)，复核范围见 [verification_summary.md](verification_summary.md)，可执行证书见 [scripts/verify_exact.py](scripts/verify_exact.py)。

## 8. 作者简介

**Russell Lyons** 是美国概率论学者，现任 Indiana University Bloomington 的 James H. Rudy Professor of Mathematics，并兼任 Adjunct Professor of Statistics。他的研究横跨图上的概率、组合学、统计力学、遍历论、调和分析和几何群论。有限 DPP 熵凹性猜想由他在 2003 年提出，并在 2014 年 ICM 邀请报告中重述。[Indiana University 官方简介](https://math.indiana.edu/about/faculty/lyons-russell.html)

他的论文 *Fourier–Stieltjes Coefficients and Asymptotic Distribution Modulo 1* 发表于 *Annals of Mathematics* 122 (1985), 155–170；提出本猜想的 2003 年 DPP 论文自身发表于 *Publications Mathématiques de l'IHÉS*。他还是 2014 年 ICM 邀请报告人，并自 2022 年起担任 *Annals of Mathematics* 副编辑。[Annals 1985 年卷页](https://annals.math.princeton.edu/1985/122-1)；[Lyons 官方履历](https://rdlyons.pages.iu.edu/pdf/cv-web.pdf)

**Jeffrey E. Steif** 是 Chalmers University of Technology 数学科学系分析与概率论教授，也是瑞典皇家科学院数学类成员。他的研究包括概率论、渗流、噪声敏感性、随机过程和遍历性质。他与 Lyons 共同提出的是相关的无限平稳 DPP 熵率猜想。[Chalmers 官方简介](https://www.chalmers.se/en/persons/steif/)；[瑞典皇家科学院简介](https://www.kva.se/en/contact/jeffrey-steif-2/)

**Yuzhou Gu** 的研究领域包括信息论、统计、概率和计算机科学。他于 2023 年获 MIT EECS 博士学位，曾任 Institute for Advanced Study 数学学院成员。2020 年还是 MIT 学生时，他完成了上述三个部分结果。[IAS 简介](https://www.ias.edu/scholars/yuzhou-gu)；[个人主页](https://sevenkplus.com/)

本仓库中的反例由本项目在机器辅助探索、严格有理化和相互隔离的复核流程中获得。公开稿目前保留作者占位符；该占位符不应被解释为作者名单、署名决定或优先权声明。

所有来源的集中索引和直接 PDF 链接见 [sources/README.md](sources/README.md)。
