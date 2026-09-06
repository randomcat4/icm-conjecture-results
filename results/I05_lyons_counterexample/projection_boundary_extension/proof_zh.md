# 续篇：投影边界的实凹性定理与显式连续实／复分离族

本文承接本公开包中的[五维复反例](../README.md)与[实／复机制背景](mechanism_background_zh.md)。
原始研究提交分别为 `a296f7cc555312e1dbd5a28af08836c4756659c7` 和
`9c27d3f25cb68ae8fe1d6596ffbb49d023d45892`；公开迁移没有改变定理、矩阵或证书算法。

**范围：下述结构证明与纯有理全参数证书不解决实限制版的全局凹性；本文没有实反例，
也没有成功的复到实反例转移。公开核验状态见 [README](README.md)。**

## 1. 新结果及准确范围

### 定理 A：任意维数、连通投影的严格内向邻域

设实正交投影 \(P\) 的非零非对角元图连通，\(Q=I-P\)。一维的 \(P=0,1\) 直接归结为 Bernoulli 熵；非平凡情形为 \(0<\mathrm{rank}P<n\)。设实对称 \(N\) 满足
\[
PNP\prec0\quad\text{在 }\mathrm{ran}P\text{ 上},\qquad
QNQ\succ0\quad\text{在 }\ker P\text{ 上}. \tag{A1}
\]
允许 \(PNQ\ne0\)。则存在 \(\varepsilon_0(P,N)>0\)，使
\[
0<P+\varepsilon N<I,\qquad
D^2H(P+\varepsilon N)[C,C]<0 \tag{A2}
\]
对所有 \(0<\varepsilon<\varepsilon_0\) 及所有非零实对称 \(C\) 成立。

\(\varepsilon_0\) 可在严格内向方向的任意紧集上统一选择。例如固定 \(\kappa>0\)，
对
\[
\|N\|_{\rm op}\le1,\quad PNP\preceq-\kappa P,\quad QNQ\succeq\kappa Q \tag{A3}
\]
中的全部 \(N\)，存在共同的正阈值。这排除了沿这样的**非切向内向锥**趋近一个固定
连通投影的实反例序列。它不是所有趋近投影的路径定理：内向余量相对距离趋零的
切向路径、退化为不连通投影的极限、以及不在足够小邻域的核都没有被排除。

### 推论 B：每个实投影的各向同性正则化

对任意实正交投影，包括不连通的投影，令
\[
K_\varepsilon=\varepsilon I+(1-2\varepsilon)P.
\]
当 \(\varepsilon>0\) 充分小时，其实 Hessian 半负定。零空间恰好是连接 \(P\)
不同坐标连通块的实对称方向。若 \(P\) 连通，则负定。

### 定理 C：一个显式、全参数认证的五维分离族

令
\[
V=\begin{pmatrix}
2719&-3449\\-3449&5840\\1009&-2889\\-818&-1513\\-2490&1322
\end{pmatrix},\qquad
P=V(V^TV)^{-1}V^T,\quad Q=I-P, \tag{C1}
\]
并取原冻结方向的实反对称部分
\[
B_*=\frac1{50}\begin{pmatrix}
0&-3&-13&-16&4\\3&0&-4&15&17\\13&4&0&6&-15\\
16&-15&-6&0&1\\-4&-17&15&-1&0
\end{pmatrix},\qquad
\widehat B=PB_*P+QB_*Q. \tag{C2}
\]
这是一个完全有理、固定的构造，不是浮点谱投影。

对**每个**
\[
\boxed{0<\varepsilon\le2^{-32}} \tag{C3}
\]
以及每个非零实对称 \(C\)，证书证明
\[
\boxed{D^2H(K_\varepsilon)[C,C]<-\|C\|_F^2,}\qquad
\boxed{D^2H(K_\varepsilon)[i\widehat B,i\widehat B]>\frac1{40}.} \tag{C4}
\]
所以这不是把前一轮的单个中心再算一次，而是明确参数区间内的一整族实／复
曲率分离。对每个这样的中心，连续性还给出一个实对称开邻域，其中实限制严格凹，
而固定纯虚方向保持正曲率；没有把不同中心的邻域半径声称为统一常数。

### 定理 D：独立实块的混合四阶也不能产生正曲率

对任意两个内部实核 \(K_1,K_2\)，在 \(K_1\oplus K_2\) 处，令
\(E(X)=\left(\begin{smallmatrix}0&X\\X^T&0\end{smallmatrix}\right)\)，则
\[
D^4H[E(X),E(X),E(Y),E(Y)]
=-12\,\mathbb E\!\left[\left(\mathrm{tr}(R_S X T_UY^T)\right)^2\right]\le0. \tag{D1}
\]
这里 \(S,U\) 独立服从两个块的完整事件分布，\(R_S=(K_1-I_{S^c})^{-1}\)、
\(T_U=(K_2-I_{U^c})^{-1}\)。若两个块均连通，且 \(X,Y\ne0\)，不等号严格。
这关闭了“纯四阶径向负，但混合四阶或许正”的一个具体漏洞；不是全耦合区域的凹性定理。

## 2. 预备：边界处的概率多项式，而不是对零概率求对数

固定完整事件空间 \(2^{[n]}\)，自然对数熵为
\[
p_S(K)=(-1)^{|S^c|}\det(K-I_{S^c}),\qquad
H(K)=-\sum_Sp_S(K)\log p_S(K). \tag{1}
\]
在投影处可以有零概率。下文只在投影处对**概率多项式**求导；所有熵求导都在
\(0<K<I\) 内进行，再取渐近极限。原对象与 Lyons [1] PDF 第 2 页 §2.1、
第 5 页 (2.11)--(2.14) 一致；本轮再次直接打开并视觉核对这两页。

**引理 1（实主子式映射的微分核）。** 对任意实对称矩阵 \(M\)，不要求可逆或正定，
若实对称 \(C\) 使全部主子式的一阶导数为零，则 \(C\) 在 \(M\) 每个坐标连通块内部为零；
反过来也成立。等价地，全部 \(Dp_S(M)[C]\) 为零恰好刻画同一空间。

证明沿用主报告的最短路径论证，但注意它不依赖内部性。一阶主子式给 \(C_{ii}=0\)，
二阶给 \(M_{ij}C_{ij}=0\)。沿长度 \(d\) 的最短路径，对诱导主子式求导，唯一尚未
消去的项是
\[
2(-1)^d C_{i_0i_d}\prod_{j=0}^{d-1}M_{i_ji_{j+1}}.
\]
于是端点项也为零。逐图距离归纳完成。反向由块对角行列式展开得到。
完整事件与包含概率之间的 Möbius 变换是可逆线性变换，所以微分核相同。证毕。

这个结论在内部与实 Fisher 零空间刻画相容；已有 \(L\)-参数结果见主报告引用的 [2]。
这里使用的边界版本已给出完整直接证明，不靠把 \(P\) 当成可逆的 \(L\)-核。

## 3. 每个事件的消失阶与严格内向扰动

设 \(r=\mathrm{rank}P\)、\(k=|S|\)，取 \(F\) 为 \(\mathrm{ran}P\) 的实正交列基，
并令
\[
d_S=\mathrm{rank}F_S,\qquad w_S=r+k-2d_S. \tag{2}
\]
则
\[
\ker(P-I_{S^c})
=(\mathrm{ran}P\cap\mathbb R^{S^c})
\oplus(\ker P\cap\mathbb R^S), \tag{3}
\]
两个维数分别为 \(r-d_S\)、\(k-d_S\)，故余维亏损为 \(w_S\)。

证明：写 \(x=u+v\)，\(u=Px,v=Qx\)，则
\((P-I_{S^c})x=(u_S,-v_{S^c})\)。零条件恰好是 (3)。

取满足 (A1) 的 \(N\)。其在 (3) 上的压缩，按上述两个子空间分块为
\(\left(\begin{smallmatrix}N_{--}&N_{-+}\\N_{+-}&N_{++}\end{smallmatrix}\right)\)，
其中 \(N_{--}\prec0,N_{++}\succ0\)。Schur 补证明它可逆，且正负惯性符合这两个维数。
因此沿 \(K_\varepsilon=P+\varepsilon N\)，
\[
p_S(K_\varepsilon)=\varepsilon^{w_S}
\bigl(a_S(N)+O(\varepsilon)\bigr),\qquad a_S(N)>0. \tag{4}
\]
正性也可由严格可行性及首项不为零得到。严格可行性来自 \(K_\varepsilon\) 与
\(I-K_\varepsilon\) 在 \(P\oplus Q\) 分解下的 Schur 补。

这里 \(w_S=0\) 恰为投影分布的正概率基事件；\(w_S=1\) 只可能在
\(k=r-1\) 或 \(r+1\)；若 \(k=r\) 而 \(p_S(P)=0\)，则 \(w_S\ge2\)。
概率的 \(j\) 阶方向导数至少含 \(\varepsilon^{\max(w_S-j,0)}\)，这是行列式余维亏损的直接结果。
这些陈述及余项在 (A1) 内的紧集上统一成立。

对于 full-spark 投影（所有 \(r\) 阶主子式正）和各向同性 \(N=I-2P\)，公式简化为
\[
w_S=|k-r|,\qquad
a_S=\begin{cases}\det P_S,&k\le r,\\\det Q_{S^c},&k\ge r.\end{cases} \tag{5}
\]
定理 C 的证书逐项核对了这两个表达。

## 4. Hessian 的三个主阶：完整渐近公式

任一实对称方向正交分解为
\[
C=D+E,\quad D=PCP+QCQ,\quad E=PCQ+QCP. \tag{6}
\]
在 \(P\oplus Q\) 基下写
\(D=\mathrm{diag}(U,V)\)、\(E=\left(\begin{smallmatrix}0&W\\W^T&0\end{smallmatrix}\right)\)。
记 \(L=\log(1/\varepsilon)\)。

**引理 2（含算子范数余项的展开）。** 存在仅依赖 \(P,N\) 的有限常数 \(M\)，使
\[
\boxed{
D^2H(P+\varepsilon N)[C,C]
=-\frac{G_{P,N}(D)}{\varepsilon}
+L\{\Lambda_P(D)-T_P(E)\}
+O_{P,N}(\|C\|_F^2),} \tag{7}
\]
其中
\[
G_{P,N}(D)=\sum_{w_S=1}\frac{\alpha_S(D)^2}{a_S(N)},\qquad
\alpha_S(D)=Dp_S(P)[D], \tag{8}
\]
\[
T_P(E)=\sum_{p_S(P)>0}\frac{(Dp_S(P)[E])^2}{p_S(P)},\qquad
\Lambda_P(D)=\sum_S w_S D^2p_S(P)[D,D]. \tag{9}
\]
余项绝对值不超过 \(M\|C\|_F^2\)，且在严格内向方向的紧集上可统一选择 \(M\)。

### 4.1 Fisher 主阶

内部通式为
\[
H''=-\sum_S\frac{(p'_S)^2}{p_S}-\sum_Sp''_S\log p_S. \tag{10}
\]
由 (3)，\(w_S=1\) 的零向量完全位于 \(P\) 或 \(Q\) 子空间，所以
\(Dp_S(P)[E]=0\)。由 (4) 及消失阶：\(w=1\) 的 Fisher 项恰有
\(\alpha_S(D)^2/(\varepsilon a_S)\) 主部；\(w=0\) 的项有界；\(w\ge2\) 的项为
\(O(\varepsilon^{w-2})\)，也有界。因此 Fisher 总和为 \(G(D)/\varepsilon+O(\|C\|^2)\)。

### 4.2 对数主阶与为什么没有 \(D,E\) 混合对数项

由 (4)，第二项为
\[
L\sum_Sw_S D^2p_S(P)[C,C]+O(\|C\|^2). \tag{11}
\]
先用粒子数生成函数计算 \(|k-r|\) 加权的和。在 \(P\oplus Q\) 分解下，
\[
\begin{aligned}
\sum_k\left.\frac{d^2}{dt^2}\Pr_{P+tC}(|X|=k)\right|_0 z^k
=(z-1)^2\{&z^{r-2}[(\mathrm{tr}U)^2-\mathrm{tr}U^2]\\
&+2z^{r-1}[\mathrm{tr}U\mathrm{tr}V-\|W\|_F^2]\\
&+z^r[(\mathrm{tr}V)^2-\mathrm{tr}V^2]\}.
\end{aligned} \tag{12}
\]
若出现 \(r=1\) 的形式负幂，其系数恒为零，先删去该项。乘以 \(|k-r|\) 求和给
\[
4\mathrm{tr}U\mathrm{tr}V-4\|W\|_F^2. \tag{13}
\]

\(w_S-|k-r|\) 带来的额外贡献，只可能来自 \(k=r,w_S=2\) 的零基事件：
其余更高亏损事件的二阶概率导数为零。此时 (3) 恰有一个 \(P\) 零向量和一个
\(Q\) 零向量。在这两个向量上，\(D\) 是对角块，\(E\) 是非对角块；二阶首项是
其 \(2\times2\) 行列式，只有“两个对角元之积减非对角元平方”，没有 \(DE\) 项。
故 (11) 的对数系数在 \(\!D,E\) 之间完全分裂。

### 4.3 纯横向项等于负的边界 Fisher

取实正交旋转的投影曲线 \(P(t)\)，使 \(P(0)=P,P'(0)=E\)。在固定坐标外积基中，
写实 Slater 振幅为 \(u_S(t)\)，于是 \(p_S(P(t))=u_S(t)^2\)（\(|S|=r\)）。
可选择水平正交列基，使
\[
\sum_{|S|=r}(u'_S(0))^2=\|W\|_F^2. \tag{14}
\]
这是外积微分的正交性：把每个占据向量替换为一个空向量的项两两正交，系数为 \(W\) 的条目。

对零基事件，\(Dp_S(P)=0\)，所以仿射线 \(P+tE\) 与投影曲线的二阶概率导数相同，
等于 \(2(u'_S(0))^2\)。亏损至少四的基，其行列式振幅的一阶导数也为零。
把这些零基事件对 (13) 的修正加回，纯 \(E\) 对数系数为
\[
-4\|W\|_F^2+4\sum_{u_S(0)=0}(u'_S(0))^2
=-4\sum_{u_S(0)\ne0}(u'_S(0))^2=-T_P(E). \tag{15}
\]
这证明 (7)。特别地，full-spark 时没有零基修正，
\[
\Lambda_P(D)=4\mathrm{tr}U\mathrm{tr}V,
\qquad T_P(E)=4\|W\|_F^2=2\|E\|_F^2. \tag{16}
\]

## 5. 正定性关闭实方向的全部量词

若 \(G(D)=0\)，所有 \(w=1\) 的一阶概率导数为零。沿与 \(P\) 对易的 \(D\)，
粒子数 \(r-1,r+1\) 的一阶总质量分别为
\(-\mathrm{tr}U,\mathrm{tr}V\)，所以两迹为零。
在正概率 \(r\)-基事件上，内部扰动满足
\[
Dp_S(P)[D]=(\mathrm{tr}U-\mathrm{tr}V)p_S(P). \tag{17}
\]
这是固定特征子空间的 Bernoulli 占据展开；改变一个占据与一个空模式至少是二阶。
而 \(w\ge2\) 的事件一阶导数本来就为零。故全部事件的一阶导数为零。
由引理 1 与 \(P\) 连通，\(D=0\)。因此 \(G\) 在内部方向空间正定。

同理，纯 \(E\) 的一阶事件导数只可能出现在正概率基事件上。
\(T(E)=0\) 于是使全部一阶导数为零，引理 1 给 \(E=0\)。故 \(T\) 在横向空间正定。

有限维性给常数 \(g,\tau>0\)、\(b,M<\infty\)，使
\[
H''\le[-g/\varepsilon+bL+M]\|D\|_F^2+[-\tau L+M]\|E\|_F^2. \tag{18}
\]
\(1/\varepsilon\) 压过 \(L\)，而 \(L\to\infty\)，所以两个系数最终均严格负，证明定理 A。
对于紧的内向方向集，\(a_S(N)\) 有统一正下界，系数和余项一致有界，\(G\) 有统一正定余量；
\(T\) 本身不依赖 \(N\)。这也证明紧集一致版本。

推论 B：按 \(P\) 的坐标连通分量分块，各向同性正则化保留这些独立块。
每个非平凡块应用定理 A；孤立的 0/1 坐标是参数 \(\varepsilon\) 或 \(1-\varepsilon\) 的
Bernoulli 熵，其二阶严格负。不同坐标块之间的方向在中心满足 \(p'=0,H''=0\)，
与块内方向的混合 Hessian 也为零（对所有块对角核，跨块一阶导数恒为零，再作块内求导）。
所以 Hessian 的零空间恰为跨坐标块方向。证毕。

## 6. 二、三、四阶的渐近分离，不只是一个 Hessian 数值

对于非零实内部方向 \(D\)，设 \(\alpha_S=Dp_S(P)[D]\)。上面的消失阶与主报告的完整
三、四阶通式给出
\[
\begin{aligned}
\varepsilon H''[D,D]&\longrightarrow-\sum_{w=1}\alpha_S^2/a_S<0,\\
\varepsilon^2 H'''[D,D,D]&\longrightarrow\sum_{w=1}\alpha_S^3/a_S^2,\\
\varepsilon^3 H''''[D,D,D,D]&\longrightarrow-2\sum_{w=1}\alpha_S^4/a_S^3<0.
\end{aligned} \tag{19}
\]
第二行可随方向变号；它没有共轭奇偶消失。严格号使用 \(P\) 连通。

相反，若 \(A=iB\ne0\)、\(B\) 实反对称且 \([B,P]=0\)，则沿任何实
\(P+\varepsilon N\) 的仿射纯虚方向，所有事件的一、三阶导数恒为零。
令 \(b_S=D^2p_S(P)[A,A]\)。有
\[
\sum_S w_S b_S=0. \tag{20}
\]
证明：内部块的迹为零，使 (13) 为零；亏损二的零基上，\(A\) 在一维 \(P\)、\(Q\)
零向量上的对角元为零，跨块也为零，故该事件的二阶导数为零。
所以
\[
\begin{aligned}
H''[A,A]&\longrightarrow \ell(P,N,B):=-\sum_S b_S\log a_S(N),\\
H'''[A,A,A]&=0,\\
\varepsilon^2H''''[A,A,A,A]&\longrightarrow-3\sum_{w=2}b_S^2/a_S(N)<0.
\end{aligned} \tag{21}
\]
最后的严格性可由非零内部反对称块的 \(r-2\) 或 \(r+2\) 粒子数二阶系数推出：
它等于相应 Hermitian 块的 \(-\mathrm{tr}A_{\rm block}^2\ne0\)，故至少一个亏损二事件的
\(b_S\ne0\)。

**机制至此明确：** 实内部方向在二阶就承担严格的 \(1/\varepsilon\) Fisher 损失；
纯虚内部方向的这个损失完全消失，对数发散也取消，留下一个可以为正的有限对数和。
负的概率曲率平方项到四阶才以 \(1/\varepsilon^2\) 回来。这不是把原反例的正值归因于
一个四次 Pfaffian 项。

## 7. 显式五维族的全区间有理认证

[`scripts/certify_projection_boundary.py`](scripts/certify_projection_boundary.py)
仅依赖同目录的 `certify.py` 中的有理矩阵与区间基础函数；
不读取旧 JSON，不调用数值特征分解，不扫描小 \(\varepsilon\)。

它逐项证明 \(P^2=P=P^T,\mathrm{tr}P=2\)、\(\widehat B^T=-\widehat B\)、
\([P,\widehat B]=0\)，以及 (5) 的 32 个 \(a_S>0\)。所有多项式系数由
**主子式的置换展开 + 完整事件 Möbius 反演**得到。

作为第二算法，程序在 \(\varepsilon=1/4\) 对全部 32 个完整事件，用直接有理高斯消元和
逆矩阵迹公式核对概率、全部 15 个一阶导数、全部 \(15\times15\) 混合二阶导数与纯虚二阶导数。
这检查多项式构造与事件符号，不是用一个点验证全参数不等式。

另用一个固定四维测试核核查非 full-spark 部分：\(V=\left(\begin{smallmatrix}1&0\\1&0\\1&1\\1&-1\end{smallmatrix}\right)\)
的投影连通但有零二阶主子式；取 \(N=I-2P+(PC_0Q+QC_0P)/3\)、\(C_0=\mathrm{diag}(1,0,0,0)\)。
程序精确核对所有事件消失阶、全部混合／横向对数系数，并对十维二次型 \(G+T\) 作有理 LDL 正定认证。
这是检验一般证明中“零基修正”与非各向同性扰动的单个代数例，不代替一般证明。

### 7.1 纯虚极限严格为正

80 项有理 atanh 展开、几何余项和向外取整给
\[
\frac1{30}<\ell(P,I-2P,\widehat B)<\frac1{29},\qquad
\ell=0.034036295104490266247982151599\ldots. \tag{22}
\]
承重的是有理区间，不是这个小数。\(\sum_Sw_Sb_S=0\) 是精确有理恒等式。

### 7.2 不用“足够小”：整个 \((0,2^{-32}]\) 的误差界

写 \(p_S(\varepsilon)=\varepsilon^{w_S}f_S(\varepsilon)\)，\(f_S(0)=a_S\)，
\(R_S=\sum_{j\ge1}|[\varepsilon^j]f_S|\)。这里 \(R_S\) 是多项式系数界，不是前文的事件逆矩阵。
在 \(\varepsilon\le\min_S a_S/(2R_S)\) 时，
\[
a_S/2\le f_S(\varepsilon)\le3a_S/2,
\quad |\log f_S-\log a_S|\le2\varepsilon R_S/a_S. \tag{23}
\]
对一、二阶导数的多项式系数作有限绝对值求和，得到有理常数
\(R_1,A_0,B_0,A_i,B_i,q_0\)，使
\[
F(C)\ge\frac{G(D)}{4\varepsilon}-\varepsilon R_1\|C\|_F^2,\tag{24}
\]
\[
Q(C)\le L\{4\mathrm{tr}U\mathrm{tr}V-4\|W\|_F^2\}
+\{q_0+\varepsilon(A_0L+B_0)\}\|C\|_F^2,\tag{25}
\]
\[
|H''[i\widehat B,i\widehat B]-\ell|\le\varepsilon(A_iL+B_i).\tag{26}
\]
(24) 使用 \((u+v)^2\ge u^2/2-v^2\)，不是错误地丢弃实 Fisher 项。
(25) 中常数二次型为 \(-\sum_Sp''_S(P)\log a_S\)，\(q_0\) 是其 15 维坐标矩阵
严格区间的绝对行和上界。全部 225 个对数发散矩阵条目还逐项核对 (16)。

为免隐藏特征值下界，利用该 \(P\) 所有非对角元非零，令 \(m=\min_{i<j}|P_{ij}|\)。
引理 1 的一、二阶主子式与 Cauchy--Schwarz 给显式
\[
G(D)\ge g\|D\|_F^2,\qquad
g=\left(4n^2+\frac{9n^2(n-1)}{m^2}\right)^{-1}. \tag{27}
\]
具体证明：\(\sum_{w=1}a_S=n\)，所以亏损一事件导数的绝对值总和不超过 \(\sqrt{nG}\)。
(17) 及两个粒子数导数使全部事件导数的绝对值总和不超过 \(2\sqrt{nG}\)。
包含概率的一、二阶主子式于是给 \(|D_{ii}|\le2\sqrt{nG}\)、
\(|D_{ij}|\le3\sqrt{nG}/m\)。平方求和即 (27)。

在本构造中，证书的量级显示为
\[
g\approx7.4911385\times10^{-6},\quad q_0<21,\quad
R_1\approx275127.542,\quad A_0\approx1247.032,\quad B_0\approx156439.272.
\]
\(A_i\approx92.3041,B_i\approx10864.466\)。JSON 保存全部精确有理常数。

取 \(\varepsilon_*=2^{-32}\)。用 \(1/2<\log2<1\) 及 \(\varepsilon\log(1/\varepsilon)\)
在 \((0,1/4]\) 递增，程序严格检查
\[
\delta_*:=\varepsilon_*\{32A_0+B_0+R_1\}<1,
\quad 32>q_0+\delta_*+1,\tag{28}
\]
\[
\varepsilon_*(6\cdot32+q_0+\delta_*+1)<g/4,
\qquad \varepsilon_*(32A_i+B_i)<1/120. \tag{29}
\]
实际 \(\delta_*\) 小于 \(0.000110\)，纯虚误差小于 \(0.000003218\)。
由 \(4\mathrm{tr}U\mathrm{tr}V\le6\|D\|_F^2\)、
\(4\|W\|_F^2=2\|E\|_F^2\)，(24)--(29) 对整个参数区间给出 (C4)。

这是一张**全参数系数与余项证书**。\(2^{-32}\) 由固定构造的有限误差常数选择，
不伪称事前登记的神奇步长，也没有靠在小 \(\varepsilon\) 处浮点求逆判号。

## 8. 独立块的混合四阶：完整正化

令 \(S\sim\mathrm{DPP}(K_1),U\sim\mathrm{DPP}(K_2)\) 独立，事件逆矩阵为 \(R_S,T_U\)，
并记
\[
q(X,Y)=\mathrm{tr}(R_S X T_UY^T).
\]
块边缘固定，Schur 行列式给
\[
p_{S,U}(X)=p_Sp_U\{1-q(X,X)+O(\|X\|^4)\}.
\]
因此互信息的四次齐次首项为
\[
\Phi(X)=\tfrac12\mathbb E q(X,X)^2. \tag{30}
\]
不能仅凭“平方和”断言它凸：被平方的二次型逐事件可以不定。这里还需如下恒等式。

完整事件行列式的归一化 \(\sum_Sp_S(K)=1\) 作为**一般矩阵条目**的多项式恒等式也成立。
对任意两条目求二阶导，在实对称基点得到
\[
\mathbb E[R_{ij}R_{kl}]=\mathbb E[R_{il}R_{kj}]. \tag{31}
\]
结合对称性，这个四指标 Fisher 张量完全对称。两个块分别应用 (31)，重排指标得
\[
\mathbb E[q(X,X)q(Y,Y)]=\mathbb E q(X,Y)^2. \tag{32}
\]
直接微分 (30)，
\[
D^2\Phi(X)[Y,Y]
=2\mathbb E[q(X,X)q(Y,Y)]+4\mathbb E q(X,Y)^2
=6\mathbb E q(X,Y)^2\ge0. \tag{33}
\]
极化给 (D1)。若两块连通，实 Fisher 正定意味着各块事件逆矩阵分别张成全部实对称矩阵。
若右侧为零，则对所有实 \(u,v\)，\((u^TXv)(u^TYv)=0\)。这是两个实多项式乘积恒为零，
所以 \(X=0\) 或 \(Y=0\)，得到严格号。证毕。

程序另用两个固定有理 \(2\times2\) 内部核及两个固定非零交叉矩阵，精确检查了 (31)--(32)。
这些小例只是排错，任意维数结论由上述证明承担。

## 9. 本轮尝试、复现与尚缺的部分

本轮先读取了远端分支、上一轮报告与 Canonical Issue 的实际评论；分支头当时仍为
`a296f7c...`。旧数学没有改动。数值尝试按 `continuation_plan.md` 保存的顺序进行：
16 个固定局部相位实化（首相位 0，其余在 \(\{0,\pi/4\}\)）没有给出正候选；标准实化那项
的约 \(10^{-14}\) 数值残差由已知恒等式判为零，不能称命中。其余 15 项仅有浮点负值，
不构成全相位定理。另在旧 \(K_*\) 的 26 个坐标、常向量、二坐标和／差向量上检查了
351 个混合 rank-one 对，未找到反例；这个更强归纳捷径仍未被证明。

随后按事前固定的唯一构造，取旧整数核第 0、1 列形成 (C1)，没有试别的投影列或随机中心。
最初的边界极限诊断用符号有理行列式与高精度对数；承重程序重新只用标准库有理数及严格对数区间。
未运行随机扫描、重型搜索或远端计算。不同算法、自检重跑、以及本段说明都不构成独立验证。

复现：
```sh
python scripts/certify_projection_boundary.py --output certificates/projection_boundary_certificate.json
```
该命令重建 `certificates/projection_boundary_certificate.json`。背景证书不被改写。

**现在的缺口比上一轮窄，但没有消失。** 一般实核仍需要证明或反驳 \(Q_K(C)\le F_K(C)\)。
定理 A 关闭固定连通投影附近、内向余量与接近速度可比较的全实方向区域；推论 B 关闭所有
固定投影的充分小各向同性正则化；定理 D 关闭独立实块耦合的混合四阶正曲率尝试。
这些结果没有覆盖一般内部核、非平衡切向边界路径，或不连通极限上同时改变块间耦合与噪声的
多尺度路径。一般正交混合、其他块嵌入和条件标签熵转移也仍未建立有效转移定理。

不将这些剩余集合说成“必有反例”，不将局部／锥形结论说成全局定理，也不主张优先权或新颖性认证。

### 来源

[1] Russell Lyons, *Determinantal Probability: Basic Properties and Conjectures*,
arXiv:1406.2707v1, 10 June 2014；PDF 第 2 页 §2.1 及第 5 页 (2.11)--(2.14)。
`https://arxiv.org/pdf/1406.2707`

[2] Victor-Emmanuel Brunel, Ankur Moitra, Philippe Rigollet, John Urschel,
*Rates of estimation for determinantal point processes*, PMLR 65 (2017), Theorem 2。
已有 \(L\)-参数 Fisher 零空间结果；本篇边界主子式微分论证不假设 \(L\) 可逆。
`https://proceedings.mlr.press/v65/brunel17a/brunel17a.pdf`

原反例与前一轮机制报告现随本公开包提供。本篇的新渐近展开、混合四阶论证与
全区间证书的核验范围见 [verification_summary.md](verification_summary.md)。
文献定向检索不承担实版本当前状态或优先权认证。
