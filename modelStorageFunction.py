#!/usr/bin/env python
# coding: utf-8

# # 木村の貯留関数法

# ## 流出高の計算

# 木村の貯留関数法では、連続式、運動方程式を次のようにモデル化する。
# $$
# \begin{align}
# \dfrac{ds}{dt} &= r_e - q_l \\
# s &= k {q_l}^p
# \end{align}
# $$
# 
# $s$：貯留高[mm/hr]、$q_l$：遅れ時間 $T_l$[hr]を考慮した流出量[mm/hr]、$r_e$：有効降雨強度[mm/h]、$k,p$：パラメータ
# 
# 流域を流出域と浸透域に分けてそれぞれに本式を適用する。
# 流出域と浸透域の面積率は$f_1$、$1-f_1$とする。なお、$f_1$はパラメータとする。
# 
# 以降、流出域は添字1と浸透域は添字2を付ける。
# 
# 流出域、浸透域の有効降雨強度$r_e$は次のように与える。
# 
# $$
# \begin{align}
# {r_e}_1 &= r \\
# {r_e}_2 &= \begin{cases}
#     0  & (R_{sa} <= \sum r) \\
#     r  & (R_{sa} > \sum r ) \\
#   \end{cases}
# \end{align}
# $$
# 
# $r$：流域平均降雨強度[mm/hr]，$R_{sa}$：飽和雨量[mm]
# 
# これらを用いて、流出域、浸透域の流出高、${q_l}_1,{q_l}_2$を計算する。

# ###  離散化

# ### ソースコード

# In[1]:


import numpy as np
import pandas as pd
from scipy import interpolate


# In[4]:


def calql(RainT, Rain, dt, k, p, Rsa):
    # 降雨は慣習的に階段上に補間する。
    funcRain = interpolate.interp1d(RainT, Rain, kind='previous')
    
    nmax = int(RainT[-1]/dt)
    ql1 = np.zeros(nmax)
    ql2 = np.zeros(nmax)
    sumR = float(0)
    isSaturated = False
    for n in range(nmax-1):
    # 流出域
        re1 = funcRain(dt*n)
        ql1[n+1] = (ql1[n]**p + dt/k*(re1-ql1[n]))**(float(1)/p)
    # 浸透域
        sumR += dt*re1
        if isSaturated:
            re2 = re1
        else: 
            if Rsa < sumR:
                re2 = (sumR-Rsa)/dt
                isSaturated = True
            else:
                re2 = float(0)
                
        ql2[n+1] = (ql2[n]**p + dt/k*(re2-ql2[n]))**(float(1)/p)
    
    return np.arange(nmax)*dt, ql1, ql2


# ## 流出量の計算

# 遅れ時間$T_l$[hr]を考慮した流域からの流出量$Q_l$[$\rm{m^3/s}$]は流出高${q_l}_1,{q_l}_2$を用いて以下のように示す。
# $$
# \begin{align}
# Q_l = \dfrac{1}{3.6}f_1 \cdot A \cdot {q_l}_1 + \dfrac{1}{3.6}(1-f_1) \cdot A \cdot {q_l}_2 + Q_b
# \end{align}
# $$
# 
# $A$：流域面積[$\rm{km}^2$]、$Q_b$：基底流量[$\rm{m^3/s}$]、$f_1$：流出域の面積率
# 
# $f_1$はパラメータである。

# ### ソースコード

# In[3]:


def calQall(ql1, ql2, f1, A, Qb):
    return (f1*ql1+(1.0-f1)*ql2)*A/3.6 + Qb


# ## 遅れ時間のシフト

# 流出量$Q$[$\rm{m^3/s}$]は次式で示される。
# 
# $$
# \begin{align}
# Q_l(t) = Q(t+T_l)
# \end{align}
# $$
# 
# $T_l$[hr]：遅れ時間
# 
# $T_l$[hr]もパラメータとなる。
# 
# コードは時間をシフトさせるだけなので省略する。
