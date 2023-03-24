# -*- coding: utf-8 -*-
# @Author: Zhixuan Li
# @Date:   2023-02-10 08:51:16
# @Last Modified by:   Zhixuan Li
# @Last Modified time: 2023-03-09 15:18:50
from ADC2Energy import *
import numpy as np
# from scipy import interpolate
import matplotlib.pyplot as plt
# from astropy.modeling import models, fitting
# from astropy.stats import sigma_clip
#from astropy.stats import histogram
from burst_specific import *
from burst_fig import *


def E_split(fitsname):
    """用于将数据按照能量值分为5-15keV, 15-50keV, 50-300keV, 300keV以上四段
    parameters: 
        fitsname: str
            数据文件名称
    returns: 
        [data_5_15: list,
         data_15_50: list,
         data_50_300: list,
         data_300_: list]
        分别为能量值处于四个区间内的数据,包含在一个列表内
    """
    fitspath='/home/bnugrid/test_Grid02/GRID/find_burst/'+fitsname  
    event=ADC2Energy(fitspath)
    #分能段和频道
    data_5_15  = []
    data_15_50 = []
    data_50_300= []
    data_300_  = []
    for i in range(4):
        for j in range(len(event[str(i)])):
            if 15>event[str(i)][j][2]>5:
                data_5_15.append(event[str(i)][j][0])
            elif event[str(i)][j][2]<50:
                data_15_50.append(event[str(i)][j][0])
            elif event[str(i)][j][2]<300:
                data_50_300.append(event[str(i)][j][0])
            else:
                data_300_.append(event[str(i)][j][0])
    return [data_5_15,data_15_50,data_50_300,data_300_]

def trigger(data):
    """用于数据的处理
    parameters: 
        data: list
        分能段过后的数据, 可直接用直方图分析
    returns: 
        nburst: list
        包含几个触发数据列表纵坐标的列表
        tburst: list
        包含几个触发数据列表横坐标的列表
        n: list
        总数据列表纵坐标值
        t: list
        总数据列表横坐标值
    """
    t_bins=0.5
    t_span=(max(data)-min(data))*50/10**9
    n, bins,patches =plt.hist(data,bins=int(t_span/t_bins)+1)
    plt.close()
    t=bins[0:len(bins)-1]*50/10**9
    #整理数据(用了插值，不太靠谱)
    index=np.ndarray.flatten(np.argwhere(n==0))
    diff_index=np.diff(index)
    pindex=np.ndarray.flatten(np.argwhere(diff_index>3))
    split_index=np.split(index,pindex+1)
    
    if 0<index.size<len(n)/5:       #假设异常段小于每个fits文件的1/5，这个数据还需要确定
        for i in range(len(split_index)):
            if split_index[i][0]==0:
                n[split_index[i][0:-1]]=n[split_index[i][-1]+1]
            elif split_index[i][-1]==len(n)-1:
                n[split_index[i][0:-1]]=n[split_index[i][0]-1]         #直接采用连直线的方式手动插值并排除可能存在于首尾的异常值(0)情况
            else:
                n[split_index[i][0]-1:split_index[i][-1]+1]=np.linspace(n[split_index[i][0]-1],n[split_index[i][-1]+1],split_index[i][-1]-split_index[i][0]+2) 
            
        # n_remove0=np.delete(n,index)
        # t_remove0=np.delete(t,index)
        # f=interpolate.interp1d(t_remove0,n_remove0,kind='cubic')
        # n[index]=f(t[index])

    #拟合基线
    # Order = 1
    # p1 = models.Polynomial1D(Order)
    # pfit = fitting.LinearLSQFitter()                   #认真考虑一下需不需要这一段
    # sigma_clip_fit = fitting.FittingWithOutlierRemoval(pfit, sigma_clip, niter=3, sigma=3.0)
    # model_pI, mask_I = sigma_clip_fit(p1, t, n)
    # n_new  = n - model_pI(n)

    sigma,aver,tfix,nfix=fit(t,n,5)
    burst_index=burst_specific(n,sigma,aver)
    nburst=[]
    tburst=[]
    #提取爆发段
    if len(burst_index)>0:
        for i in range(len(burst_index)):
            if len(burst_index[i])>80:            #需调研gamma暴的持续时间对这个参数以及最开始的bins做调整
                nburst.append(n[burst_index[i]])
                tburst.append(t[burst_index[i]])

    return nburst,tburst,n,t
