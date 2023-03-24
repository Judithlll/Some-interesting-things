# -*- coding: utf-8 -*-
# @Author: Zhixuan Li
# @Date:   2022-11-17 23:44:46
# @Last Modified by:   Zhixuan Li
# @Last Modified time: 2023-03-09 15:22:23

import matplotlib.pyplot as plt
import numpy as np
from astropy.stats import histogram
import copy

def burst_fig(data_bins,E_span):
    """用于画筛选出来的疑似爆发的图片并保存
        parameters:
                data_bins: ndarray
                暴区域分段数据, 一般包括若干个array
                E_span: str
                用于输入能段名称
        returns: 
                直接保存图片到guest-lzx/burstfigure路径下, 包括数据段中的每个疑似爆发
        """
    for i in range(len(data_bins)):
        plt.figure(figsize=(25,20))
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        figname=E_span+'_'+'burst'+str(i)     #burst name后续需要精确命名，应包括burst到达时间和其他基本信息
        bins_bayesian=histogram(data_bins[i],bins='blocks')
        bins_ori=1
        h_ori=histogram(data_bins[i],bins=bins_ori)
        
        while np.mean(h_ori[0])>np.mean(bins_bayesian[0]):
            bins_ori+=1                                    #自动选择bins
            h_ori=histogram(data_bins[i],bins=bins_ori)

        plt.hist(data_bins[i],bins=bins_ori,color='blue',alpha=0.4,label='hist')    
        plt.hist(data_bins[i],histtype='step',color='red',bins=bins_bayesian[1],label='beyesian blocks')
        plt.legend(fontsize=30,loc='best')
        plt.savefig('burstfigure/'+figname+'.png')
        plt.close()
    print('成功输出分能段trigger图至burstfigure文件夹下')
    return 

def fit(t,n,count):
    """用于生成基线数据序列和相关参数
        parameter:
                t: ndarray
                对从fits中提取出来的原始时间序列取bins以后得出的初始时间序列
                n: ndarray
                对从fits中提取出来的原始时间序列取bins以后得出的初始到达光子数
                count: int
                迭代次数, 迭代次数越多, sigma距离基线越近
        return: 
                sigma: float
                基线部分的标准差
                tfix: ndarray
                经迭代过后得出的基线部分时间序列
                nfix: ndarray
                经迭代过后得出的基线部分高度序列
        """
    tfix=copy.deepcopy(t)
    nfix=copy.deepcopy(n)
    for i in range(count):
        sigma=np.nanstd(nfix)
        aver=np.nanmean(nfix)
        for i in range(len(t)):
            if nfix[i]-aver>3*sigma:
                tfix[i]=np.NaN
                nfix[i]=np.NaN
    return sigma,aver,tfix,nfix

def light_curve(t,n,count,fitsname,E_span):
    """用于画出光变曲线
        parameter: 
                t: ndarray
                对从fits中提取出来的原始时间序列取bins以后得出的初始时间序列
                n: ndarray
                对从fits中提取出来的原始时间序列取bins以后得出的初始到达光子数
                count: int
                迭代次数, 迭代次数越多, sigma距离基线越近
                fitsname: str
                用于生成光变曲线的fits文件的名称
                E_span: str
                用于输入能段名称
        return: 
                sigma,aver: float
                分别为基线部分的标准差和平均值
            """
    sigma,aver,tfix,nfix=fit(t,n,count)
    sig=np.linspace(sigma,sigma,len(t))
    nr=copy.deepcopy(n)
    nr[nr<3*sigma+aver]=np.NaN

    plt.figure(figsize=(20,15))
    plt.subplot(211)
    plt.title(fitsname,fontsize=25,color='black')
    plt.plot(t,n,'blue',label='light_curve')
    plt.plot(tfix,3*sig+aver,'red',linewidth=3,label='3sigma+average')
    plt.plot(tfix,-3*sig-aver,'red',linewidth=3,label='-3sigma-average')
    plt.plot(t,nr,'yellow',linewidth=3,label='burst_candidate')
    plt.legend(fontsize=15)
    plt.savefig('lightcurve/'+fitsname+'_'+E_span+'.png')
    plt.close()
    print('lightcurve已成功生成至lightcurve文件夹下')
    return sigma , aver 
