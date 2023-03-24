# -*- coding: utf-8 -*-
# @Author: lizh
# @Contributor: Zhixuan Li
# @Date:   2022-11-16 11:08:24
# @Last Modified by:   Zhixuan Li
# @Last Modified time: 2023-02-10 17:27:58
import numpy as np
def burst_specific(residual,sigma,aver):
        """该函数用于在拟合后的残差里面,利用3sigma找出burst的位置
        parameters:
                residual: ndarray
                拟合后的残差
                sigma: number
                残差的标准差
                aver: number
                数据均值
        Returns:  out: tuple
                index: ndarray
                返回一个元组,第一个元素是burst的值data_burst,第二个元素是burst的在residual中的位置index
                """
        burst_index = []
        index=[]
        for i, k in enumerate(residual):
                if k > 3*sigma+aver:
                        burst_index.append(i)
        if len(burst_index)>0:
                diff_i=np.diff(burst_index)
                g=np.where(diff_i>1)[0]
                index=np.split(burst_index,list(g+1))
                #index表征爆发在residual的序号
                around_burst=20
                
                for i in range(len(index)):
                        if min(index[i])-around_burst<0:
                                index[i]=np.arange(0,max(index[i])+around_burst)
                        elif max(index[i])+around_burst>len(residual)-1:
                                #这里是index的比较，所以要在len(t)后减去1
                                index[i]=np.arange(min(index[i])+around_burst,len(residual)-1)
                        else:
                                index[i]=np.arange(min(index[i])-around_burst,max(index[i])+around_burst)                
        #        data_burst = [residual[index[i]] for i in range(len(index))]
        return index
