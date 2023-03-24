# -*- coding: utf-8 -*-
# @Author: Zhixuan Li
# @Date:   2023-02-10 15:09:12
# @Last Modified by:   Zhixuan Li
# @Last Modified time: 2023-03-09 20:57:30
from main import *

def gain_trigdata(fitsname):
    """该函数用于得到处理后的分能段光变数据和提取出的疑似爆发数据
    parameters:
        fitsname: str
            数据文件名称(包含后缀)
    returns: 
        trig_data: list 
            一个包含两个字典元素的列表, 两个字典元素分别是触发数据的时间和强度值
        data: list
            一个包含两个字典元素的列表, 两个字典元素分别是输入数据的时间和强度值
            
    """
    # fitsname='20210125_1033_unpack-202101231033_S0.fits'
    datasum=E_split(fitsname)
    E_span={'5-15keV':0,'15-50keV':1,'50-300keV':2,'300keV_&_above':3}
    nburst={'5-15keV':[],'15-50keV':[],'50-300keV':[],'300keV_&_above':[]}
    tburst={'5-15keV':[],'15-50keV':[],'50-300keV':[],'300keV_&_above':[]}
    n={'5-15keV':[],'15-50keV':[],'50-300keV':[],'300keV_&_above':[]}
    t={'5-15keV':[],'15-50keV':[],'50-300keV':[],'300keV_&_above':[]}
    for keys,values in E_span.items(): #四个能段分别得出爆发区域
        nburst[keys],tburst[keys],n[keys],t[keys]=trigger(datasum[values])
    trig_data=[tburst,nburst]
    data=[t,n]
    
    
    return trig_data,data

trig_data,data=gain_trigdata('20210125_1033_unpack-202101231033_S0.fits')
