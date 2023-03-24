# -*- coding: utf-8 -*-
# @Author: Zhixuan Li
# @Date:   2023-01-18 14:44:00
# @Last Modified by:   Zhixuan Li
# @Last Modified time: 2023-01-18 23:09:08
def ADC2Energy(fitspath):
    import numpy as np
    from astropy.io import fits
    from scipy.interpolate import interp1d
    import json

    def get_low_Channel_factor(ec_coef_low):
        """该函数提取出json文件里的低能段ADC和energy的函数系数a,b,c
            parameters: ec_coef_low: dict
                    表示json文件里的"low"key下的含四个channel的函数系数的字典
            return: low_Channel_factor: dict
                    key代表channel,value代表[a,b,c]的列表
                    """
        low_Channel_factor={}
        for i in range(0,4):
            low_Channel_factor['{}'.format(i)] = [list(ec_coef_low['ch{}'.format(i)].values())[k] for k in [0,2,4]]
        return low_Channel_factor

    def get_high_Channel_factor(ec_coef_high):
        """该函数提取出json文件里的高能段ADC和energy的函数系数a,b,c
            parameters: ec_coef_high: dict
                    表示json文件里的"high"key下的含四个channel的字典
            return: high_Channel_factor: dict
                    key代表channel,value代表[a,b,c]的列表
                    """
        high_Channel_factor={}
        for i in range(0,4):
            high_Channel_factor['{}'.format(i)] = [list(ec_coef_high['ch{}'.format(i)].values())[k] for k in [0,2,4]]
        return high_Channel_factor

    def get_ADC_bound(low_Channel_factor, high_Channel_factor,high_bound,low_bound):
        """该函数解出不同Channel对应的高低能段界限high_boun,low_bound对应ADC的bound
            parameters: 
                    low_Channel_factor: dict
                        通过函数get_low_Channel_factor得到
                    high_Channel_factor: dict
                        通过函数get_high_Channel_factor得到
                    high_bound: number
                        json文件里的high_bound
                    low_bound: number
                        json文件里的low_bound      
            return: 
                    ADC_bound: dict
                    key代表channel,value代表[low_ADC_bound,high_ADC_bound]的列表
                    """
        ## 解出低能段和高能段不同channel对应的ADC的bound
        ## energy的单位是kev
        ADC_bound={}
        for i in range(0,4):
            low_func = low_Channel_factor['{}'.format(i)][0:2]+[low_Channel_factor['{}'.format(i)][2]-low_bound]
            high_func = high_Channel_factor['{}'.format(
                i)][0:2]+[high_Channel_factor['{}'.format(i)][2]-high_bound]
            low_root_interval = np.roots(low_func)
            high_root_interval = np.roots(high_func)
            ## 正常根是绝对值小的那两个，c/a=x1*x2;b/a=x1+x2
            bound=np.sort(abs(np.hstack((low_root_interval,high_root_interval))))[0:2]  
            ADC_bound['{}'.format(i)] = bound
        return ADC_bound

    def get_ADC2Energy_func(i,low_Channel_factor,high_Channel_factor,ADC_bound):
        """该函数求出第i个channel的整个能段的ADC和energy的响应函数
            parameters:
                    i: int
                        i,代表第i个channel
                    low_Channel_factor: dict
                        通过函数get_low_Channel_factor得到
                    high_Channel_factor: dict
                        通过函数get_high_Channel_factor得到
                    ADC_bound: dict
                    通过函数get_ADC_bound得到 
            return:
                    ADC2Energy_func: list
                        代表第i个channel整个能段的ADC和energy的响应函数
                    """
        ADC2Energy_func = {'0': [], '1': [], '2': [], '3': []}
        ADC2Energy_func['{}'.format(i)].append(
            lambda x: x**2*low_Channel_factor['{}'.format(i)][0]+x*low_Channel_factor['{}'.format(i)][1]+low_Channel_factor['{}'.format(i)][2])
        ADC2Energy_func['{}'.format(i)].append(
            lambda x: x**2*high_Channel_factor['{}'.format(i)][0]+x*high_Channel_factor['{}'.format(i)][1]+high_Channel_factor['{}'.format(i)][2])
    
        ##插值得中能段响应曲线
        x1 = np.linspace(ADC_bound['{}'.format(i)][0], ADC_bound['{}'.format(i)][0], 1)
        x2 = np.linspace(ADC_bound['{}'.format(i)][1], ADC_bound['{}'.format(i)][1], 1)
        y1 = ADC2Energy_func['{}'.format(i)][0](x1)
        y2 = ADC2Energy_func['{}'.format(i)][1](x2)
        mid_func = interp1d(np.hstack((x1, x2)),
                            np.hstack((y1, y2)), kind='linear')
        ADC2Energy_func['{}'.format(i)].insert(1, mid_func)

        return ADC2Energy_func['{}'.format(i)]


    def get_Event_Energy(hdul, ec_coef_low, ec_coef_high, low_bound, high_bound):
        """该函数分channel,并求出每个event的能量
            parameters: hdul: fits file
                            通过函数fits.open得到
                        ec_coef_low: dict
                            表示json文件里的"low"key下的含四个channel的函数系数的字典
                        ec_coef_high: dict
                            表示json文件里的"high"key下的含四个channel的函数系数的字典
                        low_bound: number
                            json文件里的low_bound
                        high_bound: number
                            json文件里的high_bound
            return: dict: dict
                    key代表channel,value代表该channel的event的[[uscount,ADC,Energy],[uscount,ADC,Energy],....]的二维列表"""
        low_Channel_factor = get_low_Channel_factor(ec_coef_low)
        high_Channel_factor = get_high_Channel_factor(ec_coef_high)
        ADC_bound = get_ADC_bound(
        low_Channel_factor, high_Channel_factor, low_bound, high_bound)
        Channel = hdul[1].data['Channel']
        usCount = hdul[1].data['usCount']
        ADC = hdul[1].data['ADC']
        ## 合并相同通道
        dict = {'0': [], '1': [], '2': [], '3': []}
        for i,j in enumerate(Channel):
            dict['{}'.format(j)].append([usCount[i], ADC[i]])
        ## 计算能量
        for i in range(4):
            ADC2Energy_func = get_ADC2Energy_func(
                i, low_Channel_factor, high_Channel_factor, ADC_bound)
            for j in range(len(dict['{}'.format(i)])):
                if dict['{}'.format(i)][j][1] < ADC_bound['{}'.format(i)][0]:
                    dict['{}'.format(i)][j].append(
                        ADC2Energy_func[0](dict['{}'.format(i)][j][1]))
                elif dict['{}'.format(i)][j][1] < ADC_bound['{}'.format(i)][1]:
                    dict['{}'.format(i)][j].append(
                        ADC2Energy_func[1](dict['{}'.format(i)][j][1]))
                else:
                    dict['{}'.format(i)][j].append(
                        ADC2Energy_func[2](dict['{}'.format(i)][j][1]))
        return dict

    ## 主程序
    ## 读取json文件的能量和ADC标定系数
    with open('/home/bnugrid/test_Grid02/GRID/source_plan/ec_coef.json',encoding='utf-8') as json_file:
        ec_coef = json.load(json_file)
        ec_coef_low=ec_coef['low']
        ec_coef_high=ec_coef['high']
        low_bound=ec_coef['low_bound']
        high_bound=ec_coef['high_bound']


    ## 读取fits文件
    with fits.open(fitspath) as hdul:
        ## 读取fits文件的ADC和channel,时间数据
        event=get_Event_Energy(hdul, ec_coef_low, ec_coef_high, low_bound, high_bound)
        #最后event这个字典，key代表channel，value是一个二维列表，每个子列表里面是时间，ADC，能量
    return event
