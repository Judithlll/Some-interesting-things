# Lifetime of Ions of Li 

## Cut out the part of that image

![green1](D:\桌面\green1.jpg)

## Thoughts(By Python)

### 1. Extract the data of RGB of this image, and use the formula below to calculate the luminance.

$$
L=R×0.299+G×0.587+B×0.114
$$

### 2. Calculate the average vertically.

### 3. Then, take x as the horizontal axis, take the average luminance as the vertical axis, and make a diagram.

### 4. Through discussion, we(my classmates and I) use the model below to fit the curve.

$$
y = c_1 \exp{-\frac{x}{v\cdot \tau_1}}+c_2 \exp{-\frac{x}{v\cdot \tau_2}}
$$

## Program

```python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit

img=plt.imread("D:/Desktop/green1.jpg")
#Read the pictures

R=img[:,:,0]/255
G=img[:,:,1]/255
B=img[:,:,2]/255
AR=sum(R[0:,:])/78
AG=sum(G[0:,:])/78
AB=sum(B[0:,:])/78
#Read RGB color values and average them
#Zhixuan Li provided this passage

L=AR
for i in range(1183):
    L[i]=AR[i]*0.299*AG[i]*0.587+AB[i]*0.114
#Standard RGB conversion to brightness calculations

q=L[91]
for i in range(1183):
    L[i]=L[i]/q
#normalization

x=np.arange(0.00,0.3500,0.35/1000,"float")
def func(x, a1, b1, a2, b2):
    return a1 * np.exp(-b1 * x) + a2*np.exp(-b2*x)
popt, pcov = curve_fit(func, x, L[91:-92])
print(popt)
#fitting

plt.figure(dpi=600)
plt.plot(x,L[91:-92],'b-')
LF=[func(i, popt[0],popt[1],popt[2],popt[3],) for i in np.arange(0.00,0.3500,0.35/1000,'float')]
plt.plot(x,LF,'g-')
#plot
```

## Results

![image-20220707225143506](C:\Users\LZX\AppData\Roaming\Typora\typora-user-images\image-20220707225143506.png)
$$
[a1,b1,a2,b2]=[0.40887566   ， 4.54058842，  0.59840332 ， 36.54179452]
$$

### Calculate the lifetime:

$$
\begin{equation}
\left\{
\begin{aligned}
\frac{1}{v\cdot \tau_1}&=4.541\\
\frac{1}{v\cdot \tau_2}&=36.542
\end{aligned}
\right.
\end{equation}
$$

$$
\begin{equation}
\left\{
\begin{aligned}
\tau_1&=1.46912\times10^{-8}\\
\tau_2&=1.82565\times10^{-9}
\end{aligned}
\right.
\end{equation}
$$

