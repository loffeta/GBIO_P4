# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 11:24:35 2020

@author: Dany
"""  
import numpy as np 
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
    
def omeg1(f1,f2,m):  
    
    qa = np.linspace(-np.pi-1,np.pi+1,2000)
    exp_i_q_a_demi = np.exp(1j*qa/2)
    exp_i_q_a = np.exp(1j*qa)
    exp_i_q_a_trois_demi = np.exp(3*1j*qa/2)
           
    #y = np.sqrt((2*m))
    
    omega_1 = np.zeros((len(qa),))
    omega_2 = np.zeros((len(qa),))
    
    for i in range(len(qa)):
        
        a_0 = exp_i_q_a_demi[i]*m**2
        a_1 = -(f1+f2)*exp_i_q_a_demi[i]*m*(1+exp_i_q_a[i])
        a_2 = exp_i_q_a_trois_demi[i]*((f1+f2)**2)-(f1+f2*exp_i_q_a[i])*(f2+f1*exp_i_q_a[i])
    
        p_x = np.poly1d([a_0,a_1,a_2])
        X = np.abs(p_x.roots) 
        
        if len(X)==1:
            omega_1[i] = np.sqrt(X[0]) 
            omega_2[i] = np.sqrt(X[0]) 
        else:
            omega_1[i] = np.sqrt(X[0]) 
            omega_2[i] = np.sqrt(X[1]) 
            
    plt.close() 
    plt.plot(qa,omega_2) 
    plt.plot(qa,omega_1)
    plt.xlabel('q*a')
    plt.ylabel('w')
    plt.vlines(np.pi,0,np.max(omega_1))
    plt.vlines(-np.pi,0,np.max(omega_1))
    plt.show()
    
def omeg2(f1,f2s): 
    n = int(len(f2s)/2)
    qa = np.linspace(-np.pi-1,np.pi+1,2000)
    omega = 0 ; domega = [] ; y = []
    for f2 in f2s:
         
        fun = [lambda x:f1*np.sin(x/2)**2,
               lambda x:f2*np.sin(x)**2]
        
        omega = 2*np.sqrt(fun[0](qa)+fun[1](qa)) 
        d,_ = find_peaks(omega) 
        domega.append(d) ; y.append(omega[d])
        
        #plt.plot(qa[d],omega[d],'ro') 
        plt.plot(qa,omega,label='f2 = '+str(f2)) 
        plt.title('f1 = 1')
        plt.xlabel('q*a')
        plt.ylabel('w*m**0.5')
    
    domega = np.sort(domega) ; y = np.sort(y) 
    plt.plot(qa[domega],y,':r')  
    plt.legend()    
    plt.vlines(np.pi,0,np.max(omega))
    plt.vlines(-np.pi,0,np.max(omega))     
    plt.show()    
    
omeg2(1,[0.1,0.5,1,2,4,6,8,10]) 
   
#poly =[121,40,-242,404,121,40]
#p = np.poly1d(poly[-1:-6:-1])  
#qa = np.linspace(-np.pi,np.pi,1000)
#exp_i_q_a_quart = np.exp(1j*qa/4)
#omega_1 = lambda x: np.sqrt(11+11*x**4+(1/x)*np.sqrt(p(x**2)))
#omega_2 = lambda x: np.sqrt(11+11*x**4-(1/x)*np.sqrt(p(x**2)))
# 
#plt.plot(qa,omega_1(exp_i_q_a_quart))
#plt.plot(qa,omega_2(exp_i_q_a_quart)) 
#plt.show()