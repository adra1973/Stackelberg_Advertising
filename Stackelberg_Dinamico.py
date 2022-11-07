#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 17:53:12 2019

@author: alandanielroblesaguilar
"""

# Program for modeling a stochastic and dynamical game among two enterprises 
# competing for market share, with the Solution Concept of Stackelberg.

from random import uniform

from scipy import *
from numpy import *

from scipy import linalg
import numpy as np

import scipy.optimize as optimize

ns=101
nstg=102
nest=-1
ntry=1

nf=11
fi=0.5
ff=1.5

ng=11
gi=0.5
gf=1.5

file= open("Stackelberg_Dinamico_100_17MAR22.txt","w")   

V0=zeros(ns)
V1=zeros(ns)
V2=zeros(ns)
W0=zeros(ns)
W1=zeros(ns)
W2=zeros(ns)
E1=zeros(2)

PFG=zeros((max(nf,ng),9,2))

PF1=zeros(nf)
PG1=zeros(ng)    

PF=zeros(nf)
PG=zeros(ng)

F=zeros(nf)
G=zeros(ng)

p=700.0
beta=(1.0/1.15)
xi=0.03
zeta=0.03
alfa=0.5
a=0.0

wb1=0.0
wb0=0.0

# pf & pg are the values of probability for the Binomial Distributions.

#pf=0.48489
#pg=0.51511

pf=0.4
pg=0.4

def max_f1(y):

    s=0.0    
    
    for i in range(nf):
        for j in range(ng):                        
            r=min(1.0,max(0.0,x+F[i]*(1.0-x)*xi*a**alfa-G[j]*x*zeta*y**alfa))
            s+=(wb1*r+wb0)*PF[i]*PG[j]
        
    fx=-(p*(1.0-x)-y + beta * s)

    return fx


def max_l1(u):

    s1=0.0    
    a=u
    res=optimize.minimize_scalar(max_f1,bounds=(0.0,200.0),method='bounded')
    b=res.x

    
    for i in range(nf):
        for j in range(ng):          
            r1=min(1.0,max(0.0,x+F[i]*(1.0-x)*xi*u**alfa-G[j]*x*zeta*b**alfa))    
            t=min((ns-2),int((ns-1)*r1))
            s1+=(((V0[t+1]-V0[t])/(1.0/float(ns-1)))*(r1-t/float(ns-1))+V0[t])*PF[i]*PG[j]
    
    fx1=-(p*x-u + beta * s1)
        
    return fx1


def reg_lin(A):
    
    E=zeros(2)
    wsx=0.0
    wsy=0.0
    wsx2=0.0
    wsxy=0.0  
    
    for i in range(ns):

        wsx+=i/float(ns-1)        
        wsy+=A[i]
        wsx2+=(i/float(ns-1))**2
        wsxy+=(i/float(ns-1))*A[i]        
        

    E[1]=(wsxy-wsx*wsy/float(ns))/(wsx2-(1.0/float(ns))*wsx**2)
    E[0]=wsy/float(ns)- E[1]*wsx/float(ns)
    
    return E

def est_emp(d):
    A=zeros((max(nf,ng),2))
    B=zeros((max(nf,ng),9,2))
    c=0
    v=0
    while v<=d:
        
        u=c
        while u<10**v:
        
            u1 = uniform(0,1)
            c1 = float(pf)/float((1-pf))
            pr1 = ((1-pf)**(nf-1)) 
            F1 = pr1
            for l in range(nf):
                if u1<F1:
                    A[l,0]+=1
                    break
                pr1=pr1*c1*float(nf-1-l)/float(l+1) 
                F1=F1+pr1               

            u2 = uniform(0,1)
            c2 = float(pg)/float((1-pg))
            pr2 = ((1-pg)**(ng-1)) 
            F2 = pr2
            for i in range(ng):
                if u2<F2:
                    A[i,1]+=1
                    break
                pr2=pr2*c2*float(ng-1-i)/float(i+1) 
                F2=F2+pr2  
            u+=1
            c+=1
            
        for m in range(nf):
            B[m,v,0]=float(A[m,0])/float((c))
        for j in range(ng):
            B[j,v,1]=float(A[j,1])/float((c))
        print (v)
        v+=1
      
    return B




PFG=est_emp(nest)

incf=(ff-fi)/(nf-1)
for i in range(nf):
    F[i]=fi+incf*i

incg=(gf-gi)/(ng-1)
for i in range(ng):
    G[i]=gi+incg*i

# PF[] & PG[] are vectors of probabilities from the Binimial Distributions.                  

t=1

while t<=ntry:
    
    w=0
    while w<=(nest+1):
        
        if w==(nest+1):
            for l in range(nf):
                PF[l]=(math.factorial(nf-1)/(math.factorial(nf-1-l)* \
                math.factorial(l)))*(pf**l)*((1-pf)**(nf-1-l)) 

            for l in range(ng):
                PG[l]=(math.factorial(ng-1)/(math.factorial(ng-1-l)* \
                math.factorial(l)))*(pg**l)*((1-pg)**(ng-1-l))    

        if w<=nest:
            for m in range(nf):
                PF[m]=PFG[m,w,0]
            for j in range(ng):
                PG[j]=PFG[j,w,1]
                
                
        print (PF,PG )  


        for j in range(ns):    
            V0[j]=p*j/float(ns-1)
            W0[j]=p*(1.0-j/float(ns-1))

        x=0.0
        z=1
        while z<=nstg:
            E1=reg_lin(W0)
            wb0=E1[0]
            wb1=E1[1]    
            print ("")
            print ("Trajectory",t,"Estimation",w, "Stage", z)
            print ("")
    
            x=0.0
            k=0
            mv=0.0
            while k<ns:
                x=k/float(ns-1)
        
                res=optimize.minimize_scalar(max_l1,bounds=(0.0,200.0),method='bounded')
                V2[k]=res.x
                V1[k]=-max_l1(V2[k])  
        
                a=V2[k]
                res=optimize.minimize_scalar(max_f1,bounds=(0.0,200.0),method='bounded')
                W2[k]=res.x
                W1[k]=-max_f1(W2[k])  

                #print (round(x,2),round(V1[k],2),round(W1[k],2),round(V2[k],2),round(W2[k],2))
        
                if z==101:
                    file.write('% s' %w+ ' , ')
                    file.write('% s' %x+ ' , ')                    
                    for i in range(201):
                        a=i
                        res=optimize.minimize_scalar(max_f1,bounds=(0.0,200.0),method='bounded')
                        b1=res.x
                        file.write('% s' %b1+ ' , ')
                    file.write('\n')
                    
                if z==100:
                    file.write('% s' %w+ ' , ')
                    file.write('% s' %x+ ' , ')
                    file.write('% s' %V1[k]+ ' , ')
                    file.write('% s' %W1[k]+ ' , ')
                    file.write('% s' %V2[k]+ ' , ')
                    file.write('% s' %W2[k]+ ' , ')
                    file.write('\n')
                    
        
        
                #if max(abs(V0[k]-V1[k]),abs(W0[k]-W1[k]))>mv:
                #    mv=max(abs(V0[k]-V1[k]),abs(W0[k]-W1[k]))
        

                k+=1
            #print mv
            #if mv<.001:
            #    break
    
    
            for i in range(ns):
                V0[i]=V1[i]
                W0[i]=W1[i]
                
                
          
                
    
    
            z+=1    
    
        w+=1
    
    t+=1
    

file.close()