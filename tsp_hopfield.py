#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 22:11:32 2018

@author: yi-chun
"""

# 數學科學院 1801210118 馮逸群 Yi-Chun, Feng

import numpy as np



location=np.matrix([[0.4000,0.4439],[0.2439,0.1416],[0.1707,0.2293],[0.2293,0.7610],\
                    [0.5171,0.9414],[0.8732,0.6536],[0.6878,0.5219],[0.8488,0.3609],\
                    [0.6683,0.2536],[0.6195,0.2634]])
def distance(location):
    N=len(location)
    Dis=np.zeros([N,N])#創立距離矩陣
    for i in range(N):
        for j in range(N):
            Dis[i,j]=np.sqrt((location[i,0]-location[j,0])**2+(location[i,1]-location[j,1])**2)
            #算出各城市相互距離
    return Dis




def delta(i,j):
    if i==j:
        return 1
    else:
        return 0
    #用在weight的建立中
    #weight=-Q*delta(1-delta)-S*delta(1-delta)-T-P*D*(delta+delta)
    
    


def Weight(A,B,C,D,Dis):
    N=len(Dis)
    W=np.zeros([N,N,N,N]);
    for x in range(N):
        for i in range(N):
            for y in range(N):
                for j in range(N):
                    #W[x,i,y,j]=-A*delta(x,y)*(1-delta(i,j))-D*Dis[x,y]*delta(j,x+1)
                    W[x,i,y,j]=-A*delta(x,y)*(1-delta(i,j))-B*delta(i,j)*(1-delta(x,y))-C-D*Dis[x,y]*(delta(j,i+1)+delta(j,i-1))
    return W


def CalcuDeltaU(U,V,A,B,C,D,dis):
    N=len(V)
    deltaU=np.zeros([N,N])
    for x in range(N):
        for i in range(N):
            t1=0
            for j in range(N):
                if j!=i:
                    t1+=V[x,j]
            t1=A*t1
            t2=0
            for y in range(N):
                if y!=x:
                    t2+=V[y,i]
                    #t2 += V[x, y]
            t2=B*t2
            
            t3=0
            
            
            for j in range(N):
                t3+=V[x,j]
            for y in range(N):
                t3 += V[y, i]
            t3 = t3 - 2*N
            t3 = C * t3

            t4=0
            for y in range(N):
                if y!=x:
                    i1=i+1
                    i2=i-1
                    if i1==N:
                        i1=0
                    if i2==-1:
                        i2=N-1
                    t4+=dis[x,y]*(V[y,i1]+V[y,i2])
            t4=D*t4
            deltaU[x,i]=-t1-t2-t3-t4#-U[x,i]
    return deltaU



def tanh(x):
    
    if x>5:
        return 0.99995
    if x<-5:
        return -0.99995
    res=(np.exp(x)-np.exp(-x))/(np.exp(x)+np.exp(-x))
    return res

def CalcuV(U,U0):
    N=len(U)
    V=np.zeros([N,N])
    for x in range(N):
        for i in range(N):
            
            V[x,i]=0.5*(1+tanh(U[x,i]/U0))
            if np.isnan(V[x,i]):
                
                break
    return V
def CalcuU(U,deltaU,step):
    for x in range(N):
        for i in range(N):
            U[x,i]+=deltaU[x,i]*step
    return U

def CacuEnergy( V,Dis,A,B,C,D,N):
    V=np.matrix(V)
    t1 = np.sum(np.power(np.sum(V, 1) - 1,2))
    t2 = np.sum(np.power(np.sum(V, 0) - 1,2))
    idx = [i for i in range(1, N)]
    idx = idx + [0]
    Vt = V[:, idx]
    t3 = Dis* Vt
    t3=np.sum(np.sum(np.multiply(V,t3)))
    E=0.5*(A*t1+B*t2+D*t3)
    return E

def isOK(V):
    N=len(V)
    NewV=np.zeros([N,N])
    Route=[]
    for i in range(N):
        mm=np.max(V[:,i])
        for j in range(N):
            if V[j,i]==mm:
                NewV[j,i]=1
                Route+=[j]
                break
    return Route,NewV

def CheckValue(V):
    N=len(V)
    for i in range(N):
        for j in range(N):
            if V[i,j]>0.3 and V[i,j]<0.7:
                return 0
    return 1

def CheckPosition(V):
    N=len(V)
    V1=np.zeros([N,N])
    for i in range(N):
        row_sum=0;
        clomn_sum=0;
        for j in range(N):
            if V[i,j]>=0.9:
                V1[i,j]=1
            if V[j,i]>=0.9:
                V1[j,i]=1
            row_sum+=V1[i,j]
            clomn_sum+=V1[j,i]
    if row_sum!=N:
        return 0
    if clomn_sum!=N:
        return 0
    return 1
def CalcuDis(Route,Dis):
    distance=0
    for i in range(N-1):
        distance+=Dis[Route[i],Route[i+1]]
    distance+=Dis[Route[N-1],Route[0]]
    return distance
def Plot(Route,Data):
    Route+=[Route[0]]
    Points=Data[Route,:]
    



Dis=distance(location)
N=len(Dis)
A=500
D=200
U0=0.02
step=0.0001
B=500
C=500
#P=200,Q=S=T=500,u_0=0.02
maxecho=10

W=Weight(A,B,C,D,Dis)
U=0.002*(2*np.random.rand(N,N)-1)

V=1/N*(2*np.random.rand(N,N)-1)
Dis_best=np.sum(Dis)
Route_final="No answear"
E_all=[]
Dis_all=[]
Route_all=[]
for r in range(maxecho):
    count=0
    flag=1
    while count<500:
       
        E_all+=[CacuEnergy( V,Dis,A,B,C,D,N)]
        count+=1
        deltaU=CalcuDeltaU(U,V,A,B,C,D,Dis)
        U=CalcuU(U,deltaU,step)
        V=CalcuV(U,U0)
    Route, NewV = isOK(V)
    
    if len(np.unique(Route))==8:
        print(r,"Route found:",Route)
        Dis_r=CalcuDis(Route,Dis)
        Dis_all+=[Dis_r]
        Route_all+=[Route]
        print("Distance is",Dis_r)
        
        if Dis_r<Dis_best:
            Dis_best=Dis_r
            Route_final=Route
    #V=V+2*(2*np.random.rand(N,N)-1)
    U=0.002*(2*np.random.rand(N,N)-1)
    V= 1 / N * (2 * np.random.rand(N, N) - 1)
    A+=10*(2*np.random.rand(1,1)-1)
    B=A
    C+=10*(2*np.random.rand(1,1)-1)
    D+=10*(2*np.random.rand(1,1)-1)


if Route_final!="No answear":
    Plot(Route_final,Data)
print(V)
Route,NewV=isOK(V)
print(NewV)

print("The final best Route is :",Route_final)
print("The distance of the Route is:",Dis_best) 







