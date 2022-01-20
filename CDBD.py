import numpy as np
import cmath

#Giro                       %Rango (a partir de Home)
q1 = (   0+  0)*np.pi/180   #  -50<q1<180
q2 = ( 115-115)*np.pi/180   # -140<q2<10
q3 = (  90- 90)*np.pi/180   # -110<q3<140
q4 = (   0+  0)*np.pi/180   #  -10<q4<130
q5 = ( -90+ 90)*np.pi/180   # -180<q5<180
q6 = (   0+  0)*np.pi/180   #  -70<q6<90

d1 = 183.73 #dist. del pecho al hombro (q2)
a2 = 58.28
d3 = 244.4
a3 = 13.76 + 5
a4 = 18.53
d5 = 78.53 + 44.5
d7 = 150 #150 mm is dist. from q6 to ee; 71 mm, from q7.
k = 65*np.pi/180; #brazo derecho. Izquierdo would be -65 deg.

def DH(q, d, a, alpha):
    Cq = np.cos(q)
    Sq = np.sin(q)
    Calpha = np.cos(alpha)
    Salpha = np.sin(alpha)
    Tzrx = np.array([[Cq,-Calpha*Sq,Salpha*Sq,a*Cq], [Sq,Calpha*Cq,-Cq*Salpha,a*Sq], [0,Salpha,Calpha,d], [0,0,0,1]])
    return Tzrx

def WristTransformation(Q1,Q2,Q3,Q4,Q5):
    T0B = DH(0,         0,  0,  k)
    TBH = DH(np.pi/2,   0,  0,  0)
    T01 = DH(Q1,        d1, 0,  np.pi/2)
    T12 = DH(Q2,        0,  a2, -np.pi/2)
    T23 = DH(Q3,        d3, a3, np.pi/2)
    T34 = DH(Q4,        0,  a4, -np.pi/2)
    T45 = DH(Q5,        d5, 0,  np.pi/2)
    return np.matmul(np.matmul(np.matmul(np.matmul(np.matmul(np.matmul(T0B,TBH),T01),T12),T23),T34),T45)

def EndEfectorTransformation(Q6):
    T56 = DH(Q6,    0,  0,  -np.pi/2)
    T67 = DH(0,     d7, 0,  0)
    return np.matmul(T56,T67)

def FinalTransformation(Q1,Q2,Q3,Q4,Q5,Q6):
    r = np.matmul(WristTransformation(Q1,Q2,Q3,Q4,Q5),EndEfectorTransformation(Q6))
    return r

#--- Matriz de orientación-posición del efector final (end effector) ---#
Tee = FinalTransformation(q1,q2,q3,q4,q5,q6)
np.set_printoptions(precision=4, suppress=True)
print("Tee:\n{}".format(Tee))