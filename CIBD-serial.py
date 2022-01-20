import numpy as np
import cmath
import serial
from serial.serialutil import SerialException

#Giro                       %Rango (a partir de Home)
q1 = (   0+  0)*np.pi/180   #  -50<q1<180
q2 = ( 115-  0)*np.pi/180   # -140<q2<10
q3 = (  90-  0)*np.pi/180   # -110<q3<140
q4 = (   0+  0)*np.pi/180   #  -10<q4<130
q5 = ( -90+  0)*np.pi/180   # -180<q5<180
q6 = ( -90+  0)*np.pi/180   #  -70<q6<90

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


#--- CINEMATICA INVERSA ---#

#--- Pocición del Shoulder y del Wrist
Xs = np.array([[0], [-d1*np.cos(5*np.pi/36)], [d1*np.sin(5*np.pi/36)]])
#print("Xs:\n{}".format(Xs))

#--- Matriz de orietación posción conocida o del rastreador (tracker) del end effector
#R70=[Tee(1,1) Tee(1,2) Tee(1,3) 0; 
#     Tee(2,1) Tee(2,2) Tee(2,3) 0;
#     Tee(3,1) Tee(3,2) Tee(3,3) 0;
#           0        0        0  1];
R70 = np.array([[Tee[0,0], Tee[0,1], Tee[0,2]], [Tee[1,0], Tee[1,1], Tee[1,2]], [Tee[2,0], Tee[2,1], Tee[2,2]]])
#X70=[Tee(1,4);
#     Tee(2,4);
#     Tee(3,4);
#           1];
X70 = np.array([[Tee[0,3]], [Tee[1,3]], [Tee[2,3]]])

Xw = X70-np.matmul(R70,np.array([[0],[0],[d7]]))
#print("Xw\n{}".format(Xw))

#--- Vector Shoulder->Wrist
Xws = (Xw - Xs)

#--- Vector Wrist->Efector Final 
X7w = (X70 - Xw)

#--- Variables para el calculo de la CI
N2 = pow(np.linalg.norm(Xws), 2)


#--- Cálculo de q4
# Solución de la ecuación A4sin(q4) + B4sin(q4) = C4
A4 = 2*(a4*d3 - a3*d5 - a2*d5*np.cos(q3))
B4 = 2*(a3*a4 + d3*d5 + a2*a4*np.cos(q3))
C4 = N2 - (pow(a2,2) + pow(a3,2) + pow(a4,2) + pow(d3,2) + pow(d5,2) + 2*a2*a3*np.cos(q3))
D4 = pow(C4,2) - (pow(B4,2) + pow(A4,2))
#print("A4, B4, C4, D4 = {}, {}, {}, {}".format(A4, B4, C4, D4))

# Discernir entre q4's
q4a = np.angle((C4+cmath.sqrt(D4))/(complex(B4,-A4)))
q4b = np.angle((C4-cmath.sqrt(D4))/(complex(B4,-A4)))
if (q4a*180/np.pi<-10) or (q4a*180/np.pi>130):
    q4 = q4b
elif (q4b*180/np.pi<-10) or (q4b*180/np.pi>130):
    q4 = q4a
else:
    q4 = min(q4a, q4b)

#--- Cálculo de q2
# Solución de la ecuación A2sin(q2) + B2sin(q2) = C2
A2 = -(a2 + np.cos(q3)*(a3 + a4*np.cos(q4) - d5*np.sin(q4)))
B2 = -(d3 + d5*np.cos(q4) + a4*np.sin(q4))
C2 = Xws[1]*np.cos(5*np.pi/36) - Xws[2]*np.sin(5*np.pi/36)
D2 = pow(C2,2) - (pow(B2,2) + pow(A2,2))
#print("A2, B2, C2, D2 = {}, {}, {}, {}".format(A2, B2, C2, D2))
q2 = np.angle((C2-cmath.sqrt(D2))/complex(B2,-A2))

#--- Cálculo de q1
# Solución de la ecuación A1sin(q1) + B1sin(q1) = C1
A1 = np.sin(q2)*(d3 + d5*np.cos(q4) + a4*np.sin(q4)) - np.cos(q2)*(a2 + np.cos(q3) * (a3 + a4*np.cos(q4) - d5*np.sin(q4)))
B1 = -np.sin(q3)*(a3 + a4*np.cos(q4) - d5*np.sin(q4))
C1 = Xws[0]
D1 = pow(C1,2) - (pow(B1,2) + pow(A1,2))
#print("A1, B1, C1, D1 = {}, {}, {}, {}".format(A1, B1, C1, D1))

# Discernir entre q1's
q1a=np.angle((C1+cmath.sqrt(D1))/complex(B1,-A1))
q1b=np.angle((C1-cmath.sqrt(D1))/complex(B1,-A1))
#print("q1a, q1b = {}, {}".format(q1a, q1b))

# Posición a la que llega con q1a y q1b
# Twa and Twb se definen como la transformación de la muñeca con q1a y q1b
Twa = WristTransformation(q1a[0],q2[0],q3,q4,0)
Twb = WristTransformation(q1b[0],q2[0],q3,q4,0)
# Pwa and Pwb se definen como la posición de la muñeca con q1a y q1b
Pwa = np.array([[Twa[0,3]], [Twa[1,3]], [Twa[2,3]]])
Pwb = np.array([[Twb[0,3]], [Twb[1,3]], [Twb[2,3]]])

if (np.linalg.norm(Xw-Pwa)<=np.linalg.norm(Xw-Pwb)):
    q1=q1a
    A6temp = Twa
else:
    q1=q1b
    A6temp = Twb
    
#--- Cálculo de q6    
A6 = np.matmul(A6temp[0:3,0:3],np.array([[1,0,0],[0,0,1],[0,-1,0]]))
B6 = np.matmul(A6.transpose(),X7w)
q6 = np.arccos(B6[2]/d7)

#--- Cálculo de q5
# Solución de la ecuación A1sin(q1) + B1sin(q1) = C1
A5 = -d7*np.sin(q6)
B5 = -d7*np.sin(q6)
C5 = B6[0]+B6[1]
D5 = pow(C5,2) - (pow(B5,2) + pow(A5,2))
# avoid division by zero by...
#q5a=np.angle((C5+cmath.sqrt(D5))/complex(B5,-A5))
#q5b=np.angle((C5-cmath.sqrt(D5))/complex(B5,-A5))
q5a=np.angle(C5+cmath.sqrt(D5)) - np.angle(complex(B5,-A5))
q5b=np.angle(C5-cmath.sqrt(D5)) - np.angle(complex(B5,-A5))
        
# Posición a la que llega con q5a y q5b
# Teea and Teeb se definen como la transformación hasta el efector final
Teea = FinalTransformation(q1[0],q2[0],q3,q4,q5a[0],q6[0])
Teeb = FinalTransformation(q1[0],q2[0],q3,q4,q5b[0],q6[0])

# Pwa and Pwb se definen como la posición del efector final
Peea = np.array([[Teea[0,3]], [Teea[1,3]], [Teea[2,3]]])
Peeb = np.array([[Teeb[0,3]], [Teeb[1,3]], [Teeb[2,3]]])

if (np.linalg.norm(X70-Peea)<=np.linalg.norm(X70-Peeb)):
    q5=q5a
else:
    q5=q5b    

q1deg =int(q1*180/np.pi)
q2deg =int(q2*180/np.pi)
q3deg =int(q3*180/np.pi)
q4deg =int(q4*180/np.pi)
q5deg =int(q5*180/np.pi)
q6deg =int(q6*180/np.pi)
print("q1:{:4}    q2:{:4}    q3:{:4}    q4:{:4}    q5:{:4}    q6:{:4}".format(q1deg, q2deg, q3deg, q4deg, q5deg, q6deg))

command = "<1.2.1.{}.{}.{}.{}.{}.{}>".format(q1deg, q2deg, q3deg, q4deg, q5deg, q6deg)
print("Command: {}".format(command))