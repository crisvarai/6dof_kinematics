import numpy as np
import cmath

#Giro                       %Rango (a partir de Home)
q1 = (   0-  0)*np.pi/180   #  -50<q1<180
q2 = ( 115-  0)*np.pi/180   # -140<q2<10
q3 = (  90-  0)*np.pi/180   # -110<q3<140
q4 = (   0+  0)*np.pi/180   #  -10<q4<130
q5 = ( -90+  0)*np.pi/180   # -180<q5<180
q6 = (   0+  1)*np.pi/180   #  -70<q6<90

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

np.set_printoptions(precision=4, suppress=True)

#--- Matriz de orientación-posición del efector final (end effector) ---#
print("Input Config: <{},{},{},{},{},{}>".format(int(q1*180/np.pi), int(q2*180/np.pi), int(q3*180/np.pi), int(q4*180/np.pi), int(q5*180/np.pi), int(q6*180/np.pi)))

Tee = FinalTransformation(q1,q2,q3,q4,q5,q6)
print("Input Tee:\n{}".format(Tee))


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
q4a = np.angle((C4+cmath.sqrt(D4))/complex(B4,-A4))
q4b = np.angle((C4-cmath.sqrt(D4))/complex(B4,-A4))
q4adeg = q4a*180/np.pi
q4bdeg = q4b*180/np.pi
print("q4a,q4b:{},{}".format(q4adeg,q4bdeg))
if q4adeg>(-10) and q4adeg<130:
    if q4bdeg>(-10) and q4bdeg<130:
        print("Both q4's are viable.")
        if abs(q4a)<abs(q4b):
            q4 = q4a
        else:
            q4 = q4b
    else:
        print("q4a then!")
        q4 = q4a
elif (q4b*180/np.pi)>-10 and (q4b*180/np.pi)<130:
    print("q4b then!")
    q4 = q4b
else:
    print("None in range. Are we outside the viable workspace? Defaulting to positive sqrt.")
    q4 = q4a


#--- Cálculo de q2
# Solución de la ecuación A2sin(q2) + B2sin(q2) = C2
A2 = -(a2 + np.cos(q3)*(a3 + a4*np.cos(q4) - d5*np.sin(q4)))
B2 = -(d3 + d5*np.cos(q4) + a4*np.sin(q4))
C2 = Xws[1]*np.cos(5*np.pi/36) - Xws[2]*np.sin(5*np.pi/36)
D2 = pow(C2,2) - (pow(B2,2) + pow(A2,2))
#print("A2, B2, C2, D2 = {}, {}, {}, {}".format(A2, B2, C2, D2))
# Discernir entre q2's
q2a = np.angle((C2-cmath.sqrt(D2))/complex(B2,-A2))
q2b = np.angle((C2+cmath.sqrt(D2))/complex(B2,-A2))
q2adeg = q2a*180/np.pi
q2bdeg = q2b*180/np.pi
print("q2a,q2b:{},{}".format(q2adeg,q2bdeg))
if q2adeg>(-25) and q2adeg<115:
    if q2bdeg>(-25) and q2bdeg<115:
        print("Both q2's are viable.")
        if abs(q2a)<abs(q2b):
            q2 = q2a
        else:
            q2 = q2b
    else:
        print("q2a then!")
        q2 = q2a
elif (q2b*180/np.pi)>(-25) and (q2b*180/np.pi)<(115):
    print("q2b then!")
    q2 = q2b
else:
    print("None in range. Are we outside the viable workspace? Defaulting to positive sqrt.")
    q2 = q2a


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
q6a = np.arccos(B6[2]/d7)
q6b = -np.arccos(B6[2]/d7)


#--- Cálculo de q5
# Solución de la ecuación A1sin(q1) + B1sin(q1) = C1
A5a = -d7*np.sin(q6a)
B5a = -d7*np.sin(q6a)
C5a = B6[0]+B6[1]
D5a = pow(C5a,2) - (pow(B5a,2) + pow(A5a,2))
# avoid division by zero by...
#q5a=np.angle((C5+cmath.sqrt(D5))/complex(B5,-A5))
#q5b=np.angle((C5-cmath.sqrt(D5))/complex(B5,-A5))
q5aa=np.angle(C5a+cmath.sqrt(D5a)) - np.angle(complex(B5a,-A5a))
q5ab=np.angle(C5a-cmath.sqrt(D5a)) - np.angle(complex(B5a,-A5a))

A5b = -d7*np.sin(q6b)
B5b = -d7*np.sin(q6b)
C5b = B6[0]+B6[1]
D5b = pow(C5b,2) - (pow(B5b,2) + pow(A5b,2))
# avoid division by zero by...
#q5a=np.angle((C5+cmath.sqrt(D5))/complex(B5,-A5))
#q5b=np.angle((C5-cmath.sqrt(D5))/complex(B5,-A5))
q5ba=np.angle(C5b+cmath.sqrt(D5b)) - np.angle(complex(B5b,-A5b))
q5bb=np.angle(C5b-cmath.sqrt(D5b)) - np.angle(complex(B5b,-A5b))

# Posición a la que llega con q5a y q5b
# Teea and Teeb se definen como la transformación hasta el efector final
Teeaa = FinalTransformation(q1[0],q2[0],q3,q4,q5aa[0],q6a[0])
Teeab = FinalTransformation(q1[0],q2[0],q3,q4,q5ab[0],q6a[0])
Teeba = FinalTransformation(q1[0],q2[0],q3,q4,q5ba[0],q6b[0])
Teebb = FinalTransformation(q1[0],q2[0],q3,q4,q5bb[0],q6b[0])

# Pwa and Pwb se definen como la posición del efector final
Peeaa = np.array([[Teeaa[0,3]], [Teeaa[1,3]], [Teeaa[2,3]]])
Peeab = np.array([[Teeab[0,3]], [Teeab[1,3]], [Teeab[2,3]]])
Peeba = np.array([[Teeba[0,3]], [Teeba[1,3]], [Teeba[2,3]]])
Peebb = np.array([[Teebb[0,3]], [Teebb[1,3]], [Teebb[2,3]]])

if (np.linalg.norm(X70-Peeaa)<=np.linalg.norm(X70-Peeab)):
    q5a=q5aa
    Teea=Teeaa
else:
    q5a=q5ab
    Teea=Teeab

if (np.linalg.norm(X70-Peeba)<=np.linalg.norm(X70-Peebb)):
    q5b=q5ba
    Teeb=Teeba
else:
    q5b=q5bb
    Teeb=Teebb

if (abs(3-np.trace(np.matmul(Teea[0:3,0:3],R70.transpose())))<=abs(3-np.trace(np.matmul(Teeb[0:3,0:3],R70.transpose())))):
    q5=q5a
    q6=q6a
else:
    q5=q5b
    q6=q6b

q1deg =int(q1*180/np.pi)
q2deg =int(q2*180/np.pi)
q3deg =int(q3*180/np.pi)
q4deg =int(q4*180/np.pi)
q5deg =int(q5*180/np.pi)
q6deg =int(q6*180/np.pi)
#print("Output Config:<{:4},{:4},{:4},{:4},{:4},{:4}".format(q1deg, q2deg, q3deg, q4deg, q5deg, q6deg))
print("Output Config: <{},{},{},{},{},{}>".format(q1deg, q2deg, q3deg, q4deg, q5deg, q6deg))
Tee = FinalTransformation(q1[0],q2[0],q3,q4,q5[0],q6[0])
print("Output Tee:\n{}".format(Tee))