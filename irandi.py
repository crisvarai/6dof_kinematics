import numpy as np
import cmath
import serial
from serial.serialutil import SerialException
#import matplotlib.pyplot as plt
#from mpl_toolkits import mplot3d

#Giro               %Rango
q1=   0*np.pi/180;  # -20<q1<180
q2= -115*np.pi/180;	# -45<q2<125
q3=  90*np.pi/180;	# -45<q3<180
q4=   0*np.pi/180;  #  0<q4<135
q5= -90*np.pi/180;  # -180<q5<180
q6=   0*np.pi/180;  # -90<q6<90
#q7=?

d1=183.73;
a2=58.28;
d3=244.4;
a3=13.76+5;
a4=18.53;
d5=78.53+44.5;
d7=71;  #150 is d from q6
#k=((-65*np.pi)/180); #brazo izq
#k=((25*np.pi)/180)
k=((65*np.pi)/180); #brazo derecho

def DH(q, d, a, alpha):
    Rz = np.array([[np.cos(q), -np.sin(q), 0, 0], [np.sin(q), np.cos(q), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    Tr = np.array([[1, 0, 0, a], [0, 1, 0, 0], [0, 0, 1, d], [0, 0, 0, 1]])
    Rx = np.array([[1, 0, 0, 0], [0, np.cos(alpha), -np.sin(alpha), 0], [0, np.sin(alpha), np.cos(alpha), 0], [0, 0, 0, 1]])
    return np.matmul(np.matmul(Rz, Tr), Rx)

#----- Tabla de Denavit-Hartenberg
T0B=DH(       0,    0,  0,  k);
T0H=DH(np.pi/2,    0,  0,  0);
T01=DH(      q1,    d1, 0,  -np.pi/2);
T12=DH(      q2,    0,  a2, np.pi/2);
T23=DH(      q3,    d3, a3, -np.pi/2);
T34=DH(      q4,    0,  a4, np.pi/2);
T45=DH(      q5,    d5, 0,  -np.pi/2);
T56=DH(      q6,    0,  0,  np.pi/2);
T67=DH(0,           d7, 0,  0);

#----- Matriz de orientación-posición del end effector
Tee=np.matmul(np.matmul(np.matmul(np.matmul(np.matmul(np.matmul(np.matmul(np.matmul(T0B,T0H),T01),T12),T23),T34),T45),T56),T67)
#print(Tee)
#print("Tee\n")

#----- CINEMATICA INVERSA -------------------
#---- Matriz de orietación posción conocida o del Tracker del end effector
#forward extension
# curve
t = np.linspace(-np.pi,np.pi,60)
x = 0*t + 100
y = -250 - 100*np.cos(t)
z = 80 + 100*np.sin(t)

for i in range(len(t)):
    print("Iter: {}".format(i))
    #----- Matriz de rotación y posicion del end effector
    #R70=[Tee(1,1) Tee(1,2) Tee(1,3) 0; 
    #     Tee(2,1) Tee(2,2) Tee(2,3) 0;
    #     Tee(3,1) Tee(3,2) Tee(3,3) 0;
    #           0        0        0  1];
    R70 = np.array([[Tee[0,0], Tee[0,1], Tee[0,2]], [Tee[1,0], Tee[1,1], Tee[1,2]], [Tee[2,0], Tee[2,1], Tee[2,2]]])
    #print(R70)
    #print("R70\n")

    #X70=[Tee(1,4);
    #     Tee(2,4);
    #     Tee(3,4);
    #           1];
    X70 = np.array([[x[i]], [y[i]], [z[i]]])
    print(X70)
    #print("X70\n")

    #----- Pocición del Shoulder y del Wrist
    #Xs = transp([0, d1*cos(5*pi/36), d1*sin(5*pi/36), 1]);
    Xs = np.array([[0], [d1*np.cos(5*np.pi/36)], [d1*np.sin(5*np.pi/36)]])
    #print(Xs)
    #print("Xs\n")
    #Xw = X70 - (R70*transp([0,0,d7,0]));
    Xw = X70#>to be able to ignore q5 & q6 - np.matmul(R70, np.array([[0],[0],[d7]])) 
    #print(Xw)
    #print("Xw\n")

    #----- Vector del Shoulder al Wrist
    Xws=(Xw-Xs)
    #print(Xws)
    #print("Xws\n")

    #----- Variables para el calculo de la CI
    #N2 = Xws(1)^2 + Xws(2)^2 + Xws(3)^2 + Xws(4)^2  %-- Norma al cuadrado ||xsw||^2
    #N2I =(2*(-a4*d3+a3*d5+a2*d5*cos(q3)))*sin(q4)+(2*(a3*a4+d3*d5+a2*a4*cos(q3)))*cos(q4)+(a2^2 + a3^2 + a4^2 + d3^2 + d5^2 +2*a2*a3*cos(q3))
    N2 = pow(np.linalg.norm(Xws), 2)
    #print(N2)
    #print("N2\n")
    #hws = Xws(2)*cos(5*pi/36) + Xws(3)*sin(5*pi/36)
    #hswI=(-(a2+cos(q3)*(a3+a4*cos(q4))+d5*sin(q4)))*sin(q2)+(d3+d5*cos(q4)-a4*sin(q4))*cos(q2)
    hws=Xws[1]*np.cos(5*np.pi/36) + Xws[2]*np.sin(5*np.pi/36)
    #print(hws)
    #print("hws\n")
    #Xswy = Xws(1)#the matlab reference was pulling the first element on the table, that is the x component, right?
    #XwsyI=(sin(q2)*(d3+d5*cos(q4)-a4*sin(q4))+cos(q2)*(a2+cos(q3)*(a3+a4*cos(q4)+d5*sin(q4))))*sin(q1)+(sin(q3)*(a3+a4*cos(q4)+d5*sin(q4)))*cos(q1)
    #Xwsy=Xws[1]
    #print(Xwsy)
    #print("Xwsy\n")
    Xwsx=Xws[0]#are we sure we don't actually need the x component of this vector.
    #print(Xwsx)
    #print("Xwsx\n")

    #-- Calculo de q4
    # Solución de la ecuación A4sin(q4)+B4sin(q4)=C4
    #A4 = (2*(-a4*d3+a3*d5+a2*d5*cos(q3)));
    A4 = 2*(-a4*d3 + a3*d5 + a2*d5*np.cos(q3))
    #print(A4)
    #print("A4\n")
    #B4 = (2*(a3*a4+d3*d5+a2*a4*cos(q3)));
    B4 = 2*(a3*a4 + d3*d5 + a2*a4*np.cos(q3))
    #print(B4)
    #print("B4\n")
    #C4 = N2 -(a2^2 + a3^2 + a4^2 + d3^2 + d5^2 +2*a2*a3*cos(q3));
    C4 = N2 - (pow(a2,2) + pow(a3,2) + pow(a4,2) + pow(d3,2) + pow(d5,2) + 2*a2*a3*np.cos(q3))
    #print(C4)
    #print("C4\n")
    #D4 = C4^2 -(B4^2 + A4^2);
    D4 = pow(C4,2) - (pow(B4,2) + pow(A4,2))
    #print(D4)
    #print("D4\n")
    #q4=angle((C4+sqrt(D4))/(B4-i*A4))%*180/pi %check for 
    q4 = np.angle((C4+cmath.sqrt(D4))/(complex(B4,-A4)))#*180/np.pi
    #print(q4*180/np.pi)
    #print("q4\n")

    #-- Calculo de q2
    # Solución de la ecuación A2sin(q2)+B2sin(q2)=C2
    #A2 = (-(a2+cos(q3)*(a3+a4*cos(q4))+d5*sin(q4)));
    A2 = -(a2 + np.cos(q3)*(a3 + a4*np.cos(q4)) + d5*np.sin(q4))
    #print(A2)
    #print("A2\n")
    #B2 = (d3+d5*cos(q4)-a4*sin(q4));
    B2 = d3 + d5*np.cos(q4) - a4*np.sin(q4)
    #print(B2)
    #print("B2\n")
    #C2 = (Xws(2)*cos(5*pi/36) + Xws(3)*sin(5*pi/36));
    #C2 = hws + a3 * np.cos(q3) * np.sin(q4)
    C2 = (Xws[1]*np.cos(5*np.pi/36) + Xws[2]*np.sin(5*np.pi/36));
    #print(C2)
    #print("C2\n")
    #D2 = C2^2 - (B2^2 + A2^2);
    D2 = pow(C2,2) - (pow(B2,2) + pow(A2,2))
    #print(D2)
    #print("D2\n")
    #q2=angle((C2-sqrt(D2))/(B2-i*A2))%*180/pi;
    q2 = np.angle((C2-cmath.sqrt(D2))/complex(B2,-A2))
    #print(q2)
    #print("q2\n")

    #%-- Calculo de q1
    #% Solución de la ecuación A1sin(q1)+B1sin(q1)=C1
    #A1=(sin(q2)*(d3+d5*cos(q4)-a4*sin(q4))+cos(q2)*(a2+cos(q3)*(a3+a4*cos(q4)+d5*sin(q4))));
    #A1 = np.sin(q2)*(d3 + d5*np.cos(q4) - a4*np.sin(q4)) + np.cos(q3)*(d5*np.sin(q4)*np.cos(q4) + d5*a4*np.cos(q4)*np.cos(q2) + a3*np.cos(q2)) - a2*np.cos(q2)
    A1 = (np.sin(q2)*(d3 + d5*np.cos(q4) - a4*np.sin(q4)) + np.cos(q2)*(a2 + np.cos(q3) * (a3 + a4*np.cos(q4) + d5*np.sin(q4))));
    #print(A1)
    #print("A1\n")
    #B1=(sin(q3)*(a3+a4*cos(q4)+d5*sin(q4)));
    B1 = np.sin(q3)*(a3 + a4*np.cos(q4) + d5*np.sin(q4))
    #print(B1)
    #print("B1\n")
    #C1=Xwsx(1);
    C1 = Xwsx
    #print(C1)
    #print("C1\n")
    #D1=C1^2 - (B1^2 + A1^2);
    D1 = pow(C1,2) - (pow(B1,2) + pow(A1,2))
    #print(D1)
    #print("D1\n")
    #q1=angle((C1+sqrt(D1))/(B1-i*A1))%*180/pi;
    #q1 = (np.angle(C1-cmath.sqrt(D1)) - np.angle(complex(B1,-A1)))
    q1=np.angle((C1+cmath.sqrt(D1))/complex(B1,-A1));
    #print(q1)
    #print("q1\n")

    print("<{},{},{},{}>".format(int(q1*180/np.pi), int(q2*180/np.pi), int(q3*180/np.pi), int(q4*180/np.pi)))