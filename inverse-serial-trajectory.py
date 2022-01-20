import numpy as np
import cmath
import serial
from serial.serialutil import SerialException
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

#q1=0*np.pi/180
#q2=115*np.pi/180
q3=90*np.pi/180
#q4=90*np.pi/180
q5=-90*np.pi/180
q6=0*np.pi/180

#forward extension
t = np.linspace(100,500,400)
x = 0*t+225
y = t
z = 0*t

ax = plt.axes(projection = '3d')
ax.view_init(25,165)
ax.plot(x,y,z)
plt.show()

d1=183.73
a2=58.28
d3=244.4
a3=13.76+5
a4=18.53 
d5=78.53+44.5
d7=150
#k=((65*np.pi)/180) #not used when rotation matrix is given

for i in range(len(t)-1):
    #%Tee=TBH*T01*T12*T23*T34*T45*T56*T67
    #mano palma abajo; probar c/palma hacia adentro, pulgar y dedos rotando al estirar
    Tee=np.array([[0,0,-1,x[i]],[-1,0,0,y[i]],[0,1,0,z[i]],[0,0,0,1]])
    #print(Tee)
    #print("Tee\n")

    R70 = np.array([[Tee[0,0], Tee[0,1], Tee[0,2], 0], [Tee[1,0], Tee[1,1], Tee[1,2], 0], [Tee[2,0], Tee[2,1], Tee[2,2], 0], [0, 0, 0, 1]])
    #print(R70)
    #print("R70\n")
    X70 = np.array([[Tee[0,3]], [Tee[1,3]], [Tee[2,3]], [1]])
    #print(X70)
    #print("X70\n")
    Xs = np.array([[d1*np.cos(5*np.pi/36)], [0], [d1*np.sin(5*np.pi/36)], [1]])
    #print(Xs)
    #print("Xs\n")
    Xw = X70 - np.matmul(R70, np.array([[0],[0],[d7],[0]]))
    #print(Xw)
    #print("Xw\n")
    Xws= Xw - Xs
    #print(Xws)
    #print("Xws\n")
    N2 = pow(np.linalg.norm(Xws), 2)
    #print(N2)
    #print("N2\n")
    hws=Xws[0]*np.cos(5*np.pi/36) + Xws[2]*np.sin(5*np.pi/36)
    #print(hws)
    #print("hws\n")
    Xwsy=Xws[1]
    #print(Xwsy)
    #print("Xwsy\n")

    A4 = 2*(-a4*d3 + a3*d5 + a2*d5*np.cos(q3))
    #print(A4)
    #print("A4\n")
    B4 = 2*(a3*a4 + d3*d5 + a2*a4*np.cos(q3))
    #print(B4)
    #print("B4\n")
    C4 = N2 - (pow(a2,2) + pow(a3,2) + pow(a4,2) + pow(d3,2) + pow(d5,2) + 2*a2*a3*np.cos(q3))
    #print(C4)
    #print("C4\n")
    D4 = pow(C4,2) - (pow(B4,2) + pow(A4,2))
    #print(D4)
    #print("D4\n")
    q4 = np.angle((C4+cmath.sqrt(D4))/(complex(B4,-A4)))*180/np.pi
    #print(q4)
    #print("q4\n")

    A2 = a2 - a4 * np.cos(q3) * np.cos(q4) - d5 * np.cos(q3)*np.sin(q4)
    #print(A2)
    #print("A2\n")
    B2 = d3 + d5 * np.cos(q4) - a4 * np.sin(q4)
    #print(B2)
    #print("B2\n")
    C2 = hws + a3 * np.cos(q3) * np.sin(q4)
    #print(C2)
    #print("C2\n")
    D2 = pow(C2,2) - (pow(B2,2) + pow(A2,2))
    #print(D2)
    #print("D2\n")
    q2 = -(np.angle(C2-cmath.sqrt(D2))+np.angle(complex(B2,-A2)))*180/np.pi
    #print(q2)
    #print("q2\n")

    A1 = np.sin(q2)*(d3 + d5*np.cos(q4) - a4*np.sin(q4)) + np.cos(q3)*(d5*np.sin(q4)*np.cos(q4) + d5*a4*np.cos(q4)*np.cos(q2) + a3*np.cos(q2)) - a2*np.cos(q2)
    #print(A1)
    #print("A1\n")
    B1 = np.sin(q3)*(a3 + a4*np.cos(q4) + d5*np.sin(q4))
    #print(B1)
    #print("B1\n")
    C1 = Xwsy
    #print(C1)
    #print("C1\n")
    D1 = pow(C2,2)-(pow(B2,2) + pow(A2,2))
    #print(D1)
    #print("D1\n")
    q1 = (np.angle(C1-cmath.sqrt(D1)) - np.angle(complex(B1,-A1)))*180/np.pi
    #print(q1)
    #print("q1\n")

    #return q5 to degrees; remove later
    q5 = q5*180/np.pi
    #print(q5)
    #print("q5\n")
    #return q6 to degrees; remove later
    q6 = q6*180/np.pi
    #print(q6)
    #print("q6\n")

    #return q3 to degrees
    q3 = q3*180/np.pi
    #print(q3)
    #print("q3\n")

    set_thetas = "<1.2.1.{}.{}.{}.{}.{}.{}>".format(int(q1),int(q2),int(q3),int(q4),int(q5),int(q6))
    print("{} {}".format(Tee,set_thetas))

    try:
        serial_port = serial.Serial(
            port="/dev/ttyTHS1",
            baudrate=2000000,
            bytesize=serial.EIGHTBITS,
            parity=serial.PARITY_NONE,
            stopbits=serial.STOPBITS_ONE,
        )
        serial_port.write((f'{set_thetas}').encode())
        serial_port.close()
    except SerialException:
        print("Warning: Couldn't open serial port; command not sent.\n")
    finally:
        pass

    #go back to rad for next iter; remove later
    q5 = q5*np.pi/180
    #print(q5)
    #print("q5\n")
    #go back to rad for next iter; remove later
    q6 = q6*np.pi/180
    #print(q6)
    #print("q6\n")

    #go back to rad for next iter; remove later
    q3 = q3*np.pi/180
    #print(q3)
    #print("q3\n")