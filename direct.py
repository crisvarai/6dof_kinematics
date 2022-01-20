import numpy

q1=0*numpy.pi/180
q2=115*numpy.pi/180
q3=90*numpy.pi/180
q4=0*numpy.pi/180
q5=-90*numpy.pi/180
q6=0*numpy.pi/180

d1=183.73
a2=58.28
d3=244.4
a3=13.76+5
a4=18.53 
d5=78.53+44.5
d7=150
k=((65*numpy.pi)/180)

def Denavit_Hartenberg(q, d, a, alpha):
    Rz = numpy.array([[numpy.cos(q), -numpy.sin(q), 0, 0], [numpy.sin(q), numpy.cos(q), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    Tr = numpy.array([[1, 0, 0, a], [0, 1, 0, 0], [0, 0, 1, d], [0, 0, 0, 1]])
    Rx = numpy.array([[1, 0, 0, 0], [0, numpy.cos(alpha), -numpy.sin(alpha), 0], [0, numpy.sin(alpha), numpy.cos(alpha), 0], [0, 0, 0, 1]])
    return numpy.matmul(numpy.matmul(Rz, Tr), Rx)

T0B=Denavit_Hartenberg((numpy.pi/2), 0, 0, k);
T0H=Denavit_Hartenberg(-numpy.pi/2, 0, 0, 0);
T01=Denavit_Hartenberg(q1, d1, 0, -numpy.pi/2);
T12=Denavit_Hartenberg(q2, 0, -a2, numpy.pi/2);
T23=Denavit_Hartenberg(q3, d3, a3, -numpy.pi/2);
T34=Denavit_Hartenberg(q4, 0, a4, numpy.pi/2);
T45=Denavit_Hartenberg(q5, d5, 0, -numpy.pi/2);
T56=Denavit_Hartenberg(q6, 0, 0, numpy.pi/2);
T67=Denavit_Hartenberg(0, d7, 0, 0);

Tee=numpy.matmul(numpy.matmul(numpy.matmul(numpy.matmul(numpy.matmul(numpy.matmul(numpy.matmul(numpy.matmul(T0B,T0H),T01),T12),T23),T34),T45),T56),T67)
print(Tee)