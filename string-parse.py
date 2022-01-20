import numpy as np
text = "[-2.96174414467395e-17 8.1373251596936e-17 1 254.27 -4.93879172355211e-16 1 -8.1373251596936e-17 -224.795929712244 -1 -4.93879172355211e-16 -2.96174414467394e-17 -185.282346770381 0]"
text = text[1:len(text)-1]
print("message:", text)
mylist = text.split()
print("list:", mylist)
R70 = np.array([[mylist[0],mylist[1],mylist[2]],[mylist[4],mylist[5],mylist[6]],[mylist[8],mylist[9],mylist[10]]]).astype(np.float)
X70 = np.array([[mylist[3]],[mylist[7]],[mylist[11]]]).astype(np.float)
print(R70)
print(X70)
