import socket
import numpy as np

UDP_IP = "192.168.1.128"
UDP_PORT = 20777
sock = socket.socket(	socket.AF_INET, # Internet
			socket.SOCK_DGRAM) # UDP
sock.bind((UDP_IP, UDP_PORT))

i = 0
while True:
	data, addr = sock.recvfrom(65536)
	text = data.decode("ascii")
	print("received message:", text)
	if i > 1:
		Rtxt = text.split("\n")
		print("list:", Rtxt)
		R70 = np.array([[Rtxt[0],Rtxt[1],Rtxt[2]],[Rtxt[4],Rtxt[5],Rtxt[6]],[Rtxt[8],Rtxt[9],Rtxt[10]]]).astype(np.float)
		X70 = np.array([[Rtxt[3]],[Rtxt[7]],[Rtxt[11]]]).astype(np.float)
		print(R70)
		print(X70)
	i += 1
