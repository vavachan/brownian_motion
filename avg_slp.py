import numpy as np
infile=open("slopes.txt","r")
list=infile.readlines()
sum=np.zeros(5)
for i in range (0,5,1):
	for j in range (i,200,5):
		sum[i]=sum[i]+float(list[j])
for i in range (0,5):
	print i+1,sum[i]/40
