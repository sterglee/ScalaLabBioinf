from math import *

def softkmeans(data_input,center,k,dim,beta,iter):
	import numpy as np

	hm=np.zeros((k,len(data_input)))
	data=np.array(data_input)
	# center=data[0:k]
	center=np.array(center)
	
	for n in range(0, iter):
		for i in range(0,len(center)):
			for j in range(0,len(data)):
				diff=data[j]-center[i]
				dist=np.sqrt(np.dot(diff,diff))
				# dist=np.dot(diff,diff)	
				hm[i][j]=exp(-beta*dist)	
		hm=hm/np.sum(hm,0)		
		# update center
		for i in range(0,len(center)):
			center[i]=np.dot(hm[i],data)/np.sum(hm[i])

	return hm, center


# mt is a distance matrix. row_n means distances between the element n relative to all other element. D(i,j) means distance between element i and j
def HierarchCluster(mt):
	import numpy as np
	mtrx=np.array(mt)	# distance matrix
	'''
	grouping={}
	for i in range(1,len(mt)+1): 		
		grouping.setdefault(i,[])
	'''
	mtrx_init=np.copy(mtrx)				# original distance matrix
	gene_list=range(1,len(mtrx_init)+1)
		
	ans=[]
	while(len(gene_list)>1):
		# if(len(n_mtrx)>0): mtrx=n_mtrx
		mtrx_st=np.transpose(np.sort(mtrx,1))		#sorted, transposed
		minv=min(mtrx_st[1])		
		# print minv

		# finding minimal distance value
		for r in range(0, len(mtrx)):
			r_break=0
			for c in range(0,len(mtrx[r])):
				if mtrx[r][c]==minv:
					if r<c: 
						row=r;col=c
						r_break=1
						break
			if(r_break==1): break	
	
		# joining genes
		# print "row:%d,col:%d" %(row,col)
		if isinstance(gene_list[row],list): 
			if isinstance(gene_list[col], list):
				gene_list[row].extend(gene_list[col])
			else:
				gene_list[row].append(gene_list[col])
		else: 
			if isinstance(gene_list[col],list):
				gene_list[col].insert(0,gene_list[row])
				gene_list[row]=gene_list[col]
			else: gene_list[row]=[gene_list[row],gene_list[col]]
		del gene_list[col]

		# merge two clusters
		# calclulating D(Cnew, C) for each C in cluster	
		
		# print gene_list
		mtrx=np.delete(mtrx,col,1) 		# delete the column
		mtrx=np.delete(mtrx,col,0)  		# delete the corresponding row
		# using Davg(C1,C2)
		for cl in range(0, len(gene_list)):
			C1=gene_list[row]

			if cl!=row: C2=gene_list[cl]
			else: continue	
			# print "C1:", C1, ", cl:", cl, ", C2: ", C2
			mtrx[row][cl]=clusterDist(C1,C2,mtrx_init)
			# mtrx[col][cl]=clusterDist(C1,C2,mtrx_init)
			mtrx[cl][row]=mtrx[row][cl]
			# mtrx[cl][col]=mtrx[col][cl]
			mtrx[cl][cl]=0
		
		mtrx[row][row]=0
		# mtrx[col][col]=0

		# print mtrx

		# answer output
		print gene_list[row]
		for gl in gene_list[row]:
			ans.append(str(gl))
			ans.append(' ')
		del ans[-1]
		ans.append('\n')
	

	return ''.join(ans)
	
def flatten(matrx):
	m=[0]*len(matrx)*len(matrx[0])
	for i in range(0,len(matrx)):
		for j in range(0,len(matrx[i])):
			m[i*len(matrx)+j]=matrx[i][j]

	return m

def clusterDist(C1,C2,Dmtrx):
	if not isinstance(C1,list): C1=[C1]
	if not isinstance(C2,list): C2=[C2]
	# print "C1: ", C1,", C2: ", C2
	dist=0
	for i in range(0,len(C1)):
		pt1=C1[i]-1				# point in C1
		for j in range(0,len(C2)):
			pt2=C2[j]-1			# point in C2
			# print "pt1: ", pt1, ", pt2: ", pt2
			dist+=Dmtrx[pt1][pt2]
			
	dist=dist/len(C1)/len(C2)
	return dist
