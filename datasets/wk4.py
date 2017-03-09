def FarthestFirstTraversal(data,k,dim):
	# convert str to data
	for i in range(0,len(data)):
		data[i]=data[i].split()
		for j in range(0, dim):
			data[i][j]=float(data[i][j])

	centers=[data[0]]

	while(len(centers)<k):
		datapoint=0
		dist_data2center=[0]*len(data)
			
		# max d(data, center)
		for d in range(0,len(data)):
			dist_tmp=[0]*len(centers)
			for c in range(0,len(centers)):
				for m in range(0,dim):
					dist_tmp[c]+=(data[d][m]-centers[c][m])**2
			# find the min dist first
			dist_data2center[d]=min(dist_tmp)	

		max_dist=max(dist_data2center)
		datapoint=dist_data2center.index(max_dist)
		# print data[datapoint]
		centers.append(data[datapoint])

	return centers


def SquaredError(data, centers, dim):
	sqErr=0
	for d in range(0,len(data)):

		dist_tmp=[0]*len(centers)
		for c in range(0,len(centers)):
			for m in range(0,dim): 
				dist_tmp[c]+=(data[d][m]-centers[c][m])**2		
		sqErr+=min(dist_tmp)
	
	sqErr=sqErr/len(data)
	return sqErr

# to test:
# from scipy import cluster
# cluster.vq.kmeans2(numpy.array(data),k,iter)

import numpy
def kmeanClustering(data,k,dim,niter):
	centers=data[0:k]
	centers_tmp=centers[:]
	centers2data={}
	for c in range(0,k): centers2data.setdefault(c, [])
	data2centers=[0]*len(data)

	#sqerr=SquaredError(data,centers,dim)
	#sqerr_tmp=sqerr

	# centers to clusters initialize
	for iter in range(0,niter):
		for c in range(0,k): centers2data[c]=[]
		sqerr=0
		for d in range(0,len(data)):
			dist_tmp=[0]*k
			for c in range(0,k):			
				for m in range(0,dim): dist_tmp[c]+=(data[d][m]-centers[c][m])**2
			min_dist_ind=dist_tmp.index(min(dist_tmp))
			sqerr+=min(dist_tmp)	
			centers2data[min_dist_ind].append(data[d])
			data2centers[d]=min_dist_ind

		for c in range(0,k):
			centers_tmp[c]=numpy.mean(centers2data[c],0)
		
		sqerr_tmp=0		
		for d in range(0,len(data)):
			for m in range(0,dim): sqerr_tmp+=(data[d][m]-centers_tmp[data2centers[d]][m])**2
	 	# print sqerr, sqerr_tmp
		# sqerr_tmp=SquaredError(data,centers_tmp,dim)
	
		if sqerr_tmp<sqerr: 
			centers=centers_tmp[:]
			sqerr=sqerr_tmp
	return centers #, centers2data
	


def list2str(data_list):
	s=''
	for i in data_list:
		if(len(i)>1):
			for j in i: s+=str(j)+' '
		else: s+=str(i)+' '
		s+='\n'

def str2list(str_list):
	for i in range(0, len(str_list)):
		str_list[i]=str_list[i].split()
		for j in range(0, len(str_list[i])):
			str_list[i][j]=float(str_list[i][j])
	return str_list


'''
s=''
for i in ans[0]:
	for j in i: s+='%.2f' %j + ' '
	s+='\n'


'''
# https://datasciencelab.wordpress.com/2013/12/12/clustering-with-k-means-in-python/
def cluster_points(X,mu):
	clusters={}
	for x in X:
		bestmukey=min([(i[0], numpy.linalg.norm(numpy.array(x)-numpy.array(mu[i[0]]))) \
			for i in enumerate(mu)], key=lambda t:t[1])[0]
		try:
			clusters[bestmukey].append(x)
		except KeyError:
			clusters[bestmukey]=[x]
	return clusters

def reevaluate_centers(mu, clusters):
	newmu=[]
	keys=sorted(clusters.keys())
	for k in keys:
		newmu.append(numpy.mean(clusters[k], axis=0))
	return newmu

def has_converged(mu,oldmu):
	return set([tuple(a) for a in mu]) == set([tuple(a) for a in oldmu])

import random
def find_centers(X,K):
	#oldmu=random.sample(X,K)
	mu=random.sample(X,K)
	oldmu=X[0:K]
	# mu=X[0:K]
	while not has_converged(mu,oldmu):
		oldmu=mu
		clusters=cluster_points(X,mu)
		mu=reevaluate_centers(oldmu,clusters)
	return mu
 
