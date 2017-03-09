from wk1 import *
from wk2 import *

def Composition(text,k):
	for i in range(0, len(text)-k+1): print text[i:i+k]
	# return sorted(compos_str)
'''
def GenomePath(text_list):
	ans=''

	k_heads=['x']*len(text_list)
	k_tails=['x']*len(text_list)
	
	for(i=0;i<len(text_list);i++):
		k_heads[i]=text_list[i][0:k]
		k_tails[i]=text_list[i][-k:]	
	
	for(i=0;i<len(text_list);i++):
		try: 
			j=k_tails.index(k_heads[i])
			ans+=
		except: 
		# this should be the sequence to start with
			
'''

def GenomePath(textlist):
	ans=textlist[0]
	k=len(textlist[0])-1

	for i in range(1, len(textlist)):
		ans+=textlist[i][-1]

	return ans


def ReadTextList(filename):
	text_list=[]
	with open(filename, 'r') as f:
		for line in f:
			if line[-1]=='\r' or line[-1]=='\n': 
				line=line[:-2]	
			text_list.append(line)
	
	f.close()
	return text_list

# to generate kmer_list
# kmer_list=[]
# for i in range(0,16): kmer_list.append(bin[i][2:].zfill(4))

def OverlapGraph(kmer_list):
	K=len(kmer_list)
	prefix=['x']*K
	suffix=['x']*K

	for i in range(0, K):
		prefix[i]=kmer_list[i][0:-1]
		suffix[i]=kmer_list[i][1:]

	overlap_list=dict(zip(range(0,K),[0]*K))
	for i in range(0, K):
		overlap_list[i]=[]
		for j in range(0, K):
			if (prefix[j]==suffix[i]) & (j!=i):
				overlap_list[i].append(j)	
		
	return overlap_list

# Hamiltonian graph
def HGraph(dictionary):
	#i=dictionary.keys()[-1]  # start from 0 will enter infinite loop
	k=0
	while(k<len(dictionary.keys())):
		break_flag=False
		i=dictionary.keys()[k]
		gph=[i]
			
		while(len(gph)<len(dictionary.keys())):
			i=gph[-1]
			new_added=False
			for j in dictionary[i]:
				if j not in gph:
					gph.append(j)
					new_added=True
					break
				else: 
					new_added=False
						
				if j==dictionary[-1] & new_added==False:
					break_flag=True
					break
		
		if break_flag: 
			k+=1
			break
			
		print gph

	return		


# reload(wk3)

def DeBruijin(text,k):
	kmer_list=[]
	prefix=[]
	suffix=[]

	for i in range(0, len(text)-k+1):
		kmer_list.append(text[i:i+k])
		prefix.append(text[i:i+k-1])
		suffix.append(text[i+1:i+k])
		
	gph=dict(zip(set(prefix),[[]]*len(set(prefix))))
	
	for j in range(0,len(prefix)): 		
		if len(gph[prefix[j]])==0: gph[prefix[j]]=[suffix[j]]
		else: gph[prefix[j]].append(suffix[j])	
	ans=''
	for p in gph.keys():
		ans+=p + ' ->  ' + ','.join(gph[p]) + '\n'	

	return gph, ans

# kl for kmer list
def DeBruijin1(kl):
	prefix=[]
	suffix=[]
	
	for i in range(0,len(kl)):
		prefix.append(kl[i][0:-1])
		suffix.append(kl[i][1:])

	gph=dict(zip(set(prefix),[[]]*len(set(prefix))))
	
	for j in range(0,len(prefix)):
		if len(gph[prefix[j]])==0: gph[prefix[j]]=[suffix[j]]
		else: gph[prefix[j]].append(suffix[j])


	return gph

def kmerlst2str(kl):
	ans=kl[0]
	for i in range(1,len(kl)-1):
		ans+=kl[i][-1]

	return ans

def binarystrs(text,k):
	strlist=[]
	for i in range(0,len(text)-k+1):
		strlist.append(text[i:i+k])
	
	return strlist

# bioinformatictsII Wk2, with input as an adjacenylist
def EulerianCycle(adjlist,startnode):
	import random
	nodes=adjlist.keys()
	
	unexploredEdg=[]
	for n in nodes:
		if isinstance(adjlist[n], list):
			for i in adjlist[n]:
				unexploredEdg.append([n,i])
		else: unexploredEdg.append([n,adjlist[n]])
					
	n_ind=nodes.index(startnode)
	n=nodes[n_ind]
	edges=unexploredEdg[:]
	edge_trash=[]
	edge_recycle=[]
	cycle=[]
	cycle1=[]
		
	while(len(set(cycle1))<len(nodes)):

		if len(cycle)>0:
			cycle1.extend(cycle)
			# print cycle1
			edges.extend(edge_recycle)
			# print edge_recycle
			edge_recycle=[]
			cycle=[]
		# print len(edge_trash)
		while(len(edges)>0):
			if n in nodes:				# if the input is a unbalanced path
				nextNlist=adjlist[n]
			else: 
				if len(cycle)>1:
					cycle.pop() 
					n=cycle[-1]
					if len(edge_trash)>0: edge_recycle.append(edge_trash.pop())
				else: 
					n_ind=random.randint(0,len(nodes)-1)
					while(nodes[n_ind]==startnode):	
						n_ind=random.randint(0,len(nodes)-1)
					n=nodes[n_ind]
				nextNlist=adjlist[n]
			# print "line 197 n=",n, ", nextNlist:",  nextNlist
			for k in nextNlist:
				if [n,k] in edges:
					ind=edges.index([n,k])
					edge_trash.append([n,k])
					del edges[ind]
					if len(cycle)==0 and  len(cycle1)==0: 
						cycle.append(n)
					cycle.append(k)
					# print "line 205, n=",n, ", k=",k
					n=k
					break
				else:
					if k==nextNlist[-1]:
						# if len(cycle)==0: print n, cycle1
						if len(cycle)>1:
							cycle.pop()
							n=cycle[-1]
							# print "line 213 n:", n
							edge_recycle.append(edge_trash.pop())
							while(len(adjlist[n])==1):
								if len(cycle)>1:
									cycle.pop() 
									n=cycle[-1]
									# print "line 220, cycle,pop(),  n=", n
									edge_recycle.append(edge_trash.pop())
								else:	
									n_ind=random.randint(0,len(nodes)-1)
									# print "line 223", n_ind, nodes[n_ind]
									while(nodes[n_ind]==startnode):
										n_ind=random.randint(0,len(nodes)-1)
									n=nodes[n_ind]
						else:	
							n_ind=random.randint(0,len(nodes)-1)
							# print "line 229", n_ind, nodes[n_ind]
							while(nodes[n_ind]==startnode):
								n_ind=random.randint(0,len(nodes)-1)
								# print "line 232", n_ind, nodes[n_ind]
							n=nodes[n_ind]		
						break
					else: continue	
	# print len(unexploredEdg)				
	return cycle1	
	
def list2connection(clist):
	congraph=''
	for i in range(0, len(clist)-1):
		congraph+=str(clist[i])+'->'
	congraph+=str(clist[-1])
	return congraph

def connection2adjl(filename):
	keys=[]
	kvalues=[]

	con=ReadTextList(filename)
	for i in range(0,len(con)):
		c=con[i].split(' -> ')
		keys.append(int(c[0]))
		val=c[1].split(',')
		if i==len(con)-1: print c, val	
		for j in range(0,len(val)):
			val[j]=int(val[j])
			print val[j]	
		kvalues.append(val)

	adjlist=dict(zip(keys,kvalues))
	return adjlist
					

def EulerianPath(adjlist):
	# find unbalanced nodes
	nodes=adjlist.keys()
	in_deg=[0]*len(nodes)
	out_deg=[0]*len(nodes)
	values_flat=[]
	for n in nodes: values_flat.extend(adjlist[n])
	for n_ind in range(0,len(nodes)):
		out_deg[n_ind]=len(adjlist[nodes[n_ind]])
		in_deg[n_ind]=values_flat.count(nodes[n_ind])
	
	out_nodes=[]
	in_nodes=[]	
	for n_ind in range(0,len(nodes)):
		if out_deg[n_ind]>in_deg[n_ind]:
			out_nodes.append(nodes[n_ind])
		if in_deg[n_ind]>out_deg[n_ind]:
			in_nodes.append(nodes[n_ind])		

	for i in set(values_flat):
		if i not in nodes:
			adjlist[i]=[out_nodes[0]]
			endnode=i	
			print "endnode:", endnode
	return adjlist,out_nodes			
