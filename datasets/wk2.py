from wk1 import *

def MinimumSkew(str):
	skew=[0]*(len(str)+1)

	for k in range(0,len(str)):
		if (str[k]== 'G'): skew[k+1]=skew[k]+1
		elif (str[k] == 'C'): skew[k+1]=skew[k]-1
		else: skew[k+1]=skew[k]

	return skew

def min_index(a):
	minv=min(a)
	ind=[]
	for i in range(0,len(a)):
		if a[i]==minv: ind.append(i)	

			
	return ind

def hammingdist(str1,str2):
	ind=[]
	for i in range(0,len(str1)):
		if(str1[i]!=str2[i]): ind.append(i)
	
	return len(ind)

def ApproxPattern(pattern,text,d):
	ind=[]
	for i in range(0,len(text)-len(pattern)+1):
		str=text[i:i+len(pattern)]
		if hammingdist(str,pattern)<=d:
			ind.append(i)
	return ind

def ApproximatePatternCount(text,pattern,d):
	ans=ApproxPattern(pattern,text,d)
	count=[]
	for i in ans:
		count.append(text[i:len(pattern)+i])

	return count

def ImmediateNeighbors(pattern):
	neighborhood=[]
	nucleotides=['A','T','C','G']

	for i in range(0,len(pattern)):
		symbol=pattern[i]
		for j in nucleotides:
			if j!=symbol: 
				neighbor=pattern[0:i]+j+pattern[i+1:]
				neighborhood.append(neighbor)
	
	return neighborhood		

def Neighbors(pattern,d):
	nucleotides=['C','T','G','A']

	if d==0: return {pattern}
	if len(pattern)==1:
		return {'A','C','G','T'}
	Neighborhood=[]
	SuffixNeighbors=Neighbors(pattern[1:],d)
	for suff_text in SuffixNeighbors:
		if hammingdist(pattern[1:],suff_text)<d:
			for x in nucleotides:
				Neighborhood.append(x+suff_text)
		else:
			Neighborhood.append(pattern[0]+suff_text) 	
	
	return Neighborhood

#problem may be with ApproximatePatternCount()
def FrequentWordsWithMismatches(text,k,d):
	frequent_patterns=[]
	close=[0]*4**k
	frequency_arr=[0]*4**k

	for i in range(0,len(text)-k+1):
		neighborhood=Neighbors(text[i:i+k],d)
		
		for p in neighborhood:
			index=PatternToNumber(p)
			close[index]=1
	
	for i in range(0,4**k):
		if close[i]==1:
			pattern=NumberToPattern(i,k)
			frequency_arr[i]=ApproximatePatternCount(text,p,d)
	maxCount=max(frequency_arr)
	print maxCount

	for i in range(0,4**k):		
		if frequency_arr[i]==maxCount:
			pattern=NumberToPattern(i,k)
			frequent_patterns.append(pattern)

	return frequency_arr, frequent_patterns 	
