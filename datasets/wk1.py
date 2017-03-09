def readtext(filename):
	file=open(filename)
	str1=file.readline()
	str2=file.readline()
	return str1[0:-1], str2[0:-1]	

def PatternCount(text, pattern):
	count=[]
	for i in range(0,len(text)-len(pattern)+1):
		if(text[i:i+len(pattern)]==pattern):
			count.append(i)
	return count

def FrequentWords(text,k):
	#frequent patterns
	fp=[]
	count=range(0,len(text)-k+1);

	for i in range(0,len(count)):
		pattern=text[i:i+k]
		count[i]=len(PatternCount(text,pattern))
	maxCount=max(count)
	for i in range(0,len(count)):
		if count[i]==maxCount:
			fp.append(text[i:i+k])
	
	return list(set(fp))
	
def RevComp(text):
	rev={'A':'T', 'C':'G', 'T':'A', 'G':'C'}
	revcomp=range(len(text))

	for i in range(0,len(text)):
		revcomp[-(i+1)]=rev[text[i]]

	return ''.join(revcomp) 	

def PatternToNumber(Pattern):
	symbol_num={'A':0, 'C':1, 'G':2, 'T':3}
	num=0
	ind=0

	for i in Pattern:
		num+=symbol_num[i]*4**(len(Pattern)-ind-1)
		ind+=1
	
	return num	

def NumberToPattern(num,k):
	num_symbol={0:'A', 1:'C', 2:'G', 3:'T'}

	exponent=0
	while(num/4**exponent!=0):
		exponent+=1
	if(k>=exponent): pattern=['A']*k
	else: pattern=['A']*exponent
	for i in range(exponent-1,-1,-1):
		n=num/4**i
		m=num%4**i
		pattern[-1-i]=num_symbol[n]
		num=m
	return ''.join(pattern)

def ComputingFrequencies(text, k):
	frequencyarr=range(0,4**k)

	for i in range(0,4**k):
		frequencyarr[i]=0
	
	for i in range(0,len(text)-k+1):
		pattern=text[i:i+k]
		j=PatternToNumber(pattern)
		frequencyarr[j]+=1
	
	return frequencyarr



def ClumpFinding(genome,k,L,t):
	frequencypatterns=[]
	clump=range(0,4**k)

	for i in range(0,4**k):
		clump[i]=0
	
	text=genome[0:L]
	frequencyarr=ComputingFrequencies(text,k)
		
	for i in range(0,4**k):
		if frequencyarr[i]>=t:
			clump[i]=1

	for i in range(1, len(genome)-L+1):
		first_pattern=genome[i-1:i-1+k]
		index=PatternToNumber(first_pattern)
		frequencyarr[index]=frequencyarr[index]-1
		last_pattern=genome[i+L-k:i+L]
		index=PatternToNumber(last_pattern)
		frequencyarr[index]=frequencyarr[index]+1
		if frequencyarr[index]>=t:
			clump[index]=1
	
	for i in range(0,4**k):
		if clump[i]==1:
			pattern=NumberToPattern(i,k)
			frequencypatterns.append(pattern)
	return frequencypatterns	


# make array into string
# s=[]
# for i in farr: s.append(str(i))

# with open(filename, 'w') as f:
#	f.write(str)
# f.close()
