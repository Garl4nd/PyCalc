# prime factorization using Pollard's rho algorithm
import math
from operator import mul
import functools as ft
primes=[]
nonprimes=set()

def sieve(n):
    #lb=int(math.sqrt(n))
#    cdef int lb,x,i
    lb=n
    #print(lb)
    for x in range(2,lb):
        if x in nonprimes:
            continue
        else:
            primes.append(x)
            for i in range(x**2,lb,x):
                

                nonprimes.add(i)
    return primes
def prime_fact(n):
    if n<=0 or isinstance(n,float):
        raise ValueError("The number must be a positive integer!")
    #candidates=(num for num in sieve(n//2+1))
    lb=int(math.sqrt(n))+1
    candidates=(num for num in sieve(lb))
    primes=[]
    for candidate in candidates:
        #print(candidate)
        while not n%candidate:
           # print("match",candidate)
            n=n//candidate
            primes.append(candidate)
            if n==1:
                return primes
    return primes+[n]

def combs(nums,k):
	if k==0:
		yield [],nums
	else:
		for i,el in enumerate(nums):
			rest=nums[i+1:]
			old=nums[:i]
			for res,compl in combs(rest,k-1):
				yield [el]+res,old+compl

def groupings(nums):
	def _groups(nums,start=1,forb=tuple()):
		if not nums:
			yield []
		for gn,_ in enumerate(nums,start):
			for chosen,rest in combs(nums,gn):
				if chosen in forb:
					continue
				for res in _groups(rest,gn,forb):
					yield [chosen]+res
				forb=forb+(chosen,)

	yield from _groups(nums)

def mulgroups(nums):
	for grouping in groupings(nums):
		yield [ft.reduce(mul,group) for group in grouping]

def sortedmulgroups(nums):
	return sorted(mulgroups(nums),key=len,reverse=True)
#print(sieve(546545))
#print(prime_fact(54654546789944))