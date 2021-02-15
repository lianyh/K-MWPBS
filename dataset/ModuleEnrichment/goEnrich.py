from argparse import ArgumentParser
from signal import signal, SIGPIPE, SIG_DFL
from gurobipy import *
import networkx as nx
from networkx.algorithms import bipartite
from networkx.classes.function import *
import collections
import math 
import random
import subprocess
import time
import os

def first2(s):
    return s[:2]

def removeEdges(G,new_G):
	for i,j in new_G.edges():
		if G.has_edge(i,j):
			G.remove_edge(i, j)
		if G.has_edge(j,i):
			G.remove_edge(j, i)

	return G


def mkArgParser():
	parser = ArgumentParser()

	# Positional arguments, i.e. arguments that HAVE to be provided
	parser.add_argument("out", help="out folder")
	parser.add_argument("orig", help="original bpm results,.i.e topbottomnodes")
	parser.add_argument("s", help="s minimum number of genes in a module")
	parser.add_argument("b", help="b maximum number of genes in a module")

	return parser

if __name__ == '__main__':
	signal(SIGPIPE,SIG_DFL)
	args = mkArgParser().parse_args()
	out=str(args.out)
	s=int(args.s)
	b=int(args.b)
	
	module1=collections.OrderedDict()
	module2=collections.OrderedDict()
	module1GO=collections.OrderedDict()
	module2GO=collections.OrderedDict()
	BPMEnrichedSameFunc=collections.OrderedDict()
	BPMEnrichedDiffFunc=collections.OrderedDict()

	totalModuleSize=0
	idx1=0
	idx2=0
	with open(args.orig, "r") as f:
		for line in f:
			line_ = line.strip()
			if "Mod1" in line_:
				hhh=line_.split('\t')
				h2=hhh[1:]
				sizeh = len(h2)
				if sizeh >= s and sizeh <= b:
					module1[idx1]=sizeh
				idx1 +=1

			elif "Mod2" in line_:
				hhh=line_.split('\t')
				h2=hhh[1:]
				sizeh = len(h2)
				if sizeh >= s and sizeh <= b:
					module2[idx2]=sizeh
				idx2 +=1

	totalModuleSize=len(module1) + len(module2)
	print("Size Module 1="+str(len(module1)))
	print("Size Module 2="+str(len(module2)))
	print("Mod1="+str(module1.keys()))
	print("Mod2="+str(module2.keys()))
	print("##########################################################")
	print("")

	fileList = os.listdir(out)
	for i in fileList:
		if "module1" in i:
			with open(out+i, "r") as f:
				idx=i.split('_')[1]
				if int(idx) in module1:
					for line in f:
						line_ = line.strip()
						hhh=line_.split(',')
						go = hhh[6].strip()
						if not int(idx) in module1GO:
							module1GO[int(idx)]=[]

						if not go in module1GO[int(idx)]:
							module1GO[int(idx)].append(go)

		
		elif "module2" in i:
			with open(out+i, "r") as f:
				idx=i.split('_')[1]
				if int(idx) in module2:
					for line in f:
						line_ = line.strip()
						hhh=line_.split(',')
						go = hhh[6].strip()

						if not int(idx) in module2GO:
							module2GO[int(idx)]=[]

						if not go in module2GO[int(idx)]:
							module2GO[int(idx)].append(go)

	totalBPM=0
	for k in module1:
		if k in module1GO and k in module2GO:
			list1=module1GO[k]
			list2=module2GO[k]
			totalBPM +=1

			matchSize=len([w for w in list2 if w in list1])
			if matchSize > 0:
				if not k in BPMEnrichedSameFunc:
					BPMEnrichedSameFunc[k]=[w for w in list2 if w in list1]
			else:
				if not k in BPMEnrichedDiffFunc:
					BPMEnrichedDiffFunc[k]=1

	totalGO=0
	tlist=[]
	for a in module1GO:
		h=module1GO[a]
		for j in h:
			tlist.append(j)

	for b in module2GO:
		h=module2GO[b]
		for j in h:
			tlist.append(j)

	myset = set(tlist)

	print("Total Module Size (a)="+str(totalModuleSize))

	totalEnrichedModuleSize=len(module1GO) + len(module2GO)
	print("Go Annotation (b), Unique ENRICHED in Modules="+str(len(myset)))
	print("Total ENRICHED Module Size (c)="+str(totalEnrichedModuleSize))
	print("Total BPM="+str(totalBPM))
	print("")
