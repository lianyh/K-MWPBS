from argparse import ArgumentParser
from signal import signal, SIGPIPE, SIG_DFL
from gurobipy import *
import networkx as nx
from networkx.algorithms import bipartite
from networkx.classes.function import *
import collections
import math 
import random

def first2(s):
    return s[:2]

def removeEdges(G,new_G):
	for i,j in new_G.edges():
		if G.has_edge(i,j):
			G.remove_edge(i, j)
		if G.has_edge(j,i):
			G.remove_edge(j, i)

	return G

def runILP(G,pathwayDict,pathwayIDName):
	try:
		# Create a new model
		model = Model("mip1")
	
		# Constant D
		D=4
		x={}
		y={}
		z={}
		w={}
		e={}
		r={}
		l={}
		B=len(G)
		m=len(pathwayDict)
		q={}
		t={}

		# Create variables
		nodes=G.nodes()
		for n in nodes:
			x[n] = model.addVar(vtype=GRB.INTEGER, lb=-1,ub=1,name="x_%d"%n)
			r[n] = model.addVar(vtype=GRB.BINARY,name="r_%d"%(n)) 
			l[n] = model.addVar(vtype=GRB.BINARY,name="l_%d"%(n)) 

			for o in nodes:
				if n !=o:
					#add var yij
					y[n,o]=model.addVar(vtype=GRB.INTEGER, lb=-2,ub=2,name="y_%d_%d"%(n,o)) 
					
					#add var zij
					z[n,o]= model.addVar(vtype=GRB.BINARY,name="z_%d_%d"%(n,o)) 
					
					#add var eij
					e[n,o] = model.addVar(vtype=GRB.BINARY,name="e_%d_%d"%(n,o))

					w[n,o]=0

					if G.has_edge(n,o):
						w[n,o]=G[n][o]['weight']
						w[o,n]=G[n][o]['weight']

		for pathway in pathwayDict:
			q[pathway]= model.addVar(vtype=GRB.BINARY,name="q_%s"%(pathway)) 
			t[pathway]= model.addVar(vtype=GRB.BINARY,name="t_%s"%(pathway)) 
		
		#add constraint ==
		for n in nodes:
			model.addConstr( r[n] - l[n], GRB.EQUAL, x[n],name="r_min_r"+str(n)+"_l"+str(n)+"_constr")
			model.addConstr( r[n] + l[n] <=1,name="r_add_r"+str(n)+"_l"+str(n)+"_constr")
			for o in nodes:
				if n != o:
					model.addConstr( x[n] - x[o],GRB.EQUAL,y[n,o],name="x_"+str(n)+"_"+str(o)+"_constr")
					model.addConstr(y[n,o] <= (D*(1-z[n,o])) + 1,name="y_lt_"+str(n)+"_"+str(o)+"_constr")
					model.addConstr(e[n,o] <= D*(1-z[n,o]),name="e_"+str(n)+"_"+str(o)+"_lt_constr")
					model.addConstr(e[n,o] >= 1 - D*(z[n,o]),name="e_"+str(n)+"_"+str(o)+"_gt_constr")
					model.addConstr(y[n,o] >= 2 - D*(z[n,o]),name="y_gt_"+str(n)+"_"+str(o)+"_constr")

		for pathway in pathwayDict:
			model.addConstr((quicksum(l[i] for i in x) - quicksum(l[i]  for i in x if i in pathwayDict[pathway])) <= B*q[pathway],name="l_%s_pathway_sum_constr"%(pathway))
			model.addConstr((quicksum(r[i] for i in x) - quicksum(r[i]  for i in x if i in pathwayDict[pathway])) <= B*t[pathway],name="r_%s_pathway_sum_constr"%(pathway))
			#make sure pathway exist in either of the partite, not in bipartite
			model.addConstr(q[pathway] + t[pathway] >= 1,name="q_t_not_equal_%s_pathway"%(pathway))

		model.addConstr(quicksum(q[i] for i in pathwayDict) <= m - 1,name="q_quicksum_lt_m_min1_constr")
		model.addConstr(quicksum(t[i] for i in pathwayDict) <= m - 1,name="t_quicksum_lt_m_min1_constr")

		#set objective
		objective =  quicksum(w[i,j] * e[i,j] for i,j in y)
		model.ModelSense = GRB.MAXIMIZE
		model.setObjective(objective)

		model.optimize()
		model.write("work_schedule.lp")
		#initialize empty graph to store results
		B = nx.Graph()
		leftn=[]
		rightn=[]
		l_edges=[]
		lcolor=""
		rcolor=""
		for v in model.getVars():
			if first2(v.varName) == "e_" and v.x==1:
				ss=v.varName.split("_")
				n=int(ss[1])
				rn=r[n].x
				o=int(ss[2])
				lo=l[o].x
				result=w[n,o] * v.x
				if result >= 1:
					if rn > 0.0:
						rightn.append(n)
					if lo > 0.0:
						leftn.append(o)
					if rn > 0.0 and lo > 0.0:
						l_edges.append(tuple([n, o]))
						#print("e="+v.varName+",v.x="+str(v.x)+";weight="+str(w[n,o])+",r="+str(r))
			if (first2(v.varName) == "q_" and v.x==0) or (first2(v.varName) == "t_" and v.x==0):
				print("pathway ID="+str(v.varName)+",pathway value="+str(v.x))
				idh = v.varName[2:]
				print("pathway name="+str(pathwayIDName[idh]))
				if first2(v.varName) == "q_":
					lcolor=str(v.varName)
				if first2(v.varName) == "t_":
					rcolor=str(v.varName)
	

		B.add_nodes_from(rightn, bipartite=0,color=rcolor)
		B.add_nodes_from(leftn, bipartite=1,color=lcolor)
		B.add_edges_from(l_edges)
			
		for (u, v) in B.edges:
			if B.has_edge(u,v):
				B[u][v]['weight'] = 1

		#print("----------------------- Bipartite Score ----------------------")
		#print('Score:', model.objVal)
		#print("--------------------------------------------------------------------")
		return B
	except GurobiError:
		print('Error reported')


def mkArgParser():
	parser = ArgumentParser()

	# Positional arguments, i.e. arguments that HAVE to be provided
	parser.add_argument("n", help="n The number of nodes in the first bipartite set")
	parser.add_argument("m", help="m The number of nodes in the second bipartite set.")
	parser.add_argument("p", help="p Probability for edge creation, i.e. 0.2")
	parser.add_argument("k", help="k add k new edges to the graph G, i.e.  2")
	parser.add_argument("j", help="j new nodes")

	return parser

def subRoutine(n,m,p,k,x):

	weight_total_G=0	
	#simulated small known bipartite graph with node attribute as color
	G = bipartite.random_graph(n, m, p)
	G_top = {n for n, d in G.nodes(data=True) if d['bipartite']==0}
	G_bottom = set(G) - G_top
	#print("************G_top="+str(G_top))
	#print("************G_top="+str(G_to))

	for (u, v) in G.edges:
		if G.has_edge(u,v):
			G[u][v]['weight'] = 1

	number_of_nodes_G=len(G)
	edges_total_G =G.size(weight='weight')

	H =  nx.Graph()
	H.add_nodes_from(G)
	for i,j in G.edges():
		if G.has_edge(i,j):
			H.add_edge(i, j)
			H[i][j]['weight'] = 1

	#add j new nodes color
	colors=['red','yellow','black']
	t_nm=n+m
	for c in colors:
		for i in range(0,x):
			H.add_node(t_nm,color=c)
			#print("*********** add node t_m="+str(t_nm))
			t_nm += 1
	
	for i in H.nodes():
		if i in G_top:
			H.node[i]['color']='red'
		if i in G_bottom:
			H.node[i]['color']='yellow'


	#add new edges k
	k_counter=0
	while k_counter < k:
		u = random.choice(range(m))
		v = random.choice(range(m))
		if u != v and not H.has_edge(u,v):
			H.add_edge(u, v,weight=1.0)
			k_counter += 1

	number_of_nodes_H=len(H)
	edges_total_H =H.size(weight='weight')

	#get 3 color pathway dict list
	pathwayDict={}
	pathwayIDName={}
	pathwayIDName={'pathway_red': 'pathway_red', 'pathway_yellow': 'pathway_yellow', 'pathway_black': 'pathway_black'}
	pathwayDict['pathway_red']=set()
	pathwayDict['pathway_yellow']=set()
	pathwayDict['pathway_black']=set()

	for p in H.nodes():
		if H.node[p]['color']=='red':
			pathwayDict['pathway_red'].add(int(p))
		if H.node[p]['color']=='yellow':
			pathwayDict['pathway_yellow'].add(int(p))
		if H.node[p]['color']=='black':
			pathwayDict['pathway_black'].add(int(p))
	

	isbipartite=0
	#run ILP
	if  H.number_of_edges() > 0:
		#RUN ILP
		B=runILP(H,pathwayDict,pathwayIDName)
		number_of_nodes_B=len(B)
		number_of_edges_B=B.size()
		bottom_nodes = {n for n, d in B.nodes(data=True) if d['bipartite']==0}
		top_nodes = set(B) - bottom_nodes
		lcolor=[]
		rcolor=[]

		for i in bottom_nodes:
			rcolor.append(B.node[i]['color'])

		for i in top_nodes:
			lcolor.append(B.node[i]['color'])

		if is_empty(B):
			print("B is empty, no results!")
		if bipartite.is_bipartite(B):
			isbipartite=1

	return number_of_nodes_G, edges_total_G, number_of_nodes_H, edges_total_H, number_of_nodes_B, number_of_edges_B,lcolor,rcolor,isbipartite

if __name__ == '__main__':
	signal(SIGPIPE,SIG_DFL)
	args = mkArgParser().parse_args()
	n=int(args.n)
	m=int(args.m)
	p=float(args.p)
	k=int(args.k)
	j=int(args.j)
	x=int(j/3)

	file1 = open("simulateres","a")
	for i in range(0,10):
		num_node_G,num_edge_G,num_node_H,num_edge_H,num_node_ILP,num_edge_ILP,lcolor,rcolor,isbipartite=subRoutine(n,m,p,k,x)
		file1.write(str(n)+"\t"+str(m)+"\t"+str(p)+"\t"+str(k)+"\t"+str(j)+"\t"+str(num_node_G)+"\t"+str(num_edge_G)+"\t"+str(num_node_H)+"\t"+str(num_edge_H)+"\t"+str(num_node_ILP)+"\t"+str(num_edge_ILP)+"\t"+str(lcolor)+"\t"+str(rcolor)+"\t"+str(isbipartite)+"\n")

	file1.close() 
