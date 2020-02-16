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

def runILP(G):
	try:
		# Create a new model
		model = Model("mip1")
	
		# Constant M
		M=4
		x={}
		y={}
		z={}
		w={}
		e={}
		# Create variables
		nodes=G.nodes()
		for n in nodes:
			x[n] = model.addVar(vtype=GRB.INTEGER, lb=-1,ub=1,name="x_%d"%n)
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

		
		#add constraint ==
		for n in nodes:
			for o in nodes:
				if n != o:
					model.addConstr( x[n] - x[o],GRB.EQUAL,y[n,o],name="x_"+str(n)+"_"+str(o)+"_constr")
					model.addConstr(y[n,o] <= (M*(1-z[n,o])) + 1,name="y_lt_"+str(n)+"_"+str(o)+"_constr")
					model.addConstr(e[n,o] <= M*(1-z[n,o]),name="e_"+str(n)+"_"+str(o)+"_lt_constr")
					model.addConstr(e[n,o] >= 1 - M*(z[n,o]),name="e_"+str(n)+"_"+str(o)+"_gt_constr")
					model.addConstr(y[n,o] >= 2 - M*(z[n,o]),name="y_gt_"+str(n)+"_"+str(o)+"_constr")
		
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
		for v in model.getVars():
			if first2(v.varName) == "e_" and v.x==1:
				ss=v.varName.split("_")
				n=int(ss[1])
				o=int(ss[2])
				r=w[n,o] * v.x
				if r >= 1:
					leftn.append(n)
					rightn.append(o)
					l_edges.append(tuple([n, o]))
					#print("e="+v.varName+",v.x="+str(v.x)+";weight="+str(w[n,o])+",r="+str(r))
		
		B.add_nodes_from(leftn, bipartite=0)
		B.add_nodes_from(rightn, bipartite=1)
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

	return parser

def subRoutine(n,m,p,k):

	weight_total_G=0	
	#simulated small known bipartite graph
	G = bipartite.random_graph(n, m, p)
	G_top = {n for n, d in G.nodes(data=True) if d['bipartite']==0}
	G_bottom = set(G) - G_top
	for (u, v) in G.edges:
		if G.has_edge(u,v):
			#print('(%d,%d)' % (u,v),end =",")
			G[u][v]['weight'] = 1
	
	number_of_nodes_G=len(G)
	edges_total_G =G.size(weight='weight')

	H =  nx.Graph()
	H=G

	k_counter=0
	while k_counter < k:
		u = random.choice(range(m))
		v = random.choice(range(m))
		if u != v and not H.has_edge(u,v):
			H.add_edge(u, v,weight=1.0)
			k_counter += 1
			#print("u="+str(u)+",v="+str(v)+",weight="+str(H[u][v]['weight']))
	#print("*********************k_counter="+str(k_counter))

	number_of_nodes_H=len(H)
	edges_total_H =H.size(weight='weight')
	isbipartite=0
	#run ILP
	if  H.number_of_edges() > 0:
		#RUN ILP
		B=runILP(H)
		number_of_nodes_B=len(B)
		number_of_edges_B=B.size()

		if is_empty(B):
			print("B is empty, no results!")
		if bipartite.is_bipartite(B):
			isbipartite=1

	return number_of_nodes_G, edges_total_G, number_of_nodes_H, edges_total_H, number_of_nodes_B, number_of_edges_B,isbipartite

if __name__ == '__main__':
	signal(SIGPIPE,SIG_DFL)
	args = mkArgParser().parse_args()
	n=int(args.n)
	m=int(args.m)
	p=float(args.p)
	k=int(args.k)

	file1 = open("simulateres","a")
	for i in range(0,10):
		num_node_G,num_edge_G,num_node_H,num_edge_H,num_node_ILP,num_edge_ILP,isbipartite=subRoutine(n,m,p,k)
		file1.write(str(n)+"\t"+str(m)+"\t"+str(p)+"\t"+str(k)+"\t"+str(num_node_G)+"\t"+str(num_edge_G)+"\t"+str(num_node_H)+"\t"+str(num_edge_H)+"\t"+str(num_node_ILP)+"\t"+str(num_edge_ILP)+"\t"+str(isbipartite)+"\n")

	file1.close() 
