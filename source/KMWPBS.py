from argparse import ArgumentParser
from signal import signal, SIGPIPE, SIG_DFL
from gurobipy import *
import networkx as nx
from networkx.algorithms import bipartite
from networkx.classes.function import *
import collections
import math 
import random
import matplotlib.pyplot as plt

plt.switch_backend('agg') 

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

		B.add_nodes_from(rightn, bipartite=0)
		B.add_nodes_from(leftn, bipartite=1)
		B.add_edges_from(l_edges)
			
		for (u, v) in B.edges:
			if B.has_edge(u,v):
				B[u][v]['weight'] = 1

		print("----------------------- Bipartite Score ----------------------")
		print('Score:', model.objVal)
		print("--------------------------------------------------------------------")
		return B
	except GurobiError:
		print('Error reported')

def mkArgParser():
  parser = ArgumentParser()

  # Positional arguments, i.e. arguments that HAVE to be provided
  parser.add_argument("allgenes", help="input i.e 143genes_discover_me_unique file")
  parser.add_argument("me_results", help="input i.e. BRCA_me_DISCOVER file")
  parser.add_argument("pathway", help="input i.e. kegg.gmt file")
  parser.add_argument("plotName", help="input i.e. plotfigure name")

  return parser

if __name__ == '__main__':
	signal(SIGPIPE,SIG_DFL)
	args = mkArgParser().parse_args()

	#Read gene list and add the node (gene) to graph G
	G = nx.Graph()
	node_index = {}
	index_node = {}
	pval_ori={}
	index1=1
	with open(args.allgenes, "r") as f:
		for line in f:
			line_ = line.split('\t')

			node1= line_[0].strip()
			node_index[node1]=index1
			index_node[index1]=node1
			G.add_node(index1)
			index1 += 1


	pathwayDict={}
	pathwayIDName={}
	#READ IN PATHWAY
	with open(args.pathway, "r") as f:
		for line in f:
			line=line.strip()
			line_ = line.split('\t')
			pathwayid= line_[0].strip()
			pathwayname= line_[1].strip()
			line_ = line_[2:]
			if not pathwayid in pathwayDict:
				pathwayDict[pathwayid]=set()
				pathwayIDName[pathwayid]=pathwayname

			for a in line_:
				if a.strip() in node_index:
					pathwayDict[pathwayid].add(node_index[a.strip()])


	#Read ME results and fill the edges of node (select those with pvalue < 0.05 only)
	#convert p-value to weight and store as edge's weight
	with open(args.me_results, "r") as f:
		for line in f:
			line_ = line.split('\t')
			gene1=line_[0].strip()
			gene2=line_[1].strip()
			pval=line_[2].strip()

			edge1=node_index[gene1]
			edge2=node_index[gene2]
			
			w = float(pval)
			if float(pval) < 0.05 and float(pval) > 0:
				w = (-1*(math.log2(float(pval)))) / (math.log2(100))
				pval_ori[edge1,edge2]=pval
				pval_ori[edge2,edge1]=pval
				G.add_edge(edge1,edge2,weight=w)
			if float(pval) == 0:
				pval_ori[edge1,edge2]=pval
				pval_ori[edge2,edge1]=pval
				G.add_edge(edge1,edge2,weight=1)


	#run ILP
	runidx=1
	while  G.number_of_edges() > 0:
		number_of_nodes_G=len(G)
		number_of_edges_G=G.size()
		weight_total_G =G.size(weight='weight')
		
		#RUN ILP
		B=runILP(G,pathwayDict,pathwayIDName)
		number_of_nodes_B=len(B)
		number_of_edges_B= B.number_of_edges()
		weight_total_B = B.size(weight='weight')

		if is_empty(B):
			print("B is empty, no results! End Here!")
			break
		print("----------------------- Bipartite Extracted Results ----------------------")

		if bipartite.is_bipartite(B):
			print("B is Bipartite graph")
		else:
			print("B is not Bipartite graph")

		print("")
		print("weight_total_G="+str(weight_total_G))
		print("weight_total_B="+str(weight_total_B))
		print("number_of_nodes_G="+str(number_of_nodes_G))
		print("number_of_nodes_B="+str(number_of_nodes_B))
		print("number_of_edges_G="+str(number_of_edges_G))
		print("number_of_edges_B="+str(number_of_edges_B))

		print("----------------------- Bipartite End Results ----------------------")

		print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		try:
			bottom_nodes = {n for n, d in B.nodes(data=True) if d['bipartite']==0}
			top_nodes = set(B) - bottom_nodes
			
			mapping={}
			pos = {}
			print("Bottom Nodes:")
			print(bottom_nodes)
			for s in bottom_nodes:
				print(str(index_node[s]),end =";")
				mapping[s]=str(index_node[s])
			print("")
				
			print("Top Nodes:")
			print(top_nodes)
			for t in top_nodes:
				print(str(index_node[t]),end =";")
				mapping[t]=str(index_node[t])
			print("")
			
			plt.figure()
			H=nx.relabel_nodes(B,mapping)

			r= {n for n, d in H.nodes(data=True) if d['bipartite']==0}
			l = set(H) - r
			color = bipartite.color(H)
			color_dict = {0:'r',1:'lightcyan'}
			color_list = [color_dict[i[1]] for i in H.nodes.data('bipartite')]
			nodesize_list = [1000] * len(H)

			# Update position for node from each group
			pos.update((node, (1, index)) for index, node in enumerate(l))
			pos.update((node, (2, index)) for index, node in enumerate(r))
			
			nx.draw(H,pos,with_labels=True, node_color = color_list,node_size=nodesize_list)
			plotname="solution_"+args.plotName+"_"+str(runidx)+".png"

			plt.savefig(plotname,dpi=350, format="PNG")
			runidx += 1
		except nx.exception.AmbiguousSolution:
			print("Ambiguous Solution From NetworkX Bipartite, Just Used the Above Returned Bipartite Graph Results.")
		
		G=removeEdges(G,B)
		weight_total_G=0
		print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print("New G After Edges Removal*************************************************")
		#for (u, v) in G.edges:
		#	print("u="+str(u)+",v="+str(v))
