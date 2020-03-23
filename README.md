# K-MWPBS (Top-K Maximum Weight Pathway Constrained Bipartite Subgraph)
--------------------------------
*Prerequisite*
--------------------------------
Please install:<br/>
python >= 3.6 version<br/>
Gurobi version 9<br/>
networkx version 2.4<br/>

Command:
------------------------------
python KMWPBS.py uniq_gene_list discover_me_gene_pairs kegg_pathway.gmt kwbps_plot_res > kwbps_results <br/>

Input parameters:<br/>
uniq_gene_list: file contains unique gene list found in discover_me_gene_pairs<br/>
discover_me_gene_pairs: discover mutual exclusivity gene pairs with third column as p-value (cut-off < 0.05)<br/>
kegg_pathway.gmt: kegg pathway<br/>
kwbps_plot_res: output KMWPBS plot name<br/>
kwbps_results: kwbps print out results, bipartite graph KMWPBS results<br/>

Output:<br/>
.png plot (with kwbps_plot_res name)<br/>
KMWPBS bipartite results <br/>


Simulation Command:
--------------------------------
python simulate01.py n m p k<br/>

Input parameters:<br/>
n (int) – The number of nodes in the first bipartite set<br/>
m (int) – The number of nodes in the second bipartite set<br/>
p (float) – Probability for edge creation, i.e. 0.5<br/>
k (int) - The number of random edges to be added<br/>
