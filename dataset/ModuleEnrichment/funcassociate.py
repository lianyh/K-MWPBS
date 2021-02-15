#!/usr/bin/env python
import sys
import json
import httplib
import numpy as np
import collections
from argparse import ArgumentParser
from signal import signal, SIGPIPE, SIG_DFL

def mkArgParser():
	parser = ArgumentParser()
	parser.add_argument("genelist", help="predicted bipartite gene set from partitions")
	parser.add_argument("m", help="module 1 or 1")
	parser.add_argument("out", help="output folder")
	return parser

class FuncassociateClient(object):
    """Query funcassociate to find out enriched GO terms."""
    
    host = 'llama.mshri.on.ca'
    query_url = '/cgi/funcassociate/serv'

    def __init__(self):
        self.c = httplib.HTTPConnection(self.host)

    def close_conn(self):
        self.c.close()

    def jsonify(self, data):
        return (json.dumps(data)).encode('utf-8')

    def request(self, payload):
        self.c.request('POST', self.query_url, self.jsonify(payload), headers={'Content-type': 'application/json'})
        response = self.c.getresponse()
        if response.status == httplib.OK:
            return response.read()

    def available_species(self):
        payload = {'method': 'available_species','id': 0,}
        return self.request(payload)

    def available_namespaces(self, species=['Homo sapiens']):
        payload = {'method': 'available_namespaces','params': species,'id': 123123123}
        return self.request(payload)

    def go_associations(self, 
        params=['Homo sapiens', 'uniprot_swissprot_accession', ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP']]):
        payload = {'method': 'go_associations','params': params,'id': 1,}
        return self.request(payload)


    def functionate(self, query, species="Homo sapiens", namespace="entrezgene", genespace=None, mode="unordered", reps=2000, support=None, associations=None):
        params_dict = { "query": query,"species": species,"namespace": namespace,"mode": mode,"reps": reps}
        if associations is not None:
            params_dict["support"] = support
        if associations is not None:
            params_dict["associations"] = associations
            del params_dict["species"]
            del params_dict["namespace"]
        if genespace is not None:
            params_dict["genespace"] = genespace
        payload = {'id': 1,'method': 'functionate','params': [ params_dict ],'jsonrpc': '2.0'}
        response = self.request(payload=payload)
		#print response
        response = json.loads(response)
		#self.close_conn()
				
        if "error" in response and response["error"] is not None:
            raise ValueError("Error %s" % response['error']['message'])

        response = response['result']
        return response


if __name__ == '__main__':
    signal(SIGPIPE,SIG_DFL)
    args = mkArgParser().parse_args()
    inputall=collections.OrderedDict()
    idx=0
    with open(args.genelist, "r") as f:
        for line in f:
            # ignore header lines
            if not line.startswith("@"):
                line_ = line.strip()
                hhh=line_.split('\t')
                h2 = hhh[1:]
                inputall[idx]=h2
                idx += 1

    fa = FuncassociateClient()
    for i in range(idx):
        file1=open(args.out+"BPM_"+str(i)+"_module"+str(args.m),"w+")
        inputlist=inputall[i]
        payload = {'id': 1,'method': 'functionate','params': [ { "query": inputlist,"species":"Saccharomyces cerevisiae","namespace": "sgd_systematic","mode":"unordered"}], 'jsonrpc': '2.0'}
        response=fa.request(payload=payload)
        response = json.loads(response)
        response = response['result']
        print(response)
        reslist=[]
        for a in response:
            if a=="over":
                reslist=response[a]
	
        for b in reslist:
            file1.write(str(b))

        file1.close()
	#print response
    fa.close_conn()
