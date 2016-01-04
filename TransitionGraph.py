from igraph import *

P_ABS = 1.0
P_D = 0.5
P_QA = 0.3
P_QB = 0.1
P_QB2 = 0.1
P_PQ = 0.0

P_QA_r = 0.0
P_QB_r = 0.0
P_QB2_r = 0.0

GraphPSIITransitions = Graph([(0,1), (1,2), (2,3), (2,4), (3,5), (4,5), (5,6), (6,7), (6,8), (7,9), (8,9), (9,10), (10,11),
							  (1,0), (2,1), (3,2), (4,2), (5,3), (5,4), (6,5), (7,6), (8,6), (9,7), (9,8), (10,9), (11,10),
							  (10,2), (8,0)], directed = True)

GraphPSIITransitions.es["Transition"] = ["Absorption", "Pheo to QA", "Absorption", "QA- to QB", "QA- to QB", "Absorption", "Pheo to QA", "Absorption", "QA- to QB-", "QA- to QB-", "Absorption", "Pheo to QA", "Absorption",
													   "Excitation decay", "QA- to Pheo", "Excitation decay", "QB- to QA", "QB- to QA", "Excitation decay", "QA- to Pheo", "Excitation decay", "QB2- to QA", "QB2- to QA", "Excitation decay", "QA- to Pheo", "Excitation decay",
													   "PQH2 diffusion", "PQH2 diffusion"]

GraphPSIITransitions.es["Probability"] = [P_ABS, P_QA, P_ABS, P_QB, P_QB, P_ABS, P_QA, P_ABS, P_QB2, P_QB2, P_ABS, P_QA, P_ABS,
													   P_D, P_QA_r, P_D, P_QB_r, P_QB_r, P_D, P_QA_r, P_D, P_QB2_r, P_QB2_r, P_D, P_QA_r, P_D,
													   P_PQ, P_PQ]

GraphPSIITransitions.vs["label"] = ["A QA QB", "A* QA QB", "A QA- QB", "A* QA- QB", "A QA QB-", "A* QA QB-", "A QA- QB-", "A* QA- QB-", "A QA QB2-", "A* QA QB2-", "A QA- QB2-", "A* QA- QB2-"]
GraphPSIITransitions.vs["A"] = ["ground", "excited", "ground", "excited", "ground", "excited", "ground", "excited", "ground", "excited", "ground", "excited"]
GraphPSIITransitions.vs["QA"] = [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1]
GraphPSIITransitions.vs["QB"] = [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2]

#print GraphPSIITransitions.es[0].target

#print GraphPSIITransitions.degree()
#print GraphPSIITransitions.degree(type = "in")
#print GraphPSIITransitions.degree(type = "out")

#length of the list
vertex = 1

transitionList = GraphPSIITransitions.es.select(_source=vertex)[:]["Transition"]

transitionListProbabilities = []
for transition in range(0, len(transitionList)):
	if transitionList[transition] != 'Absorption' and transitionList[transition] != 'Excitation decay':
		transitionListProbabilities.append([transitionList[transition], GraphPSIITransitions.es.select(_source=vertex)[:]["Probability"][transition]])

transitionListProbabilities.sort(key=lambda x: x[1])

print transitionListProbabilities
for transfer in transitionListProbabilities:
	print transfer

print GraphPSIITransitions.es.select(_source=vertex, Transition = transfer[0])[0]

print GraphPSIITransitions.vs[0]["A"]

#print GraphPSIITransitions.es.select(_source=vertex)[:]["Transition"]
#print GraphPSIITransitions.es.select(_source=vertex)[:]["Probability"]

#print GraphPSIITransitions.es.select(_source=vertex, Transition = "Absorption")[:]

#for nr in range(0, nrPossibleTransitions):
#	print 'Target: %i Probability: %i' % (GraphPSIITransitions.es.select(_source=vertex)[nr].target, GraphPSIITransitions.es.select(_source=vertex)[nr]["Probability"])
#	print  GraphPSIITransitions.es.select(_source=vertex)[nr]["Transition"]

#print GraphPSIITransitions.vs[vertex]["A"]

#transitionList = GraphPSIITransitions.es.select(_source=vertex)[:]["Transition"]
#print transitionList


#layout = GraphPSIITransitions.layout("rt")
#layout = GraphPSIITransitions.layout("kk")
#layout = GraphPSIITransitions.layout("star")
#layout = GraphPSIITransitions.layout("sugiyama")
layout = GraphPSIITransitions.layout("tree", root = 0)
#layout = GraphPSIITransitions.layout("grid_fr")
#plot(GraphPSIITransitions, layout = layout, bbox = (700, 700), margin = (50, 50, 50, 50))#, vertex.shape = "Rectangle")

#print GraphPSIITransitions.vs[0]