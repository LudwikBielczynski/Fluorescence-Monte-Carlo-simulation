import random
import matplotlib.pyplot as plt
import numpy as np
from igraph import *

def createPSIITransitionGraph(P_ABS, P_D, P_QA, P_QB, P_QB2, P_PQ, P_QA_r, P_QB_r, P_QB2_r):
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

    return GraphPSIITransitions

def chunks(l, numberOfGroups):
    """
    Yield n successive chunks from l.
    """
    chunkLength = len(l)/numberOfGroups
    for i in xrange(0, len(l), chunkLength):
        yield l[i:i + chunkLength]

class PSII(object):
    """
    Representation of a PSII particle with a blocked QA cite
    """   
    def __init__(self, layer = 1, size = 1, photonFlux = 1000, leafArea = 1000):
        """
        Initialize a PSII instance, saves all parameters as attributes of the instance.
        
        Input values:

        Values calculated:

        """
        self.layer = layer
        self.photonFlux = photonFlux
        self.leafArea = leafArea
        self.size = size

        self.P_ABS = self.photonFlux/float(self.leafArea) * self.size
        self.P_D = 1.0
        self.P_F = 0.3
        self.P_QA = 1.0                  #Based on Lazar and Schansker, 2009
        self.P_QB = 0.35              
        self.P_QB2 = 0.175

        self.P_QA_r = 0.0
        self.P_QB_r = 0.0175
        self.P_QB2_r = 0.0035

        self.P_PQ_initial = 0.08
        self.P_PQ = self.P_PQ_initial
        self.P_PQ_r = 0.005


        self.graphPSIITransitions = createPSIITransitionGraph(self.P_ABS, self.P_D, self.P_QA, self.P_QB, self.P_QB2, self.P_PQ, self.P_QA_r, self.P_QB_r, self.P_QB2_r)
        self.state = 0      #Node in the GraphPSIITransitions that the PSII is currently occupying
        self.stateA = self.graphPSIITransitions.vs[self.state]["A"]
        self.stateQA = self.graphPSIITransitions.vs[self.state]["QA"]
        self.stateQB = self.graphPSIITransitions.vs[self.state]["QB"]

    def updatePhotonFlux(self, PhotonFlux):
        """
        Updates the photonFlux, dependent variables and the transition graph: P_ABS, graphPSIITransitions

        input: int
        """    
        self.photonFlux = PhotonFlux
        self.P_ABS = self.photonFlux /float(self.leafArea) * self.size
        self.graphPSIITransitions = createPSIITransitionGraph(self.P_ABS, self.P_D, self.P_QA, self.P_QB, self.P_QB2, self.P_PQ, self.P_QA_r, self.P_QB_r, self.P_QB2_r)

    def updateState(self):
        self.stateA = self.graphPSIITransitions.vs[self.state]["A"]
        self.stateQA = self.graphPSIITransitions.vs[self.state]["QA"]
        self.stateQB = self.graphPSIITransitions.vs[self.state]["QB"]        

    def update(self, light = "on"):
        """
        Function describing the space of possible transition without or under illumination. The PSII is assumed have a QB site that is able to accept an electron.
        The reopening of RCs is added.
        When the light is off only the closed, excited RCs can fluoresce or decay non-radiatively.
        When the light is on following scheme is assumed:

        Input:
            light: str "on" or "off" representing if the photon flux will be hitting the RCs during a timestep

        returns a tuple of three booleans: True/False if a photon is absorbed, True/False if a photon is fluoresced, True/False if a PQ is diffused
        """            
        transitionList = self.graphPSIITransitions.es.select(_source=self.state)[:]["Transition"]
        transitionListProbabilities = self.graphPSIITransitions.es.select(_source=self.state)[:]["Probability"]
        Absorbed = False
        Fluoresced = False
        PQdiffuse = False

        if light == "off":
            if "Excitation decay" in transitionList:
                print 'not implemented'

        if light == "on":
            #Excitation
            randomizedNumber = random.random()
            if 'Absorption' in transitionList and randomizedNumber <= self.P_ABS:
                Absorbed = True
                self.state = self.graphPSIITransitions.es.select(_source=self.state, Transition = "Absorption")[0].target
                #if the excited PSII ends up with an empty QA in the time scale of 100us it makes always a transition to QA
                if self.state == 1:
                    self.state = 2
                if self.state == 5:
                    self.state = 6
                if self.state == 9:
                    self.state = 10                    
                self.updateState()

                transitionList = self.graphPSIITransitions.es.select(_source=self.state)[:]["Transition"]
                transitionListProbabilities = self.graphPSIITransitions.es.select(_source=self.state)[:]["Probability"]

            #All Electron transfer events
            transitionsElectronTransferProbabilities = self.graphPSIITransitions.es.select(_source=self.state)[:]["Probability"]
            transitionsElectronTransfer = []
            P_ETr = 0 #Probability of Electron transfer
            for transition in range(0, len(transitionList)):
                if transitionList[transition] != 'Absorption' and transitionList[transition] != 'Excitation decay' and transitionList[transition] != 'PQH2 diffusion':
                    transitionsElectronTransfer.append([transitionList[transition], transitionsElectronTransferProbabilities[transition]])
                    P_ETr += transitionsElectronTransferProbabilities[transition]
            if len(transitionsElectronTransfer) != 0:
                transitionsElectronTransfer.sort(key=lambda x: x[1])
                randomizedNumber = random.random()

                if len(transitionsElectronTransfer) == 1:
                    transfer = transitionsElectronTransfer[0]
                    if randomizedNumber <= transfer[1]:
                        self.state = self.graphPSIITransitions.es.select(_source=self.state, Transition = transfer[0])[0].target
                        self.updateState()
                        transitionList = self.graphPSIITransitions.es.select(_source=self.state)[:]["Transition"]

                if len(transitionsElectronTransfer) == 2:
                    if randomizedNumber <= P_ETr:
                        transfer = transitionsElectronTransfer[0]
                        if randomizedNumber <= transfer[1]/P_ETr:
                            self.state = self.graphPSIITransitions.es.select(_source=self.state, Transition = transfer[0])[0].target
                            self.updateState()
                            transitionList = self.graphPSIITransitions.es.select(_source=self.state)[:]["Transition"]
                        transferBack = transitionsElectronTransfer[1]
                        if randomizedNumber >= transfer[1]/P_ETr:
                            self.state = self.graphPSIITransitions.es.select(_source=self.state, Transition = transferBack[0])[0].target
                            self.updateState()
                            transitionList = self.graphPSIITransitions.es.select(_source=self.state)[:]["Transition"]

            #Decay
            if 'Excitation decay' in transitionList and randomizedNumber <= self.P_D:
                self.state = self.graphPSIITransitions.es.select(_source=self.state, Transition = "Excitation decay")[0].target
                self.updateState()
                transitionList = self.graphPSIITransitions.es.select(_source=self.state)[:]["Transition"]
                transitionListProbabilities = self.graphPSIITransitions.es.select(_source=self.state)[:]["Probability"]                
                randomizedNumber = random.random()
                if randomizedNumber <= self.P_F:
                    Fluoresced = True

            #PQ diffusion
            if 'PQH2 diffusion' in transitionList and randomizedNumber <= self.P_PQ:
                PQdiffuse = True
                self.state = self.graphPSIITransitions.es.select(_source=self.state, Transition = "PQH2 diffusion")[0].target
                self.updateState()

        return Absorbed, Fluoresced, PQdiffuse

    def updatePQ(self, actualPQPool, maxPQPool):
#        if actualPQPool <= 1:
#            self.P_PQ = 0
#        else:
#            self.P_PQ = self.P_PQ_initial
        self.P_PQ = self.P_PQ_initial * (float(actualPQPool)/maxPQPool)

class Layer(object):
    """
    Representation of layers in the leaf.
    """
    def __init__(self, PSIIs, layersNumber, PQnumber):
        """
        Initialization function, saves the PSIIs.

        Input:
            PSIIs: list of PSII objects (representing a fraction of total PSIIs assigned to this layer)
            layersNumber: int representing the number of layers into which the RCs will be destributed

        Fcount: int representing the number of photons fluoresced from a layer
        AbsorbedCount: int representing the number of photons absorbed by the layer        
        """
        self.maxPQPool = PQnumber
        self.actualPQPool = PQnumber
        self.PSIIs = PSIIs
        self.layersNumber = layersNumber
        self.FCount = 0
        self.AbsorbedCount = 0

    def updateCounters(self, Absorbed, Fluoresced, PQ):
        """
        Changes the count of the absorbed and fluoresced photons by a specific layer      
        """
        if Absorbed == True:
            self.AbsorbedCount += 1
        if Fluoresced == True:
            self.FCount += 1
        if PQ == True:
            if self.actualPQPool > 0:
                self.actualPQPool -= 1

    def updatePSIIs(self, light, PhotonFlux, PSIIs = None):
        """
        Updates the PSIIs and adjusts the PhotonFlux parameters

        returns: a pair of int representing the number of Fluoresced and Absorbed photons by the layer 
        """       
        if PSIIs == None:
            PSIIs = self.PSIIs
        self.FCount = 0
        self.AbsorbedCount = 0
        self.QA = 0
        self.QAminus = 0        
        self.QB = 0
        self.QBminus = 0
        self.QB2minus = 0
        for psii in PSIIs:
            psii.updatePhotonFlux(PhotonFlux)
            if light == "on":
                if psii.P_ABS < 1:                            #case of single excitations
                    Absorbed, Fluoresced, PQ = psii.update(light)
                    self.updateCounters(Absorbed, Fluoresced, PQ)
                    if random.random() <= psii.P_PQ_r and self.actualPQPool < self.maxPQPool:
                        self.actualPQPool += 1
                    psii.updatePQ(self.actualPQPool, self.maxPQPool)
                else:
                    for excitation in range(0,int(psii.P_ABS)):   #case if multiple excitations
                        Absorbed, Fluoresced, PQ = psii.update(light)
                        self.updateCounters(Absorbed, Fluoresced, PQ)
                        if random.random() <= psii.P_PQ_r and self.actualPQPool < self.maxPQPool:
                            self.actualPQPool += 1                        
                        psii.updatePQ(self.actualPQPool, self.maxPQPool)

            if light == "off":
                Absorbed, Fluoresced, PQ = psii.update(light)
                self.updateCounters(Absorbed, Fluoresced, PQ)
                psii.updatePQ(self.actualPQPool, self.maxPQPool)                

            if psii.stateQA == 0:
                self.QA += 1
            if psii.stateQA == 1:
                self.QAminus += 1                    
            if psii.stateQB == 0:
                self.QB += 1
            if psii.stateQB == 1:
                self.QBminus += 1
            if psii.stateQB == 2:
                self.QB2minus += 1

        return self.FCount, self.AbsorbedCount, self.QA, self.QAminus, self.QB, self.QBminus, self.QB2minus

class Leaf(object):
    """
    Representation of a simplified leaf.
    """    
    def __init__(self, PSIIs, layersNumber, PQnumber):
        """
        Initialization function, saves the PSIIs

        Input:
            PSIIs: list representing all PSIIs objects in the leaf
            LayersNumber: int representing the number of layers

        Layers: list representing the Layers of PSIIs present in the leaf
        totalFluoresced: int representing the number of photons Fluoresced by the whole leaf
        totalAbsorbed: int representing the number of photons Absorbed by the whole leaf
        """
        self.PQnumber = PQnumber
        self.QA = 0
        self.QAminus = 0
        self.QB = 0
        self.QBminus = 0
        self.QB2minus = 0

        self.PSIIs = PSIIs
        self.layersNumber = layersNumber

        self.Layers = []

        self.totalFluoresced = 0
        self.totalAbsorbed = 0

    def assignPSIIToLayers(self):
        """
        Change the PSII layer number
        """
        layerPSIIs = chunks(self.PSIIs, self.layersNumber)
        try:
            layer = 0
            while True:
                for psii in layerPSIIs.next():
                    psii.layer = layer
                layer += 1
        except StopIteration:
            pass

    def getPSIIsInLayer(self, layer):
        """
        Returns a list of PSIIs from a specific layer
        """
        layersPSIIs = list(chunks(self.PSIIs, self.layersNumber))
        return layersPSIIs[layer]

    def createLayers(self):
        """
        Creates the layer objects and assigns PSIIs to selected ones.
        """        
        for layer in range(0, self.layersNumber):
            PSIIs = self.getPSIIsInLayer(layer)
            self.Layers.append(Layer(PSIIs, layer, self.PQnumber))

    def updateLayers(self, light):
        """
        Calculates how much light is fluoresced and absorbed by passing through all the layers.

        returns: a pair of int representing the total amount of Fluoresced and Absorbed lught by the leaf
        """            
        self.totalFluoresced = 0
        self.totalAbsorbed = 0
        self.PQ = 0
        self.QA = 0
        self.QAminus = 0        
        self.QB = 0
        self.QBminus = 0
        self.QB2minus = 0
        PSIIs = self.PSIIs 
        PhotonFlux = PSIIs[0].photonFlux
        for layer in self.Layers:
            Fluoresced, Absorbed, QA, QAminus, QB, QBminus, QB2minus = layer.updatePSIIs(light, PhotonFlux)
            self.totalFluoresced += Fluoresced
            self.totalAbsorbed += Absorbed
            self.PQ += layer.actualPQPool
            self.QA += QA
            self.QAminus += QAminus            
            self.QB += QB
            self.QBminus += QBminus
            self.QB2minus += QB2minus
            PhotonFlux = PhotonFlux - Absorbed + Fluoresced
          
        return self.totalFluoresced, self.totalAbsorbed, self.PQ, self.QA, self.QAminus, self.QB, self.QBminus, self.QB2minus

def simulatingLeaf(numPSIIs = 1000, timeSteps = 100, trialsNum = 1, size = 1, photonFlux = 1000, layers = 1, PQnumber = 1000):
    """
    Runs simulations and plots graphs for PSIIs in the leaf.
    """
    timepoint = 1000
    trialsSum = [0]
    PQsum = [100]
    QAsum = [100] 
    QAminusSum = [0]        
    QBsum = [100]
    QBminusSum = [0]
    QB2minusSum = [0]            
    
    for time in range(1, timeSteps + 1):
        trialsSum.append(0)
        PQsum.append(0)
        QAsum.append(0)
        QAminusSum.append(0)        
        QBsum.append(0)
        QBminusSum.append(0)
        QB2minusSum.append(0)

    for trial in range(0, trialsNum):
        PSIIs = []                                      #Creating PSIIs
        for nr in range(0, numPSIIs):
            #PSIIs.append(PSII_DCMU(size = size, state = "open", photonFlux = photonFlux, leafArea = 10000, probabilityFluorescence = 0.3, probablilityAnihilation = 0.01))
            PSIIs.append(PSII(size = size, photonFlux = photonFlux, leafArea = 10000))
        simulatedLeaf = Leaf(PSIIs, layers, PQnumber)             #Creating the leaf
        simulatedLeaf.assignPSIIToLayers()              #Creating layers in the leaf
        simulatedLeaf.createLayers()

        Fluorescence = [0]
        PQlist = [PQnumber]

        for time in range(1, timeSteps+1):
            if time%250 == 0:
                print time
            Fluoresced, Absorbed, PQ, QA, QAminus, QB, QBminus, QB2minus = simulatedLeaf.updateLayers(light = "on")
            #print "Fluoresced: %i Absorbed: %i" % (Fluoresced, Absorbed)
            trialsSum[time] += Fluoresced
            PQsum[time] += PQ
            QAsum[time] += QA
            QAminusSum[time] += QAminus
            QBsum[time] += QB
            QBminusSum[time] += QBminus
            QB2minusSum[time] += QB2minus

        if trial%1 == 0:
            print 'Trial nr: %i' % trial

    timeReal = range(0,timeSteps + 1)
    for time in range(1, timeSteps + 1):
        trialsSum[time] /= float(trialsNum)
        PQsum[time] /= (float(trialsNum) * PQnumber)/100
        QAsum[time] /= (float(trialsNum) * numPSIIs)/100
        QAminusSum[time] /= (float(trialsNum) * numPSIIs)/100        
        QBsum[time] /= (float(trialsNum) * numPSIIs)/100
        QBminusSum[time] /= (float(trialsNum) * numPSIIs)/100
        QB2minusSum[time] /= (float(trialsNum) * numPSIIs)/100
        timeReal[time] /= float(timeSteps + 1)

    fig, ax1 = plt.subplots()

    ax1.set_xlim(xmin = timeReal[1], xmax = timeReal[timeSteps])
    Fluorescence = ax1.plot(timeReal, trialsSum, 'b-', label = "Size: " + str(size) + " PhotonFlux: " + str(photonFlux) + " Layers: " + str(layers) )
    ax1.set_xlabel("Time [s]")
    ax1.set_xscale('log')
    ax1.set_ylabel("Ft [counts]")

    ax2 = ax1.twinx()
    ax2.set_xlim(xmin = timeReal[1], xmax = timeReal[timeSteps])
    ax2.set_ylim(ymin = 0, ymax = 100)
    ax2.set_ylabel("[%]")        
    QA = ax2.plot(timeReal, QAsum, 'm-', label = "QA pool" )
    QAminus = ax2.plot(timeReal, QAminusSum, 'm--', label = "QA- pool" )
    QB = ax2.plot(timeReal, QBsum, 'g-', label = "QB pool" )
    QBminus = ax2.plot(timeReal, QBminusSum, 'g--', label = "QB- pool" )
    QB2minus = ax2.plot(timeReal, QB2minusSum, 'g:', label = "QB2- pool")
    PQ = ax2.plot(timeReal, PQsum, 'r-', label = "PQ pool" )

    lns = Fluorescence + QA + QAminus + QB + QBminus + QB2minus + PQ
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc = 'best', fontsize = 'small')
    return trialsSum, PQsum

#####################################################
###################SIMULATION########################
#####################################################
numPSIIs = 1000
timeSteps = 10000
trialsNum = 10
size = 1
layer = 1
PQnumber = 3750 # based on Oja et al., 2011
photonFluxList = [1000]

projectPath = '/home/ludwik/Documents/python/Monte-Carlo/Layers and size/'
def Simulate(numPSIIs, timeSteps, trialsNum, photonFluxList, size, layer):
    for light in photonFluxList:
        print 'Light: %i' % light
        simulatingLeaf(numPSIIs = numPSIIs, timeSteps = timeSteps, trialsNum = trialsNum, size = size, photonFlux = light, layers = layer, PQnumber = PQnumber)
    fileName = str('numPSIIs%i timeSteps%i trialsNum%i size%.2f layers%i lightDependency.svg' % (numPSIIs, timeSteps, trialsNum, size, layer))
    plt.savefig(projectPath + "QA" + fileName, width = 30, height = 8)
    plt.close()
Simulate(numPSIIs, timeSteps, trialsNum, photonFluxList, size, layer)
#
#size = 1.0
#photonFlux = 1000
#a = PSII(size = size, photonFlux = photonFlux, leafArea = 10000)
#for nr in range(0,100):
#    print a.update()
#    print a.graphPSIITransitions.vs[a.state]["label"]