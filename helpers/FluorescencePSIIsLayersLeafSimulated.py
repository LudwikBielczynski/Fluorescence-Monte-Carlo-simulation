import random
import matplotlib.pyplot as plt
import numpy as np

def chunks(l, numberOfGroups):
    """
    Yield n successive chunks from l.
    """
    chunkLength = len(l)/numberOfGroups
    for i in xrange(0, len(l), chunkLength):
        yield l[i:i + chunkLength]

class PSII_DCMU(object):
    """
    Representation of a PSII particle with a blocked QA cite
    """   

    def __init__(self, layer = 1, size = 1, state = "open", photonFlux = 1000, leafArea = 1000, probabilityFluorescence = 0.3, probablilityAnihilation = 0.01):
        """

        Initialize a PSII instance, saves all parameters as attributes of the instance.
        
        Input values:
            layer: int representing the number of layers between which the PSIIs are distributed
            size: int representing size of the PSII
            state: str representing the state of PSII accepted values "open", "closed ground" or "closed excited"
            photonFlux: int representing the number of photons that are appearing in one light: "on" event
            probabilityFluorescence: float in the range 0-1 representing the probability of a closed excited RC to fluoresce a photon
            probablilityAnihilation: float in the range 0-1 representing the probability of a closed excited RC to decay to the closed ground state during a double excitation event

        Values calculated:
            probabilityAbsorbed: float from 0-1 representing the probability of a RC to absorb a photon. If above 1 represents multiple photon excitation. 
            probabilityDecay
        """
        self.layer = layer
        self.photonFlux = photonFlux
        self.leafArea = leafArea
        self.size = size
        self.state = state
        self.probabilityAbsorbed = self.photonFlux/float(self.leafArea) * self.size

        self.probabilityFluorescence = probabilityFluorescence 
        self.probabilityDecay = 1 - np.exp(-1/3.5)                            #probability of decay based on the fluorescence lifetime and selected timestep


        self.probablilityAnihilation = probablilityAnihilation

    def updatePhotonFlux(self, PhotonFlux):
        """
        Updates the photonFlux and dependent variables: probabilityAbsorbed

        input: int
        """    
        self.photonFlux = PhotonFlux
        self.probabilityAbsorbed = self.photonFlux /float(self.leafArea) * self.size

    def doesFluoresce(self):
        """
        Checks if the photon will be fluoresced from a PSII.
        The excited Reaction Center of a PSII can decay in two manners fluorescence or non-radiative decay (ISC).
            If an RC fluoresced it decays to a ground state.
            Otherwise it decays non-radiatevely

        returns boolean: True if photon is fluoresced False otherwise
        """          
        if self.state == "closed excited":
            if random.random() <= self.probabilityDecay:
                self.state = "closed ground"
                if random.random() <= self.probabilityFluorescence:	
                    return True
                else:
                    return False                                    #radiationless decay
            else:
                return False
        else:
            return False

    def update(self, light):
        """
        Function describing the space of possible transition without or under illumination. The PSII is assumed to be simplified and DCMU treated.
        When the light is off only the closed, excited RCs can fluoresce or decay non-radiatively.
        When the light is on following scheme is assumed:

        open -> closed ground <-> closed excited

        Input:
            light: str "on" or "off" representing if the photon flux will be hitting the RCs during a timestep

        returns a pair of booleans: True/False if a photon is absorbed and True/False if a photon is fluoresced 
        """            
        if light == "off":
            if self.state == "closed excited":
                return self.doesFluoresce()

        if light == "on":         
            if random.random() <= self.probabilityAbsorbed:
                Absorbed = True
                if self.state == "open":
                    self.state = "closed ground"
                    return Absorbed, False
                if self.state == "closed ground":
                    self.state = "closed excited"
                    return Absorbed, self.doesFluoresce()
                if self.state == "closed excited":
                    if random.random() <= self.probablilityAnihilation: #Check what is really happening during double excitation
                        return Absorbed, False
                    else:
                        return Absorbed, self.doesFluoresce()
            else:
                Absorbed = False
                if self.state == "closed excited":                 #Radiationless decay if the light is not absorbed   
                    return Absorbed, self.doesFluoresce()
                else:
                    return Absorbed, False  

class PSII_QB(PSII_DCMU):
    """
    Representation of a PSII particle that can transfer the electron from the QA to QB site without a limitation on the PlastoQuinone pool
    """ 
    def __init__(self, layer = 1, size = 1, state = "open", photonFlux = 1000, leafArea = 1000, probabilityFluorescence = 0.3, probablilityAnihilation = 0.01):
        """

        Initialize a PSII instance, saves all parameters as attributes of the instance.
        
        Input values:
            layer: int representing the number of layers between which the PSIIs are distributed
            size: int representing size of the PSII
            state: str representing the state of PSII accepted values "open", "closed ground" or "closed excited"
            photonFlux: int representing the number of photons that are appearing in one light: "on" event
            probabilityFluorescence: float in the range 0-1 representing the probability of a closed excited RC to fluoresce a photon
            probablilityAnihilation: float in the range 0-1 representing the probability of a closed excited RC to decay to the closed ground state during a double excitation event

        Values calculated:
            probabilityAbsorbed: float from 0-1 representing the probability of a RC to absorb a photon. If above 1 represents multiple photon excitation. 
            probabilityDecay
        """
        PSII_DCMU.__init__(self, layer, size, state, photonFlux, leafArea, probabilityFluorescence, probablilityAnihilation)
        self.stateQAReduction = 0
        self.stateQBReduction = 0
        self.probabilityPQdiffusion = 0.3

        self.probablityPheotoQA = 0.1
        self.probabilityQAtoQB = 0.1

    def update(self, light):
        """
        Function describing the space of possible transition without or under illumination. The PSII is assumed have a QB site that is able to accept an electron.
        The reopening of RCs is added.
        When the light is off only the closed, excited RCs can fluoresce or decay non-radiatively.
        When the light is on following scheme is assumed:

        open <-> closed ground <-> closed excited

        Input:
            light: str "on" or "off" representing if the photon flux will be hitting the RCs during a timestep

        returns a pair of booleans: True/False if a photon is absorbed and True/False if a photon is fluoresced 
        """            
        if light == "off":
            if self.state == "closed excited":
                return self.doesFluoresce()

        if light == "on": 
            PQdiffuse = False
            if self.stateQBReduction == 2 and random.random() <= self.probabilityPQdiffusion:
                self.stateQBReduction = 0
                PQdiffuse = True

            if self.stateQAReduction == 0 and random.random() <= self.probablityPheotoQA:
                self.state = "open"                 
                self.stateQAReduction += 1

            if self.state == "closed ground" and self.stateQAReduction == 1 and self.stateQBReduction != 2 and random.random() <= self.probabilityQAtoQB: #Added reopening of RCs
                self.state = "open"  
                self.stateQBReduction += 1
                self.stateQAReduction -= 1

            if random.random() <= self.probabilityAbsorbed:
                Absorbed = True
                if self.state == "open":
                    self.state = "closed ground"
                    return Absorbed, False, PQdiffuse
                if self.state == "closed ground":
                    self.state = "closed excited"
                    return Absorbed, self.doesFluoresce(), PQdiffuse
                if self.state == "closed excited":
                    if random.random() <= self.probablilityAnihilation: #Check what is really happening during double excitation
                        return Absorbed, False, PQdiffuse
                    else:
                        return Absorbed, self.doesFluoresce(), PQdiffuse
            else:
                Absorbed = False
                if self.state == "closed excited":                 #Radiationless decay if the light is not absorbed   
                    return Absorbed, self.doesFluoresce(), PQdiffuse
                else:
                    return Absorbed, False, PQdiffuse

    def updatePQ(self, actualPQPool, maxPQPool):
        self.probabilityPQdiffusion = float(actualPQPool)/maxPQPool 

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
                if psii.probabilityAbsorbed < 1:                            #case of single excitations
                    Absorbed, Fluoresced, PQ = psii.update(light)
                    self.updateCounters(Absorbed, Fluoresced, PQ)
                    psii.updatePQ(self.actualPQPool, self.maxPQPool)
                else:
                    for excitation in range(0,int(psii.probabilityAbsorbed)):   #case if multiple excitations
                        Absorbed, Fluoresced, PQ = psii.update(light)
                        self.updateCounters(Absorbed, Fluoresced, PQ)
                        psii.updatePQ(self.actualPQPool, self.maxPQPool)

                if psii.stateQAReduction == 0:
                    self.QA += 1
                if psii.stateQAReduction == 1:
                    self.QAminus += 1                    
                if psii.stateQBReduction == 0:
                    self.QB += 1
                if psii.stateQBReduction == 1:
                    self.QBminus += 1
                if psii.stateQBReduction == 2:
                    self.QB2minus += 1


            if light == "off":
                if psii.update(light) == True:
                    self.FCount += 1

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

selectedTimepoint = []

def simulatingLeaf(numPSIIs = 1000, timeSteps = 100, trialsNum = 1, size = 1, photonFlux = 1000, layers = 1, PQnumber = 1000):
    """
    Runs simulations and plots graphs for PSIIs in the leaf.
    """
    global selectedTimepoint
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
            PSIIs.append(PSII_QB(size = size, state = "open", photonFlux = photonFlux, leafArea = 10000, probabilityFluorescence = 0.3, probablilityAnihilation = 0.01))
        simulatedLeaf = Leaf(PSIIs, layers, PQnumber)             #Creating the leaf
        simulatedLeaf.assignPSIIToLayers()              #Creating layers in the leaf
        simulatedLeaf.createLayers()

        Fluorescence = [0]
        PQlist = [PQnumber]

        for time in range(1, timeSteps+1):
            Fluoresced, Absorbed, PQ, QA, QAminus, QB, QBminus, QB2minus = simulatedLeaf.updateLayers(light = "on")
            #print "Fluoresced: %i Absorbed: %i" % (Fluoresced, Absorbed)
            trialsSum[time] += Fluoresced
            PQsum[time] += PQ
            QAsum[time] += QA
            QAminusSum[time] += QAminus
            QBsum[time] += QB
            QBminusSum[time] += QBminus
            QB2minusSum[time] += QB2minus

        #for time in range(timeSteps + 1, timeSteps*2 + 2):
        #    Fluoresced, Absorbed = simulatedLeaf.updateLayers(light = "off")
        #    #print "Fluoresced: %i Absorbed: %i" % (Fluoresced, Absorbed)
        #    Fluorescence.append(Fluoresced)
        #    trialsSum[time] += Fluorescence[time]
        if trial%10 == 0:
            print 'Trial nr: %i' % trial

    for time in range(1, timeSteps+1):
        trialsSum[time] /= float(trialsNum)
        PQsum[time] /= (float(trialsNum) * PQnumber)/100
        QAsum[time] /= (float(trialsNum) * numPSIIs)/100
        QAminusSum[time] /= (float(trialsNum) * numPSIIs)/100        
        QBsum[time] /= (float(trialsNum) * numPSIIs)/100
        QBminusSum[time] /= (float(trialsNum) * numPSIIs)/100
        QB2minusSum[time] /= (float(trialsNum) * numPSIIs)/100

    fig, ax1 = plt.subplots()
    selectedTimepoint.append(trialsSum[timepoint])
    
    Fluorescence = ax1.plot(range(0,timeSteps + 1), trialsSum, 'b-', label = "Size: " + str(size) + " PhotonFlux: " + str(photonFlux) + " Layers: " + str(layers) )
    ax1.set_xlabel("Time [ms]")
    ax1.set_xscale('log')
    ax1.set_ylabel("Ft [counts]")
    ax1.set_xlim(xmin = 2,xmax = timeSteps + 1)

    ax2 = ax1.twinx()
    PQ = ax2.plot(range(0,timeSteps + 1), PQsum, 'r-', label = "PQ pool" )
    ax2.set_ylabel("[%]")
    QA = ax2.plot(range(0,timeSteps + 1), QAsum, 'm-', label = "QA pool" )
    QAminus = ax2.plot(range(0,timeSteps + 1), QAminusSum, 'm--', label = "QA- pool" )    
    QB = ax2.plot(range(0,timeSteps + 1), QBsum, 'g-', label = "QB pool" )
    QBminus = ax2.plot(range(0,timeSteps + 1), QBminusSum, 'g--', label = "QB- pool" )
    QB2minus = ax2.plot(range(0,timeSteps + 1), QB2minusSum, 'g:', label = "QB2- pool" )        
    lns = Fluorescence + PQ + QA + QAminus + QB + QBminus + QB2minus
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc = (.05, .5), fontsize = 'small')
    return trialsSum, PQsum

#####################################################
###################SIMULATION########################
#####################################################

#numPSIIs = 10000
#numPSIIsList = [10000,10000,10000]
#timeSteps = 5
#trialsNum = 1000
#size = 0.5
#photonFluxList = range(500,1001,500)
#layersList = range(1,5)
#numPSIIsList = [10000,10000,10000,10000]
#photonFlux = 1000
#for layer in layersList:
#    #print 'Light: %i' % light
#    simulatingLeaf(numPSIIs = numPSIIsList[layer-1], timeSteps = timeSteps, trialsNum = trialsNum, size = size, photonFlux = photonFlux, layers = layer)
#plt.legend(loc = "best")
#projectPath = '/home/ludwik/Documents/python/Monte-Carlo/'
#fileName = str('numPSIIs%i timeSteps%i trialsNum%i size%i layers.svg' % (numPSIIs, timeSteps, trialsNum, size))
#plt.savefig(projectPath + fileName, width = 30, height = 8)
#plt.close()
#plt.show()

numPSIIs = 1000
timeSteps = 1000
trialsNum = 100
size = 1
layer = 1
#photonFluxList = range(200,1001,200)
photonFluxList = [500]

projectPath = '/home/ludwik/Documents/python/Monte-Carlo/Layers and size/'
def Simulate(numPSIIs, timeSteps, trialsNum, photonFluxList, size, layer):
    for light in photonFluxList:
        print 'Light: %i' % light
        simulatingLeaf(numPSIIs = numPSIIs, timeSteps = timeSteps, trialsNum = trialsNum, size = size, photonFlux = light, layers = layer)
    fileName = str('numPSIIs%i timeSteps%i trialsNum%i size%.2f layers%i lightDependency.svg' % (numPSIIs, timeSteps, trialsNum, size, layer))
    plt.savefig(projectPath + "QA" + fileName, width = 30, height = 8)
    #plt.savefig(projectPath + fileName, width = 30, height = 8)    
    plt.close()

Simulate(numPSIIs, timeSteps, trialsNum, photonFluxList, size, layer)
selectedTimepoint1 = selectedTimepoint 
for time in range(0, len(selectedTimepoint1)):
    selectedTimepoint1[time] /= float(photonFluxList[time])
    selectedTimepoint1[time] /= trialsNum


selectedTimepoint = []
size = 0.5
layer = 1
Simulate(numPSIIs, timeSteps, trialsNum, photonFluxList, size, layer)
selectedTimepoint2 = selectedTimepoint 
for time in range(0, len(selectedTimepoint1)):
    selectedTimepoint2[time] /= float(photonFluxList[time])
    selectedTimepoint2[time] /= trialsNum


#plt.plot(photonFluxList, selectedTimepoint1, label = "PSII size = 1 layers = 1")
#plt.plot(photonFluxList, selectedTimepoint2, label = "PSII size = 0.5 layers = 1")#

#plt.legend(loc = "best", fontsize = 'x-small')
#plt.xlabel("PPFD")
#plt.ylabel("F[t=10]/PPFD")
#fileName = str('Ft to PPFD Normalised LightDependency trialsNum = %i.svg' % trialsNum)
#plt.savefig(projectPath + fileName, width = 20, height = 7)
#plt.close()












#################THROUBLESHOOTING#########################

#timeSteps = 10
#trialsNum = 10
#size = 1
#photonFluxList = range(500,1001,500)
#photonFlux = 1000
#numPSIIs = 10
#PSIIs = []
#for nr in range(0, numPSIIs):
#    PSIIs.append(PSII(size = size, state = "open", photonFlux = photonFlux, leafArea = 1000, probabilityFluorescence = 0.2, probablilityAnihilation = 0.01))
#    
