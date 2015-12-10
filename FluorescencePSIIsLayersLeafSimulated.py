import random
import matplotlib.pyplot as plt

def chunks(l, numberOfGroups):
    """
    Yield n successive chunks from l.
    """
    chunkLength = len(l)/numberOfGroups
    for i in xrange(0, len(l), chunkLength):
        yield l[i:i + chunkLength]

class PSII(object):
    """
    Representation of a PSII particle
    """   

    def __init__(self, layer = 1, size = 1, state = "ground", photonFlux = 1000, leafArea = 1000, probabilityExcited = 0.1, probabilityFluorescence = 0.05, probabilityRadiationless = 0.1, probablilityAnihilation = 0.01):
        """

        Initialize a PSII instance, saves all parameters as attributes
        of the instance.

        size: int representing size of the PSII

        state: str representing the state of PSII accepted values "ground", "closed ground" or "closed excited"
        
        absCrossection: float calculated from the the PSII size normalized to the leafArea

        probabilityAbsorbed: float 

        """
        self.layer = layer
        self.photonFlux = photonFlux
        self.leafArea = leafArea
        self.size = size
        self.state = state
        self.absCrossection = self.size/float(self.leafArea)
        self.probabilityAbsorbed = self.photonFlux * self.absCrossection

        self.probabilityExcited = probabilityExcited
        self.probabilityFluorescence = probabilityFluorescence
        self.probabilityRadiationless = probabilityRadiationless       
        self.probabilityDecay = probabilityFluorescence + probabilityRadiationless

        self.probablilityAnihilation = probablilityAnihilation


    def updatePhotonFlux(self, PhotonFlux):
        """
        Updates the photonFlux and dependent variables: absCrossection, probabilityAbsorbed
        """    
        self.photonFlux = PhotonFlux
        self.absCrossection = self.size/float(self.leafArea)
        self.probabilityAbsorbed = self.photonFlux * self.absCrossection

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
        if light == "off":
            if self.state == "closed excited":
                return self.doesFluoresce()

        if light == "on":         
            Absorbed = False
            if self.state == "closed excited":                 #Radiationless decay if the light is not absorbed   
                return Absorbed, self.doesFluoresce()

            if random.random() <= self.probabilityAbsorbed:
                Absorbed = True
                if self.state == "ground":
                    self.state = "closed ground"
                    return Absorbed, False
                if self.state == "closed ground":
                    if random.random() <= self.probabilityExcited:
                        self.state = "closed excited"
                        return Absorbed, self.doesFluoresce()
                if self.state == "closed excited":
                    if random.random() <= self.probablilityAnihilation:
                        self.state = "closed ground"
                        return Absorbed, False
                else:
                    return Absorbed, False
            else:
                return False, False  

class Layer(object):
    """
    Representation of layers in the leaf.
    """    
    def __init__(self, PSIIs, layersNumber):
        """
        Initialization function, saves the PSIIs 

        """
        self.PSIIs = PSIIs
        self.layersNumber = layersNumber
        self.FCount = 0
        self.AbsorbedCount = 0

    def updateAbsorbedFluorescedCount(self, Absorbed, Fluoresced):
        if Absorbed == True:
            self.AbsorbedCount += 1
        if Fluoresced == True:
            self.FCount += 1

    def updatePSIIs(self, light, PhotonFlux, PSIIs = None):
        if PSIIs == None:
            PSIIs = self.PSIIs
        self.FCount = 0
        self.AbsorbedCount = 0
        for psii in PSIIs:
            psii.updatePhotonFlux(PhotonFlux)
            if light == "on":
                if psii.probabilityAbsorbed < 1:
                    Absorbed, Fluoresced = psii.update(light)
                    self.updateAbsorbedFluorescedCount(Absorbed, Fluoresced)

                else:
                    for excitation in range(0,int(psii.probabilityAbsorbed)):
                        Absorbed, Fluoresced = psii.update(light)
                        self.updateAbsorbedFluorescedCount(Absorbed, Fluoresced)

            if light == "off":
                if psii.update(light) == True:
                    self.FCount += 1

        return self.FCount, self.AbsorbedCount

class Leaf(object):
    """
    Representation of a simplified leaf.
    """    
    def __init__(self, PSIIs, layersNumber):
        """
        Initialization function, saves the PSIIs 

        """
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
        for layer in range(0, self.layersNumber):
            PSIIs = self.getPSIIsInLayer(layer)
            self.Layers.append(Layer(PSIIs, layer))

    def updateLayers(self, light):
        self.totalFluoresced = 0
        self.totalAbsorbed = 0
        PSIIs = self.PSIIs 
        PhotonFlux = PSIIs[0].photonFlux
        for layer in self.Layers:
            Fluoresced, Absorbed = layer.updatePSIIs(light, PhotonFlux)
            self.totalFluoresced += Fluoresced
            self.totalAbsorbed += Absorbed
            PhotonFlux = PhotonFlux - Absorbed + Fluoresced
        return self.totalFluoresced, self.totalAbsorbed

    def howMuchAbsorbedByLayer(self, layer):
        print 'not implemented yet'

selectedTimepoint = []

def simulatingLeaf(numPSIIs = 1000, timeSteps = 100, trialsNum = 1, size = 1, photonFlux = 1000, layers = 1):
    """
    Runs simulations and plots graphs for PSIIs in the leaf.

    """
    global selectedTimepoint
    timepoint = 10
    trialsSum = [0]
    for time in range(1, timeSteps + 1):
        trialsSum.append(0)

    for trial in range(0, trialsNum):
        PSIIs = []                                      #Creating PSIIs
        for nr in range(0, numPSIIs):
            PSIIs.append(PSII(size = size, state = "ground", photonFlux = photonFlux, leafArea = 10000, probabilityExcited = 0.8, probabilityFluorescence = 0.05, probabilityRadiationless = 0.15, probablilityAnihilation = 0.01))
        simulatedLeaf = Leaf(PSIIs, layers)             #Creating the leaf
        simulatedLeaf.assignPSIIToLayers()              #Creating layers in the leaf
        simulatedLeaf.createLayers()

        Fluorescence = [0]

        for time in range(1, timeSteps+1):
            Fluoresced, Absorbed = simulatedLeaf.updateLayers(light = "on")
            #print "Fluoresced: %i Absorbed: %i" % (Fluoresced, Absorbed)
            Fluorescence.append(Fluoresced)
            trialsSum[time] += Fluorescence[time]

        for time in range(1, timeSteps+1):
            Fluorescence[time] /= photonFlux

        #for time in range(timeSteps + 1, timeSteps*2 + 2):
        #    Fluoresced, Absorbed = simulatedLeaf.updateLayers(light = "off")
        #    #print "Fluoresced: %i Absorbed: %i" % (Fluoresced, Absorbed)
        #    Fluorescence.append(Fluoresced)
        #    trialsSum[time] += Fluorescence[time]
        if trial%10 == 0:
            print 'Trial nr: %i' % trial

    selectedTimepoint.append(trialsSum[timepoint])
    plt.plot(range(0,timeSteps + 1), trialsSum, label = "Size: " + str(size) + " PhotonFlux: " + str(photonFlux) + " Layers: " + str(layers) )
    plt.xlim(xmin = 0,xmax = timeSteps + 1)
    return trialsSum

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

numPSIIs = 10000
timeSteps = 10
trialsNum = 1
size = 1
layer = 1
photonFluxList = range(200,1001,200)

def Simulate(numPSIIs, timeSteps, trialsNum, photonFluxList, size, layer):
    for light in photonFluxList:
        print 'Light: %i' % light
        simulatingLeaf(numPSIIs = numPSIIs, timeSteps = timeSteps, trialsNum = trialsNum, size = size, photonFlux = light, layers = layer)
    plt.legend(loc = "best", fontsize = 'small')
    plt.xlabel("Time [a.u.]")
    plt.ylabel("Ft [counts]")
    projectPath = '/home/ludwik/Documents/python/Monte-Carlo/Layers and size/'
    fileName = str('numPSIIs%i timeSteps%i trialsNum%i size%.2f layers%i lightDependency.svg' % (numPSIIs, timeSteps, trialsNum, size, layer))
    plt.savefig(projectPath + fileName, width = 30, height = 8)
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

selectedTimepoint = []


size = 1
layer = 2
Simulate(numPSIIs, timeSteps, trialsNum, photonFluxList, size, layer)


selectedTimepoint3 = selectedTimepoint 
for time in range(0, len(selectedTimepoint1)):
    selectedTimepoint3[time] /= float(photonFluxList[time])
    selectedTimepoint3[time] /= trialsNum

selectedTimepoint = []
size = 0.5
layer = 2
Simulate(numPSIIs, timeSteps, trialsNum, photonFluxList, size, layer)

selectedTimepoint4 = selectedTimepoint 
for time in range(0, len(selectedTimepoint1)):
    selectedTimepoint4[time] /= float(photonFluxList[time])
    selectedTimepoint4[time] /= trialsNum

selectedTimepoint = []

size = 1
layer = 3
Simulate(numPSIIs, timeSteps, trialsNum, photonFluxList, size, layer)

selectedTimepoint5 = selectedTimepoint 
for time in range(0, len(selectedTimepoint1)):
    selectedTimepoint5[time] /= float(photonFluxList[time])
    selectedTimepoint5[time] /= trialsNum

selectedTimepoint = []
size = 0.5
layer = 3
Simulate(numPSIIs, timeSteps, trialsNum, photonFluxList, size, layer)

selectedTimepoint6 = selectedTimepoint 
for time in range(0, len(selectedTimepoint1)):
    selectedTimepoint6[time] /= float(photonFluxList[time])
    selectedTimepoint6[time] /= trialsNum    

plt.plot(photonFluxList, selectedTimepoint1, label = "PSII size = 1 layers = 1")
plt.plot(photonFluxList, selectedTimepoint2, label = "PSII size = 0.5 layers = 1")

plt.plot(photonFluxList, selectedTimepoint3, label = "PSII size = 1 layers = 2")
plt.plot(photonFluxList, selectedTimepoint4, label = "PSII size = 0.5 layers = 2")

plt.plot(photonFluxList, selectedTimepoint5, label = "PSII size = 1 layers = 3")
plt.plot(photonFluxList, selectedTimepoint6, label = "PSII size = 0.5 layers = 3")
plt.legend(loc = "best", fontsize = 'x-small')
plt.xlabel("PPFD")
plt.ylabel("F[t=10]/PPFD")
fileName = str('Ft to PPFD Normalised LightDependency trialsNum = %i.svg' % trialsNum)
plt.savefig(projectPath + fileName, width = 20, height = 7)
plt.close()












#################THROUBLESHOOTING#########################

#timeSteps = 10
#trialsNum = 10
#size = 1
#photonFluxList = range(500,1001,500)
#photonFlux = 1000
#numPSIIs = 10
#PSIIs = []
#for nr in range(0, numPSIIs):
#    PSIIs.append(PSII(size = size, state = "ground", photonFlux = photonFlux, leafArea = 1000, probabilityExcited = 0.1, probabilityFluorescence = 0.2, probabilityRadiationless = 0.05, probablilityAnihilation = 0.01))
#    
