'''
Created on Sep 5, 2012

@author: dominic
'''


##( (  (One:0.2,Two:0.3):0.3, (Three:0.5,Four:0.3):0.2  ):0.3 ,Five:0.7):0.0;
##
##           +-+ One
##        +--+
##        |  +--+ Two
##     +--+
##     |  | +----+ Three
##     |  +-+
##     |    +--+ Four
##     +
##     +------+ Five

class chObj:
    def __init__(self, i, realLabel):
        self.i = i
        self.reallabel = realLabel

class child:
    def __init__(self, obj, distanceToNearestNode = 1.0):
        self.id = obj.index ##this is the SM index #This means that you have an 
                            ###object that has an interface that yields index and reallabel
        self.label = obj.reallabel ##We want this in newick
        self.obj = obj
        self.dist = distanceToNearestNode
        
    def export(self):
        elStr = str(self.label) + ':'
        return elStr
            
class Node:
    def __init__(self, idNum, terminal = 0, elemLst = []):
        self.id = idNum
        self.Proximal = 1
        self.Terminal = terminal
        self.childTypes = []
        self.children = []
        for ch in elemLst:
            self.children.append(ch)
        self.upstream = None
        self.upstreamDist = 0
        self.downstream = None
        
    def childIDLst(self):
        childid = []
        for e in self.children:
            if 'child' in str(e):
                childid.append(e.id)
            else:
                childid.extend(e.childIDLst())
        if self.downstream:
            childid.extend(self.downstrean.childIDLst())
        return childid
    def display(self):
        for elem in self.children:
            print elem.id
    def export(self):
        ##( (  (One:0.2,Two:0.3):0.3, (Three:0.5,Four:0.3):0.2  ):0.3 ,Five:0.7):0.0;
        typeLst = []
        for node in self.children:
            if not('child' in str(node)):
                typeLst.append(node.Terminal)
            else:
                typeLst.append(1)
        if len(typeLst) == 1 and typeLst[0] ==1:
            outStr = self.children[0].export() + str(self.upstreamDist)
            return outStr
        if len(typeLst) ==2:
            nodeA, nodeB = self.children
            outStr = '(' + nodeA.export() + ',' + nodeB.export() + '):' + str(self.upstreamDist)
            return outStr
            
        if 'child' in str(self.children[0]):
            outStr = '('
            for child in self.children:
                child.dist = 1.0
                outStr = child.export() + ','
            outStr = outStr[:-1] + '):' + str(self.upstreamDist) + ','
            return outStr
        else:
            outStr = '('
            for node in self.children:
                outStr = outStr + node.export() + ','
        outStr = outStr[:-1] + '):' + str(self.upstreamDist)
                    
        return outStr


class ClusterObject:
    def __init__(self, SM, elementList = [], root = 'heirtest'):
        self.SM = SM
        self.leafs = []
        self.Root = None
        self.NodeLst = [] ##mutable as you combine Nodes
        self.NodeDict = {} ##Just Tells you where a child is
        self.elementList = elementList  ##You need an element list with objects that have an interface that yields .index and .reallabel
                                        ##and the .index has to be the same indeax as is referenced within SM
        self.NodeSM = SimilarityMatrix()
        self.n = len(elementList)
        self.ProximalNodes = self.n
        self.PairDB = []
        self.root = root
        self.CompareFunction = lambda x: x
        self.GenerateNodeLst()
        
    def GenerateNodeLst(self):
        i = 0
        for element in self.elementList:
            self.leafs.append(child(element))
            self.leafs[-1].dist = 0.1
            i +=1
        i = 0
        for ch in self.leafs:
            self.NodeLst.append(Node(i, terminal = 1, elemLst=[ch]))
            i +=1
        self.GenerateNodeDict()
    def GenerateNodeDict(self):
        self.NodeDict = {}
        for n in range(len(self.NodeLst)):
            node = self.NodeLst[n]
            node.id = n ##This renumbers the nodes
            elements = node.childIDLst()
            for e in elements:
                self.NodeDict[e] = n
    def MakeClustersOFN(self, clustList = None, clthreshold = None, SM = None):
        if clustList == None:
            clustList = self.MakeClusters(clthreshold = clthreshold, SM = SM)
        branches = len(clustList)
        
        if branches == 1:##This is if there is only a single cluster containing all children
            p = clustList[0]
            ofn = 'Clustered_' + self.root + '.dnd'
            OFN = FA.path([UD,ofn], fileType = 'raw')
            OFN.outDB(p)
            return OFN
        OFNLst = US.FilePathLstObj()
        notclustered = []
        for i in range(branches):
            p = clustList[i]
            if not(p): continue
            if '(' in p[:1]:
                ofn = 'Clustered_' + str(i) + '_' + self.root + '.dnd'
                OFN = FA.path([UD,ofn], fileType = 'raw')
                OFN.outDB(p)
                OFNLst.add([OFN])
            else:
                notclustered.append(p)
        if notclustered:
            ofn = 'NotClustered_' + '_' + self.root + '.txt'
            OFN = FA.path([UD,ofn], fileType = 'tab')
            OFN.outDB(notclustered)
            OFNLst.add([OFN])       
        return OFNLst
    def MakeClusters(self, clthreshold = None, SM = None):
        if SM <> None:
            raise 'HeirarchicalClusterObject:MakeClusters: trying to use a deprecated pathway that was removed 111108'
##        if SM == None:
##            SM = self.SM
##            ##Allows you to use the same settings for multiple SM runs
        if self.SM.Similarity: ##Must be distance Not similarity Don't rely on this
            print 'HeirarchicalClusterObject:MakeClusters: You have chosen to let me decide how to normalize your SM. BAD IDEA!@!'
            self.SM.Normalize(flipOrder =1)
        #else:
            #SM.Normalize()
        self.SM.display()
        
        numClusters = self.n
        print 'HeirarchicalClusterObject:MakeClusters: processing %i elements' %numClusters
        ## This takes care of cases in which there is an infinite distanace between two clusters
        ## You have two choices: a. set threshold to default None, they get combind with a large distance 10 max
        ## b. set threshold to 0.0 (or any other float) and if the distance threshold is not met, then
        ##clusters will be split off into separate files
        if clthreshold == 0.0:
            clthreshold == self.SM.maxVal + 1 ##This means you will not combine infinite distanced elements
        if clthreshold == None:
            if self.SM.maxVal > 100:
                print 'HeirarchicalClusterObject:MakeClusters: Normalizing all data between 0 and 100'
                self.SM.Normalize(newMin = 0.0, newMax = 100)##I had one SM that had numbers much larger than 100
            clthreshold  =  self.SM.maxVal + self.SM.maxVal/100##101 ##This means you will combine anything as long as it is at the top(note the extra 1)
        cycle = 0
        self.ProximalNodes = self.n
        while self.ProximalNodes > 1:
            cycle +=1
            move = None
            print 'Cycle %i num Proximal Nodes %i' %(cycle, self.ProximalNodes)
            pairDB = self.GenerateNodeSimilarityMatrix()
            NodeA, NodeB, dist = pairDB[0]
            print '%s, %s have this dist %s which must be less than %s' %(str(NodeA), str(NodeB), str(dist), str(clthreshold))
            if dist < clthreshold: #clthreshold is the MaxDist To allow combination
                print 'Combining Node %s and Node %s distance is %2.2f which is < clthreshold %0.2f' %(str(NodeA.id), str(NodeB.id), dist, clthreshold)
                self.BranchTwoNodes(NodeA, NodeB)
                move = 1
            else:
                print 'REJECT Node combination %s and Node %s distance is %2.2f which is > clthreshold %0.2f' %(str(NodeA.id), str(NodeB.id), dist, clthreshold)
            if not(move):
                self.clusterstats()
                break  
        clustList = self.NewickExport()
        return clustList
    
    def GenerateNodeSimilarityMatrix(self):
        #This has to be in the same direction as input SM
        self.NodeSM = SimilarityMatrix(Similarity = self.SM.Similarity)
        self.pairDB = []
        for qnode in self.NodeLst:
            for snode in self.NodeLst:
                if qnode.id == snode.id: continue
                if qnode.Proximal and snode.Proximal:
                    qLst = qnode.childIDLst()
                    sLst = snode.childIDLst()
                    dist = self.SM.AvgGroupDistance(qLst, sLst)
                    if dist ==None:
                        print '\n\n\nAdjust your PM/MM ratio, you have simply run out of useable probes'
                        print 'Two of your nodes are at an infinite distance to each other'
                        print 'I dont handle such situations yet, standby for the dump\n\n'
                        raise 'GenerateNodeSimilarityMatrix: You\'ve hit a deprecated feature of this proc.'
##                    if dist == None:
##                        if self.SM.Similarity:##Note the node SM but the actual child SM!!
##                            dist = max(0.01, self.SM.minVal-1)
##                        else:
##                            dist = self.SM.maxVal + 1
                    self.NodeSM[qnode.id][snode.id] = dist
                    self.NodeSM[snode.id][qnode.id] = dist
                    self.pairDB.append([qnode, snode, dist])
                    self.pairDB.append([snode, qnode, dist])
        
        self.pairDB.sort(lambda x,y:cmp(x[-1], y[-1]))

        if self.NodeSM.Similarity:
            raise 'GenerateNodeSimilarityMatrix: You\'ve hit a deprecated feature of this proc.'
            #self.pairDB.reverse()##Thus placing the closest sequences at the top of the proc list
            ##On second thought no SM is allowed into the function that isn't a distance function
        print 'GenerateNodeSimilarityMatrix(self):'
        self.NodeSM.display()
        return self.pairDB
    def clusterstats(self):
        print 'This is the distance distribution:'
        self.SM.ShowHistogram()
        self.CountProximalNodes()
        print 'This clusterObj has %i elements in %i Proximal Nodes' %(self.n, self.ProximalNodes)
    def __getitem__(self,j):
        if 'Node' in str(j):
            return j
        else:
            return self.NodeLst[j.idNum]
    def NewickExport(self):
        self.CountProximalNodes()
        outStr = ''
        for node in self.NodeLst:
            nodeexport = None
            if node.Proximal:
                nodeexport = str(node.export())
                outStr = outStr + nodeexport
            if self.ProximalNodes > 1 and nodeexport:
                outStr = outStr + '<separation>'
        clustList = outStr.split('<separation>')
        for c in range(len(clustList)):
            p = clustList[c]
            if p:
                clustList[c] = clustList[c] + ';'
        return clustList
    def BranchTwoNodes(self, NodeA, NodeB):
        if not('Node' in str(NodeA)):
               NodeA = self[NodeA]
               NodeB = self[NodeB]
        #NodeC gets written on top of Node A and Node B
        self.n +=1
        NodeC = Node(self.n, elemLst = [NodeA, NodeB])
        NodeC.childTypes = [NodeA.Terminal, NodeB.Terminal]
        #Node A and B are no longer terminal they lose this status
        NodeA.upstream = NodeC; NodeA.Proximal = 0
        NodeB.upstream = NodeC; NodeB.Proximal = 0
        #CalculateInterNode Distance
        NodeA_IDLst = NodeA.childIDLst()
        NodeB_IDLst = NodeB.childIDLst()
        NodeAB_Dist = self.SM.AvgGroupDistance(NodeA_IDLst, NodeB_IDLst)# or 1.0
        NodeA_WorldDist = self.SM.AvgGroupDistance(NodeA_IDLst)# or 1.0
        NodeB_WorldDist = self.SM.AvgGroupDistance(NodeB_IDLst)# or 1.0 ##deleted this 111108
        NodeA.upstreamDist = NodeAB_Dist*(NodeA_WorldDist/NodeB_WorldDist)
        NodeB.upstreamDist = NodeAB_Dist*(NodeB_WorldDist/NodeA_WorldDist)
        self.NodeLst.append(NodeC)
        self.CountProximalNodes()
    def CountProximalNodes(self):
        ##Nodes that have yet to be merged
        n = 0
        for node in self.NodeLst:
            if node.Proximal:
                n +=1
        self.ProximalNodes = n
    def Brothers(self, qID):
        node = self.NodeLst[self.NodeDict[qID]]
        while 1:
            brothers = node.childIDLst()
            if qID in brothers:
                return brothers
            node = node.downstream
            if not(node):
                print 'Error query Not found Dictionary is wrong'
                return None
            
    def Cousins(self, qID):
        node = self.NodeLst[self.NodeDict[qID]]
        cousins = []
        found = None
        while 1:
            brothers = node.childIDLst()
            cousins.extend(brothers)
            node = node.downstream
            if not(node):
                break
        if qID in cousins:
            return cousins
        else:
            print 'NodeDict is wromng Error Cannot find ' + str(qID)
            return None
    def Display(self):
        for node in self.NodeLst:
            print 'id %s Proximal %i' %(str(node.id), node.Proximal)
            print str(node.childIDLst())
    
    
if __name__ == "__main__":
    print 'Newick'
    SM = [[0,1,2,2,3],
          [1,0,2,2,3],
          [2,1,0,2,3],
          [4,1,2,0,3],
          [1,1,2,1,0],]
    i = 0
    elemLst = ['One', 'Two', 'Three', 'Four', 'Five']

    
    CO = ClusterObject(SM, elementList = elemLst, root = 'heirtest')
        
