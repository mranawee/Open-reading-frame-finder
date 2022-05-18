# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 13:02:08 2020

@author: mano_
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 20:12:13 2020

@author: mano_
"""

#!/usr/bin/env python3
# Name: Mano (mranawee)
# Group Members: 

class ProteinParam :
    '''
    These tables will be used in methods for calculating:
     - molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
     - absorbance at 280 nm (aa2abs280)
     - pKa of positively charged Amino Acids (aa2chargePos)
     - pKa of negatively charged Amino acids (aa2chargeNeg)
     - constants aaNterm and aaCterm for pKa of the respective termini
     
     Example)
     Input(protein string): VLSPADKTNVKAAW
     
     Expected Output: Number of Amino Acids: 14
                      Molecular Weight: 1499.7
                      molar Extinction coefficient: 5500.00
                      mass Extinction coefficient: 3.67
                      Theoretical pI: 9.88
                      Amino acid composition:
                      A = 21.43%
                      C = 0.00%
                      D = 7.14%
                      E = 0.00%
                      F = 0.00%
                      G = 0.00%
                      H = 0.00%
                      I = 0.00%
                      K = 14.29%
                      L = 7.14%
                      M = 0.00%
                      N = 7.14%
                      P = 7.14%
                      Q = 0.00%
                      R = 0.00%
                      S = 7.14%
                      T = 7.14%
                      V = 14.29%
                      W = 7.14%
                      Y = 0.00%

    '''

 
    aa2mw = { #Molecular weight of each AA
             'A': 89.093, 'G': 75.067, 'M': 149.211, 'S': 105.093,
            'C': 121.158,
             'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103,
            'I': 131.173,
             'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188,
            'Q': 146.145,
             'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201,
            'Y': 181.189
             }
 
    mwH2O = 18.015 
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}
    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34
 
    def __init__ (self, protein):
        '''
        This method initialized the protein attribute as well as the Amino Acid Composition dictionary that will be used for
        calculating the AA composition in the aaComposition method.
        '''
        self.protein = protein.upper() #Making the inputted protein string upper case
        
        self.aaCompDict = {key:0 for key in ProteinParam.aa2mw.keys()} #Taking keys from the aa2mw dictionary
        
        '''
        This for loop will iterate through the inputted protein string and add to the composition dictionary the number of each
        AA.  The main method will first call the aaComp method, which will just return the dictionary with the counts of each
        AA.  The main method will then do the appropriate calculations for the AA composition.
        '''
        for aa in self.protein:
            if aa in self.aaCompDict:
                self.aaCompDict[aa] +=1
    
    def aaCount (self):
        '''
        This method will count the number of amino acids from the inputted string. What's very important here is that
        this counts only valid amino acids from the dictionary already made and no spaces.  To make sure this happens,
        I added an if statement within the for loop that only allows counting if the AA is already found in the dictionary.
        '''
        count = 0
        for aa in self.protein:
            if aa in ProteinParam.aa2mw.keys(): #Count won't include invalid characters or spaces
                count += 1
        
        return (count)
 
    def pI (self): 
        '''
        This method calculates the theoretical pI by finding the pH that yields the closest net charge to 0.  This 
        for loop iterates over all pH values to find the charge closest to 0.  To find the best pH accurate to 2 decimal
        places, I made a range of 1401 followed by dividing the pH value by 100.  The value closest to 0 as a pH was the 
        best value.  Each pH that is iterated then calls the charge method.  
        '''
        bestCharge = 100000
        #for pH in range(14+1):
        for pH in range(1401):
            pH = pH / 100 #Makes the calculation more precise by increments of 0.01(Help by Dennis Mulligan)
            thisCharge = self._charge_(pH)
            if abs(thisCharge) < abs(bestCharge):
                bestCharge = thisCharge
                bestPH = pH
        
        return(bestPH)
           
    def aaComposition (self) :
        '''
        The dictionary is already made for calculating AA composition.  The for loop is already made in the init
        method for iteration and the proper calculations are done in the main method.  
        '''
        return (self.aaCompDict)
        
    def _charge_ (self,pH):
        '''
        This method is never used directly by the main method.  I did a for loop for positive and negative charges,
        followed by an if statement for each separately.  I then wrote equations for calculating the charges
        of the termini.
        '''
        pos = 0
        neg = 0
        for aa in self.protein:
            
            if aa in ProteinParam.aa2chargePos.keys():
                posTerm1 = float(10**ProteinParam.aa2chargePos.get(aa)) #numerator
                posTerm2 = float(posTerm1 + 10**pH) #denominator
                pos += float(posTerm1 / posTerm2)
            
            
            elif aa in ProteinParam.aa2chargeNeg.keys():
                negTerm1 = float(10**pH) #numerator
                negTerm2 = float(10**ProteinParam.aa2chargeNeg.get(aa) + 10**pH) #denominator
                neg += float(negTerm1 / negTerm2)
            
            else:
                netCharge = 0
               
        NTerminusTerm1 = float(10**ProteinParam.aaNterm)
        NTerminusTerm2 = float((10**ProteinParam.aaNterm) + 10**pH)
        NTerm = NTerminusTerm1 / NTerminusTerm2
        
                              
        CTerminusTerm1 = float((10**pH)) 
        CTerminusTerm2 = float(ProteinParam.aaCterm + 10**pH)
        CTerm = CTerminusTerm1 / CTerminusTerm2
        
        
        netCharge = float((pos + NTerm) - (neg + CTerm))
        
        return(netCharge)
        
    def molarExtinction (self):
        '''
        From the AA composition dictionary, I can use that to calculate the molar extinction coefficient of the 
        inputted protein.  From that dictionary, we will know the composition of Ys, Ws, and Cs.  The coefficients
        for each of those are already provided, and are calculated.  
        
        I first set this method up by taking the count of each of the three AAs in the protein.  
        
        In the next three lines, I multiplied the counts by each molar coefficient and added it all
        to get the total.
        '''
        Y = self.aaCompDict['Y'] #Count of each of the AAs
        W = self.aaCompDict['W']
        C = self.aaCompDict['C']
        
        molarY = Y * self.aa2abs280['Y']#Count multiplied by the molar coefficient of each
        molarW = W * self.aa2abs280['W']
        molarC = C * self.aa2abs280['C']
        
        totalMolar = molarY + molarW + molarC #Total molar extinction coefficient calculated
        return(totalMolar)

    def massExtinction (self):
        '''
        Mass extinction coefficient is calculated by dividing the molar extinction coefficient calculated by the 
        molecular weight. Molecular weight and molar extinction is called in this method.
        '''
        myMW = self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0
 
    def molecularWeight (self):
        '''
        This method used to calculate the molecular weight included a for loop that started by iterating through the 
        inputted string. If each amino acid inputted was found in the aa2mw dictionary, the value of the AAs mw was calculated
        and added to the next AA mw found.  The molecular weight of water was then subtracted.
        '''
        mw = 0
        for aa in self.protein:
            if aa in ProteinParam.aa2mw.keys():
                mw += ProteinParam.aa2mw.get(aa) - ProteinParam.mwH2O 
        newMW = sum([mw]) + ProteinParam.mwH2O
        return(newMW)
        
class NucParams:
    '''
    This nucParams class will contain 3 dictionaries that will store the counts for making appropriate calculations.  The 
    rnaCodonTable maps each AA to a codon, and the other dictionaries can pull data from this table for its own initialization.
    '''

    rnaCodonTable = {

    # RNA codon table

    # U

    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU

    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC

    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA

    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG

    # C

    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU

    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC

    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA

    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG

    # A

    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU

    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC

    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA

    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG

    # G

    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU

    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC

    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA

    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG

    }

    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()} 

    
    def __init__ (self, inString = ''):
      '''
        This method initialized the new sequence attribute.  All three dictionaries were included here as well, which include 
        codon composition(codonCompDict), nucleotide composition(nucleoCompDict), and amino acid composition(aaCompDict).
        codonCompDict took the keys from rnaCodonTable, iterating through it.  nucloCompDict will hold the counts of each 
        nucleotide.  aaCompDict will hold the counts of each amino acid. 

        The methods returning the codon, nucleotide, and amino acid compositions just need to return the appropriate
        dictionaries as all the iterations are done in the init and addSequence methods.
      '''

      self.codonCompDict = {codon:0 for codon in NucParams.rnaCodonTable.keys()} #Count of each codon

      self.nucleoCompDict = {'A':0, 'C':0, 'G':0, 'T':0, 'U':0, 'N':0} #Count of each nucleotide

      self.aaCompDict = { #Count of each amino acid
                    'A': 0, 'G': 0, 'M': 0, 'S': 0, 'C': 0,
                    'H': 0, 'N': 0, 'T': 0, 'D': 0, 'I': 0,
                    'P': 0, 'V': 0, 'E': 0, 'K': 0, 'Q': 0,
                    'W': 0, 'F': 0, 'L': 0, 'R': 0, 'Y': 0,
                    '-':0
                    } 

      self.addSequence(inString) #Initializing the addSequence method.


    def addSequence (self, thisSequence):
      '''
      This method takes a sequence and iterate through it.  I used the range method from the 0th position to the 
      end of the sequence, and iterating by threes to count codons and add the count to the codon dictionary.  The
      if statement makes sure that the codon is counted only if found in the codon dictionary.  

      I then used the RNA codon table for translating each codon to amino acid and adding to the amino acid count. 
      '''
          
      for nuc in thisSequence.upper():
        if nuc in self.nucleoCompDict:
          self.nucleoCompDict[nuc] += 1

     
      rnaAA = thisSequence.replace('T', 'U')
      for p in range(0, len(rnaAA), 3):
        codon = rnaAA[p:p+3]
        if codon in self.codonCompDict:
          self.codonCompDict[codon] += 1 #Counting each codon and adding each count to codonCompDict
          aa = NucParams.rnaCodonTable[codon]
          self.aaCompDict[aa] += 1 #Counting each amino acid and addind each count to aaCompDict
    

    def aaComposition(self):
      '''
      This method returns the composition of the amino acids in the inputted sequence
      '''
      return (self.aaCompDict)

        
    def nucComposition(self):
      '''
      This method returns the nucleotide composition of the amino acids in the inputted sequence.
      '''
      return (self.nucleoCompDict)
       

    def codonComposition(self):
      '''
      This method returns the codon composition of the amino acids in the inputted sequence.
      '''
      return (self.codonCompDict)


    def nucCount(self):
      '''
      This method returns the number of nucleotides in the inputted sequence.  I just took the values from the nucleotide
      dictionary and added them all.
      '''
      return sum(self.nucleoCompDict.values())

    
    def gcContent(self):
      '''
      This method calculates the gc content by finding the sum of Cs and Gs in the nucleotide dictionary and dividing by the sum,
      which nucCount already calculates.
      '''
      gCount = self.nucleoCompDict['G']
      cCount = self.nucleoCompDict['C']

      gcContent = ((gCount + cCount) / sum(self.nucleoCompDict.values())) * 100

      return (gcContent)

import sys #module is used in doOpen function

 

class FastAreader :

    

    def __init__ (self, fname=''):

        '''contructor: saves attribute fname '''

        self.fname = fname

            

    def doOpen (self):

        if self.fname is '':

            return sys.stdin

        else:

            return open(self.fname)

 

    def readFasta (self):

        

        header = ''

        sequence = ''

        

        with self.doOpen() as fileH:

            

            header = ''

            sequence = ''

 

            # skip to first fasta header

            line = fileH.readline()

            while not line.startswith('>') :

                line = fileH.readline()

            header = line[1:].rstrip()

 

            for line in fileH:

                if line.startswith ('>'):

                    yield header,sequence

                    header = line[1:].rstrip()

                    sequence = ''

                else :

                    sequence += ''.join(line.rstrip().split()).upper()

 

                 

        yield header,sequence
        
class OrfFinder:
    '''
    This class will take an inputted FASTA file and return a list of lists.  
    '''
    
    def __init__(self, sequence):

        '''
        The inputted sequence, list of start codons, stop codons, the dictionary of the 
        complement strand and the range of frames were instantiated below.
        '''
        self.sequence = sequence
        
        self.startCodons = ['ATG','GTG','CTG'] #list
        self.stopCodons = ['TGA', 'TAA', 'TAG'] #list
        self.complementDict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
        self.frames = range(3) #frames are indexed 0,1,2
        
    def orf(self,isReverse): 
        '''
        The main algorithm for finding genes lies in this method.  It iterates for both strands, finding all 6 ORFs.

        '''

        '''
        If its the complement strand, call the complement method that matches each nucleotide to the corresponding strand.
        '''

        if isReverse:
            self.sequence = self.getComplement()


        self.openReadingFrames = [] #list 
        
        '''
        for iterating through the range of 3 frames:
          There are no starts it has come across and no stop codons
        '''
        for frame in self.frames:
            startList = [0]
            foundStop = False
            '''
            Iterate through the codons:
              add to the list of start positions if you come across a start codon.

            Iterate through the codons to so the same for stop codons. Once you find the end,
            add three to get the end position of the stop codon  
            '''
            for index in range(frame, len(self.sequence),3):
                codon = self.sequence[index:index+3]
                if codon in self.startCodons:
                    startList.append(index)

                elif codon in self.stopCodons:
                    end = index + 3
                    foundStop = True
                    '''
                    iterate through the start list and find the frame position.  If it's on the reverse 
                    strand, make the frame position negative.  
                    '''
                    for start in startList:
                        tempFrame = frame + 1
                        if isReverse:
                            tempFrame *= -1
                        '''
                        set the value of tempstart
                        '''
                        tempStart = start + 1
                        if not start == 0:
                            tempStart -= frame
                        '''
                        Set the value of length: end - start, adding one to have the proper length.
                        Wash the startlist after.
                        '''
                        length = end - tempStart + 1
                        result = [tempFrame, tempStart, end, length] #brackets indicate LIST
                        self.openReadingFrames.append(result)
                    startList = []
            '''
            If there's no stop just output the rest of the gene.
            '''
            if not foundStop:
                length = len(self.sequence)
                tempFrame = frame + 1
                if isReverse:
                    tempFrame *= -1
                
                 
                
            
            if startList:
              result = [0] 
              self.openReadingFrames.append([tempFrame, startList[0], len(self.sequence), (len(self.sequence) - startList[0])])  
        '''
        Return the list of lists: List of the frame, start position, end position, and length.
        '''    
        return (self.openReadingFrames)



    def getComplement(self):
        '''
        Returns the complementary strand by iterating through the dictionary and mapping each nucleotide to it's complementary.
        '''
        reverseSequence = self.sequence[::-1]
        comp = ''
        for nuc in reverseSequence:
            comp += self.complementDict.get(nuc)
        return (comp)
    
    def getFinal(self,forward, backward):
        '''
        Sort the data by length and then by frame
        '''
        finalList = sorted(forward + backward, key = lambda x: (-x[3],x[1])) #3rd index is length and 1st index is sort
        return(finalList)

   