# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 08:17:47 2020

@author: mano_
"""

from sequenceAnalysis import FastAreader, OrfFinder

class CommandLine() :
    '''
    Handle the command line, usage and help requests.
    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various
    argument options, a standard usage and help.
    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool 
    was set as an option using add_argument, then myCommandLine.args.requiredbool 
    will name that option.
    '''
    
    def __init__(self, inOpts=None) :
         '''
         Implement a parser to interpret the command line argv string using 
         argparse.
         '''
         import argparse
         
         self.parser = argparse.ArgumentParser(
                 description = 'Program prolog - a brief description of what this thing does',
                 epilog = 'Program epilog- some other stuff you feel compelled to say',
                 add_help = True, #default is True
                 prefix_chars = '-',
                 usage = '%(prog)s[options] -option1[default] <input >output'
                 )
         
         #self.parser.add_argument('inFile', action = 'store', 
                                  #help='input file name')
         
         #self.parser.add_argument('outFile', action = 'store', 
                                  #help='output file name')
         
         self.parser.add_argument('-lG', '--longestGene', action = 'store', 
                                  nargs='?', const=True, default=False, 
                                  help='longest Gene in an ORF')
         
         self.parser.add_argument('-mG', '--minGene', type=int, 
                                  choices= (100,200,300,500,1000), default=100, 
                                  action = 'store', help='minimumGene length')
         
         self.parser.add_argument('-s', '--start', action = 'append', default = ['ATG'],
                                  nargs='?', help='start Codon') #allows multiple list options
         
         self.parser.add_argument('-t', '--stop', action = 'append', 
                                  default = ['TAG','TGA','TAA'], nargs='?', 
                                  help='stop Codon') #allows multiple list options
         
         self.parser.add_argument('-v', '--version', action='version',
                                  version='%(prog)s 0.1')
         
         if inOpts is None :
             self.args = self.parser.parse_args()
             
         else:
             self.args = self.parser.parse_args(inOpts)
             
def main(inCL = None): #if CL = None
    '''
    This method uses the fastAReader and OrfFInder class to return the appropriate format the the 
    algorithm calculated.
    '''
    
    if inCL is None:
        myCommandLine = CommandLine()
             
    else:
        myCommandLine = CommandLine(inCL)
    

    myReader = FastAreader()

    
    for header, sequence in FastAreader().readFasta() :
        newHeader = header.strip(">").rstrip()
       
        myORF = OrfFinder(sequence) 
        
    

        #need to call the total orf list
        forwardData = myORF.orf(False)
        #complement = myORF.getComplement()
        backwardData = myORF.orf(True)
        final = myORF.getFinal(forwardData, backwardData)
    
        for framedata in final:
            #This is the proper format needed
            if framedata[0] > 0:
                print("+{}..\t{}\t{}\t{}".format(framedata[0], framedata[1], framedata[2], framedata[3]))
            else:
                print("{}..\t{}\t{}\t{}".format(framedata[0], framedata[1], framedata[2], framedata[3]))


if __name__ == "__main__":
    main()         
         
         
         
