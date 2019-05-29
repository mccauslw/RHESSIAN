import sys, string, re, codecs, operator
from itertools import product
from scipy.stats import binom
from math import exp

letterCodes = [
    'A', 'B', 'C', 'D', 'A', # 1-10
    'C', 'D', 'B', 'B', 'A',
    'A', 'C', 'D', 'B', 'C', # 11-20
    'A', 'D', 'D', 'C', 'A',
    'A', 'C', 'B', 'A', 'B', # 21-30
    'C', 'D', 'B', 'A', 'B',
    'C', 'D', 'D', 'C', 'A', # 31-40
    'B', 'D', 'D', 'B', 'A',
    'C', 'A', 'D', 'C', 'B', # 41-50
    'B', 'B', 'C', 'C', 'A',
    'D', 'D', 'A', 'B', 'C', # 51-60
    'C', 'D', 'D', 'B', 'A',
    'B', 'B', 'A', 'B', 'C', # 61-70
    'D', 'A', 'C', 'D', 'A',
    'C', 'D', 'A', 'B', 'B', # 71-80
    'D', 'C', 'C', 'D', 'D',
    'A', 'A', 'B', 'C', 'D', # 81-90
    'A', 'B', 'C', 'D', 'A',
    'B', 'D', 'C', 'A', 'B', # 91-100
    'C', 'D', 'A', 'B', 'C',
    'D', 'A', 'B', 'C', 'D', # 101-110
    'A', 'B', 'C', 'D', 'A',
    'B', 'B', 'C', 'D', 'D', # 111-120
    'C', 'A', 'C', 'A', 'B',
    'D', 'C', 'A', 'C', 'D', # 121-130
    'A', 'A', 'B', 'D', 'A',
    'C', 'D', 'A', 'A', 'B', # 131-140
    'C', 'D', 'D', 'B', 'C',
    'A' # 141 
]

# Iterate over participant's data files
allChoiceCounts = []

miscDir = 'Misc_Tables/'
distractorFile = open(miscDir + 'distractor_table.txt', 'w')
distractorFile.write('\tnBest\tnWorst\tnBoth\n')
lexicoFile = open(miscDir + 'lexico_table.txt', 'w')
lexicoFile.write('\tnABCDE\tnEDCBA\n')
levelFile = open(miscDir + 'level_table.txt', 'w')
levelFile.write('\tn0\tn1\tn2\tn3\tn4\tn5\tn6\n')
objectFile = open(miscDir + 'object_table.txt', 'w')
objectFile.write('\tA\tB\tC\tD\tE\n')
testFile = open(miscDir + 'test_table.txt', 'w')
testFile.write('\tH\tt\tplo\tphi\n')
MIFile = open(miscDir + 'MI_table.txt', 'w')
#MIFile.write('\t%s\n' % string.join((('x%i\ty%i\tz%i' % (i, i, i)) for i in range(10)), '\t'))
MIFile.write('par\tset\tNAx\tNAy\tNAz\tNxy\tNyz\tNxz\tMIx\tMIy\tMIz\tpixyz\tpixzy\tpiyxz\tpiyzx\tpizxy\tpizyx\n')
tripleFile = open(miscDir + 'triple_table.txt', 'w')
tripleFile.write('\t%s\n' % string.join((('Ax%i\tAy%i\tAz%i\txy%i\tyz%i\txz%i' % (i, i, i, i, i, i)) for i in range(1)), '\t'))
binaryA = [3,5,6,9,10,12,17,18,20,24]

for index, letter in enumerate(letterCodes):
    
    # Construct file name of this partipant's file
    indexAsString = str(index+1)
    fileName = 'Raw_Data/' + indexAsString + letter + '.txt'
    distractorFile.write(indexAsString + '\t')
    lexicoFile.write(indexAsString + '\t')
    levelFile.write(indexAsString + '\t')
    objectFile.write(indexAsString + '\t')
    testFile.write(indexAsString + '\t')
    
    # Read in all lines of this participant's file
    pFile = codecs.open(fileName, 'r', 'utf-16')
    lines = pFile.readlines()

    # First line in participant file is original filename of raw data, for index=0,...,80
    #
    if index<81:
        fName = lines[0]
        # Second line gives variable names
        varNames = lines[1].split('\t')
        startLine = 2
    else:
        varNames = lines[0].split('\t')
        startLine = 1
    
    # Subsequent lines in participant file give choice data for given participant
    # x goes from 0 to 31 and gives subsets of master set of size five
    # y goes from 0 to 4 and gives object indices
    # set matrix to zero
    participantChoiceCounts = [[0 for y in range(5)] for x in range(32)]
    participantChoiceSequences = [[] for x in range(32)]
    nvec = [0 for x in range(4)]
    cvec = [0 for x in range(4)]
    
    # nCountLevels gives, for each i=0,1,2,3,4,5,6, the number of times an object is chosen
    # that number of times in the same choice set
    nCountLevels = [0 for i in range(7)]
    # nCountObjects gives, for each i=0,1,2,3,4, the number of times the participant chooses
    # object i
    nCountObjects = [0 for i in range(5)]

    # Reset counts pertaining to choice in distractor trials
    nWorstAvail = 0
    nWorstChosen = 0
    nBestAvail = 0
    nBestNotChosen = 0
    nBothAvail = 0
    nDoubleError = 0
    
    # Reset counts pertaining to choice in gamble trials
    ABCDE_count = 0
    EDCBA_count = 0

    # Compute choice counts for participant
    # There is one line in a participant's file for every click.
    # A click is ignored if the click is not on a choice object
    # A choice is recorded differently if it is a gamble or a distractor trial
    # Sets {A}, {B}, {C}, {D} and {E} are represented by 1, 2, 4, 8, 16
    for nline, line in enumerate(lines[startLine:]):
        varValues = line.split('\t') # Values of all variables for this log item
        lineDict = dict(zip(varNames, varValues))
        choiceSet = 0
        if lineDict['ClickedGamble'] == '': # Ignored click, no choice object
            continue
        sheetChoice = int(lineDict['ClickedGamble']) - 1
        choice = -1
        # Build choice set, update choice counts
        for i, g in enumerate(['g1', 'g2', 'g3', 'g4', 'g5']):
            for j, gamble in enumerate(['A', 'B', 'C', 'D', 'E']):
                if lineDict[g] == ('Gamble_' + gamble + '.bmp'):
                    trialType = 'gamble'
                    # Gamble is in choice set, update choice set
                    choiceSet += pow(2, j)
                    if sheetChoice == i:
                        # Gamble is the chosen gamble, set choice to index of gamble
                        choice = j
            for j, distractor in enumerate(['1', '2', '3', '4', '5', '6', '7']):
                if lineDict[g] == ('D' + distractor + '.bmp'):
                    trialType = 'distractor'
                    choiceSet += pow(2, j)
                    if sheetChoice == i:
                        # Distractor is the chosen distractor
                        choice = j
        if trialType == 'gamble':
            # trial is a regular trial, with gambles from {A, B, C, D, E},
            # update appropriate choice count
            participantChoiceCounts[choiceSet][choice] += 1
            participantChoiceSequences[choiceSet].append(choice)
            if (choiceSet / pow(2, choice)) == 1:
                EDCBA_count += 1
            if (choiceSet % pow(2, choice)) == 0:
                ABCDE_count += 1
            nCountObjects[choice] += 1
        if trialType == 'distractor':
            # trial is a distractor trial, with gambles from D1, D2, ..., D7
            if (choiceSet & 1)>0:
                nBestAvail += 1
                if (choice != 0):
                    nBestNotChosen += 1
            if (choiceSet & 64)>0:
                nWorstAvail += 1
                if (choice == 6):
                    nWorstChosen += 1
            if ((choiceSet & 1)>0) and ((choiceSet & 64)>0):
                nBothAvail += 1
                if (choice == 6):
                    nDoubleError += 1
            # Next lines give detail verifing that distractor counts are right
            #print choice, bin(choiceSet)
            #print (choiceSet & 1)>0, choice == 0, (choiceSet & 64)>0, choice == 6
            #if ((choiceSet & 64) and (choice == 6)):
            #    print('Distractor 7 presented and chosen')
            #if ((choiceSet & 1) and (choice != 0)):
            #    print('Distractor 1 presented and not chosen')
    # Uncomment next line for more information on distractors
    # print nBestAvail, nBestNotChosen, nWorstAvail, nWorstChosen, nBothAvail, nDoubleError
    distractorFile.write('%d\t%d\t%d\n' % (nBestNotChosen, nWorstChosen, nDoubleError))
    lexicoFile.write('%d\t%d\n' % (ABCDE_count, EDCBA_count))
    for i in range(32):
        nChoices = 0
        for j in range(5):
            nCountLevels[participantChoiceCounts[i][j]] += 1
            nChoices = nChoices + participantChoiceCounts[i][j]
        if not(nChoices in [0,6]):
            print index+1, i, nChoices
    nCountLevels[0] -= 85 # 85 is number of zeros that are associated with unavailable options
    levelFile.write(string.join([str(x) for x in nCountLevels], '\t') + '\n')
    objectFile.write(string.join([str(x) for x in nCountObjects], '\t') + '\n')
    allChoiceCounts.append(participantChoiceCounts)
    participantBinaryChoiceSequences = [participantChoiceSequences[x] for x in binaryA]
    for s in participantBinaryChoiceSequences:
        d = map(operator.sub, s[1:], s[:-1])
        Ai = [i!=0 for i in d].count(True)    # Number of transition
        Ni = [i==s[0] for i in s].count(True) # count number of times first symbol in s appears in s
        pvec = [0, 2.0/6, 2.0/15, 2.0/20]
        nReps = 6
        if Ni > nReps/2:
            j = nReps - Ni
        else:
            j = Ni
        nvec[j] += 1
        if Ai == 1:
            cvec[j] += 1
    # Fill in information on multiplicative inequality violations
    #MIFile.write('%i' % index)
    tripleFile.write('%i' % index)
    for x in range(5):
        xSet = pow(2, x)
        for y in range(x+1, 5):
            ySet = pow(2, y)
            for z in range(y+1, 5):
                zSet = pow(2, z)
                A = xSet + ySet + zSet
                xy = xSet + ySet
                yz = ySet + zSet
                xz = xSet + zSet
                NA = [participantChoiceCounts[A][i] for i in [x, y, z]]
                nxy = participantChoiceCounts[xy][x]
                nyz = participantChoiceCounts[yz][y]
                nxz = participantChoiceCounts[xz][x]
                pixyz = (nyz - NA[1])
                pixzy = ((6-nyz) - NA[2])
                piyxz = (nxz - NA[0])
                piyzx = ((6-nxz) - NA[2])
                pizxy = (nxy - NA[0])
                pizyx = ((6-nxy) - NA[1])
                xDiff = 6*participantChoiceCounts[A][x] - participantChoiceCounts[xy][x] * participantChoiceCounts[xz][x]
                yDiff = 6*participantChoiceCounts[A][y] - participantChoiceCounts[xy][y] * participantChoiceCounts[yz][y]
                zDiff = 6*participantChoiceCounts[A][z] - participantChoiceCounts[xz][z] * participantChoiceCounts[yz][z]
                #MIFile.write('\t%i\t%i\t%i' % (xDiff, yDiff, zDiff))
                MIFile.write('%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n' %
                    (index+1, A, NA[0], NA[1], NA[2], nxy, nyz, nxz,
                    6 * NA[0] - nxy * nxz,
                    6 * NA[1] - (6-nxy) * nyz,
                    6 * NA[2] - (6-nxz) * (6-nyz),
                    pixyz, pixzy, piyxz, piyzx, pizxy, pizyx))
                tripleFile.write('\t%i\t%i\t%i\t%i\t%i\t%i'
                    % (participantChoiceCounts[A][x], participantChoiceCounts[A][y], participantChoiceCounts[A][z],
                    participantChoiceCounts[xy][x], participantChoiceCounts[yz][y], participantChoiceCounts[xz][x]))
    #MIFile.write('\n')
    tripleFile.write('\n')

    # Vectors of length 3. First, second and third elements correspond to patterns 1-5, 5-1; 2-4, 4-2, 3-3
    # nvec[i] gives the number of binary choice sets where the participant's choices follow pattern i
    # cvec[i] gives the number of times, out of nvec[i] where the participant's choices have a single breakpoint
    # pvec[i] gives the probability, under the iid hypothesis, that there will be a single breakpoint
    n, c, p = nvec[1:], cvec[1:], pvec[1:]
    
    # Construct iterator for {0,1,...,n_0} x {0,1,...,n_1} x {0,1,...,n_2}
    margins = [range(i+1) for i in n]
    table = product(*margins)

    # List of lists giving log binomial probability mass function for three values of (n[i],p[i])
    log_pmf_list = [binom.logpmf(range(n[i]+1), n[i], p[i]) for i in range(3)]

    # Compute probabilities for each element of iterator (all possible values of c vector)
    ptable = [sum([log_pmf_list[i][v[i]] for i in range(3)]) for v in table]
    # Compute probability for particular observed c vector
    ppart = sum([log_pmf_list[i][c[i]] for i in range(3)])

    # Compute entropy H under iid, the negative expected value of ppart under iid
    H = -sum([exp(el) * el for el in ptable])
    # Compute p values for < ppart and <= ppart
    p_hi = sum([exp(ptable[i]) for i in range(len(ptable)) if ptable[i]-10e-9 <= ppart])
    p_lo = sum([exp(ptable[i]) for i in range(len(ptable)) if ptable[i]+10e-9 <= ppart])
    testFile.write('%e\t%e\t%e\t%e\n' % (H, ppart, p_lo, p_hi))

# C boilerplate
header_string = '''#include "RCM.h"
#include "RCM_multi_lottery.h"

'''
extern_C_function_string = '''
void set_multi_lottery_data(Universe *univ, RCM *rcm, int participant)
{
	unsigned i, A;
    for (A=0; A<univ->nSets; A++)
        if (univ->cardinality[A] > 1) {
            for (i=0; i<univ->cardinality[A]; i++) {
                unsigned x = univ->element[A][i];
                rcm->N[A][x] = multi_lottery[participant][A][x];
            }
            rcm->N_total[A] = 6;
        }
}
'''

# LaTeX boilerplate

LaTeXTablePreamble = '\\begin{table}[H]\n\t\\centering\n\t\\begin{tabular}{ccccc}\n\t\tA & B & C & D & E \\\\\n\t\t\\hline\n'
LaTeXTablePostamble = '\t\t\\hline\n\t\\end{tabular}\\caption{Data of subject %i}\\label{t:%i}\n\end{table}\n\n'

cFile = open('C_Code/RCM_multi_lottery.c', 'w')
cFile.write(header_string)
cFile.write('static int multi_lottery[141][32][5] = \n')
LaTeXFile = open('Data_Tables/RCM_multi_data.tex', 'w')
allChoiceList = []
for i, participantChoiceCounts in enumerate(allChoiceCounts):
    LaTeXFile.write(LaTeXTablePreamble)
    participantChoiceList = []
    for j, setChoiceCounts in enumerate(participantChoiceCounts):
        objectChoiceList = [str(x) for x in setChoiceCounts]
        setChoiceString = '{' + string.join(objectChoiceList, ', ') + '}'
        participantChoiceList.append(setChoiceString)
        LaTeXObjectChoiceList = objectChoiceList
        doPrint = True
        for choice in range(5):
            if not (j & pow(2, choice)): # Object j not in choice set
                LaTeXObjectChoiceList[choice] = '-'
            if ((j & pow(2, choice) == j) or j==0): # Choice set is singleton or empty
                doPrint = False
        if doPrint:
            LaTeXFile.write('\t\t' + string.join(LaTeXObjectChoiceList, ' & ') + ' \\\\\n')
    participantChoiceString = '{ /* participant ' + str(i+1) + ' */\n    ' + string.join(participantChoiceList, ',\n    ') + ' }'
    allChoiceList.append(participantChoiceString)
    LaTeXFile.write(LaTeXTablePostamble % (i+1, i+1))
allString = '{ ' + string.join(allChoiceList, ',\n  ') + ' };\n'
cFile.write(allString)
cFile.write(extern_C_function_string)
