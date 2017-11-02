#I worked with Alex Orimoloye

import sys
import math

#Reading the file
try:
	read_file = sys.argv[1]
except:
	raise Exception("Did not specify a correct filename argument on the command line")
try:
	read_me = open(read_file, 'r')
	write_me = sys.stdout
except:
	raise Exception("File does not exist or is not readable")

"""
Dictionaries in use:
column_dic
	key: column #
	value: string of column characters
probability_dic
	key: column #
	value: dictionary of bases as keys and their corresponding probabilities
entropy_dic
	key: column #
	value: entropy of the column
dualProb_dic
	key: tuple of 2 columns (i,j) where i < j
	value: dictionary of 2 base tuples as keys and their corresponding probabilities
MI_dic
	key: tuple of 2 columns (i,j) where i < j
	value: mutual information of the 2 columns
"""
column_dic = {}
probability_dic = {}
entropy_dic = {}
dualProb_dic = {}
MI_dic = {}

def main():
	readSequenceFile()
	calcProb()
	calcEntropy()
	calcDualProb()
	calcMI()
	writeBottom()
	writeTop()

	#Exit and save
	write_me.close()

"""
Reads the sequence file, storing the abundance and entropy of the 4 bases by column in a dictionary.
"""
def readSequenceFile():
	for line in read_me:
		if (line[0] == '#'):
			continue
		elif (line[0] == '/'):
			break
		try:
			words = line.split()
			sequence = words[1]
			for i in range(len(sequence)):
				if i in column_dic:
					column_dic[i] = column_dic[i] + sequence[i]
				else:
					column_dic[i] = sequence[i]
		except:
			raise Exception("Incorrectly formatted sequence file")
	
"""
Store the probabilities in a dictionary from the column dictionary.
"""
def calcProb():
	size = len(column_dic[0])
	for key in column_dic:
		stringSeq = column_dic[key]
		probability_dic[key] = {'A':0, 'U':0, 'C':0, 'G':0, '.':0}
		for i in range(len(stringSeq)):
			letter = stringSeq[i]
			(probability_dic[key])[letter] = (probability_dic[key])[letter] + 1
		(probability_dic[key])['A'] = float((probability_dic[key])['A']) / size
		(probability_dic[key])['U'] = float((probability_dic[key])['U']) / size
		(probability_dic[key])['C'] = float((probability_dic[key])['C']) / size
		(probability_dic[key])['G'] = float((probability_dic[key])['G']) / size
		(probability_dic[key])['.'] = float((probability_dic[key])['.']) / size
	read_me.close()

"""
Calculates entropy from the probability dictionary.
"""
def calcEntropy():
	for key in probability_dic:
		total = 0
		for letter in probability_dic[key]:
			prob = (probability_dic[key])[letter]
			if (prob != 0):
				total = total + (prob * math.log(prob, 2))
		entropy_dic[key] = -total

"""
Calculates probability pairs from the probability dictionary.
"""
def calcDualProb():
	for first in probability_dic:
		for second in probability_dic:
			if first < second:
				dualProb_dic[(first,second)] = { \
				'AA':0, 'AC':0, 'AG':0, 'AU':0, 'A.':0, \
				'CA':0, 'CC':0, 'CG':0, 'CU':0, 'C.':0, \
				'GA':0, 'GC':0, 'GG':0, 'GU':0, 'G.':0, \
				'UA':0, 'UC':0, 'UG':0, 'UU':0, 'U.':0, \
				'.A':0, '.C':0, '.G':0, '.U':0, '..':0 \
				}
				size = len(column_dic[first])

				"""
				There are 2 ways we tried to find dual probability.
				The first way is correct, but takes a long time. You can run with the print statement to see the probability be calculated.
				The second way is incorrect but runs faster. Runs under the assumption that bases are independent.
				"""

				### Method 1: Finds the dual probability correctly but takes a long ass time
				for i in range(size):
					for j in range(size):
						letter1 = (column_dic[first])[i]
						letter2 = (column_dic[second])[i]
						(dualProb_dic[(first,second)])[letter1 + letter2] = (dualProb_dic[(first,second)])[letter1 + letter2] + 1
				for x in dualProb_dic[(first,second)]:
					(dualProb_dic[(first,second)])[x] = float((dualProb_dic[(first,second)])[x]) / (size * size)
				# print(dualProb_dic[(first,second)])
				
				### Method 2: Finds the dual_probability assuming bases are independent
				# for x in dualProb_dic[(first,second)]:
				# 	letter1 = x[0]
				# 	letter2 = x[1]
				# 	jointProb = (probability_dic[first])[letter1] * (probability_dic[second])[letter2]
				# 	(dualProb_dic[(first,second)])[x] = jointProb

"""
Calculates MI from the dual probability and probability dictionary.
"""
def calcMI():
	for first in entropy_dic:
		for second in entropy_dic:
			if first < second:
				pair_dic = dualProb_dic[(first,second)]
				total = 0
				for x in pair_dic:
					prob = pair_dic[x]
					if (prob != 0):
						total = total + (prob * math.log(prob, 2))
				mi = entropy_dic[first] + entropy_dic[second] - total
				MI_dic[(first,second)] = mi

"""
Writes the bottom 10 columns with the lowest entropy
"""
def writeBottom():
	ten_dic = {}
	largest_key = 0
	largest_val = 0
	for x in entropy_dic:
		if len(ten_dic) < 10:
			ten_dic[x] = entropy_dic[x]
			if len(ten_dic) == 10:
				for y in ten_dic:
					if ten_dic[y] > largest_val:
						largest_key = y
						largest_val = ten_dic[y]
		else:
			if entropy_dic[x] < largest_val:
				del ten_dic[largest_key]
				largest_val = 0
				ten_dic[x] = entropy_dic[x]
				for y in ten_dic:
					if ten_dic[y] > largest_val:
						largest_key = y
						largest_val = ten_dic[y]
	for x in ten_dic:
		write_me.write(str(x) + '\n')

"""
Writes the top 50 column pairs with the highest MI
"""
def writeTop():
	fifty_dic = {}
	smallest_keyPair = (0,0)
	smallest_val = float("inf")
	for x in MI_dic:
		if len(fifty_dic) < 50:
			fifty_dic[x] = MI_dic[x]
			if len(fifty_dic) == 50:
				for y in fifty_dic:
					if fifty_dic[y] < smallest_val:
						smallest_keyPair = y
						smallest_val = fifty_dic[y]
		else:
			if MI_dic[x] > smallest_val:
				del fifty_dic[smallest_keyPair]
				smallest_val = float("inf")
				fifty_dic[x] = MI_dic[x]
				for y in fifty_dic:
					if fifty_dic[y] < smallest_val:
						smallest_keyPair = y
						smallest_keyPair = fifty_dic[y]
	for x in fifty_dic:
		write_me.write(str(x[0]) + ',' + str(x[1]) + '\n')

#Run main
if __name__== "__main__":
	main()