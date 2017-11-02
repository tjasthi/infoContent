The program infoContent opens a files, whose filenames iare specified by the user as a command-line argument. The first file is an multiple sequence alignment (MSA) in Stockholm format.
Reads the MSA file and finds relationships between entropy and mutual information between columns. Uses columns of the different sequences to calculate probability.

The file with write:
The "bottom 10" columns with the lowest entropy
The "top 50" column pairs with the highest mutual information

Calculating entropy, joint entropy, and mutual information: https://people.cs.umass.edu/~elm/Teaching/Docs/mutInf.pdf

5 Tasks:
1. Read in an MSA in Stockholm format
2. Calculate p_i(x)pi(x) and H(i)H(i) for every column i
3. Calculate p_{i,j}(x, y)pi,j(x,y) for every pair of columns ii and jj such that i < j
4. Calculate the mutual information I(i,j)I(i,j) for every pair of columns ii and jj such that i < j
5. Print the following information to the screen/stdout

Ideally we can limit the amount of iterations to increase the runtime of the function when comparing probability, entropy, and mi.
