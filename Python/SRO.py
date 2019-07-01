import numpy as np

c=12		# FCC coordination number
f1=0.1		# Fraction of Fe atoms
f2=0.6		# Fraction of Ni atoms
f3=0.3		# Fraction of Cr atoms

# Read in the RDF file and skip all lines - 1 to get to very end of the file
f=open('rdf-1202.txt','r')
g=open('SRO_params.txt','w')
length=len(f.readlines())
f.seek(0)
for i in xrange(length-1):
	f.readline()
x=f.readline().split()

# Read in the effective coordination numbers of each atom interaction
c11=float(x[3])
c12=float(x[5])
c13=float(x[7])
c21=float(x[9])
c22=float(x[11])
c23=float(x[13])
c31=float(x[15])
c32=float(x[17])
c33=float(x[19])

# Calculate the Warren-Cowley Parameter for each interaction
a11=1-c11/(f1*c)
a12=1-c12/(f2*c)
a13=1-c13/(f3*c)
a21=1-c21/(f1*c)
a22=1-c22/(f2*c)
a23=1-c23/(f3*c)
a31=1-c31/(f1*c)
a32=1-c32/(f2*c)
a33=1-c33/(f3*c)

# Print out the coordination numbers and SRO values

g.write('Fe-Fe coord #: ' + str(c11) + '\n')
g.write('Fe-Ni coord #: ' + str(c12) + '\n')
g.write('Fe-Cr coord #: ' + str(c13) + '\n')
g.write('Total Fe coord #: ' + str(c11+c12+c13) + '\n\n')
g.write('Fe-Fe SRO param: ' + str(a11) + '\n')
g.write('Fe-Ni SRO param: ' + str(a12) + '\n')
g.write('Fe-Cr SRO param: ' + str(a13) + '\n\n')

g.write('Ni-Fe coord #: ' + str(c21) + '\n')
g.write('Ni-Ni coord #: ' + str(c22) + '\n')
g.write('Ni-Cr coord #: ' + str(c23) + '\n')
g.write('Total Ni coord #: ' + str(c21+c22+c23) + '\n\n')
g.write('Ni-Fe SRO param: ' + str(a21) + '\n')
g.write('Ni-Ni SRO param: ' + str(a22) + '\n')
g.write('Ni-Cr SRO param: ' + str(a23) + '\n\n')

g.write('Cr-Fe coord #: ' + str(c31) + '\n')
g.write('Cr-Ni coord #: ' + str(c32) + '\n')
g.write('Cr-Cr coord #: ' + str(c33) + '\n')
g.write('Total Cr coord #: ' + str(c31+c32+c33) + '\n\n')
g.write('Cr-Fe SRO param: ' + str(a31) + '\n')
g.write('Cr-Ni SRO param: ' + str(a32) + '\n')
g.write('Cr-Cr SRO param: ' + str(a33) + '\n\n')

#Total SRO parameter
g.write(str(a11+a12+a13+a21+a22+a23+a31+a32+a33))

f.close()
g.close()
