import sys
import numpy as np
import matplotlib.pyplot as plt
import mysql.connector as db
import h5py
import csv
from pyparsing import Word, Optional, OneOrMore, Group, ParseException

###################################
# @author Shimaa_Gamil
# KAUST-CCRC
###################################

def parseChemicalFormula(formula):
    # Paul McGuire
    # Modified by @shimaa_gamil to combine similar atoms together
    
    caps = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    lowers = caps.lower()
    digits = "0123456789"

    element = Word( caps, lowers )
    integer = Word( digits )
    elementRef = Group( element + Optional( integer, default="1" ) )
    chemicalFormula = OneOrMore( elementRef )

    # try parsing the formula
    formulaData = chemicalFormula.parseString(formula)

    # print the results before combine similar atoms
    #print (formula, "->", formulaData)

    # combine similar atoms
    for i in range(len(formulaData)):
        if (i >= len(formulaData)):
            break
        for j in range(i+1, len(formulaData)):
            if (j >= len(formulaData)):
                break
            if (formulaData[i][0] == formulaData[j][0]): # same atom
                formulaData[i][1] = int(formulaData[i][1]) + int (formulaData[j][1])
                del formulaData[j]

    return formulaData

'''
filename = "C:\\Shimaa_Work\\Therm\\Test_Data\\1M2BBEN.LST"
file = open(filename, "r")
lines = file.readlines()
dataLine = lines[2]
dataArray = dataLine.split()


speciesName = dataArray[0]
Hf = float(dataArray[1])
S = float(dataArray[2])
Cp300 = float(dataArray[3])
Cp400 = float(dataArray[4])
Cp500 = float(dataArray[5])
Cp600 = float(dataArray[6])
Cp800 = float(dataArray[7])
Cp1000 = float(dataArray[8])
Cp1500 = float(dataArray[9])
'''

speciesName = sys.argv[1]
speciesFormula = sys.argv[2]
Hf = float(sys.argv[3])
S = float(sys.argv[4])
Cp300 = float(sys.argv[5])
Cp400 = float(sys.argv[6])
Cp500 = float(sys.argv[7])
Cp600 = float(sys.argv[8])
Cp800 = float(sys.argv[9])
Cp1000 = float(sys.argv[10])
Cp1500 = float(sys.argv[11])

symmetry = float(sys.argv[12])
NRotors = float(sys.argv[13])
hdf_fileName = sys.argv[14]

#print (group_id, Hf , S, Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500, symmetry, NRotors)

#get the number of atoms
formulaData = parseChemicalFormula(speciesFormula)
atom_count = 0
for atom in formulaData:
         atom_count += int(atom[1])

# Gas Constant (R)
R = 1.987 # Gas Constant in CAL

# Calculate CPINF
N = atom_count
Cpinf = R * (3 * N - (2 + NRotors/2)) # non linear molecule
#print (N, Cp1000, Cpinf)


# get A and B
Y_300 = -1 * np.log (1- (Cp300/Cpinf)) # A + BT
Y_1000 = -1 * np.log (1- (Cp1000/Cpinf))

a = np.array([[1,300], [1,1000]])
b = np.array([Y_300, Y_1000])
ab = np.linalg.solve(a, b)
A = ab[0]
B = ab[1]
#print (Y_300, Y_1000)
##print (np.allclose(np.dot(a, ab), b)) # to check if linear solve is correct

# estimate cp1500 if doesn't exist
if (Cp1500 < 0):
    Y_1500 = A + B * 1500
    Cp1500 = (1 - np.exp(-1*Y_1500)) * Cpinf
    #print ('calculated CP 1500 is', Cp1500)


# calculate cp at high temp i.e. 2000, 2500, 3000, 3500, 4000, 4500, 5000

Y_2000 = A + B * 2000
Cp2000 = (1 - np.exp(-1*Y_2000)) * Cpinf

Y_2500 = A + B * 2500
Cp2500 = (1 - np.exp(-1*Y_2500)) * Cpinf

Y_3000 = A + B * 3000
Cp3000 = (1 - np.exp(-1*Y_3000)) * Cpinf

Y_3500 = A + B * 3500
Cp3500 = (1 - np.exp(-1*Y_3500)) * Cpinf


Y_4000 = A + B * 4000
Cp4000 = (1 - np.exp(-1*Y_4000)) * Cpinf

Y_4500 = A + B * 4500
Cp4500 = (1 - np.exp(-1*Y_4500)) * Cpinf

Y_5000 = A + B * 5000
Cp5000 = (1 - np.exp(-1*Y_5000)) * Cpinf


# Prepare the polynomial regression
lowTpoints = np.array([(int(300), Cp300/R), (int(400), Cp400/R), (int(500), Cp500/R), (int(600), Cp600/R),(int(800), Cp800/R), (int(1000), Cp1000/R), (int(1500), Cp1500/R)]) #  removed:(int(3000), Cp3000/R)
highTpoints = np.array([(int(1000), Cp1000/R), (int(1500), Cp1500/R), (int(2000), Cp2000/R), (int(2500), Cp2500/R), (int(3000), Cp3000/R), (int(3500), Cp3500/R),(int(4000), Cp4000/R), (int(4500), Cp4500/R),(int(5000), Cp5000/R)])


# get x and y vectors
lx = lowTpoints[:,0]
ly = lowTpoints[:,1]

hx = highTpoints[:,0]
hy = highTpoints[:,1]

degree = 4

# low T coeffecients
lcoeffs = np.polyfit(lx, ly, degree) # Polynomial coefficients, highest power first.
#print ("lcoeffs", lcoeffs)
#print ("lx:", lx)
#print ("ly:", ly)


#p = np.poly1d(lcoeffs)
#test_900 = p(900) * R
#print (p(300)*R , p(400)*R , p(500)*R , p(600)*R ,p(800)*R ,p(1000)*R)       
'''       
# fit values, and mean
yhat = p(lx)                      # or [p(z) for z in x]
ybar = np.sum(ly)/len(ly)          # or sum(y)/len(y)
ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
sstot = np.sum((ly - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
lR2 = ssreg / sstot
'''
#print ("R^2 for Low Temp =", lR2)




# high T coeffecients
hcoeffs = np.polyfit(hx, hy, degree) # Polynomial coefficients, highest power first.
#print ("hcoeffs", hcoeffs)

#ph = np.poly1d(hcoeffs)
#test_1200 = ph(1200)*R
#print (ph(1000)*R , ph(1500)*R )      
'''       
# fit values, and mean
hyhat = ph(hx)                      # or [p(z) for z in x]
hybar = np.sum(hy)/len(hy)          # or sum(y)/len(y)
hssreg = np.sum((hyhat-hybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
hsstot = np.sum((hy - hybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
hR2 = hssreg / hsstot
'''
#print ("R^2 for High Temp =", hR2)

### Find the mid-point
lp = np.poly1d(lcoeffs)
hp = np.poly1d(hcoeffs)
min_diff = lp(1000)
mid_point = 1000
for i in range(1000, 1500, 1):
     l = lp(i)*R
     h = hp(i)*R
     if ( np.abs(l - h)  < min_diff ):
         min_diff = np.abs(l - h)
         mid_point = i

#print (mid_point)
#print ( lp (mid_point) )
#print ( hp (mid_point) )

#print (mid_point, min_diff)

## do refitting using the mid-point
lowTpoints_refit = np.array([(int(300), Cp300/R), (int(400), Cp400/R), (int(500), Cp500/R), (int(600), Cp600/R),(int(800), Cp800/R), (int(1000), Cp1000/R), (int(mid_point), lp(mid_point)), (int(1500), Cp1500/R)])
highTpoints_refit = np.array([(int(mid_point), lp(mid_point)), (int(1500), Cp1500/R), (int(2000), Cp2000/R), (int(3000), Cp3000/R),(int(5000), Cp5000/R)]) #(int(2500), Cp2500/R) , (int(3000), Cp3000/R), (int(3500), Cp3500/R),(int(4000), Cp4000/R), (int(4500), Cp4500/R)


# get x and y vectors
lx_refit = lowTpoints_refit[:,0]
ly_refit = lowTpoints_refit[:,1]

hx_refit = highTpoints_refit[:,0]
hy_refit = highTpoints_refit[:,1]

degree = 4

# low T coeffecients
lcoeffs = np.polyfit(lx_refit, ly_refit, degree) # Polynomial coefficients, highest power first.

# high T coeffecients
hcoeffs = np.polyfit(hx_refit, hy_refit, degree) # Polynomial coefficients, highest power first.

lp = np.poly1d(lcoeffs)
hp = np.poly1d(hcoeffs)
#print ( lp (mid_point) )
#print ( hp (mid_point) )

output_file = 'D:\\Python34\\nasa_log.csv'
with open(output_file, 'at') as out_thermo:
    thermo_writer = csv.writer(out_thermo)
    thermo_writer.writerow ([speciesName, mid_point, lp (mid_point)*R, hp (mid_point)*R])
			

## get a1 to a7 of low coeffs
a5 = lcoeffs[0]
a4 = lcoeffs[1]
a3 = lcoeffs[2]
a2 = lcoeffs[3]
a1 = lcoeffs[4]

T = 298.15

rightSide_a6 = a1 + ( (a2/2) * T ) + ( (a3/3) * (T**2) ) + ( (a4/4) * (T**3) ) + ( (a5/5) * (T**4))
leftSide_a6 = (Hf*1000)/(R * T) # multiply hf by 1000 to be in CAL instead of KCAL

a6 = ( leftSide_a6 - rightSide_a6 ) * T

rightSide_a7 = (a1 * np.log(T)) + (a2 * T) + ( (a3/2) * (T**2) ) + ( (a4/3) * (T**3) ) + ( (a5/4) * (T**4) )
leftSide_a7 = S / R

a7 = leftSide_a7 - rightSide_a7

NASA_POLY_FIT_Low = ['%e' % a1, '%e' % a2, '%e' % a3, '%e' % a4, '%e' % a5, '%e' % a6, '%e' % a7]

## calculate Hf and S at mid point
T = mid_point
Hf_mid = ( a1 + ( (a2/2) * T ) + ( (a3/3) * (T**2) ) + ( (a4/4) * (T**3) ) + ( (a5/5) * (T**4) ) + (a6/T) ) * (R * T)
Hf_mid = Hf_mid/1000 # to be in KCAL

S_mid = ( (a1 * np.log(T)) + (a2 * T) + ( (a3/2) * (T**2) ) +  ( (a4/3) * (T**3) ) + ( (a5/4) * (T**4) + a7) ) * R 
#print (Hf_mid)

## high coeffs
a5_h = hcoeffs[0]
a4_h = hcoeffs[1]
a3_h = hcoeffs[2]
a2_h = hcoeffs[3]
a1_h = hcoeffs[4]

rightSide_a6_h = a1_h + ( (a2_h/2) * T ) + ( (a3_h/3) * (T**2) ) + ( (a4_h/4) * (T**3) ) + ( (a5_h/5) * (T**4) )
leftSide_a6_h = (Hf_mid*1000)/(R * T) # multiply hf by 1000 to be in CAL instead of KCAL

a6_h = ( leftSide_a6_h - rightSide_a6_h ) * T

rightSide_a7_h = ( a1_h * np.log(T) ) + ( a2_h * T ) + ( (a3_h/2) * (T**2) ) +  ( (a4_h/3) * (T**3) ) + ( (a5_h/4) * (T**4))
leftSide_a7_h = S_mid / R

a7_h = leftSide_a7_h - rightSide_a7_h

NASA_POLY_FIT_High = ['%e' % a1_h, '%e' % a2_h, '%e' % a3_h, '%e' % a4_h, '%e' % a5_h, '%e' % a6_h, '%e' % a7_h]


#print (NASA_POLY_FIT_High[0], NASA_POLY_FIT_High[1], NASA_POLY_FIT_High[2], NASA_POLY_FIT_High[3], NASA_POLY_FIT_High[4], NASA_POLY_FIT_High[5], NASA_POLY_FIT_High[6])
#print (NASA_POLY_FIT_Low[0], NASA_POLY_FIT_Low[1], NASA_POLY_FIT_Low[2], NASA_POLY_FIT_Low[3], NASA_POLY_FIT_Low[4], NASA_POLY_FIT_Low[5], NASA_POLY_FIT_Low[6])

lp_hf = np.poly1d(np.append(lcoeffs, [a6]))
hp_hf = np.poly1d(np.append(hcoeffs, [a6_h]))

lp_s = np.poly1d(np.append(lcoeffs, [a7]))
hp_s = np.poly1d(np.append(hcoeffs, [a7_h]))


x_l = []
y_l = []
x_h = []
y_h = []
y_l_hf = []
y_h_hf = []
y_l_s = []
y_h_s = []

for j in range(300, mid_point, 1):
    x_l.append(j)
    y_l.append(lp(j)*R)
    # Hf
    Hf = ( a1 + ( (a2/2) * j ) + ( (a3/3) * (j**2) ) + ( (a4/4) * (j**3) ) + ( (a5/5) * (j**4) ) + (a6/j) ) * (R * j)
    Hf = Hf/1000 # to be in KCAL
    # S
    S = ( (a1 * np.log(j)) + (a2 * j) + ( (a3/2) * (j**2) ) +  ( (a4/3) * (j**3) ) + ( (a5/4) * (j**4) + a7) ) * R
    y_l_hf.append (Hf)
    y_l_s.append (S)

for j in range(mid_point, 5000, 1):
    x_h.append(j)
    y_h.append(hp(j)*R)
    # Hf
    Hf = ( a1_h + ( (a2_h/2) * j ) + ( (a3_h/3) * (j**2) ) + ( (a4_h/4) * (j**3) ) + ( (a5_h/5) * (j**4) ) + (a6_h/j) ) * (R * j)
    Hf = Hf/1000 # to be in KCAL
    # S
    S = ( (a1_h * np.log(j)) + (a2_h * j) + ( (a3_h/2) * (j**2) ) +  ( (a4_h/3) * (j**3) ) + ( (a5_h/4) * (j**4) + a7_h) ) * R
    y_h_hf.append (Hf)
    y_h_s.append (S)

#print (y_h)   

'''
plt.plot(x_l, y_l, 'r-')
plt.plot( x_h, y_h, 'g-')
plt.plot(x_l, y_l_hf, 'b-')
plt.plot( x_h, y_h_hf, 'y-')
plt.plot(x_l, y_l_s, 'r-')
plt.plot( x_h, y_h_s, 'y-')
plt.show()
'''

#########################################
### Write Polynomial Result to HDF5 file
#########################################

if (len(hdf_fileName) > 3): # check if there is  a passing file name in args
    hf = h5py.File(hdf_fileName, 'w')


    ##write species name to hdf5
    sn_set = hf.create_dataset("thermoData/speciesName", (1,), dtype = 'S32')
    sn_set[0] = speciesName.encode("ascii", "ignore")

    ##write species formula to hdf5
    sf_set = hf.create_dataset("thermoData/speciesFormula", (1,), dtype = 'S32')
    sf_set[0] = speciesFormula.encode("ascii", "ignore")

    ##write units to hdf5
    uset = hf.create_dataset("thermoData/units", (1,), dtype = 'S32')
    uset[0] = 'KCAL'.encode("ascii", "ignore")


    ##write thermo parameters to hdf5
    #thset = hf.create_dataset("thermoData/thermoParameters", data = thermoParam, dtype = np.floating)

    ##write reference temperature to hdf5
    rset = hf.create_dataset("thermoData/referenceTemperature", (1,), dtype = np.floating)
    rset[0] = 298.15



    ##write thermo polynomials to hdf5
    low_coeffs = [300, mid_point, a1, a2, a3, a4, a5, a6, a7]
    high_coeffs = [mid_point, 5000, a1_h, a2_h, a3_h, a4_h, a5_h, a6_h, a7_h]

    thermoCoeffs = []
    thermoCoeffs.append(list(high_coeffs))
    thermoCoeffs.append(list(low_coeffs))

    coset = hf.create_dataset("thermoData/thermoPolynomials/thermoCoeff", data=thermoCoeffs, dtype = np.floating)


    ##write thermo polynomials type to hdf5
    tset = hf.create_dataset("thermoData/thermoPolynomials/thermoType", (1,), dtype='S32')
    tset[0] = np.string_("NASA7")

    hf.close()

#format species NASA fit data into chemkin format
chem_string_line_1 = ''
chem_string_line_2 = ''
chem_string_line_3 = ''
chem_string_line_4 = ''

# Line 1
chem_string_line_1 += '{0:<16}CFTHERMO'.format(speciesName)

#move the next line to the beginning so that we need the number of atoms to use rotors
#formulaData = parseChemicalFormula(speciesFormula)

if len(formulaData) <= 4:
    # Use the original Chemkin syntax for the element counts
    for atom in formulaData:
         chem_string_line_1 += '{0!s:<2}{1:<3d}'.format(atom[0], int(atom[1]))
    chem_string_line_1 += '     ' * (4 - len(formulaData))
else:
    chem_string_line_1 += '     ' * 4

chem_string_line_1 += 'g{0:<10.3f}{1:<10.3f}{2:<8.2f}     {3:.0f}1'.format(300, 5000, mid_point, NRotors)

if len(formulaData) > 4:
    chem_string_line_1 += '&\n'
    # Use the new-style Chemkin syntax for the element counts
    # This will only be recognized by Chemkin 4 or later
    for atom in formulaData:
        chem_string_line_1 += '{0!s:<2}{1:<3d}'.format(atom[0], atom[1])

# Line 2
chem_string_line_2 += '{0:< 15.8E}{1:< 15.8E}{2:< 15.8E}{3:< 15.8E}{4:< 15.8E}    2'.format(a1_h, a2_h, a3_h, a4_h, a5_h)

# Line 3
chem_string_line_3 += '{0:< 15.8E}{1:< 15.8E}{2:< 15.8E}{3:< 15.8E}{4:< 15.8E}    3'.format(a6_h, a7_h, a1, a2, a3)

# Line 4
chem_string_line_4 += '{0:< 15.8E}{1:< 15.8E}{2:< 15.8E}{3:< 15.8E}                   4'.format(a4, a5, a6, a7)

print (chem_string_line_1)
print (chem_string_line_2)
print (chem_string_line_3)
print (chem_string_line_4)

    
sys.exit()


