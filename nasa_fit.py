import sys
import numpy as np
import matplotlib.pyplot as plt

# @author Shimaa_Gamil
# KAUST-CCRC


def doNASAFit (props):

    #print (props)

    speciesName = props[0]
    Hf          = float(props[1])
    S           = float(props[2])
    Cp300       = float(props[3])
    Cp400       = float(props[4])
    Cp500       = float(props[5])
    Cp600       = float(props[6])
    Cp800       = float(props[7])
    Cp1000      = float(props[8])
    Cp1500      = float(props[9])
    NRotors     = int (props [10])
    atom_count  = int (props [11])
        
    #print (speciesName, Hf , S, Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500)

    # Gas Constant (R)
    R = 1.987 # Gas Constant in CAL

    # Calculate CPINF
    
    Cpinf = R * (3 * atom_count - (2 + NRotors/2)) # non linear molecule    

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
    lowTpoints = np.array([(int(300), Cp300/R), (int(400), Cp400/R), (int(500), Cp500/R), (int(600), Cp600/R),(int(800), Cp800/R), (int(1000), Cp1000/R), (int(1500), Cp1500/R), (int(3000), Cp3000/R)])
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

    a5_h = hcoeffs[0]
    a4_h = hcoeffs[1]
    a3_h = hcoeffs[2]
    a2_h = hcoeffs[3]
    a1_h = hcoeffs[4]

    rightSide_a6_h = a1_h + ( (a2_h/2) * T ) + ( (a3_h/3) * (T**2) ) + ( (a4_h/4) * (T**3) ) + ( (a5_h/5) * (T**4) )
    leftSide_a6_h = (Hf*1000)/(R * T) # multiply hf by 1000 to be in CAL instead of KCAL

    a6_h = ( leftSide_a6_h - rightSide_a6_h ) * T

    rightSide_a7_h = ( a1_h * np.log(T) ) + ( a2_h * T ) + ( (a3_h/2) * (T**2) ) +  ( (a4_h/3) * (T**3) ) + ( (a5_h/4) * (T**4))
    leftSide_a7_h = S / R

    a7_h = leftSide_a7_h - rightSide_a7_h

    NASA_POLY_FIT_High = ['%e' % a1_h, '%e' % a2_h, '%e' % a3_h, '%e' % a4_h, '%e' % a5_h, '%e' % a6_h, '%e' % a7_h]

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

    #print (mid_point, min_diff)


    #print (NASA_POLY_FIT_High[0], NASA_POLY_FIT_High[1], NASA_POLY_FIT_High[2], NASA_POLY_FIT_High[3], NASA_POLY_FIT_High[4], NASA_POLY_FIT_High[5], NASA_POLY_FIT_High[6])
    #print (NASA_POLY_FIT_Low[0], NASA_POLY_FIT_Low[1], NASA_POLY_FIT_Low[2], NASA_POLY_FIT_Low[3], NASA_POLY_FIT_Low[4], NASA_POLY_FIT_Low[5], NASA_POLY_FIT_Low[6])

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



    #format species NASA fit data into chemkin format; ONLY 3 lines, escape the first one
    chem_string_line_2 = ''
    chem_string_line_3 = ''
    chem_string_line_4 = ''


    # Line 2
    chem_string_line_2 += '{0:< 15.8E}{1:< 15.8E}{2:< 15.8E}{3:< 15.8E}{4:< 15.8E}    2'.format(a1_h, a2_h, a3_h, a4_h, a5_h)

    # Line 3
    chem_string_line_3 += '{0:< 15.8E}{1:< 15.8E}{2:< 15.8E}{3:< 15.8E}{4:< 15.8E}    3'.format(a6_h, a7_h, a1, a2, a3)

    # Line 4
    chem_string_line_4 += '{0:< 15.8E}{1:< 15.8E}{2:< 15.8E}{3:< 15.8E}                   4'.format(a4, a5, a6, a7)

    #print (chem_string_line_2)

    return chem_string_line_2, chem_string_line_3, chem_string_line_4


