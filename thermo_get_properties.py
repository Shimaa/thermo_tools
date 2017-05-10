from __future__ import print_function

from collections import defaultdict
import logging
import os.path
import numpy as np
import re
import itertools
import sys
import csv
import nasa_fit 

########
#
# This moduele used some functions from cantera/ck2cti.py to read Chemkin therm data
# Modified by Shimaa_Gamil / KAUST-CCRC to modify thermo parameters, then write a new thermo file
# 
########


class Species(object):
    def __init__(self, label):
        self.label = label
        self.thermo = None
        self.transport = None
        self.note = None
        self.composition = None

    def __str__(self):
        return self.label

    def to_cti(self, indent=0):
        lines = []
        atoms = ' '.join('{0}:{1}'.format(*a)
                         for a in self.composition.items())

        prefix = ' '*(indent+8)

        lines.append('species(name={0!r},'.format(self.label))
        lines.append(prefix + 'atoms={0!r},'.format(atoms))
        if self.thermo:
            lines.append(prefix +
                         'thermo={0},'.format(self.thermo.to_cti(15+indent)))
        if self.transport:
            lines.append(prefix +
                         'transport={0},'.format(self.transport.to_cti(14+indent)))
        if self.note:
            lines.append(prefix + 'note={0!r},'.format(self.note))

        lines[-1] = lines[-1][:-1] + ')'
        lines.append('')

        return '\n'.join(lines)


class ThermoModel(object):
    """
    A base class for thermodynamics models, containing several attributes
    common to all models:

    =============== =================== ========================================
    Attribute       Type                Description
    =============== =================== ========================================
    `Tmin`          ``float``           The minimum temperature at which the model is valid, or ``None`` if unknown or undefined
    `Tmax`          ``float``           The maximum temperature at which the model is valid, or ``None`` if unknown or undefined
    `comment`       ``str``             Information about the model (e.g. its source)
    =============== =================== ========================================

    """

    def __init__(self, Tmin=None, Tmax=None, comment=''):
        if Tmin is not None:
            self.Tmin = Tmin
        else:
            self.Tmin = None
        if Tmax is not None:
            self.Tmax = Tmax
        else:
            self.Tmax = None
        self.comment = comment


class NASA(ThermoModel):
    """
    A single NASA polynomial for thermodynamic data. The `coeffs` attribute
    stores the seven or nine polynomial coefficients
    :math:`\\mathbf{a} = \\left[a_{-2}\\ a_{-1}\\ a_0\\ a_1\\ a_2\\ a_3\\ a_4\\ a_5\\ a_6 \\right]`
    from which the relevant thermodynamic parameters are evaluated via the
    expressions

    .. math:: \\frac{C_\\mathrm{p}(T)}{R} = a_{-2} T^{-2} + a_{-1} T^{-1} + a_0 + a_1 T + a_2 T^2 + a_3 T^3 + a_4 T^4

    .. math:: \\frac{H(T)}{RT} = - a_{-2} T^{-2} + a_{-1} T^{-1} \\ln T + a_0 + \\frac{1}{2} a_1 T + \\frac{1}{3} a_2 T^2 + \\frac{1}{4} a_3 T^3 + \\frac{1}{5} a_4 T^4 + \\frac{a_5}{T}

    .. math:: \\frac{S(T)}{R} = -\\frac{1}{2} a_{-2} T^{-2} - a_{-1} T^{-1} + a_0 \\ln T + a_1 T + \\frac{1}{2} a_2 T^2 + \\frac{1}{3} a_3 T^3 + \\frac{1}{4} a_4 T^4 + a_6

    For the 7 coefficient form, the first two coefficients are taken to be zero.
    """

    def __init__(self, coeffs, **kwargs):
        ThermoModel.__init__(self, **kwargs)
        if len(coeffs) not in (7,9):
            raise InputParseError('Invalid number of NASA polynomial coefficients; '
                                  'should be 7 or 9.')
        self.coeffs = coeffs

    def to_cti(self, indent=0):
        prefix = ' '*indent
        vals = ['{0: 15.8E}'.format(i) for i in self.coeffs]
        if len(self.coeffs) == 7:
            lines = ['NASA([{0:.2f}, {1:.2f}],'.format(self.Tmin[0], self.Tmax[0]),
                     prefix+'     [{0}, {1}, {2},'.format(*vals[0:3]),
                     prefix+'      {0}, {1}, {2},'.format(*vals[3:6]),
                     prefix+'      {0}]),'.format(vals[6])]
        else:
            lines = ['NASA9([{0:.2f}, {1:.2f}],'.format(self.Tmin[0], self.Tmax[0]),
                     prefix+'      [{0}, {1}, {2},'.format(*vals[0:3]),
                     prefix+'       {0}, {1}, {2},'.format(*vals[3:6]),
                     prefix+'       {0}, {1}, {2}]),'.format(*vals[6:9])]

        return '\n'.join(lines)
    
class MultiNASA(ThermoModel):
    """
    A set of thermodynamic parameters given by NASA polynomials. This class
    stores a list of :class:`NASA` objects in the `polynomials`
    attribute. When evaluating a thermodynamic quantity, a polynomial that
    contains the desired temperature within its valid range will be used.
    """

    def __init__(self, polynomials=None, **kwargs):
        ThermoModel.__init__(self, **kwargs)
        self.polynomials = polynomials or []

    def to_cti(self, indent=0):
        prefix = ' '*indent
        lines = []
        for i,p in enumerate(self.polynomials):
            if i == 0:
                lines.append('({0}'.format(p.to_cti(indent+1)))
            elif i != len(self.polynomials)-1:
                lines.append(prefix + ' {0}'.format(p.to_cti(indent+1)))
            else:
                lines.append(prefix + ' {0})'.format(p.to_cti(indent+1)[:-1]))

        return '\n'.join(lines)
    
def warn(message):
        logging.warning(message)

def parseComposition(elements, nElements, width):
    """
    Parse the elemental composition from a 7 or 9 coefficient NASA polynomial
    entry.
    """
    composition = {}
    for i in range(nElements):
        symbol = elements[width*i:width*i+2].strip()
        count = elements[width*i+2:width*i+width].strip()
        if not symbol:
            continue
        try:
            count = int(float(count))
            if count:
                composition[symbol.capitalize()] = count
        except ValueError:
            pass
    return composition    

def readThermoEntry(lines, TintDefault):
        """
        Read a thermodynamics entry for one species in a Chemkin-format file
        (consisting of two 7-coefficient NASA polynomials). Returns the label of
        the species, the thermodynamics model as a :class:`MultiNASA` object, the
        elemental composition of the species, and the comment/note associated with
        the thermo entry.
        """
        identifier = lines[0][0:24].split()
        species = identifier[0].strip()

        if len(identifier) > 1:
            note = ''.join(identifier[1:]).strip()
        else:
            note = ''

        # Extract the NASA polynomial coefficients
        # Remember that the high-T polynomial comes first!
        try:
            Tmin = fortFloat(lines[0][45:55])
            Tmax = fortFloat(lines[0][55:65])
            
            try:
                Tint = fortFloat(lines[0][65:75])
                NRotors = lines[0][78]
            except ValueError:
                Tint = TintDefault
                NRotors = 0

            coeffs_high = [fortFloat(lines[i][j:k])
                           for i,j,k in [(1,0,15), (1,15,30), (1,30,45), (1,45,60),
                                         (1,60,75), (2,0,15), (2,15,30)]]
            coeffs_low = [fortFloat(lines[i][j:k])
                           for i,j,k in [(2,30,45), (2,45,60), (2,60,75), (3,0,15),
                                         (3,15,30), (3,30,45), (3,45,60)]]

        except (IndexError, ValueError) as err:
            raise Exception('Error while reading thermo entry for species {0}:\n{1}'.format(species, err))

        composition = parseComposition(lines[0][24:44], 4, 5)

        

        # Non-standard extended elemental composition data may be located beyond
        # column 80 on the first line of the thermo entry
        if len(lines[0]) > 80:
            elements = lines[0][80:]
            composition2 = parseComposition(elements, len(elements)//10, 10)
            composition.update(composition2)

        # Construct and return the thermodynamics model
        thermo = MultiNASA(
            polynomials=[
                NASA(Tmin=(Tmin,"K"), Tmax=(Tint,"K"), coeffs=coeffs_low),
                NASA(Tmin=(Tint,"K"), Tmax=(Tmax,"K"), coeffs=coeffs_high)
            ],
            Tmin=(Tmin,"K"),
            Tmax=(Tmax,"K"),
        )
       

        return species, thermo, composition, note, Tint, NRotors
    
def readline(ck_file, line_number):

    line = ck_file.readline()
   # print (line_number, line)
    if '!' in line:
        return line.split('!', 1)
    elif line:
        return line, ''
    else:
        return None, None
    
def contains(seq, value):
    if isinstance(seq, str):
        return value.lower() in seq.lower()
    else:
        return get_index(seq, value) is not None

def fortFloat(s):
    """
    Convert a string representation of a floating point value to a float,
    allowing for some of the peculiarities of allowable Fortran representations.
    """
    s = s.strip()
    s = s.replace('D', 'E').replace('d', 'e')
    s = s.replace('E ', 'E+').replace('e ', 'e+')
    return float(s)


def getThermoChemistryFromNASA (lowCoeffs, highCoeffs, Tint):

    T = 298.15
    R = 1.987 # Gas Constant in CAL

    # a values
    a1 = lowCoeffs[0]
    a2 = lowCoeffs[1]
    a3 = lowCoeffs[2]
    a4 = lowCoeffs[3]
    a5 = lowCoeffs[4]
    a6 = lowCoeffs[5]
    a7 = lowCoeffs[6]
    
    # Calculate Hf
    Hf = ( a1 + ( (a2/2) * T ) + ( (a3/3) * (T**2) ) + ( (a4/4) * (T**3) ) + ( (a5/5) * (T**4) ) + (a6/T) ) * (R * T)
    Hf = Hf/1000 # to be in KCAL
    
    # Calculate S
    S = ( (a1 * np.log(T)) + (a2 * T) + ( (a3/2) * (T**2) ) +  ( (a4/3) * (T**3) ) + ( (a5/4) * (T**4) + a7) ) * R 
    
    lowCoeffs = lowCoeffs[0:5] # from a1 to a5
    highCoeffs = highCoeffs[0:5] # from a1 to a5
    
    lowCoeffs.reverse() # reverse as poynomial coeffs are in reverse order comparing to a coeffs
    highCoeffs.reverse()
    p_low = np.poly1d(lowCoeffs)
    p_high = np.poly1d(highCoeffs)    
    
    cp300 = p_low(300) * R
    cp400 = p_low(400) * R
    cp500 = p_low(500) * R
    cp600 = p_low(600) * R
    cp800 = p_low(800) * R
    cp1000 = p_low(1000) * R

    cp1500 = p_high (1500) * R

    return cp300, cp400, cp500, cp600, cp800, cp1000, cp1500, Hf, S
	

def readThermoFile (path, output_file):	
        
	#path = sys.argv[1] # 'C:\\Shimaa_Work\\Thermo_UQ\\PECS_therm.therm'
	#output_file = sys.argv[2] # 'C:\\Shimaa_Work\\Thermo_UQ\\PECS_thermo_props.csv'
	#print (path)
	line_number = 0
	species_count = 0
	
	try:

		with open(output_file, 'w', newline='') as f_thermo:
			thermo_writer = csv.writer(f_thermo)
			thermo_writer.writerow (['Species', 'Hf', 'S', 'Cp300', 'Cp400', 'Cp500', 'Cp600', 'Cp800', 'Cp1000','Cp1500', 'NRotors', 'Atom#'])
			
			with open(path, 'rU') as ck_file:
				line_number += 1
				line, comment = readline(ck_file, line_number)
				
				while line is not None:
					#print ('here')
					line_number += 1
					line, comment = readline(ck_file, line_number)

					tokens = line.split() or ['']

					#print (tokens)

					
					# List of thermodynamics (hopefully one per species!)
					if tokens[0].upper().startswith('THER'):
						line_number += 1            
						line, comment = readline(ck_file, line_number)
						
						if line is not None and not contains(line, 'END'):
							TintDefault = float(line.split()[1])
						
						#print (TintDefault)

						thermo = []   
						while line is not None and not contains(line, 'END'):
							line_number += 1
							line, comment = readline(ck_file, line_number)

							# Grudging support for implicit end of section
							if contains(line, 'REAC') or contains(line, 'TRAN'):
								warn('"THERMO" section implicitly ended by start of '
										  'next section on line {0}.'.format(line_number))
								advance = False
								tokens.pop()
								break

							if len(line) >= 80 and line[79] in ['1', '2', '3', '4']:
								thermo.append(line)

								if line[79] == '4':

									species_count += 1
									
									label, thermo, comp, note, Tint, NRotors  = readThermoEntry(thermo, TintDefault)
									nasa_low = thermo.polynomials[0]
									nasa_high = thermo.polynomials[1]

									'''
									print ('\nSpecies # ', species_count)
									print (label)
									print ('Tint =', Tint)

									
									for coeff in nasa_low.coeffs:
										print ( float(coeff))

									for coeff in nasa_high.coeffs:
										print ( float(coeff))

								
									print  (nasa_low.coeffs)
									print  (nasa_high.coeffs)
									print (comp)
									print (note)
									print ("rotors", NRotors)
									'''									
									
									
									cp300, cp400, cp500, cp600, cp800, cp1000, cp1500, Hf, S = getThermoChemistryFromNASA (nasa_low.coeffs, nasa_high.coeffs, Tint)

								  #  print (cp300, cp400, cp500, cp600, cp800, cp1000, cp1500)
								  #  print ('Hf = ', Hf, 'S = ', S)

									
									
									# calculate atom number
									O = 'O'
									H = 'H'
									C = 'C'
									try:
										atom_num = 0
										if O in comp.keys():
											atom_num += comp[O]

										if H in comp.keys():
											atom_num += comp[H]

										if C in comp.keys():
											atom_num += comp[C]
											
										#print (label, atom_num)
									except:
										atom_num = 0
										print ('cannot calculate atom number')
										
									thermo_writer.writerow([label, Hf, S, cp300, cp400, cp500, cp600, cp800, cp1000, cp1500, NRotors, atom_num])

									'''
									try:
										species = speciesDict[label]
										# use the first set of thermo data found
										if species.thermo is not None:
											warn('Found additional thermo entry for species {0}'.format(label))
										else:
											species.thermo = thermo
											species.composition = comp
											species.note = note
									except KeyError:
										logging.info('Skipping unexpected species "{0}" while reading thermodynamics entry.'.format(label))
									'''
									
									thermo = []
												 
					  
	except:
		return


# a function to get thermo roperties of the passed species from a list of species properties
def get_species_thermo_properties (species_name, species_props_list):
	for species_prop in species_props_list:
		if species_prop[0] == species_name:							
			return species_prop


def reWriteThermoFile (butanol_species_list, species_props_list, source_thermo_file, output_thermo_file):	
        
	line_number = 0
	species_count = 0
	
	try:

		with open(output_thermo_file, 'wt') as out_thermo:
			
			with open(source_thermo_file, 'rU') as in_thermo:
				line_number += 1
				line, comment = readline(in_thermo, line_number)

				out_thermo.write (line)
				
				while line is not None:
					#print ('here')
					line_number += 1
					line, comment = readline(in_thermo, line_number)

					tokens = line.split() or ['']
					#print (tokens)
    
					out_thermo.write (line)
					
					
					# List of thermodynamics (hopefully one per species!)
					if tokens[0].upper().startswith('THER'):
						line_number += 1            
						line, comment = readline(in_thermo, line_number)

						out_thermo.write (line)
						
						if line is not None and not contains(line, 'END'):
							TintDefault = float(line.split()[1])
						
						#print (TintDefault)

						thermo = []   
						while line is not None and not contains(line, 'END'):

							#print (line)

							
							line_number += 1
							line, comment = readline(in_thermo, line_number)

							

							# Grudging support for implicit end of section
							if contains(line, 'REAC') or contains(line, 'TRAN'):
								warn('"THERMO" section implicitly ended by start of '
										  'next section on line {0}.'.format(line_number))
								advance = False
								tokens.pop()
								break
							

							if len(line) >= 80 and line[79] in ['1', '2', '3', '4']:
								thermo.append(line)

								if line[79] == '4':

									species_count += 1
									
									label, thermoData, comp, note, Tint, NRotors  = readThermoEntry(thermo, TintDefault)
									nasa_low = thermoData.polynomials[0]
									nasa_high = thermoData.polynomials[1]

									#print (butanol_species_list)


                                                                        # check if species is in target species list, then calculate thermo data using updated properties
                                                                        # you should replce the butanol species list with your species list
									if label in butanol_species_list:

										species_thermo_properties = get_species_thermo_properties (label, species_props_list)
										#print (species_thermo_properties)

										#call nasa fit										
										fit_line_2, fit_line_3, fit_line_4 =  nasa_fit.doNASAFit (species_thermo_properties)
		
										out_thermo.write (thermo [0]) # take 1st line as it is
										out_thermo.write (fit_line_2+'\n')
										out_thermo.write (fit_line_3+'\n')
										out_thermo.write (fit_line_4+'\n')
										
										
									else:  # if not, copy same thermo to the output file										
										out_thermo.write (thermo[0])
										out_thermo.write (thermo[1])
										out_thermo.write (thermo[2])
										out_thermo.write (thermo[3])
									
									thermo = []
												 
							else:
								out_thermo.write (line)
							
                                                              
	except:
		return

