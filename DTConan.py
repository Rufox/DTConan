#from chemspipy import ChemSpider
import pubchempy as pcp
#from rdkit import Chem
import os
import sys
import shutil
import glob
import numpy as np
from sympy.utilities.iterables import multiset_permutations
import configparser

# conda activate my-rdkit-env
######################################################
#cs = ChemSpider('cp0mOsKP18B1nQQmplYVUHwol5jS6rzI')

#c = cs.get_compound(236)


#print (c.csid)
#print(c.molecular_formula)
#print(c.nominal_mass)
#print(c.smiles)
#print(c.common_name)
#print(c.mol_3d)
#print(c.image)
#for result in cs.search_by_formula('C6H6'):#cs.search('Benzene'):
 #   print(result)

######################################################
def Var_init():
	global number, Eq_Base, list_of_compounds, seeDB
	global transmutation, permutation
	global gl_scaling, Big_variable
	Big_variable = {}
	gl_scaling = 1.2
	number = 0
	seeDB = 0
def is_number(d,n):
	is_number = True
	try:
		num = float(n)
		# check for "nan" floats
		is_number = num == num   # or use `math.isnan(num)`
		if 'Pcent' in d and (num<0 or num>1):
			print("Valor de ",d, " fuera de los limites [0,1]\n")
			exit(1)
	except ValueError:
		is_number = False
		print ("Informacion en archivo de variables incorrecta\n",d ,"=", n)
		exit(1)
	return is_number
def establecerVariablesDefault():
	global transmutation, permutation,Eq_Base,gl_scaling,list_of_compounds#enlace,altura,shape,Eq_Global,totalVertices,atomos_dados
	global seeDB
	#print(Big_variable)
	# Variables Opcionales
	if "base_formula" in Big_variable.keys():
		Eq_Base = Big_variable["base_formula"]
	if "base_coumpond" in Big_variable.keys():
		list_of_compounds = Big_variable["base_coumpond"].split(' ')
	if "transmutation" in Big_variable.keys():
		transmutation = Big_variable["transmutation"].split(' ')
	if "permutation" in Big_variable.keys():
		permutation = Big_variable["permutation"].split(' ')
	## OPTIONAL NOW
	if "distance_scaling" in Big_variable.keys():
		is_number("distance_scaling",Big_variable["distance_scaling"])
		gl_scaling = float(Big_variable["distance_scaling"])
	if "statedb" in Big_variable.keys():
		is_number("statedb",Big_variable["statedb"])
		seeDB = int(Big_variable["statedb"])
# Funcion simple lectura de archivo
def leerArchivoParametros(configs):
	config = configparser.ConfigParser()
	config.read(configs)

	for secciones in config.sections():
		for (variable, valor) in config.items(secciones):
			#agregarVariableGlobal(variable, valor)
			Big_variable[variable]=valor
	establecerVariablesDefault()

#########################################
#########################################
def changeAtoms(order,pila_atoms):#old_at, old_num, new_at, new_num, pila_atoms):
	i = 0
	todo = order.split('-')
	for x in range(len(pila_atoms)):
		#print(x)
		if(pila_atoms[x] == todo[0]):
			pila_atoms[x] = todo[2]
		#	print(todo[2], pila_atoms[x])
			i+=1
		if (i == int(todo[1])):
			break;
	#print("MESSAGE:",i, "of", todo[1], todo[0], "where transmuted to", todo[2])
	print("MESSAGE", str(i)+todo[0],"->",str(i)+todo[2],"("+todo[1], "asked)")
	return pila_atoms
def escribirArchivoXYZ(name, title,atomicSymbols ,atoms_not_play,coordsList,position_play,position_not_play):
	input = open(name+".xyz","a")
	input.write(str(len(coordsList))+"\n")
	input.write(title+"\n")
	#print (coordsList)
	#print(shape)
	n_atom = len(coordsList)
	posicion=0
	for yes in range(len(atomicSymbols)):
		input.write(atomicSymbols[yes]+"\t")
		input.write(' '.join(map(str, coordsList[position_play[yes]]))+"\n")
	
	for no in range(len(atoms_not_play)):#range(totalVertices):
		#for depto in range(shape[puntos]):
		input.write(atoms_not_play[no]+"\t")#+"\t"+(coordsList[puntos]+"\n"))
		input.write(' '.join(map(str, coordsList[position_not_play[no]]))+"\n")
		#	posicion+=1			

def escribirInputGaussian(name, title,atomicSymbols ,atoms_not_play,coordsList,position_play,position_not_play):
	#coordsList = np.array([[row[0],row[1], row[2], row[3]] for row in origincoords])
	global number
	input = open(title+".com","w+")
	
	#print (var.Big_variable)
	input.write("%NProc="+Big_variable["core"]+"\n")
	input.write("%mem="+Big_variable["memory"]+"GB\n")
	input.write("#"+Big_variable["header"]+"\n\n")
	input.write("Automatic Input "+name+" "+str(number)+"\n\n")
	input.write(Big_variable["charge_multi"]+"\n")
	#for line in coordsList:
	#	input.write(' '.join(map(str, line))+"\n")
	posicion=0
	for yes in range(len(atomicSymbols)):
		if (atomicSymbols[yes] != "X"):
			input.write(atomicSymbols[yes]+"\t")
			input.write(' '.join(map(str, coordsList[position_play[yes]]))+"\n")
	
	for no in range(len(atoms_not_play)):#range(totalVertices):
		#for depto in range(shape[puntos]):
		if (atoms_not_play[no] != "X"):
			#print (atoms_not_play[no])
			input.write(atoms_not_play[no]+"\t")#+"\t"+(coordsList[puntos]+"\n"))
			input.write(' '.join(map(str, coordsList[position_not_play[no]]))+"\n")
	input.write("\n")       #Porque Gaussian es espcial
	input.close()
	#exit (1)
	number+=1

def GlobalPermutation(Eq_Global):
	print ("Printing all possible permutations")
	permutation =[]
	for p in multiset_permutations(Eq_Global):
		print(p)
		permutation.append(p)
		#break
	return permutation

def savePlayingAtoms(pila_atoms, players):
	#print (players)
	atoms_play =[]
	atoms_not_play =[]
	position_not_play =[]
	position_play = []
	for x in range(len(pila_atoms)):
		if(pila_atoms[x] not in players):
			#print(pila_atoms[x])
			position_not_play.append(x)
			atoms_not_play.append(pila_atoms[x])
		else:
			atoms_play.append(pila_atoms[x])
			position_play.append(x)
	#print(position_not_play)
	return atoms_not_play,position_not_play, atoms_play, position_play
# NOTA, escoger 2d o 3D a gusto.

def lookForAllDB():
	global number
	j= 0
	number = 0
	list_of_compounds = pcp.get_compounds(Eq_Base,'formula')#,record_type='3d')#,listkey_count=1)
	for result in list_of_compounds:#,listkey_count=1):#,record_type='3d'):
		#,listkey_count=1):
		print(j,result)
		cid = result.to_dict(properties=['cid'])['cid']
		atom = result.to_dict(properties=['atoms'])['atoms']
		n_atom = len(atom)
		atom_c = []
		symbol = []
		for i in range(0,n_atom):
			#print (atom[i]['element'], atom[i]['x'],atom[i]['y'],atom[i]['z'] if 'z' in atom[i] else "0.0" )
			atom_c.append([atom[i]['x']*gl_scaling,atom[i]['y']*gl_scaling,atom[i]['z']*gl_scaling if 'z' in atom[i] else "0.0"])
			symbol.append(atom[i]['element'])
		escribirArchivoXYZ("Original","Coumpound-"+str(cid)+"\tIdentificador: "+str(j),symbol,[],atom_c,list(range(0,n_atom)),[])
		# NOTA: no se hacer esto bien, pero por el momento trabajar con array normal, luego transformar a numpy.array 
		# Cambiamos 1 H por 1 X
		for tr in transmutation:
			symbol=changeAtoms(tr, symbol)
			print(symbol)
		print("FINAL ATOMS",symbol)
		#exit(1)
		new_atom_symbols = np.asarray(symbol)
		atoms_not_play, position_not_play, atoms_play ,position_play= savePlayingAtoms(new_atom_symbols,permutation)
		### aca hacer las premutaciones SOLO los atomos que juegan
		#exit(0)
		Permu = GlobalPermutation(atoms_play)
		#IMPRESION
		#print (poly)
		#exit(0)
		for resultado in range(len(Permu)):
			escribirArchivoXYZ("Recopilacion","Pemu-"+str(j)+"-"+str(resultado),Permu[resultado],atoms_not_play,atom_c,position_play,position_not_play)
			#print ("Recopilacion","Pemu-"+str(j),Permu[resultado])#,poly)77
			#escribirInputGaussian
			#print("Permut",resultado,Permu[resultado])#,poly)
			escribirInputGaussian("Recopilacion","Pemu-"+str(j)+"-"+str(resultado),Permu[resultado],atoms_not_play,atom_c,position_play,position_not_play)
		j+=1
	print("Estructuras finales: ",j)

def lookForSpecificDB(cid,j):
	global number
	result = pcp.Compound.from_cid(cid)
	#,listkey_count=1):
	cid = result.to_dict(properties=['cid'])['cid']
	atom = result.to_dict(properties=['atoms'])['atoms']
	n_atom = len(atom)
	atom_c = []
	symbol = []
	for i in range(0,n_atom):
		#print (atom[i]['element'], atom[i]['x'],atom[i]['y'],atom[i]['z'] if 'z' in atom[i] else "0.0" )
		atom_c.append([atom[i]['x']*gl_scaling,atom[i]['y']*gl_scaling,atom[i]['z']*gl_scaling if 'z' in atom[i] else "0.0"])
		symbol.append(atom[i]['element'])
	escribirArchivoXYZ("Original","Coumpound-"+str(cid)+"\tIdentificador: "+str(j),symbol,[],atom_c,list(range(0,n_atom)),[])
	# NOTA: no se hacer esto bien, pero por el momento trabajar con array normal, luego transformar a numpy.array 
	# Cambiamos 1 H por 1 X
	for tr in transmutation:
		symbol=changeAtoms(tr, symbol)
		print(symbol)
	print("FINAL ATOMS",symbol)
	#exit(1)
	new_atom_symbols = np.asarray(symbol)
	atoms_not_play, position_not_play, atoms_play ,position_play= savePlayingAtoms(new_atom_symbols,permutation)
	### aca hacer las premutaciones SOLO los atomos que juegan
	Permu = GlobalPermutation(atoms_play)
	for resultado in range(len(Permu)):
		escribirArchivoXYZ("Recopilacion","Pemu-"+str(j)+"-"+str(resultado),Permu[resultado],atoms_not_play,atom_c,position_play,position_not_play)
		#print ("Recopilacion","Pemu-"+str(j),Permu[resultado])#,poly)77
		#escribirInputGaussian
		#print("Permut",resultado,Permu[resultado])#,poly)
		escribirInputGaussian("Recopilacion","Pemu-"+str(j)+"-"+str(resultado),Permu[resultado],atoms_not_play,atom_c,position_play,position_not_play)
	j+=1
#exit(0)

#exit(1)
if os.path.exists("Original.xyz"):
    os.remove("Original.xyz")
if os.path.exists("Recopilacion.xyz"):
    os.remove("Recopilacion.xyz")
if os.path.exists("Inputs"):
	shutil.rmtree("Inputs")
# Iniciamos las variables principales
Var_init()
# Leemos varbiables de usuario.
leerArchivoParametros(sys.argv[1])

print("Valor de cosas varias")
print (transmutation, permutation,Eq_Base,gl_scaling)

#exit(0)
if (seeDB == 1):
	print("Printing DB information for",Eq_Base)
	j = 0
	list_of_compounds=[]
	list_of_compounds = pcp.get_compounds(Eq_Base,'formula')#,listkey_count=1)
	print(list_of_compounds)
	for result in list_of_compounds:
		cid = result.to_dict(properties=['cid'])['cid']
		atom = result.to_dict(properties=['atoms'])['atoms']
		n_atom = len(atom)
		atom_c = []
		symbol = []
		for i in range(0,n_atom):
			atom_c.append([atom[i]['x'],atom[i]['y'],atom[i]['z'] if 'z' in atom[i] else "0.0"])
			symbol.append(atom[i]['element'])
		escribirArchivoXYZ("Original","Coumpound-"+str(cid)+"\tIdentificador: "+str(j),symbol,[],atom_c,list(range(0,n_atom)),[])
		j+=1
	print("Original DB coordinates shown in file Original.xyz")
	exit(0)
else:
	print("SIN DB")
	try:
		j = 0
		for cid in list_of_compounds:
			lookForSpecificDB(cid,j)
			j+=1
	except NameError:
		print("Doing a full DB search")
		lookForAllDB()

# MOVER TODO .COM A CARPETA INPUTS
if not os.path.exists("Inputs"):
    os.makedirs("Inputs")

files = glob.glob("*.com")
for file in files:
	shutil.move(file,"Inputs")

print("Original DB coordinates shown in file Original.xyz")

print("Full permutations coordinates in Recopilacion.xyz")

print("Gaussian inputs in Inputs folder")
exit(0)

#for result in pcp.get_compounds('CO37Li3','formula'):
#	print(result)

#BENCENO
m = Chem.MolFromSmiles('C1=CC=CC=C1')#('Cc1ccccc1')
print(Chem.MolToMolBlock(m))
# Busqeuda de formula




atom_c = [['2.809', '0.4594', '0.0'],
 ['2' ,'-0.1284' ,'0.0'],
 ['3.618' ,'-0.1284', '0.0'],
 ['2.309', '-1.0794' ,'0.0'],
 ['3.309' ,'-1.0794', '0.0'],
 ['2.809' ,'1.0794', '0.0'],
 ['1.4103' ,'0.0632', '0.0'],
 ['4.2077', '0.0632', '0.0'],
 ['1.9446', '-1.581' ,'0.0'],
 ['3.6734', '-1.581', '0.0']]
atom_c= [[  '1.143000000000'     '-2.004800000000'      '0.340000000000'],
[    '-0.019400000000'     '-2.331100000000'     '-0.329300000000'],
[     '1.999900000000'    '-1.141400000000'     '-0.327500000000'],
[     '0.003500000000'     '-1.386600000000'     '-1.362700000000'],
[     '1.205200000000'     '-0.681400000000'     '-1.376400000000'],
[     '2.317900000000'      '0.015400000000'      '0.344900000000'],
[     '2.004000000000'      '1.175600000000'     '-0.335800000000'],
[     '1.205400000000'      '0.714700000000'     '-1.398100000000'],
[     '1.397400000000'     '-0.011000000000'      '1.396000000000'],
[     '0.691800000000'     '-1.213800000000'      '1.398200000000'],
[    '-0.701500000000'     '-1.237200000000'      '1.383700000000'],
[     '0.712900000000'      '1.208500000000'      '1.378700000000'],
[    '-0.681600000000'      '1.183700000000'      '1.383500000000'],
[    '-1.385900000000'     '-0.018800000000'      '1.371600000000'],
[     '1.157700000000'      '2.029200000000'      '0.341800000000'],
[    '-2.015200000000'     '-1.184200000000'     '-0.332100000000'],
[    '-1.189000000000'     '-2.044600000000'      '0.350900000000'],
[    '-2.312200000000'     '-0.001000000000'      '0.334400000000'],
[    '-1.992900000000'      '1.162800000000'     '-0.335500000000'],
[    '-1.129800000000'      '2.004800000000'      '0.341500000000'],
[     '0.010800000000'      '2.327100000000'     '-0.372000000000'],
[    '-1.210300000000'     '-0.688400000000'     '-1.368600000000'],
[    '-1.206300000000'      '0.708500000000'     '-1.401100000000'],
[    '-0.004300000000'      '1.416600000000'     '-1.426900000000' ]]

symbol = ['C' ,'C', 'C', 'C', 'C','C' ,'C', 'C', 'C', 'C','C' ,'C', 'C', 'C', 'C','C' ,'C', 'C', 'C', 'C','C' ,'C', 'C', 'C']
atom_c = np.asarray(atom_c)
#symbol = np.asarray(symbol)
#symbol = symbol.astype(str)
print(atom_c)
print(symbol)
print(atom_c[0])

# NOTA: no se hacer esto bien, pero por el momento trabajar con array normal, luego transformar a numpy.array 
new_atom_symbols=changeAtoms("C", 4, "Si", 3, symbol)
print(new_atom_symbols)

new_atom_symbols = np.asarray(new_atom_symbols)
atoms_not_play, position_not_play, atoms_play ,position_play= savePlayingAtoms(new_atom_symbols,["C","Si"])
### aca hacer las premutaciones SOLO los atomos que juegan
#exit(0)
Permu = GlobalPermutation(atoms_play)
#IMPRESION
#print (poly)
#exit(0)
for resultado in range(len(Permu)):
	escribirArchivoXYZ("Recopilacion","Pemu-"+str(resultado),Permu[resultado],atoms_not_play,atom_c,position_play,position_not_play)
	print ("Recopilacion","Pemu-"+str(resultado),Permu[resultado])#,poly)77
	#escribirInputGaussian
	print("Permut",resultado,Permu[resultado])#,poly)
