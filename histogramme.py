# -*- coding: latin-1 -*-
# rene.natowicz@esiee.fr - version 2, 23 mai 2022
""" 
Le fichier des distances relatives (DR) dont vous voulez afficher l'histogramme contient 
une distance relative par ligne et rien d'autre.
Exemple : 
0.494520822486
0.0939254807197
0.365846578605
0.150577932504
...

Dans ce qui suit, ce fichier est GAINS_nruns=5000_Linf=5_Lsup=15_Cinf=5_Csup=15.csv.csv (csv = comma separated values)
mais le "DR" peut être remplacé par ce que vous voulez. 
Exemple : fichier gainsRelatifs_nruns=5000_Linf=5_Lsup=15_Cinf=5_Csup=15.csv


Ouvrir un terminal Unix ou Linux et se placer dans le répertoire qui contient
le fichier DR.CSV

Taper la ligne de commande     
   python3 histogramme.py DR   <<< attention SVP : "DR" et non pas "DR.csv"
ou 
   python histogramme.py DR
(respectivement Python version 3 et Python version 2)

Vous obtiendrez un fichier image 
   DR.png
contenant l'histogramme. 


Exemple d'exécution dans un terminal Unix : 
% python histogramme.py DR
l'histogramme est dans le fichier DR.png
% 
"""

import sys
import csv
import matplotlib.pyplot as plt

def histogramme(fileName) : 
	GR = [] # gains relatifs 
	with open(fileName+".csv") as csvfile:
		reader = csv.reader(csvfile)
		for row in reader:
			gr = row[0]
			GR.append(float(gr))
	lg = len(GR)
	print(str(lg)+" valeurs dans le fichier des gains relatifs")
	lz = len([z for z in GR if z == 0.0])
	print(str(lz)+" valeurs 0.0 dans le fichier des gains relatifs")
	# Enlever les 0 car ils écrasent l'histogramme
	GR = [g for g in GR if g != 0.0]
	lg = len(GR)
	print(str(lg)+" valeurs non nulles dans le fichier des gains relatifs")
	print("valeur non nulle minimum : " + str(min(GR)))
	print("valeur non nulle maximum : " + str(max(GR)))
	print("creation de l'histogramme des valeurs non nulles de gain relatif")
	h = plt.hist(GR,bins=int(lg/10)) #  nombre de bins à votre convenance...
	plt.savefig(fileName+".png")
	plt.close()
	csvfile.close

def main() : 
	if len(sys.argv) != 2 :
		print("Usage : python3 histogramme.py fileName")
		print("Exemple : python3 hystogramme.py DR")
		return
	fileName = sys.argv[1]
	histogramme(fileName)
	print("l'histogramme est dans le fichier " + fileName + ".png")

main()