from __future__ import division #floating point division automatic
import os
import gzip
import numpy as np
import cPickle as pickle
import copy
import sys
import subprocess32
import argparse
#from runnucmer import runnucmer, getnucsnps, which
from collections import Counter, OrderedDict
from scipy import sparse
from itertools import combinations
import pkg_resources
import logging
from shutil import rmtree


logger = logging.getLogger(__name__)

#getting package data
pkg_resources.resource_filename("quanttb", 'data/')


ref = pkg_resources.resource_filename("quanttb", 'data/GCF_000277735.2_ASM27773v2_genomic.fna')
with open(ref, "r") as reffile:
	refdata = reffile.read()
refdata = "".join(refdata.split()[6:])
refdata = "0" + refdata

# getting mutation dictionary for antibiotic resistance
mutable = {}
mutlist = pkg_resources.resource_filename("quanttb", 'data/snpmutlist.txt')
with open(mutlist, 'rb') as f:
	f.readline()
	for line in f:
		snp = line.rstrip().split(',')
		pos = int(snp[0])
		if pos in mutable:
			mutable[pos]['alt'].append(snp[2])
		else:
			mutable[pos] = {'drug': snp[3],'alt':[snp[2] ]}

ranges = pkg_resources.resource_filename("quanttb",'data/pegenes.tsv')
with open(ranges, 'r') as file:
			ranges = file.read().split('\n')

ranges = [x for x in ranges if len(x) != 0]
ranges = [y.split() for y in ranges]
highs = np.array([int(y[1]) for y in ranges])
lows = np.array([int(y[0]) for y in ranges])				

def in_range(x):
	return np.any((lows <= x) & (x <= highs))

def openFile(filename):
	if not os.path.isfile(filename):
		raise Exception(filename + " doesn't exist")
	if filename.endswith(".gz"):
		file = gzip.open(filename, "r")		
	elif filename.endswith(".vcf"):
		file = open(filename, "r")
	elif filename.endswith((".pkl", '.samp','.db')):
		with open(filename, 'rb') as f:
			file = pickle.load(f)
	elif filename.endswith("snps"):
		file = open(filename, "r")
	return file

class Snpset(object):
	'''Snpset Object.
	Makes a Snpset object from the given VCF file. A snpset object basically contains a dictionary of position allele values. 
	If Hardfilter is set to true, the resulting Snpset object does not have ambiguous sites or low coverages, only the position
	and passing allele is stored.
	Include refers to a list of sites to ensure should not be filtered from the vcf file. 
'''
	def __init__(self, fileName, hardfilter = True, include = None, data = None):

		self.fileName = fileName
		if data:
			self.sample = self.fileName
			self.getsnpsdict(data)
		elif self.fileName.endswith(".snps"):
			self.getsnpsnuc()
		elif self.fileName.endswith(('.vcf','.vcf.gz')):
			self.getsamplename()
			self.getsnps(hardfilter, include = include)
			self.hardfilter = hardfilter
		else:
			logging.error('In order to make snpsets, must input either a vcf file or a .snps file')
			sys.exit()
		self.masked =False
		self.findres()
	def getsamplename(self):
		name = os.path.basename(self.fileName)
		name = name.split('.')[0]

		self.sample = name
	def getsnps(self, hard = True, include = None):
		if not hard:
			self.refdict = dict(list(enumerate(refdata)))
		vcffile = openFile(self.fileName)
		nucs = []
		counts = {}
		positions = []
		for line in vcffile:
			if not line.startswith("#"):
				base = line.split()
				if  "Del" in base[6]:
					nucleo = "N"
				elif "LowCov" in base[6] and hard == True:
					nucleo = "N"
				elif base[4] == "." and hard == True: # skipping ref calls
					continue
				elif len(base[3]) != 1 or len(base[4]) != 1: #removing insertions
					nucleo = "N"
				elif int(base[5]) < 11: #removing low quality
					nucleo = "N"
				elif in_range(int(base[1])): #removing PPE sites
					nucleo = "N"
				elif "Amb" in base[6] and hard == True:
					nucleo = "N"				
				else:
					info = base[7].split(';')
					if len(info) > 4 : #removing reference calls based on a very high BC count
						BC = [int(i) for i in info[5].split('=')[1].split(",")]
						if sum(BC) != 0:
							thresh = 0.9 if hard else 0.99
							if any ([i/float(sum(BC)) > thresh for i in BC]):
								if base[4] == ".":
									if hard or include is None:
										continue
									elif int(base[1]) not in include:
										continue
									if int(base[1]) in include and not hard:
										nucleo = {base[3]:1}
										cov = {base[3]:max(BC)}
								else:
									nucleo = base[4] if hard else {base[4]:1} #switching to frequency
									cov = {base[4]:max(BC)}
							elif hard:
								nucleo = "N"
							else:
								cov = {nuc:count for nuc, count in zip(["A", "C", "G", "T"], BC) if count != 0 }
								nucleo = {nuc: (cov[nuc] / float(sum(cov.values()))) for nuc in cov.iterkeys()} #convert to frequency
						else:
							nucleo = "N"	
					else:
						nucleo = "N"
				positions.append(int(base[1]))
				nucs.append(nucleo)
				if not hard and nucleo != "N":
					counts[int(base[1])] = cov 
		vcffile.close()
		self.badloccount = nucs.count("N")
		self.badpos = [pos for pos,snp in zip(positions, nucs) if snp == "N"] #saving where the badlocations are
		self.snptable = {pos:snp for pos,snp in zip(positions,nucs) if snp != "N"}
		if not hard:
			self.covtable = counts
		self.hardfilter = hard
		self.orgcount = len(self.snptable)
		self.posincluded = include
		self.getsnpinfo()

	def getsnpsnuc(self):
		snpfile = openFile(self.fileName)
		nucs = []
		positions = []

		names = snpfile.readline().split()
		self.sample = os.path.splitext(os.path.basename(names[1]))[0]
		for i in [0,1,2]:
			snpfile.readline()

		for line in snpfile:
			line = line.split()
			pos = int(line[0])
			if in_range(pos):
				nucleo = "N"
			else:
				nucleo = line[2]
			nucs.append(nucleo)
			positions.append(pos)
		snpfile.close()

		self.badloccount = nucs.count("N")
		self.badpos = [pos for pos,snp in zip(positions, nucs) if snp == "N"] #saving where the badlocations are
		self.snptable = {pos:snp for pos,snp in zip(positions,nucs) if snp != "N"}
		self.hardfilter = True
		self.orgcount = len(self.snptable)
		self.getsnpinfo()

	def getsnpsdict(self, snptable):
		self.badpos = [x for x in snptable if in_range(x)]
		self.snptable = {pos:nuc for pos,nuc in snptable.items() if pos not in self.badpos}
		self.orgcount = len(self.snptable)
		if any(isinstance(nuc,dict) for nuc in self.snptable.values()):
			self.hardfilter = False
		else:
			self.hardfilter = True
		self.getsnpinfo()


	def getsnpinfo(self):
		self.snpcounts = len(self.snptable)
		#self.ambpos = [y for y in self.snptable if isinstance(self.snptable[y], dict)]
		if self.hardfilter == True:
			self.ambpos = []
		else:
			amb = []
			for snp in self.snptable:
				BC = self.snptable[snp].values()
				if all([i/float(sum(BC)) < 0.95 for i in BC]):
					amb.append(snp)
			self.ambpos = amb
			self.amcov = sum([sum(x.values()) for x in self.snptable.values()])/self.snpcounts if self.snpcounts != 0 else 0

	def getsnppos(self, amb = False):
		if amb:
			return  sorted(self.snptable.keys())
		else:
			return  [snp for snp in sorted(self.snptable.keys()) if snp not in self.ambpos]
	def getsnpseq(self, amb = False):
		if amb:
			return [self.snptable[snp] if not isinstance(self.snptable[snp], dict) else self.snptable[snp].keys() for snp in sorted(self.snptable.keys()) ]
		else:
			return [self.snptable[y] for y in sorted(self.snptable.keys()) if not isinstance(self.snptable[y], dict) ]
	def save(self, outputfile):
		with open(outputfile + ".samp", "wb") as f:
			pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)
	def posunion(self, set2, amb = False): #find snps in common
		unions = set(self.getsnppos(amb = amb) + set2.getsnppos(amb = amb))
		return list(unions - set(self.badpos))
	def intersectpos(self, set2, amb = False):
		return list(set(self.getsnppos(amb = True)) & set(set2.getsnppos(amb = True)))
	def unique(self, set2): #find snps unique to this sample
		#unique = list(set(self.getsnppos()) - set(set2.getsnppos()))
		unique = []
		for snp in self.snptable:
			if snp in set2.snptable:
				if self.snptable[snp] !=  set2.snptable[snp]:
					unique.append(snp)
				elif isinstance(self.snptable[snp], dict):
					if all(nuc in set2.snptable[snp] for nuc in self.snptable[snp].keys()):
						unique.append(snp)
			else:
				unique.append(snp)
		return unique

	def removesnps(self, snpdict, covremove = None, refremove = None):
		if covremove is None:
			covremove = 0.5
		if self.hardfilter:
			self.snptable = {snp:nuc for snp,nuc in self.snptable.items() if snp not in snpdict or nuc != snpdict[snp]}
		else:
			if refremove:
				self.refdict = {snp:nuc for snp,nuc in self.refdict.items() if snp not in snpdict}
			for pos in self.snptable.keys():
				if pos in snpdict:
					nuc = snpdict[pos]
					if nuc in self.snptable[pos]:
						if covremove is not None:
							self.snptable[pos][nuc] = max(self.snptable[pos][nuc] - covremove, 0)
							if self.snptable[pos][nuc] * sum(self.covtable[pos].values()) < self.covtable[pos][nuc]: #removing partial snps this is wrong needs to be fixed
							#if round(self.snptable[pos][nuc] * self.covtable[pos][nuc]) == 0:
								del self.snptable[pos][nuc]
							elif self.snptable[pos][nuc] < 0.05:
								del self.snptable[pos][nuc]
						else:
							del self.snptable[pos][nuc]
				elif refremove is not None: #removing influence of reference database, may help! Also removing bad positions from snpsample
					if pos not in self.refdict:
						continue
					if pos in refremove:
						self.snptable[pos] = {nuc: max(val - covremove, 0) for nuc,val in self.snptable[pos].items()} #removing bad positions from snpsample based on cov
						self.snptable[pos] = {nuc:val for nuc,val in self.snptable[pos].items() if val > 0.05}
						self.snptable[pos] = {nuc:val for nuc,val in self.snptable[pos].items() if val * sum(self.covtable[pos].values()) < self.covtable[pos][nuc]}
					else:
						rfnuc = self.refdict[pos]
						if rfnuc in self.snptable[pos]:
							self.snptable[pos][rfnuc] = max(self.snptable[pos][rfnuc] - covremove, 0)
							if self.snptable[pos][rfnuc] * sum(self.covtable[pos].values()) < self.covtable[pos][rfnuc]: #removing partial snps
							#if round(self.snptable[pos][rfnuc] * self.covtable[pos][rfnuc]) == 0:
								del self.snptable[pos][rfnuc]
							elif self.snptable[pos][rfnuc] < 0.05:
								del self.snptable[pos][rfnuc]
				if len(self.snptable[pos]) == 0:
					del self.snptable[pos]
		self.getsnpinfo()

	def buffer(self, locations):
		for pos in locations:
			if pos not in self.snptable and pos not in self.badpos:
				self.snptable[pos] = refdata[pos]
		self.getsnpinfo()

	def findres(self, mutable = mutable):
		resloc = list(set(self.snptable.keys()) & set(mutable.keys()))

		resistances = {}
		nucleotides = set(['A', 'C', 'G', 'T'])
		drugs = {}
		for pos in resloc:
			res = set(mutable[pos]['alt'])
			sus = nucleotides - res
			if self.hardfilter:
				samploc = self.snptable[pos]
			else:
				samploc = set(self.snptable[pos].keys())
			items = self.snptable[pos]

			if len(sus.intersection(samploc)) == len(samploc):
				continue #hom or het susc
			if len(res.intersection(samploc)) == len(samploc):
				if len(samploc) > 1 and max(items.values()) < 0.9 and not self.hardfilter:
					resistances[pos] = "resistant-het"
				else:
					resistances[pos] = "resistant-homo"
			elif len(sus.intersection(samploc)) > 0 and len(res.intersection(samploc)) > 0:
				top = max(items, key=items.get)
				if items[top] > 0.9:
					if top in res:
						resistances[pos] =  "resistant-homo"
				else:
					respcnt = sum([y for x,y in items.items() if x in res])
					suspcnt = 1-respcnt
					resistances[pos] = "hetres, res: " + str(round(respcnt, 2)) + '% sus: ' + str(round(suspcnt, 2)) + "%"
			else:
				logging.debug(pos, items, res)
		for x in resistances:
			d = mutable[x]['drug'] 
			if d in drugs:
				drugs[d][x] =  resistances[x]
			else:
				drugs[d] = {x:resistances[x]}

		self.resistances = drugs




def iteration(db, snpsampleset, iterationmax = 8, withref = True,  maxthresh = 0.15):
	''' Iterative classification method. Input is db as a SNPdb object. Can specify the max number of iterations, and the max threshold for iterations to stop.
	If withref is true, it will ensure the H37rv reference is included in both the sample and the db'''


	thedb = copy.deepcopy(db)
	if withref:
		if "H37rv" not in thedb.sampnames:
			thedb.makerefstrain()
		includeme = thedb.sample('H37rv').snptable.values()
	else:
		includeme = None
		if "H37rv" in thedb.sampnames:
			thedb.removeref("H37rv")

	logging.info("Comparing " + snpsampleset.sample + " against " + str(len(db.sampnames))  + " in database")
	thesample = copy.deepcopy(snpsampleset)
	if thesample.posincluded is None and withref:
		logging.info("The sample does not have reference bases included by default")

	strains = []
	maxvalue = 1
	iterations = 0
	iterationres = dict()

	if not hasattr(thedb, 'matrix'):
		logging.info("Making matrix.")
		thedb.makematrix()
		thedb.matrix = thedb.matrix.A
	else:
		thedb.matrix = thedb.matrix.A


	#getting sample matrix
	sampsnps = []
	sampfreq =  []
	sampcov = []
	for snp,d in thesample.snptable.items():
		item = [(snp,nuc) for nuc in d]
		sampfreq += d.values()
		sampcov += thesample.covtable[snp].values()
		sampsnps += item
	sampsnps = set(sampsnps)

	sampfreqs = []
	sampcovs = []
	badpos = set(thesample.badpos)
	for x,y in thedb.matcols:
		if x in thesample.snptable:
			if y in thesample.snptable[x]:
				sampfreqs.append(thesample.snptable[x][y])
				sampcovs.append(thesample.covtable[x][y])
			else:
				sampfreqs.append(0)
				sampcovs.append(0)
		else:
			sampfreqs.append(np.nan)
			sampcovs.append(np.nan)
	sampfreqs = np.array(sampfreqs)
	sampcovs = np.array(sampcovs)

	#subsetting db to only places that have in common with the sample
	knownlocs = np.where(~np.isnan(sampcovs))[0]
	thisit_cols = thedb.matcols[knownlocs]
	sampfreqs_known = sampfreqs[knownlocs]
	sampcovs_known = sampcovs[knownlocs]
	thisit_mat = thedb.matrix[:,knownlocs] #matrix subsetted by only the snps that are knwon in the sample
	thisit_loc = thedb.matloc[knownlocs]

	covlist = np.array([sum(thesample.covtable[x].values()) for x in thisit_loc])


	#removing genomes that don't have at least 15 snps in common with sample
	if 'H37rv' in thedb.matrows:
		h37 = thisit_loc[(thisit_mat[np.where(db.matrows == 'H37rv')[0],] == True)[0]]
		noth37 = np.where(~np.isin(thisit_loc,h37))[0]
		remaining =  np.count_nonzero(thisit_mat[:,noth37], axis = 1)
	else:
		remaining =  np.count_nonzero(thisit_mat, axis = 1)
	passing = np.where((remaining >= 15) + (thedb.matrows == 'H37rv'))[0]
	thisit_refs = thedb.matrows[passing]
	thisit_mat = thisit_mat[passing,]
	thedb.matrix = thedb.matrix[passing,]

	#getting the snp counts of only the genomes that will be evaluated
	snpcounts = thedb.matrix.sum(axis = 1) 

	#initializing the coverage threshold for a snp to be detected in sample to low value
	scov = sum(sampcov)/len(sampcov)
	covthresh = scov/15
	if scov > 25:
		maxthresh = 0.2

	#calculating current coverage of database in sample
	refdictsnps = set(thesample.refdict.items())
	tog=set(map(tuple, thisit_cols)) & sampsnps - refdictsnps
	ins=[x in tog for x in sampsnps]
	freqin = np.array(sampfreq)[ins]
	freqin = freqin[freqin != 0].mean()
	covin = np.array(sampcov)[ins]
	covin = covin[covin != 0].mean()

	#saving results
	iterationres['0'] = [thisit_cols, thisit_refs , thisit_mat, sampfreqs_known]	


	#conducting iterations. While the number of snps in database is greater than 15, the threshold has not been reached, 
	#iterationmax has not been reached, and the number of passing snps in sample is greater than 0.
	while(len(np.unique(thisit_loc)) > 15 and maxvalue > maxthresh and iterations < iterationmax and len(passing) > 0):


		#getting statsitics from the sample
		freqmat = thisit_mat * sampfreqs_known #frequency of SNPs for every strain in database
		samepos = freqmat.sum(axis = 1)
		fullpos = np.count_nonzero(freqmat, axis = 1) #Number of SNPs present at any frequency for every strain in database
		sampscore = samepos / len(set(thisit_loc))
		sampscore[np.where(fullpos < 16)] = 0 


		logging.debug('sample average cov: ', str(scov),  'covthresh ' + str(covthresh))
		
		#for the calculation of percent overlap, only considering SNPs that pass vcoverage threshold
		allowed = (np.where((sampfreqs_known > 0.05)  & (sampcovs_known > covthresh) )[0]) 
		top = np.count_nonzero(freqmat[:,allowed], axis = 1)
		top[np.where(fullpos < 16)] = 0 #Removing incidence of SNPs present at less than 16
		top = top.astype(float)
		pcntc = np.divide(top, snpcounts, out=np.zeros_like(top), where=snpcounts!=0)

		#calculating the average frequency and coverage per strain of SNPs in sample
		totfreqs = freqmat.sum(axis = 1)
		amcov = np.divide(totfreqs, snpcounts, out=np.zeros_like(top), where=snpcounts!=0) 
		covmat = thisit_mat * sampcovs_known
		totcovs = covmat.sum(axis = 1)
		amdepth = np.divide(totcovs, snpcounts, out=np.zeros_like(top), where=snpcounts!=0)

		#total score average between two stats
		totscore = (sampscore + pcntc) / 2


		#harsher on h37rv
		if 'H37rv' in thisit_refs:
			if top[thisit_refs == 'H37rv'][0] < snpcounts[thisit_refs == 'H37rv'][0] * 0.4:
				totscore[thisit_refs == 'H37rv'] = 0

		ind = np.argsort(totscore)[::-1][0:5] # top five score indices

		topten, topten_totscore, topten_sampscore, topten_pcntc = [x[ind] for x in [thisit_refs, totscore, sampscore, pcntc]]

		logging.debug(topten, topten_totscore)

		identified = topten[0]

		maxvalue = topten_totscore[0]	

		#checking whether the identified strain passes secondary check of it's coverage significatly influencing the total coverage of SNPs in the database to the sample
		#if it doesnt, check whether anyhing else in the top 5 does. This removes spuriously identified strains
		changelist =OrderedDict()
		for poss in ind:
			query = set(map(tuple, thisit_cols[~thisit_mat[np.where(thisit_refs == thisit_refs[poss])[0]][0]]))
			tog=(query & sampsnps) - refdictsnps
			ins=[x in tog for x in sampsnps]
			ncovin = np.array(sampcov)[ins]
			ncovin = ncovin[ncovin != 0].mean()
			changelist[poss] = ncovin
			logging.debug(thisit_refs[poss], ncovin)

		change=(covin-changelist[ind[0]])/covin

		logging.debug(change)

		if change < .015 and scov > 10 :
			for poss in changelist:
				change=(covin-changelist[poss])/covin
				if change > .015:
					logging.debug(identified + ' has been changed to ' + thisit_refs[poss] + ' with change of ',  str(change), 'and ncovin of ', str(changelist[poss]))
					ind[0] = poss
					identified = thisit_refs[poss]
					changed = True
					topten, topten_totscore, topten_sampscore, topten_pcntc = [x[ind] for x in [thisit_refs, totscore, sampscore, pcntc]]
					maxvalue = topten_totscore[0]	
					break
				else:
					changed = False
			if changed == False:
				maxvalue = 0


		iterations += 1

		#if the maxvalue of this iteration passes threshold, remove the SNPs from database and sample. 
		if maxvalue > maxthresh:
			identified = thedb.sample(identified)
			identified.totscore, identified.samepos, identified.fullpos, identified.samplescore, identified.percentcommon,identified.amcov, identified.amdepth, identified.snpcounts = [x[ind[0]] for x in [totscore, samepos, fullpos, sampscore, pcntc, amcov, amdepth, snpcounts]]

			#recalculating covthresh to make sure next identified strain has coverage at least 5% of this one
			if scov < 25:
				covthresh = (0.05 * identified.amdepth )/ 0.95
			else:
				covthresh = max((0.05 * identified.amdepth )/ 0.95, 2) 
			
			#removing snps from this strain from db
			idsnps=np.nonzero(thedb.matrix[ind[0],])[0]
			thisit_idsnps = np.nonzero(thisit_mat[ind[0]])[0]
			identified.comsnppos = map(tuple, thisit_cols[thisit_idsnps])
			clearsnps = np.where(thedb.matrix[ind[0],] == False)[0]

			thedb.matrix = thedb.matrix[:,clearsnps] 
			thedb.matloc = thedb.matloc[clearsnps]
			thedb.matcols = thedb.matcols[clearsnps]

			
			#removing snps from this strain from sample
			sampfreqs_known[thisit_idsnps] = np.maximum(sampfreqs_known[thisit_idsnps] - amcov[ind[0]] , 0)
	
			thisit_remsnps = thisit_idsnps[(sampfreqs_known[thisit_idsnps] * covlist[thisit_idsnps] < sampcovs_known[thisit_idsnps])]  
			thisit_remsnps = np.append(thisit_remsnps, thisit_idsnps[(sampfreqs_known[thisit_idsnps] < .05)])
			thisit_remsnps = np.unique(thisit_remsnps)

			logging.debug('removed snps from sample :'  , str(len(thisit_remsnps)), 'snps in id strain: ', str(len(thisit_idsnps)))
			
			sampfreqs_known[thisit_remsnps] = np.nan

	
			knownlocs = np.where(~np.isnan(sampfreqs_known))[0]
			thisit_cols = thisit_cols[knownlocs]
			sampfreqs_known = sampfreqs_known[knownlocs]
			sampcovs_known = sampcovs_known[knownlocs]
			thisit_mat = thisit_mat[:,knownlocs] #matrix subsetted by only the snps that are knwon in the sample
			thisit_loc = thisit_loc[knownlocs]
			covlist = covlist[knownlocs]

			if 'H37rv' in thisit_refs:
				h37 = thisit_loc[(thisit_mat[np.where(thisit_refs == 'H37rv')[0],] == True)[0]]
				noth37 = np.where(~np.isin(thisit_loc,h37))[0]
				remaining =  np.count_nonzero(thisit_mat[:,noth37], axis = 1)
				passing = np.where((remaining >= 15) + (thisit_refs == 'H37rv'))[0]
			else:
				remaining =  np.count_nonzero(thisit_mat, axis = 1)
				passing = np.where(remaining >= 15)[0]

			thisit_refs = thisit_refs[passing]
			thisit_mat = thisit_mat[passing,]
			thedb.matrix = thedb.matrix[passing,]
			
			#getting the snp counts of only the genomes that will be evaluated
			snpcounts = thedb.matrix.sum(axis = 1) 


			refdictsnps = set(thesample.refdict.items())
			tog=(set(map(tuple, thisit_cols)) & sampsnps) - refdictsnps

			ins=[x in tog for x in sampsnps]

			freqin = np.array(sampfreq)[ins]
			freqin = freqin[freqin != 0].mean()
			ncovin = np.array(sampcov)[ins]
			ncovin = ncovin[ncovin != 0].mean()

			change=(covin-ncovin)/covin
			logging.debug(identified.sample, change, ncovin)
			if change < 0.015 and scov > 10: #applying additional filter for high coverage samples
				maxvalue = 0
			else:
				logging.debug(identified.sample +  ' totscore: ' + str(topten_totscore[0]) +  'amcov: ' + str(identified.amcov) +  ' snpcounts: ' + str(identified.snpcounts) +  ' samepos: ' + str(identified.samepos)+  ' fullpos: ' + str(identified.fullpos))
				strains.append(identified)

			covin = ncovin
			iterationres[str(iterations)] = [thisit_cols, thisit_refs , thisit_mat, sampfreqs_known, covlist]

		iterationres[str(iterations)] = [thisit_cols, thisit_refs , thisit_mat, sampfreqs_known, covlist]



	#getting abundances
	idstrains = [db.sample(x.sample) for x in strains]
	idstrains = remsubsets(idstrains,snpsampleset)
	passed = [x.sample for x in idstrains]
	strains = [x for x in strains if x.sample in passed]

	if len(idstrains ) > 1:
		un = [uniquesnps(x, idstrains) for x in idstrains]
	else:
		un = [x.snptable.items() for x in idstrains]
	sampcovs = [getsampcov(x, snpsampleset, True) for x in un]

	for x,y in zip(strains,sampcovs):
		x.finalcov = y
		x.finalrel = y/float(sum(sampcovs))


	row1= ["refname" , "totscore",  "samplescore" , "commonpcnt" , "finalrel",  "finalcov", "amcov", "amdepth" ]
	row_format ="{:>15}" * (len(row1) + 1)
	print row_format.format("", *row1)

	for idd in strains:
		print row_format.format("", *[idd.sample , str(round(idd.totscore, 3)) , str(round(idd.samplescore, 3)) , 
			str(round(idd.percentcommon , 3))		 , str(round(idd.finalrel, 3)) , str(round(idd.finalcov, 3)), str(round(idd.amcov, 3)), str(round(idd.amdepth, 3))])

	return [strains, thedb, thesample, iterationres]

def remsubsets(idstrains, samp):
	rem = []
	for pair in combinations(idstrains,2):
		un = [uniquesnps(x, pair) for x in pair]
		compos = [len([(x,y) for x,y in snps if x in samp.covtable]) for snps in un]
		if min(compos) < 16:
			sub = compos.index(min(compos))
			rem += [pair[sub].sample]
			logging.debug('Identified strains: ' +  pair[sub].sample +  ' is subset of another. Removing')
	rem = list(set(rem))
	return [x for x in idstrains if x.sample not in rem]


#function to find uniquesnps in a snpset compared to a list of other snpsets
def uniquesnps(pov, listofsnpsets):
	povset = set(pov.snptable.items())
	unique = set()
	s = True
	for x in listofsnpsets:
		if pov == x:
			continue
		if s == True:
			unique = povset - set(x.snptable.items())
		else:
			unique = unique - set(x.snptable.items())
		s = False
	return(list(unique))

#function to get coverage of a sample based on snp tuples
def getsampcov(snps, samp, tot = True):
	compos = [(x,y) for x,y in snps if x in samp.covtable]
	#tot=  [set91.covtable[x][y] for x,y in compos if y in set91.covtable[x]]
	cov=sum([samp.covtable[x][y] for x,y in compos if y in samp.covtable[x]]) / float(len(snps)) if len(snps) > 0 else 0
	if not tot:
		den = len([samp.covtable[x][y] for x,y in compos if y in samp.covtable[x]])
		cov = sum([samp.covtable[x][y] for x,y in compos if y in samp.covtable[x]]) / float(den) if den != 0 else 0

	return(cov) 

#function to get coverage of a sample based on snp tuples
def getpcntcov(snps, samp):
	compos = [(x,y) for x,y in snps if x in samp.covtable]
	pc=[samp.covtable[x][y] for x,y in compos if y in samp.covtable[x]]
	pc = len([x for x in pc if x >= 2]) / float(len(snps)) if len(snps) > 0 else 0
	return (pc)

#function to find uniquesnps in a snpset compared to another snpset
def unsnps(pov, other):
	un=list(set(pov.snptable.items()) - set(other.snptable.items()))
	return(un)

def comsnpget(s):
	return (set([(x, y.keys()[0], y.values()[0]) for x,y in s.comsnppos.items()]))



def decide(db,otherans, sample, pov = None):
	print 'pov stats are ', pov.sample, pov.totscore
	close = [len(comsnpget(x) & comsnpget(pov)) / len(comsnpget(pov)) for x in otherans]
	closerev = [len(comsnpget(x) & comsnpget(pov)) / len(comsnpget(x)) if len(comsnpget(x)) > 0 else 0 for x in otherans ]
	#similargens  = [y for x,y in close if x > 0.95]
	closet = zip(close,closerev,otherans)
	similargens = [z for x,y,z in closet if x > 0.95 if y < 0.96] 
	print [(x,y,z.sample) for x,y,z in closet]
	if len(similargens) > 0:
		unsims = [unsnps(x,pov) for x in similargens]
		pcovs = [getpcntcov(s,sample) for s in unsims]
		pcovstr = zip(pcovs, similargens)
		print [(x, y.sample, y.totscore) for x,y in pcovstr]
		if max(pcovstr)[0] > 0.25:
			chosen = max(pcovstr)[1]
			return (chosen)
			# [otherans.remove(y) for x,y in pcovstr if y != chosen]
			# sets = otherans
		else:
			#[otherans.remove(y) for x,y in pcovstr]
			#otherans += [pov]
			#sets = otherans
			return pov
	else:
		#sets = otherans + [pov]
		return pov
	# un = [uniquesnps(x, sets) for x in sets]
	# sampcovs = [getsampcov(x, sample) for x in un]
	# # filt = [y for x,y in zip(sampcovs, [x.sample for x in sets]) if x > 1]
	# # if len(filt) == 1:
	# # 	return db.sample(filt[0])
	# # elif len(filt) != len(sets):
	# # 	sets = [db.sample(x) for x in filt]
	# # 	un = [uniquesnps(db.sample(x), sets) for x in filt]
	# # 	sampcovs = [getsampcov(x, sample) for x in un]
	# return(max(zip(sampcovs, sets))[1])



#removes snps from database object
def removepos(db, posdict, remrefs = False):
	for ref in db.samplesets:
		if not ref.masked:
			ref.removesnps(posdict)

class Snpdb(object):
	'''Database class. Basically a collection of Snpsets. Needs as input a list of Snpsets or list of vcf files. Reducedist is the minimum distance between snps in a genome
	Default is 25. Filterdist is the minimum SNP distance between genomes in database. Default is 100. If No filtering desired. Set filterdist to None. 
	'''
	def __init__(self, input, bufferr = False, reducedist = 25, filterdist = None ):
		self.input = input
		self.reducedist = reducedist
		self.filterdist = filterdist
		self.getsets()
		if bufferr:
			self.buffersamples()
		
	def getsets(self):
		if isinstance(self.input, list):
			if isinstance(self.input[0], str):
				thesets = []
				count = 0
				for sample in self.input:
					if count % 50 == 0 and count != 0:
						logging.info("Added " + str(count) + " genomes.")
					thesets.append(Snpset(sample))

					count += 1
			elif isinstance(self.input[0], Snpset):
				thesets = self.input
			else:
				raise Exception ("Please provide a list of vcf filenames or pkl filenames or SnpSets")
			self.samplesets =  thesets
			self.sampnames = [sett.sample for sett in self.samplesets]
			logging.info('Removing nearby SNPs in each genome of distance ' + str(self.reducedist) + ' to each other')
			self.reducesnpsets(snpdist = self.reducedist)
			if self.filterdist is not None:
				logging.info('Filtering genomes in the database on distancee ' + str(self.filterdist) + 'to each other')
				self.filtersnpdb(thresh = self.filterdist)
			logging.info('Making reference strain')
			self.makerefstrain()
			logging.info('Making a SNP matrix')
			self.makematrix()

	#unique positions in database
	def uniquepos(self):
		pos = []
		for sett in self.samplesets:
			pos += sett.snptable.keys()
		return set(pos)

	def uniquesnps(self):
		uniques = []
		snps = set().union(*[sett.snptable.items() for sett in self.samplesets])
		#snps = set([sett.snptable.items() for sett in self.samplesets])
		return snps
	def sample(self, name):
		return self.samplesets[self.sampnames.index(name)] 

	def snpoverlap(self, set2):
		#unsnps = dict(self.uniquesnps())
		over = []
		for snp in self.uniquesnps():
			pos =  snp[0]
			nuc =  snp[1]
			if pos in set2.snptable:
				if nuc == set2.snptable[pos]:
					over.append(pos)
				elif isinstance(set2.snptable[pos], dict):
					if nuc in set2.snptable[pos].keys():
						over.append(pos)
		return list(set(over)) #to account for situations where ambigous bases match to two different samples

	def snpoverlap2(self, set2):
		#unsnps = dict(self.uniquesnps())
		return list(self.uniquepos() & set(set2.snptable.keys()))
	def snpfreq(self):
		#unsnps = dict(self.uniquesnps())
		counts = Counter(sum([sett.snptable.keys() for sett in self.samplesets], []))
		return counts #to account for situations where ambigous bases match to two different samples

	def removeref(self, sampname):
		self.samplesets = [x for x in self.samplesets if x.sample != sampname]
		self.sampnames = [x.sample for x in self.samplesets]
		#if not self.counts:
		#	del self.counts
		#self.makematrix()

	def makerefstrain(self, cutoff = 0.75):
		counts = self.possnpcounts()
		snpamount = round(len(self.sampnames) * cutoff)
		commonsnps=[x[0] for x in counts if counts[x] > snpamount]
		refseq = {snp:refdata[snp] for snp in commonsnps}
		self.samplesets.append(Snpset(fileName="H37rv", data = refseq))
		self.sampnames = [sett.sample for sett in self.samplesets]

	def filtersnpdb(self, thresh = 150):
		logging.info('Filtering might take some time')
		sampnames = copy.deepcopy(self.sampnames)
		dists = {}
		track = 0
		for samp in self.samplesets:
			track += 1
			if track % 100 == 0 and track != 0:
				logging.info ("Processed " + str(track) + " genomes.")
			sampnames.remove(samp.sample)
			for x in sampnames:
				aset= set(samp.snptable.items())
				bset = set(self.sample(x).snptable.items())
				sd = (aset - bset).union(bset - aset)
				sdloc = [y[0] for y in sd]
				badposs = set(samp.badpos + self.sample(x).badpos)
				distt = len(set(sdloc) - badposs)
				if samp.sample not in dists:
					dists[samp.sample] = {x:distt}
				else:
					dists[samp.sample][x] = distt
	

		def most_common(lst):
			return max(set(lst), key=lst.count)

		close = [0]
		remm = []
		while (len(close) > 0):
			close=[x for x in dists if any([snpc < thresh for snpc in dists[x].values()])]
			if len(close) == 0:
				break
			allclose = sum([[x for x in dists[y] if dists[y][x] < thresh] for y in close], [])
			toremove = most_common(allclose)
			remm +=  [toremove]
			del dists[toremove]
			for gen in dists.keys():
				dists[gen] = {x:y for x,y in dists[gen].items() if x != toremove}
			close=[x for x in dists if any([snpc < thresh for snpc in dists[x].values()])]
			logging.debug(len(close))

		logging.info("Removing references " + ",".join(remm) + "references from db \n")
		for ref in remm:
			self.removeref(ref)

	def buffersamples(self):
		allpos = list(self.uniquepos())
		for sample in self.samplesets:
			for snp in allpos:
				if snp not in sample.snptable and snp not in sample.badpos:
					sample.snptable[snp] = refdata[snp]
			sample.getsnpinfo()
	def possnpcounts(self, snp = True):
		if not hasattr(self, 'counts'): 
			if snp:
				self.counts = Counter(sum([sett.snptable.items() for sett in self.samplesets], []))
			else:
				self.counts = Counter(sum([sett.snptable.keys() for sett in self.samplesets], []))
		return(self.counts)
	def snpfrequency(self, snptup):
		count = self.possnpcounts()[snptup]
		return(count)
	def weighedsnppersample(self):
		counts = self.possnpcounts()
		for sample in self.samplesets:
			sample.weighedsnp =  [countallele[snp] for snp in sample.snptable.items()]

	def reducesnpsets(self, snpdist = 25): #removes positions per sample that are within x bases of each other
		track = 0
		for ref in self.samplesets:
			track += 1
			if track % 100 == 0 and track != 0:
				logging.info ("Processed " + str(track) + " genomes.")
			closeby = np.where(np.diff(sorted(ref.snptable.keys())) < snpdist)[0]
			closeby = np.append(closeby + 1, closeby)
			toremove = [y for x,y in enumerate(sorted(ref.snptable.keys())) if x in closeby]
			ref.badpos += toremove
			[ref.snptable.pop(x, None) for x in toremove]
			ref.getsnpinfo()
	
	def makematrix(self): #adding matrix attributes to db
		del self.counts
		allsnps = set(self.possnpcounts())
		dbsets = [set(x.snptable.items()) for x in self.samplesets]
		dblist = []

		mat = np.empty((0, len(allsnps)), dtype = bool)
		num = 0
		for ref in dbsets:
			num += 1
			if num % 50 == 0 and num != 0:
				if mat is None:
					mat = np.array(dblist)
					dblist = []
				else:
					mat = np.concatenate((mat,np.array(dblist)), axis = 0)
					dblist = []
			dblist.append([True if x in ref else False for x in allsnps])
		dblist = np.array(dblist)

		mat = np.concatenate((mat, dblist), axis = 0)
		mat = sparse.csr_matrix(mat)
		self.matrix=mat
		self.matcols = np.array(list(allsnps), dtype=object)
		self.matrows = np.array(self.sampnames)
		self.matloc = np.array([x for x,y in allsnps])

	def addstrains(self, files):
		if not isinstance(files, list):
			files = [files]
		if isinstance(files[0], str):
			newsets = []
			count = 0
			for sample in files:
				if count % 50 == 0 and count != 0:
					logging.info ("Added " + str(count) + " new genomes.")
				newsets.append(Snpset(sample))
				count += 1
			#thesets = [Snpset(sample) for sample in self.input]
		elif all([isinstance(x, Snpset) for x in files]):
			newsets = files
		else:
			raise Exception ("Please provide a list of vcf filenames or pkl filenames or SnpSets")
		self.samplesets += newsets
		self.sampnames = [sett.sample for sett in self.samplesets]

	def save(self, outname):
		with open(outname + ".db" , "wb") as outfile:
			pickle.dump(self, outfile, pickle.HIGHEST_PROTOCOL)
		
def getvcf(fastas, tempdir = "temp", keepint = False):
	logging.info("Getting bam file and pileup for " + " ".join(fastas))

	if not os.path.isdir(tempdir):
		os.makedirs(tempdir)

	if len(fastas) == 0:
		logging.error("No fasta files")
		return (None)
	if len(fastas) == 1:
		logging.info ("Using just 1 fastq file ")
	fastq = " ".join(fastas)

	name = os.path.splitext(os.path.basename(fastas[0]))[0]
	bamloc = tempdir + "/" + name + ".bam"
	vcfloc = tempdir + "/" + name + ".vcf"

	aln = ""
	index = ""
	getvcf = ""
	
	if not os.path.isfile(vcfloc):
		pilon = "java -jar " + pkg_resources.resource_filename('quanttb', 'data/pilon-1.22.jar')

		getvcf = pilon + " --genome " + ref + " --bam " + bamloc +  " --output " + name + " --outdir " + tempdir + " --vcf" + " --threads 4 --fix none"
		if not os.path.isfile(bamloc):
			aln = "bwa mem -t 4 " + ref + " " + fastq +  " | samtools sort - -T sam_tempsort" + name + " > " + bamloc
			index = "samtools index " + bamloc
		else:
			logging.info("Already have bam for  " + name + " Not remaking")
	else:
		logging.info("Already have vcf for  " + name + " Not remaking")

	vcflog = open(tempdir + "/" + name + ".log", "a")
	for command in (aln,  index, getvcf):
		if command != '':
			subprocess32.call(command, shell = True, stderr=vcflog, timeout = 10000, stdout = vcflog)

	if not keepint and os.path.isfile(bamloc):
		[os.remove(x) for x in [bamloc, bamloc + ".bai"]]
	return vcfloc
