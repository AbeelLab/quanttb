from classify import *
from runnucmer import *


class SpecialFormatter(logging.Formatter):
	FORMATS = {logging.DEBUG :"DBG: %(module)s: %(lineno)d: %(message)s",
			   logging.ERROR : "ERROR: %(message)s",
			   logging.INFO : "%(message)s",
			   'DEFAULT' : "%(levelname)s: %(message)s"}

	def format(self, record):
		self._fmt = self.FORMATS.get(record.levelno, self.FORMATS['DEFAULT'])
		return logging.Formatter.format(self, record)

hdlr = logging.StreamHandler(sys.stderr)
hdlr.setFormatter(SpecialFormatter())
logging.root.addHandler(hdlr)


def makeoutput(qtobj):
		#out = self.output
		if "/" in qtobj.out:
			qtobj.outdirr = os.path.dirname(qtobj.out)
		else:
			qtobj.outdirr = "output"
			qtobj.out = qtobj.outdirr + "/" + qtobj.out

		if not os.path.isdir(qtobj.outdirr):
			os.makedirs(qtobj.outdirr)
		
		qtobj.temp = qtobj.outdirr + "/temp"
		if not os.path.isdir(qtobj.temp):
			os.makedirs(qtobj.temp)


def getvcfs(fqlist, qtobj):
	if which('samtools') is None or which('bwa') is None:
		logging.error('To run this, both samtools (v 1.7 or higher) and bwa (v. 0.7.17 or higher) need to be installed in your system and in your path')
		logging.info('samtools download: https://sourceforge.net/projects/samtools/files/samtools/1.7/' )
		logging.info('BWA download: https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2/download')
		sys.exit()

	sets = []
	for fq in fqlist: #should check if samples are fq
		if not all(os.path.isfile(readset) for readset in fq):
			logging.error ("Some readsets don't exist:")
			logging.info("\n".join ([readset for readset in fq if not os.path.isfile(readset)]))
			continue
		samplevcf = getvcf(fq, tempdir = qtobj.temp, keepint = qtobj.keepin)
		sets.append(samplevcf)
	return (sets)

class QuantTBMain(object):
	#from classify import Snpdb, Snpset
	

	def __init__(self):
		usage = '\r{}\nUsage: %(prog)s <command> [options]\nCommand: makesnpdb\t Make a reference SNP database\n\t quant\t\t Quantify sample with a ref SNP db\n\t variants\t Generate a vcf from sequencing readsets\n'.format(
			'Program: QuantTB (Detection of mixed infections)\nVersion: 1.01\n'.ljust(len('usage:')))

		parser = argparse.ArgumentParser(
			prog='quantTB', usage=usage, add_help=False)

		parser.add_argument('command', help='Subcommand to run')

		if not sys.argv[1:2]:
			parser.print_usage()
			exit(1)

		args = parser.parse_args(sys.argv[1:2])

		if not hasattr(self, args.command):
			print 'Unrecognized command: {}'.format(args.command)
			# parser.print_help()
			exit(1)

		getattr(self, args.command)()

	def makesnpdb(self):
		usage = '\r{}\nUsage: %(prog)s makesnpdb [options] <-g>\n'.format(
			''.ljust(len('usage:')))

		parser = argparse.ArgumentParser(
			description='', usage=usage, add_help=False)


		required = parser.add_argument_group('Required arguments')
		required.add_argument('-g', dest='dbfiles', nargs='+', 
		help='Files you want to use to make the reference database, Can be either .fa, .fna, .fasta, .snps .vcf(.gz), or .samp', type = str, required = True)
	
		optional = parser.add_argument_group('Optional arguments')
		optional.add_argument('-reducedist', dest = 'redist', type = int, help = 'When making database, what is the minimum distance between SNPs in a genome', default = 25)
		#optional.add_argument('-filterdist', dest = 'filterdist', type = int, help = 'When making database, what is the minimum distance between SNPs?', default = None)
		optional.add_argument('-o', dest = 'output', type = str, help = 'Directory/File where  you want the snpdb to be written to', default = 'output/snpdb')

		optional.add_argument('-k', dest = 'keep', action = 'store_true', help = 'Keep temp files?')
		optional.add_argument('-l', dest='log_level', help='Set the logging level',
				   choices=['DEBUG', 'ERROR', 'CRITICAL'], default = "INFO")


		if not sys.argv[2:]:
			parser.print_help()
			exit(1)

		args = parser.parse_args(sys.argv[2:])

		if args.log_level:
			logging.getLogger().setLevel(getattr(logging, args.log_level))

		self.out = args.output
		makeoutput(self)
		if args.keep:
			self.keepin = True
		else:
			self.keepin = False

			
		dbgenomes = args.dbfiles
		if not all(os.path.isfile(genome) for genome in dbgenomes):
			logging.error ("Some db files don't exist:")
			logging.info("\n".join ([genome for genome in dbgenomes if not os.path.isfile(genome)]))
			sys.exit()
		if all([x.endswith(('.fna', '.fasta', '.fa')) for x in dbgenomes]):
			logging.info('Running nucmer on ' + str(len(dbgenomes)) + ' fasta files to obtain SNPs relative to H37Rv reference')
			if which('nucmer') is None or which('show-snps') is None:
   				logging.error('To run this, Mummer(nucmer & show-snps) needs to be installed in your system and in your path')
   				sys.exit
			count = 0
			newloc = []
			for genome in dbgenomes:
				if count % 50 == 0 and count != 0:
					logging.info("Processed " + str(count) + " genomes.")
				count += 1
				genomename = os.path.splitext(os.path.basename(genome))[0]
				outname = os.path.join(self.temp, genomename)
				runnucmer(genome, outname)
				newloc.append(getnucsnps(outname))
				os.remove(outname + ".delta")
			dbgenomes = newloc
		logging.info("Making reference db out of " + str(len(dbgenomes)) + " ref genomes, saving in " + self.out + '.db')
		db = Snpdb(dbgenomes, reducedist = args.redist)
		db.save(self.out)
		if not self.keepin:
			rmtree(self.temp)
	
	def quant(self):
		usage = '\r{}\nUsage: %(prog)s quant [options] <-db> <-s>\n'.format(
			''.ljust(len('usage:')))
		parser = argparse.ArgumentParser(
			description='', usage=usage, add_help=False)

		#required = parser.add_argument_group('Required arguments')
		optional = parser.add_argument_group('Optional arguments')
		optional.add_argument('-v', dest='vcfsamples', nargs='+', help='VCF(s) or snp profiles that you want tested against the refdb, can either be .vcf(.gz) or .samp')
		optional.add_argument('-f', dest='fastq', nargs='+', action = 'append' ,
		help='Fastq Readset(s) that you want tested against refdb, can specify multiple times for multiple pairs, (.fq, .fastq)')
		optional.add_argument('-db', dest='db', help='Location of reference SNP DB file. (.db) If not supplied a default TB db will be used', type = str)
		optional.add_argument('-o', dest = 'output', type = str, help = 'Directory/File where  you want results written to', default = 'output/results.txt')
		optional.add_argument('-resout', dest = 'resprint', action = "store_true", help = 'Should stats from each run be output')
		optional.add_argument('-i', dest = 'iters', type = int, help = 'Number of iterations for classificaiton', default = 8)
		optional.add_argument('-abres', dest = 'abres', action='store_true', help = 'Should resistances from each sample be output?')
		optional.add_argument('-k', dest = 'keep', action = 'store_true', help = 'Keep temp files?')
		optional.add_argument('-l', dest='log_level', help='Set the logging level',
                   choices=['DEBUG', 'ERROR', 'CRITICAL'], default = "INFO")


		if not sys.argv[2:]:
			parser.print_help()
			exit(1)

		args = parser.parse_args(sys.argv[2:])

		if args.log_level:
			logging.getLogger().setLevel(getattr(logging, args.log_level))

		self.out = args.output
		makeoutput(self)

		if args.vcfsamples is None and args.fastq is None:
			logging.error ("Must supply VCF or fastq samples to classiy")
			sys.exit()
		samples = args.vcfsamples


		if args.keep:
			self.keepin = True
		else:
			self.keepin = False

		if args.abres:
			abres =self.outdirr + "/" + 'antibioticres.txt'
			if os.path.isfile(abres):
				os.remove(abres)

		if args.fastq is not None:
			sets = getvcfs(args.fastq, self)
			if samples is None:
				samples = sets
			else:
				samples += sets


		if not all(os.path.isfile(genome) for genome in samples):
			logging.error ("Some sample files don't exist:")
			logging.info("\n".join ([genome for genome in samples if not os.path.isfile(genome)]))
			sys.exit()

		if not args.db:
			self.db = resource_filename(Requirement.parse('quanttb'), 'quanttb/data/snpdb_100snpdist.db')
		else:
			self.db = args.db
			if not os.path.isfile(self.db):
				logging.error('Snpdb file does not exist')
				logging.info(self.db)
				sys.exit()
			if not self.db.endswith('.db'):
				logging.info('Snpdb needs to be a db file (.db)')
				sys.exit()

		logging.info("Using reference SNP db from " + os.path.basename(self.db) )
		self.db = openFile(self.db)


		samplenames = []

		for sample in samples:
			if sample.endswith((".samp",".pkl")):
				logging.info("Using snps from " + sample )
				samplesnps = openFile(sample)
				samplenames.append(samplesnps.sample)
			elif sample.endswith(('.vcf', '.vcf.gz')):
				logging.info("Extracting snps from " + sample)
				if "H37rv" in self.db.sampnames:
					includeme = self.db.sample("H37rv").snptable.keys()
				else:
					includeme = None
				samplesnps = Snpset(sample, False, include = includeme)

				num = 1
				if samplesnps.sample in samplenames:
					while(samplesnps.sample + '_' + str(num) in samplenames):
						num += 1
					samplesnps.sample += '_' + str(num)

				if self.keepin:
					samplesnps.save(self.outdirr + "/" + samplesnps.sample)
				samplenames.append(samplesnps.sample)
			else:
				logging.error('Improper file type')
				sys.exit()

			if args.abres:
				samplesnps.findres()
				with open(abres, 'ab') as f:
					if not os.path.isfile(abres):
						f.write("sample,drug,restype,dist\n")
					for drug,resis in samplesnps.resistances.items():
						for base in resis.values():
							f.write( ",".join([samplesnps.sample, drug,base]) + "\n")
			
			logging.info('Starting iterations')
			res = iteration(self.db, samplesnps, iterationmax = args.iters)

			if not os.path.isfile(self.out):
				if args.resprint:
					header = "sample, refname, totscore, samplescore, commonpcnt, relabundance, depth, amcov, amdepth,  fullpositions, samepos, snpcounts\n"
				else:
					header = "sample, refname, totscore, relabundance, depth\n"
			else:
				header = None


			with open(self.out, "ab") as outputfile:
				if header: outputfile.write(header)
				if res[0] is not None:
					for idd in res[0]:
						if args.resprint:
							outputfile.write(",".join([samplesnps.sample, idd.sample , str(round(idd.totscore, 3)) , str(round(idd.samplescore, 3)), str(round(idd.percentcommon, 3)), 
								str(round(idd.finalrel, 3)), str(round(idd.finalcov, 3)), str(round(idd.amcov, 3)), str(round(idd.amdepth, 3)), 
								str(idd.fullpos), str(round(idd.samepos, 3)), str(idd.snpcounts)])  + "\n")
						else:
							outputfile.write(",".join([samplesnps.sample, idd.sample , str(round(idd.totscore, 3)) ,str(round(idd.finalrel, 3)), str(round(idd.finalcov, 3))])  + "\n")
	
		if not self.keepin and os.path.isdir(self.temp):
			logging.debug('Removing temporary files')
			rmtree(self.temp)

	def variants(self):
		usage = '\r{}\nUsage: %(prog)s variants [options] <-f> \n'.format(
			''.ljust(len('usage:')))
		parser = argparse.ArgumentParser(
			description='', usage=usage, add_help=False)

		required = parser.add_argument_group('Required arguments')
		required.add_argument('-f', dest='fastq', nargs='+', action = 'append' ,
		help='Fastq Readset(s) that you want converted to a vcf (using samtools,bwa and pilon), can specify multiple times for multiple pairs, (.fq, .fastq)')
		optional = parser.add_argument_group('Optional arguments')
		optional.add_argument('-o', dest = 'output', type = str, help = 'Directory you want vcf written to', default = 'output/')
		optional.add_argument('-k', dest = 'keep', action = 'store_true', help = 'Keep temp files?')
		optional.add_argument('-l', dest='log_level', help='Set the logging level',
                   choices=['DEBUG', 'ERROR', 'CRITICAL'], default = "INFO")
		if not sys.argv[2:]:
			parser.print_help()
			exit(1)

		args = parser.parse_args(sys.argv[2:])

		self.out = args.output

		makeoutput(self)

		if args.keep:
			self.keepin = True
		else:
			self.keepin = False

		sets = getvcfs(args.fastq, self)
		
	



def main():
	QuantTBMain()

if __name__ == '__main__':
	main()


