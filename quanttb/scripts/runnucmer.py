import subprocess32 #
import os
import sys
import argparse
import logging
import pkg_resources
#run nucmer


ref = pkg_resources.resource_filename('quanttb', 'data/GCF_000277735.2_ASM27773v2_genomic.fna')
#ref = "/tudelft.net/staff-bulk/ewi/insy/tbdiagnostics/aligner/index/GCF_000277735.2_ASM27773v2_genomic.fna"

def runnucmer(fastafile, outputfile):
	if not os.path.isfile(fastafile):
		logging.error("File doesn't exist")
		return None
	elif not fastafile.endswith((".fna", '.fa', '.fasta')):
		logging.error(fastafile + " is not an fna file \n")
		return None
	logging.debug ("Running nucmer for " +  fastafile )

	cmd = "nucmer --maxmatch -c 100 -p " + outputfile + " "  + ref + " " + fastafile 
	try:
		subprocess32.call(cmd, shell = True)
	except Exception as e:
		logging.error ("Nucmer failed for following reason for " +  fastafile + "\n")
		logging.info(e)

def getnucsnps(df):
	deltafile = df + ".delta"
	if not os.path.isfile(deltafile):
		logging.error(deltafile + " doesn't exist \n")
		return None
	logging.debug ("Getting nucmer snps for " +  df)
	snpout = df +  ".snps"
	cmd = "show-snps -CrIT " + deltafile + " > "  + snpout
	try:
		subprocess32.call(cmd, shell = True)
	except Exception as e:
		logging.error ("getting snps failed for following reason for " +  deltafile + "\n")
		logging.info(e)
	return (snpout)

# dependencies
def which(pgm):
    path=os.getenv('PATH')
    for p in path.split(os.path.pathsep):
        p=os.path.join(p,pgm)
        if os.path.exists(p) and os.access(p,os.X_OK):
            return p


if __name__ == '__main__':

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



	parser = argparse.ArgumentParser()
	parser.add_argument('-f', '--list', dest='filenames', nargs='+', help='<Required> Files you want to get compared to ref genome', required=True)
	parser.add_argument('-o', dest = 'output', type = str, help = 'Directory where  you want results written to', default = '.')
	#parser.add_argument('-amb', dest='ambig', nargs='+', help='Include ambiguous sites?')
	parser.add_argument('-l', dest='log_level', help='Set the logging level',
                   choices=['DEBUG', 'ERROR', 'CRITICAL'], default = "INFO")
	args = parser.parse_args()
	if args.log_level:
   		logging.getLogger().setLevel(getattr(logging, args.log_level))

   	if which('nucmer') is None or which('show-snps') is None:
   		logging.error('To run this, Mummer(nucmer & show-snps) needs to be installed in your system and in your path')
   		sys.exit

	if not all(os.path.isfile(genome) for genome in args.filenames):
		logging.error ("Some files don't exist:\n")
		logging.info("\n".join ([genome for genome in args.filenames if not os.path.isfile(genome)]))
		sys.exit()

	outdirr = args.output

	if not os.path.isdir(outdirr):
		os.makedirs(outdirr)

	logging.info ('Running nucmer and getting snps for '+ str(len(args.filenames)) + ' files. Saving files to ' + os.path.realpath(outdirr))
	count = 0
	for genome in args.filenames:
		if count % 50 == 0 and count != 0:
			logging.info("Processed " + str(count) + " genomes.")
		count += 1
		genomename = os.path.splitext(os.path.basename(genome))[0]
		outname = os.path.join(outdirr, genomename)
		runnucmer(genome, outname)
		getnucsnps(outname)
		os.remove(outname + ".delta")

