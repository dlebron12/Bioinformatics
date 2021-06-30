from __future__ import print_function
import pyfasta
import pybedtools
import argparse
import logging
import os
from subprocess import Popen
from subprocess import PIPE

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()

MAX_ID_LEN = 100

def parse_commandline():
    parser = argparse.ArgumentParser(
        description="A utility to validate a bed file is compatibile with a fasta file/genome")

    parser.add_argument("--fasta", "-f", type=argparse.FileType('r'), help="A fasta file")
    parser.add_argument("--bed", "-b", type=argparse.FileType('r'), help="A bed file")
    parser.add_argument("--shorten", action='store_true', default=False,
                        help="Truncate chromosome names to the first white-space character. Creates .shrt files with shortened names.")
    parser.add_argument("--compressed","-c",action='store_true',default=False,help="If the fasta file is in a compressed format, do not unzip.")

    return parser.parse_args()

def validate_chromosome_names(fasta, bed):
    logger.info("Found the following contig names:")
    [logger.info("\t{}".format(chrm)) for chrm in fasta.keys()]
    n_errors = 0
    for entry in bed:
        if entry.chrom not in fasta:
            logger.info("Unrecognized chromosome: {}".format(entry.chrom))
            n_errors += 1
    logger.info("Identified {} invalid entries.".format(n_errors))

if __name__ == "__main__":
    args = parse_commandline()

    if args.shorten:
        fasta_name_base, fasta_ext = os.path.splitext(args.fasta.name)
        ## Add shrt to the filename to indicate short headers in both fasta and bed
        fasta_filename = ".".join([fasta_name_base, "shrt"]) + fasta_ext
        bed_name_base, bed_ext = os.path.splitext(args.bed.name)
        bed_filename = ".".join([bed_name_base, "shrt"]) + bed_ext
        
        with open(fasta_filename, 'wb') as out_fasta:
            fasta = pyfasta.Fasta(args.fasta.name)
            for k in fasta:
                out_fasta.write(">" + k.split()[0] + "\n")
                [out_fasta.write(ss) for ss in fasta[k]]
                out_fasta.write("\n")
            original_fasta_n_entries = len(set(fasta.keys()))
        with open(bed_filename, 'wb') as out_bed:
            bed = pybedtools.BedTool(args.bed.name)
            for l in bed:
                l.chrom = l.chrom.split()[0]
                out_bed.write(str(l) + "\n")

    else:
        fasta_filename = args.fasta.name
        bed_filename = args.bed.name
        
    ## If the fasta file is compressed then:
        
    if not args.compressed:
        logger.info("Loading fasta file: {}".format(fasta_filename))
        fasta = pyfasta.Fasta(fasta_filename)

        for k in fasta.keys():
            assert len(k) <= 100

        if args.shorten:
            assert len(set(fasta.keys())) == original_fasta_n_entries

        logger.info("Loading bed file {}".format(bed_filename))
        bed = pybedtools.BedTool(bed_filename)
        logger.info("Checking bed file's chromosome names are contained in the fasta file")
        validate_chromosome_names(fasta, bed)
    
    else: 
        if not (fasta_filename.endswith(".gz") or fasta_filename.endswith(".zip")):
            logger.info("Compressed fasta needs to be in .gz or .zip format")
            ## Raise Error
        else:
            #output=check_output(f"zcat {fasta_filename} | grep '>'", shell=True)
            p1 = Popen(["zcat" , fasta_filename], stdout=PIPE)
            p2 = Popen(["grep", ">"], stdin=p1.stdout, stdout=PIPE)
            p1.stdout.close()
            output = p2.communicate()[0]
            
            fasta_chrms=str(output).replace("b'>","").split("\\n>")
            bed = pybedtools.BedTool(bed_filename)

            ## check which in bed are not in fasta_chrms
            for entry in bed:
                if entry.chrom not in fasta_chrms:
                    logger.info("Bed file Chromosome {} is not found in fasta".format(entry.chrom))
                else:
                    logger.info("{} found in fasta".format(entry.chrom))
