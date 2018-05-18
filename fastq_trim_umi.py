from Bio import SeqIO, bgzf
from Bio.Seq import Seq
from Bio.Alphabet import SingleLetterAlphabet
from Bio.SeqRecord import SeqRecord
import gzip
import argparse

parser = argparse.ArgumentParser(prog='fastq_trim_umi')
parser.add_argument('-i', help='input fastq.gz file', dest='input_fastq_gz')
parser.add_argument('-o', help='output fastq.gz file', dest='output_fastq_gz')
parser.add_argument('-l', help='length of the UMI barcode', dest='umi_len', default=12, type=int)
args = parser.parse_args()

ifilename = args.input_fastq_gz
ofilename = args.output_fastq_gz
umilen = args.umi_len

with gzip.open(ifilename, "r") as handle, bgzf.BgzfWriter(ofilename, "wb") as fout:
    for rec in SeqIO.parse(handle, "fastq"):
	umi = rec.seq[0:umilen]

        rid = rec.id + ":" + str(umi)
        rseq = Seq(str(rec.seq)[umilen:], SingleLetterAlphabet)
        rq = rec.letter_annotations["phred_quality"][umilen:]

        nrec = SeqRecord(rseq, id=rid, description="")
        nrec.letter_annotations["phred_quality"] = rq

        SeqIO.write(sequences=nrec, handle=fout, format="fastq")
