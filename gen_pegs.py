"""
Author: Will Chen
Date: 11/03/2020
Input: 
1. target sequence
2. start end 0-based location of expected deltion OR
   deletion size 
3. 
Output: designed txt file
Useage: 
"""
import os,argparse
import numpy as np
import primedel.design as pde

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_seq",dest="sequence", metavar='fasta_file',help="A fasta file with a single sequence")
parser.add_argument("-f", "--file", dest="scored_file", metavar='scored_gRNA_file',help="Scored file from flashfry, CRISPOR or MIT gpp")
parser.add_argument("-s", "--dsize", dest="size_range",nargs=2, metavar=('dmin','dmax'), type=int, default=(100,500),help="deletion size range")
parser.add_argument("-x", "--dpos", dest="pos_range",nargs=2, metavar=('start','end'),type=int, help="0-based position of start and end of deletion, including any peg-pair candidates within ± 50bp of position provided")
parser.add_argument("-p", "--precise", dest="precise", action="store_true",default=False, help="Use with -x option, will return a few peg-pair candidates that generate exact deletion expected")
parser.add_argument("-t", "--threshold", dest="threshold",metavar=('1-99'),default= 50, help="threshold to filter out poor gRNAs, default set to 50 as CRISPROR recommended. Not applied to MIT gpp")
parser.add_argument("-l", "--len_homo", dest="homology_len",metavar=('10-50'),default= 30, help="length of homology in pegRNA for recombination")
parser.add_argument("-o", "--outdir", dest="outdir",metavar=('output_dir'),default= os.getcwd(),help="path for output files, default current working directory")
parser.add_argument("-n", "--name", dest="fname",metavar=('fname'),  default='target',help="file name for out put files")
args = parser.parse_args()

seq = pde.read_fasta(args.sequence) # only support .fa files for now
args.outdir += '/' if args.outdir[-1]  != '/' else ''
if not os.path.isdir(args.outdir):
    os.makedirs(args.outdir)

if args.scored_file:
    if '.scored' in args.scored_file:
        guides = pde.read_flashfry(args.scored_file, args.threshold)
    elif '.txt' in args.scored_file:
        guides = pde.read_gpp_designer(args.scored_file)
    elif '.xls' in args.scored_file:
        guides = pde.read_crispor(args.scored_file, args.threshold)
else:
    guides = pde.gen_guides(seq)

model=np.load('/'.join([os.path.dirname(pde.__file__),'indel_ratio.npz']))
fds = [(s[0], s[1],  # get guide and nick site
        13 + (pde.softmax(np.dot(pde.onehotencoder(s[0]),model['weights'])+model['bias'])[0]<0.75)//1 + # if ratio nick at 16 is >0.25
        ((s[0].upper().count('G')+s[0].upper().count('C'))/20<0.4)//1,s[-1])  # if the GC content is lower than 0.4 get PBS length
       for s in guides if s[2]=='FWD' and s[1]>30 ] 
rvs = [(s[0],s[1],  # get guide and nick site
        13 + (pde.softmax(np.dot(pde.onehotencoder(s[0]),model['weights'])+model['bias'])[0]<0.75)//1 + 
        ((s[0].upper().count('G')+s[0].upper().count('C'))/20<0.4)//1,s[-1])
       for s in guides if s[2]=='REV' ]

if args.pos_range:
    if args.precise:
        pairs = pde.peg_design_by_start_end(fds,rvs,seq,args.pos_range,args.homology_len,p=args.precise)
    else:
        pairs = pde.peg_design_by_start_end(fds,rvs,seq,args.pos_range,args.homology_len)
else:
    pairs = pde.peg_design_by_size(fds,rvs,seq,args.size_range,args.homology_len)

# Writing stats files and fastq files 
fname = args.outdir + args.fname + '_pegpairs_design.txt'

np.savetxt(fname,pairs,delimiter='\t',fmt='%s',
          header='RNA_1\thomology+PBS_1\tnick_1\tgRNA_2\thomology+PBS_2\tnick_2\tComments\tDeletion size\tExpected deletion result',comments='')