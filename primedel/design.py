import numpy as np
import pandas as pd
import regex as re

def reverse_complement(seq):
    """This function returns the reverse_complement sequence of the input sequence
    from 3' to 5' """
    complement = {'A':'T', 'C':'G','G':'C', 'T':'A', 'N':'N'}
    rcomp_seq = ''
    for base in seq:
        rcomp_seq = complement[base] + rcomp_seq   
    return rcomp_seq

def read_fasta(filename):
    """Read the  fasta file"""
    with open(filename,'r') as FASTA:
        for line in FASTA:
            if line[0] != '>':
                seq = line.rstrip()
    FASTA.closed
    return seq.upper()

def gen_guides(seq):
    fwd = [(seq[m.start():m.start()+20],m.start()+17,'FWD',0) for m in re.finditer(r"(?=[AGCT]{21}GG)",seq,overlapped=True)]
    rev = [(seq[m.start()+3:m.start()+23],m.start()+6,'REV',0) for m in re.finditer("CC[AGCT]{21}",seq,overlapped=True)]
    return (fwd+rev)

def read_gpp_designer(designer_output):
    score = pd.read_table(designer_output)
    guides = score[( score['sgRNA Sequence'].str.contains('TTTT' )==False) & (score['sgRNA Sequence'].str.contains('AAAA')==False)].loc[:,('sgRNA Sequence','Target Cut Length','Orientation','Combined Rank')]
    guides.replace('antisense','REV',inplace=True)
    guides.replace('sense','FWD',inplace=True)
    return (np.array(guides))

def read_flashfry(flashfry_file,threshold=50):
    score = pd.read_table(flashfry_file)
    guides = score[(score['target'].str.contains('TTTT' )==False) &
      (score['target'].str.contains('AAAA')==False ) &
      (score['dangerous_GC']=="NONE") & 
      (score['dangerous_polyT']=="NONE") &
      (score['Hsu2013']>=threshold)].sort_values(by='AggregateRankedScore_medianRank')
    guides.replace('RVS','REV',inplace=True)
    guides['seq']=guides['target'].str[:20]
    guides['start']+=np.where(guides['orientation']== 'FWD', 17, 6)
    return np.array(guides.loc[:,('seq','start','orientation','AggregateRankedScore_medianRank')])

def read_crispor(crispor_output,threshold=50):
    score = pd.read_excel(crispor_output,header=8)
    condition = (score['targetSeq'].str.contains('TTTT' )==False) &\
    (score['targetSeq'].str.contains('AAAA')==False ) &  (score['mitSpecScore']>=threshold)
    guides = pd.concat((score[condition]['targetSeq'],score[condition]['#guideId'].str.extract('(\d+)([A-Za-z]+)', expand=True)),axis=1)
    guides['start'] = guides[0].astype(int)
    guides['rank'] = guides.index+1
    guides.replace('rev','REV',inplace=True)
    guides.replace('forw','FWD',inplace=True)
    guides['start']+=np.where(guides[1]== 'FWD', -4, 5)
    guides['seq']=guides['targetSeq'].str[:20]
    return np.array(guides.loc[:,('seq','start',1,'rank')])

def onehotencoder(seq):
    '''convert to single and di-nucleotide hotencode'''
    nt= ['A','T','C','G']
    head = []
    l = len(seq)
    for k in range(l):
        for i in range(4):
            head.append(nt[i]+str(k))

    for k in range(l-1):
        for i in range(4):
            for j in range(4):
                head.append(nt[i]+nt[j]+str(k))
    head_idx = {}
    for idx,key in enumerate(head):
        head_idx[key] = idx
    encode = np.zeros(len(head_idx))
    for j in range(l):
        encode[head_idx[seq[j]+str(j)]] =1.
    for k in range(l-1):
        encode[head_idx[seq[k:k+2]+str(k)]] =1.
    return encode

def softmax(weights):
    return (np.exp(weights)/sum(np.exp(weights)))

def gen_pegpair(g1,g2,seq,len_homo,p=False,nick_start=None,nick_end=None):
    "g1 on fwd, g2 on rev. Assuming g2 returned from Flashfry and others is RC to seq"
    homology1 = reverse_complement(seq[g2[1]:g2[1]+len_homo])
    homology2 = seq[g1[1]-(len_homo):g1[1]]
    pbs1,pbs2 = reverse_complement(g1[0][17-g1[2]:17]), reverse_complement(g2[0][17-g2[2]:17])
    note = 'PolyT in homology, pass' if ('TTTT' in homology1) or ('TTTT' in homology2) else ''
    note += 'GC rich in homology,pass' if (homology1[0:5].count('C')+homology1[0:5].count('G')>=4) or (homology2[0:5].count('C')+homology2[0:5].count('G')>=4)
    if p:
        ins1,ins2 = reverse_complement(seq[g1[1]:nick_start]), seq[nick_end:g2[1]]
        note += 'PolyT in peg_insertion, pass' if ('TTTT' in ins1) or ('TTTT' in ins2) else ''
        return (g1[0],homology1+ins1+pbs1,g1[1],g2[0],homology2+ins2+pbs2,g2[1],
            note,nick_end-nick_start,seq[:nick_start]+'-'*(nick_end-nick_start)+seq[nick_end:])
    else:
        return (g1[0],homology1+pbs1,g1[1],g2[0],homology2+pbs2,g2[1],
            note,g2[1]-g1[1],seq[:g1[1]]+'-'*(g2[1]-g1[1])+seq[g2[1]:])

def peg_design_by_size(fwd,rev,seq,size_range,homology_length):
    min_dsize,max_dsize = size_range
    pairs = []
    if min_dsize<10:
        return ('Minmial deletion size is >10')
    if max_dsize>len(seq)-100:
        return ('Please reduce your deletion size')
    for a in fwd:
        for b in rev:
            if a[1] + min_dsize < b[1] and a[1] + max_dsize > b[1]:
                pairs.append(gen_pegpair(a,b,seq,homology_length))
    return np.array(pairs)

def peg_design_by_start_end(fwd,rev,seq,pos_range,homology_length,p=False):
    start,end = pos_range
    pairs = []
    for a in fwd:
        for b in rev:
            if p:
                if a[1] in range(start,start+50) and b[1] in range(end-50,end):
                    pairs.append(gen_pegpair(a,b,seq,homology_length,p=True,nick_start=start,nick_end=end))
            else:
                if a[1] in range(start-50,start+50) and b[1] in range(end-50,end+50):
                    pairs.append(gen_pegpair(a,b,seq,homology_length,p=False,nick_start=start,nick_end=end))
    if not pairs and not p:
        print ('Cannot find peg pairs match both window above, try to return pairs match either')
        for a in fwd:
            for b in rev:
                if a[1] in range(start-50,start+50) and b[1] > a[1]+20:
                    pairs.append(gen_pegpair(a,b,seq,homology_length)) 
                elif b[1] in range(end-50,end+50) and a[1] < b[1] -20:
                    pairs.append(gen_pegpair(a,b,seq,homology_length))
    return np.array(pairs)