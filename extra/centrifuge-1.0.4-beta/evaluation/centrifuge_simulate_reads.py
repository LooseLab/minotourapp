#!/usr/bin/env python

#
# Copyright 2015, Daehwan Kim <infphilo@gmail.com>
#
# This file is part of HISAT 2.
#
# HISAT 2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HISAT 2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT 2.  If not, see <http://www.gnu.org/licenses/>.
#

import sys, os, subprocess, inspect
import math, random, re
from collections import defaultdict, Counter
from argparse import ArgumentParser, FileType


"""
"""
def reverse_complement(seq):
    result = ""
    for nt in seq:
        base = nt
        if nt == 'A':
            base = 'T'
        elif nt == 'a':
            base = 't'
        elif nt == 'C':
            base = 'G'
        elif nt == 'c':
            base = 'g'
        elif nt == 'G':
            base = 'C'
        elif nt == 'g':
            base = 'c'
        elif nt == 'T':
            base = 'A'
        elif nt == 't':
            base = 'a'
        
        result = base + result
    
    return result


"""
"""
def get_genome_seq_id(genome_name):
    genome_seq_id = genome_name.split()[0]
    if len(genome_seq_id.split('|')) >= 2:
        genome_seq_id = '|'.join(genome_seq_id.split('|')[:2])
    return genome_seq_id
    

"""
Random source for sequencing errors
"""
class ErrRandomSource:
    def __init__(self, prob = 0.0, size = 1 << 20):
        self.size = size
        self.rands = []
        for i in range(self.size):
            if random.random() < prob:
                self.rands.append(1)
            else:
                self.rands.append(0)
        self.cur = 0
        
    def getRand(self):
        assert self.cur < len(self.rands)
        rand = self.rands[self.cur]
        self.cur = (self.cur + 1) % len(self.rands)
        return rand


"""
"""
def read_genomes(genomes_file, seq2taxID):
    genome_dic = {}    
    tax_id, sequence = "", ""
    for line in genomes_file:
        if line[0] == ">":
            if tax_id and sequence:
                if genome_seq_id in genome_dic:
                    genome_dic[tax_id] += sequence
                else:
                    genome_dic[tax_id] = sequence
            
            genome_name = line[1:-1]
            genome_seq_id = get_genome_seq_id(genome_name)
            assert genome_seq_id in seq2taxID
            tax_id = seq2taxID[genome_seq_id]
            sequence = ""
        else:
            sequence += line[:-1]

    if tax_id and sequence:
        if tax_id in genome_dic:
            genome_dic[tax_id] += sequence
        else:
            genome_dic[tax_id] = sequence
    
    return genome_dic


"""
"""
def read_transcript(genomes_seq, gtf_file, frag_len):
    genes = defaultdict(list)
    transcripts = {}

    # Parse valid exon lines from the GTF file into a dict by transcript_id
    for line in gtf_file:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if '#' in line:
            line = line.split('#')[0].strip()
        try:
            chrom, source, feature, left, right, score, \
                strand, frame, values = line.split('\t')
        except ValueError:
            continue
        if not chrom in genome_seq:
            continue
        
        # Zero-based offset
        left, right = int(left) - 1, int(right) - 1
        if feature != 'exon' or left >= right:
            continue

        values_dict = {}
        for attr in values.split(';')[:-1]:
            attr, _, val = attr.strip().partition(' ')
            values_dict[attr] = val.strip('"')

        if 'gene_id' not in values_dict or \
                'transcript_id' not in values_dict:
            continue

        transcript_id = values_dict['transcript_id']
        if transcript_id not in transcripts:
            transcripts[transcript_id] = [chrom, strand, [[left, right]]]
            genes[values_dict['gene_id']].append(transcript_id)
        else:
            transcripts[transcript_id][2].append([left, right])

    # Sort exons and merge where separating introns are <=5 bps
    for tran, [chr, strand, exons] in transcripts.items():
            exons.sort()
            tmp_exons = [exons[0]]
            for i in range(1, len(exons)):
                if exons[i][0] - tmp_exons[-1][1] <= 5:
                    tmp_exons[-1][1] = exons[i][1]
                else:
                    tmp_exons.append(exons[i])
            transcripts[tran] = [chr, strand, tmp_exons]

    tmp_transcripts = {}
    for tran, [chr, strand, exons] in transcripts.items():
        exon_lens = [e[1] - e[0] + 1 for e in exons]
        transcript_len = sum(exon_lens)
        if transcript_len >= frag_len:
            tmp_transcripts[tran] = [chr, strand, transcript_len, exons]

    transcripts = tmp_transcripts

    return genes, transcripts
    

"""
"""
def generate_rna_expr_profile(expr_profile_type, num_transcripts = 10000):
    # Modelling and simulating generic RNA-Seq experiments with the flux simulator
    # http://nar.oxfordjournals.org/content/suppl/2012/06/29/gks666.DC1/nar-02667-n-2011-File002.pdf
    def calc_expr(x, a):
        x, a, b = float(x), 9500.0, 9500.0
        k = -0.6
        return (x**k) * math.exp(x/a * (x/b)**2)
    
    expr_profile = [0.0] * num_transcripts
    for i in range(len(expr_profile)):
        if expr_profile_type == "flux":
            expr_profile[i] = calc_expr(i + 1, num_transcripts)
        elif expr_profile_type == "constant":
            expr_profile[i] = 1.0
        else:
            assert False

    expr_sum = sum(expr_profile)
    expr_profile = [expr_profile[i] / expr_sum for i in range(len(expr_profile))]
    assert abs(sum(expr_profile) - 1.0) < 0.001
    return expr_profile


"""
"""
def generate_dna_expr_profile(expr_profile_type, num_genomes):
    # Modelling and simulating generic RNA-Seq experiments with the flux simulator
    # http://nar.oxfordjournals.org/content/suppl/2012/06/29/gks666.DC1/nar-02667-n-2011-File002.pdf
    def calc_expr(x, a):
        x, a, b = float(x), 9500.0, 9500.0
        k = -0.6
        return (x**k) * math.exp(x/a * (x/b)**2)
    
    expr_profile = [0.0] * num_genomes
    for i in range(len(expr_profile)):
        if expr_profile_type == "flux":
            expr_profile[i] = calc_expr(i + 1, num_genomes)
        elif expr_profile_type == "constant":
            expr_profile[i] = 1.0
        else:
            assert False

    expr_sum = sum(expr_profile)
    expr_profile = [expr_profile[i] / expr_sum for i in range(len(expr_profile))]
    assert abs(sum(expr_profile) - 1.0) < 0.001
    return expr_profile


"""
"""
def getSamAlignment(dna, exons, genome_seq, trans_seq, frag_pos, read_len, err_rand_src, max_mismatch):
    # Find the genomic position for frag_pos and exon number
    tmp_frag_pos, tmp_read_len = frag_pos, read_len
    pos, cigars, cigar_descs = exons[0][0], [], []
    e_pos = 0
    prev_e = None
    for e_i in range(len(exons)):
        e = exons[e_i]
        if prev_e:
            i_len = e[0] - prev_e[1] - 1
            pos += i_len
        e_len = e[1] - e[0] + 1
        if e_len <= tmp_frag_pos:
            tmp_frag_pos -= e_len
            pos += e_len
        else:
            pos += tmp_frag_pos
            e_pos = tmp_frag_pos
            break                        
        prev_e = e

    # Define Cigar and its descriptions
    assert e_i < len(exons)
    e_len = exons[e_i][1] - exons[e_i][0] + 1
    assert e_pos < e_len
    cur_pos = pos
    match_len = 0
    prev_e = None
    mismatch, remain_trans_len = 0, len(trans_seq) - (frag_pos + read_len)
    assert remain_trans_len >= 0
    for e_i in range(e_i, len(exons)):
        e = exons[e_i]
        if prev_e:
            i_len = e[0] - prev_e[1] - 1
            cur_pos += i_len
            cigars.append(("{}N".format(i_len)))
            cigar_descs.append([])
        tmp_e_left = e_left = e[0] + e_pos
        e_pos = 0

        # Simulate mismatches due to sequencing errors
        mms = []
        for i in range(e_left, min(e[1], e_left + tmp_read_len - 1)):
            if err_rand_src.getRand() == 1:
                assert i < len(genome_seq)
                err_base = "A"
                rand = random.randint(0, 2)
                if genome_seq[i] == "A":
                    err_base = "GCT"[rand]
                elif genome_seq[i] == "C":
                    err_base = "AGT"[rand]
                elif genome_seq[i] == "G":
                    err_base = "ACT"[rand]
                else:
                    err_base = "ACG"[rand]                    
                mms.append(["", "single", i, err_base])

        tmp_diffs = mms
        def diff_sort(a , b):
            return a[2] - b[2]

        tmp_diffs = sorted(tmp_diffs, cmp=diff_sort)
        diffs = []
        if len(tmp_diffs) > 0:
            diffs = tmp_diffs[:1]
            for diff in tmp_diffs[1:]:
                _, tmp_type, tmp_pos, tmp_data = diff
                _, prev_type, prev_pos, prev_data = diffs[-1]
                if prev_type == "deletion":
                    prev_pos += prev_data
                if tmp_pos <= prev_pos:
                    continue
                diffs.append(diff)

        cigar_descs.append([])
        prev_diff = None
        for diff in diffs:
            diff_id, diff_type, diff_pos, diff_data = diff
            if prev_diff:
                prev_diff_id, prev_diff_type, prev_diff_pos, prev_diff_data = prev_diff
                if prev_diff_type == "deletion":
                    prev_diff_pos += prev_diff_data
                assert prev_diff_pos < diff_pos
            diff_pos2 = diff_pos
            if diff_type == "deletion":
                diff_pos2 += diff_data
            if e_left + tmp_read_len - 1 < diff_pos2 or e[1] < diff_pos2:
                break            
            if diff_type == "single":
                if diff_id == "" and mismatch >= max_mismatch:
                    continue                
                cigar_descs[-1].append([diff_pos - tmp_e_left, diff_data, diff_id])
                tmp_e_left = diff_pos + 1
                if diff_id == "":
                    mismatch += 1
            elif diff_type == "deletion":
                if len(cigars) <= 0:
                    continue
                del_len = diff_data
                if remain_trans_len < del_len:
                    continue
                remain_trans_len -= del_len
                if diff_pos - e_left > 0:
                    cigars.append("{}M".format(diff_pos - e_left))
                    cigar_descs[-1].append([diff_pos - tmp_e_left, "", ""])
                    cigar_descs.append([])
                cigars.append("{}D".format(del_len))
                cigar_descs[-1].append([0, del_len, diff_id])
                cigar_descs.append([])
                tmp_read_len -= (diff_pos - e_left)
                e_left = tmp_e_left = diff_pos + del_len
            elif diff_type == "insertion":
                if len(cigars) > 0:
                    ins_len = len(diff_data)
                    if e_left + tmp_read_len - 1 < diff_pos + ins_len:
                        break
                    if diff_pos - e_left > 0:
                        cigars.append("{}M".format(diff_pos - e_left))
                        cigar_descs[-1].append([diff_pos - tmp_e_left, "", ""])
                        cigar_descs.append([])
                    cigars.append("{}I".format(ins_len))
                    cigar_descs[-1].append([0, diff_data, diff_id])
                    cigar_descs.append([])
                    tmp_read_len -= (diff_pos - e_left)
                    tmp_read_len -= ins_len
                    e_left = tmp_e_left = diff_pos
            else:
                assert False
            prev_diff = diff

        e_right = min(e[1], e_left + tmp_read_len - 1)
        e_len = e_right - e_left + 1
        remain_e_len = e_right - tmp_e_left + 1
        if remain_e_len > 0:
            cigar_descs[-1].append([remain_e_len, "", ""])
        if e_len < tmp_read_len:
            tmp_read_len -= e_len
            cigars.append(("{}M".format(e_len)))
        else:
            assert e_len == tmp_read_len
            cigars.append(("{}M".format(tmp_read_len)))
            tmp_read_len = 0
            break
        prev_e = e

    # Define MD, XM, NM, Zs, read_seq
    MD, XM, NM, Zs, read_seq = "", 0, 0, "", ""
    assert len(cigars) == len(cigar_descs)
    MD_match_len, Zs_match_len = 0, 0
    cur_trans_pos = frag_pos
    for c in range(len(cigars)):
        cigar = cigars[c]
        cigar_len, cigar_op = int(cigar[:-1]), cigar[-1]
        cigar_desc = cigar_descs[c]
        if cigar_op == 'N':
            continue
        if cigar_op == 'M':
            for add_match_len, alt_base, snp_id in cigar_desc:
                MD_match_len += add_match_len
                Zs_match_len += add_match_len
                assert cur_trans_pos + add_match_len <= len(trans_seq)
                read_seq += trans_seq[cur_trans_pos:cur_trans_pos+add_match_len]
                cur_trans_pos += add_match_len
                if alt_base != "":
                    if MD_match_len > 0:
                        MD += ("{}".format(MD_match_len))
                        MD_match_len = 0
                    MD += trans_seq[cur_trans_pos]
                    if snp_id != "":
                        if Zs != "":
                            Zs += ","
                        Zs += ("{}|S|{}".format(Zs_match_len, snp_id))
                        Zs_match_len = 0
                    else:
                        Zs_match_len += 1
                    if snp_id == "":
                        XM += 1
                        NM += 1
                    read_seq += alt_base
                    cur_trans_pos += 1
        elif cigar_op == 'D':
            assert len(cigar_desc) == 1
            add_match_len, del_len, snp_id = cigar_desc[0]
            MD_match_len += add_match_len
            Zs_match_len += add_match_len
            if MD_match_len > 0:
                MD += ("{}".format(MD_match_len))
                MD_match_len = 0
            MD += ("^{}".format(trans_seq[cur_trans_pos:cur_trans_pos+cigar_len]))
            read_seq += trans_seq[cur_trans_pos:cur_trans_pos+add_match_len]
            if Zs != "":
                Zs += ","
            Zs += ("{}|D|{}".format(Zs_match_len, cigar_desc[0][-1]))
            Zs_match_len = 0
            cur_trans_pos += cigar_len
        elif cigar_op == 'I':
            assert len(cigar_desc) == 1
            add_match_len, ins_seq, snp_id = cigar_desc[0]
            ins_len = len(ins_seq)
            MD_match_len += add_match_len
            Zs_match_len += add_match_len
            read_seq += trans_seq[cur_trans_pos:cur_trans_pos+add_match_len]
            read_seq += ins_seq
            if Zs != "":
                Zs += ","
            Zs += ("{}|I|{}".format(Zs_match_len, cigar_desc[0][-1]))
            Zs_match_len = 0
        else:
            assert False

    if MD_match_len > 0:
        MD += ("{}".format(MD_match_len))

    if len(read_seq) != read_len:
        print >> sys.stderr, "read length differs:", len(read_seq), "vs.", read_len
        print >> sys.stderr, pos, "".join(cigars), cigar_descs, MD, XM, NM, Zs
        assert False

    return pos, cigars, cigar_descs, MD, XM, NM, Zs, read_seq


"""
"""
cigar_re = re.compile('\d+\w')
def samRepOk(genome_seq, read_seq, chr, pos, cigar, XM, NM, MD, Zs, max_mismatch):
    assert chr in genome_seq
    chr_seq = genome_seq[chr]
    assert pos < len(chr_seq)

    # Calculate XM and NM based on Cigar and Zs
    cigars = cigar_re.findall(cigar)
    cigars = [[int(cigars[i][:-1]), cigars[i][-1]] for i in range(len(cigars))]
    ref_pos, read_pos = pos, 0
    ann_ref_seq, ann_ref_rel, ann_read_seq, ann_read_rel = [], [], [], []
    for i in range(len(cigars)):
        cigar_len, cigar_op = cigars[i]
        if cigar_op == "M":
            partial_ref_seq = chr_seq[ref_pos:ref_pos+cigar_len]
            partial_read_seq = read_seq[read_pos:read_pos+cigar_len]
            assert len(partial_ref_seq) == len(partial_read_seq)
            ann_ref_seq += list(partial_ref_seq)
            ann_read_seq += list(partial_read_seq)
            for j in range(len(partial_ref_seq)):
                if partial_ref_seq[j] == partial_read_seq[j]:
                    ann_ref_rel.append("=")
                    ann_read_rel.append("=")
                else:
                    ann_ref_rel.append("X")
                    ann_read_rel.append("X")
            ref_pos += cigar_len
            read_pos += cigar_len
        elif cigar_op == "D":
            partial_ref_seq = chr_seq[ref_pos:ref_pos+cigar_len]
            ann_ref_rel += list(partial_ref_seq)
            ann_ref_seq += list(partial_ref_seq)
            ann_read_rel += (["-"] * cigar_len)
            ann_read_seq += (["-"] * cigar_len)
            ref_pos += cigar_len
        elif cigar_op == "I":
            partial_read_seq = read_seq[read_pos:read_pos+cigar_len]
            ann_ref_rel += (["-"] * cigar_len)
            ann_ref_seq += (["-"] * cigar_len)
            ann_read_rel += list(partial_read_seq)
            ann_read_seq += list(partial_read_seq) 
            read_pos += cigar_len
        elif cigar_op == "N":
            ref_pos += cigar_len
        else:
            assert False
    
    assert len(ann_ref_seq) == len(ann_read_seq)
    assert len(ann_ref_seq) == len(ann_ref_rel)
    assert len(ann_ref_seq) == len(ann_read_rel)
    ann_Zs_seq = ["0" for i in range(len(ann_ref_seq))]

    Zss, Zs_i, snp_pos_add = [], 0, 0
    if Zs != "":
        Zss = Zs.split(',')
        Zss = [zs.split('|') for zs in Zss]

    ann_read_pos = 0
    for zs in Zss:
        zs_pos, zs_type, zs_id = zs
        zs_pos = int(zs_pos)
        for i in range(zs_pos):
            while ann_read_rel[ann_read_pos] == '-':
                ann_read_pos += 1
            ann_read_pos += 1
        if zs_type == "S":
            ann_Zs_seq[ann_read_pos] = "1"
            ann_read_pos += 1
        elif zs_type == "D":
            while ann_read_rel[ann_read_pos] == '-':
                ann_Zs_seq[ann_read_pos] = "1"
                ann_read_pos += 1
        elif zs_type == "I":
            while ann_ref_rel[ann_read_pos] == '-':
                ann_Zs_seq[ann_read_pos] = "1"
                ann_read_pos += 1
        else:
            assert False

    tMD, tXM, tNM = "", 0, 0
    match_len = 0
    i = 0
    while i < len(ann_ref_seq):
        if ann_ref_rel[i] == "=":
            assert ann_read_rel[i] == "="
            match_len += 1
            i += 1
            continue
        assert ann_read_rel[i] != "="
        if ann_ref_rel[i] == "X" and ann_read_rel[i] == "X":
            if match_len > 0:
                tMD += ("{}".format(match_len))
                match_len = 0
            tMD += ann_ref_seq[i]
            if ann_Zs_seq[i] == "0":
                tXM += 1
                tNM += 1
            i += 1
        else:
            assert ann_ref_rel[i] == "-" or ann_read_rel[i] == "-"
            if ann_ref_rel[i] == '-':
                while ann_ref_rel[i] == '-':
                    if ann_Zs_seq[i] == "0":
                        tNM += 1
                    i += 1
            else:
                assert ann_read_rel[i] == '-'
                del_seq = ""
                while  ann_read_rel[i] == '-':
                    del_seq += ann_ref_seq[i]
                    if ann_Zs_seq[i] == "0":
                        tNM += 1
                    i += 1
                if match_len > 0:
                    tMD += ("{}".format(match_len))
                    match_len = 0
                tMD += ("^{}".format(del_seq))

    if match_len > 0:
        tMD += ("{}".format(match_len))

    if tMD != MD or tXM != XM or tNM != NM or XM > max_mismatch or XM != NM:
        print >> sys.stderr, chr, pos, cigar, MD, XM, NM, Zs
        print >> sys.stderr, tMD, tXM, tNM
        assert False
        
        
"""
"""
def simulate_reads(index_fname, base_fname, \
                       dna, paired_end, read_len, frag_len, \
                       num_frag, expr_profile_type, error_rate, max_mismatch, \
                       random_seed, sanity_check, verbose):
    random.seed(random_seed)
    
    # Current script directory
    curr_script = os.path.realpath(inspect.getsourcefile(simulate_reads))
    ex_path = os.path.dirname(curr_script)
    centrifuge_inspect = os.path.join(ex_path, "../centrifuge-inspect")

    err_rand_src = ErrRandomSource(error_rate / 100.0)
    
    if read_len > frag_len:
        frag_len = read_len

    # Read taxonomic IDs
    seq2texID = {}
    tax_cmd = [centrifuge_inspect,
               "--conversion-table",
               index_fname]
    tax_proc = subprocess.Popen(tax_cmd, stdout=subprocess.PIPE)
    for line in tax_proc.stdout:
        seq_id, tax_id = line.strip().split()
        seq2texID[seq_id] = tax_id

    # Read names
    names = {}
    name_cmd = [centrifuge_inspect,
                "--name-table",
                index_fname]
    name_proc = subprocess.Popen(name_cmd, stdout=subprocess.PIPE)
    for line in name_proc.stdout:
        tax_id, name = line.strip().split('\t')
        names[tax_id] = name

    # Genome sizes
    sizes = {}
    size_cmd = [centrifuge_inspect,
                "--size-table",
                index_fname]
    size_proc = subprocess.Popen(size_cmd, stdout=subprocess.PIPE)
    for line in size_proc.stdout:
        tax_id, size = line.strip().split('\t')
        sizes[tax_id] = int(size)

    # Read genome sequences into memory
    genomes_fname = index_fname + ".fa"
    if not os.path.exists(genomes_fname):
        print >> sys.stderr, "Extracting genomes from Centrifuge index to %s, which may take a few hours ..."  % (genomes_fname)
        extract_cmd = [centrifuge_inspect,
                       index_fname]
        extract_proc = subprocess.Popen(extract_cmd, stdout=open(genomes_fname, 'w'))
        extract_proc.communicate()
    genome_seqs = read_genomes(open(genomes_fname), seq2texID)

    if dna:
        genes, transcripts = {}, {}
    else:
        genes, transcripts = read_transcript(genome_seqs, gtf_file, frag_len)
        
    if sanity_check:
        sanity_check_input(genomes_seq, genes, transcripts, frag_len)

    if dna:
        expr_profile = generate_dna_expr_profile(expr_profile_type, min(len(genome_seqs), 100))
    else:
        num_transcripts = min(len(transcripts), 10000)
        expr_profile = generate_rna_expr_profile(expr_profile_type, num_transcripts)

    expr_profile = [int(expr_profile[i] * num_frag) for i in range(len(expr_profile))]
    assert num_frag >= sum(expr_profile)
    while sum(expr_profile) < num_frag:
        for i in range(min(num_frag - sum(expr_profile), len(expr_profile))):
            expr_profile[i] += 1
    assert num_frag == sum(expr_profile)

    if dna:
        genome_ids = genome_seqs.keys()
    else:
        transcript_ids = transcripts.keys()
        random.shuffle(transcript_ids)
        assert len(transcript_ids) >= len(expr_profile)

    # Truth table
    truth_file = open(base_fname + ".truth", "w")
    print >> truth_file, "taxID\tgenomeLen\tnumReads\tabundance\tname"
    truth_list = []
    normalized_sum = 0.0
    debug_num_frag = 0
    for t in range(len(expr_profile)):
        t_num_frags = expr_profile[t]
        if dna:
            tax_id = genome_ids[t]
        else:
            transcript_id = transcript_ids[t]
            chr, strand, transcript_len, exons = transcripts[transcript_id]
        assert tax_id in genome_seqs and tax_id in sizes
        genome_len = sizes[tax_id]
        raw_abundance = float(t_num_frags)/num_frag
        normalized_sum += (raw_abundance / genome_len)
        truth_list.append([tax_id, genome_len, t_num_frags, raw_abundance])
        debug_num_frag += t_num_frags
    assert debug_num_frag == num_frag
    for truth in truth_list:
        tax_id, genome_len, t_num_frags, raw_abundance = truth
        can_tax_id = tax_id
        if '.' in can_tax_id:
            can_tax_id = can_tax_id.split('.')[0]
        name = "N/A"        
        if can_tax_id in names:
            name = names[can_tax_id]
        abundance = raw_abundance / genome_len / normalized_sum
        print >> truth_file, "{}\t{}\t{}\t{:.6}\t{}".format(tax_id, genome_len, t_num_frags, abundance, name)
    truth_file.close()

    # Sequence Classification Map (SCM) - something I made up ;-)
    scm_file = open(base_fname + ".scm", "w")

    # Write SCM header
    print >> scm_file, "@HD\tVN:1.0\tSO:unsorted"
    for tax_id in genome_seqs.keys():
        name = ""
        if tax_id in names:
            name = names[tax_id]
        print >> scm_file, "@SQ\tTID:%s\tSN:%s\tLN:%d" % (tax_id, name, len(genome_seqs[tax_id]))

    read_file = open(base_fname + "_1.fa", "w")
    if paired_end:
        read2_file = open(base_fname + "_2.fa", "w")

    cur_read_id = 1
    for t in range(len(expr_profile)):
        t_num_frags = expr_profile[t]
        if dna:
            tax_id = genome_ids[t]
            print >> sys.stderr, "TaxID: %s, num fragments: %d" % (tax_id, t_num_frags)
        else:
            transcript_id = transcript_ids[t]
            chr, strand, transcript_len, exons = transcripts[transcript_id]
            print >> sys.stderr, transcript_id, t_num_frags

        genome_seq = genome_seqs[tax_id]
        genome_len = len(genome_seq)
        if dna:
            t_seq = genome_seq
            exons = [[0, genome_len - 1]]
        else:            
            t_seq = ""
            for e in exons:
                assert e[0] < e[1]
                t_seq += genome_seq[e[0]:e[1]+1]
            assert len(t_seq) == transcript_len
            
        for f in range(t_num_frags):
            if dna:
                while True:
                    frag_pos = random.randint(0, genome_len - frag_len)
                    if 'N' not in genome_seq[frag_pos:frag_pos + frag_len]:
                        break
            else:
                frag_pos = random.randint(0, transcript_len - frag_len)

            pos, cigars, cigar_descs, MD, XM, NM, Zs, read_seq = getSamAlignment(dna, exons, genome_seq, t_seq, frag_pos, read_len, err_rand_src, max_mismatch)
            pos2, cigars2, cigar2_descs, MD2, XM2, NM2, Zs2, read2_seq = getSamAlignment(dna, exons, genome_seq, t_seq, frag_pos+frag_len-read_len, read_len, err_rand_src, max_mismatch)
            cigar_str, cigar2_str = "".join(cigars), "".join(cigars2)
            if sanity_check:
                samRepOk(genome_seq, read_seq, chr, pos, cigar_str, XM, NM, MD, Zs, max_mismatch)
                samRepOk(genome_seq, read2_seq, chr, pos2, cigar2_str, XM2, NM2, MD2, Zs2, max_mismatch)

            if Zs != "":
                Zs = ("\tZs:Z:{}".format(Zs))
            if Zs2 != "":
                Zs2 = ("\tZs:Z:{}".format(Zs2))
            
            if dna:
                XS, TI = "", ""                
            else:
                XS = "\tXS:A:{}".format(strand)
                TI = "\tTI:Z:{}".format(transcript_id)                

            print >> read_file, ">{}".format(cur_read_id)
            print >> read_file, read_seq
            output = "{}\t{}\t{}\t{}\tNM:i:{}\tMD:Z:{}".format(cur_read_id, tax_id, pos + 1, cigar_str, NM, MD)
            if paired_end:
                print >> read2_file, ">{}".format(cur_read_id)
                print >> read2_file, reverse_complement(read2_seq)
                output += "\t{}\t{}\tNM2:i:{}\tMD2:Z:{}".format(pos2 + 1, cigar2_str, NM2, MD2)
            print >> scm_file, output
                
            cur_read_id += 1
            
    scm_file.close()
    read_file.close()
    if paired_end:
        read2_file.close()


if __name__ == '__main__':
    parser = ArgumentParser(
        description='Simulate reads from Centrifuge index')
    parser.add_argument('index_fname',
                        nargs='?',
                        type=str,
                        help='Centrifuge index')
    """
    parser.add_argument('gtf_file',
                        nargs='?',
                        type=FileType('r'),
                        help='input GTF file')
    """
    parser.add_argument('base_fname',
                        nargs='?',
                        type=str,
                        help='output base filename')
    parser.add_argument('--rna',
                        dest='dna',
                        action='store_false',
                        default=True,
                        help='RNA-seq reads (default: DNA-seq reads)')
    parser.add_argument('--single-end',
                        dest='paired_end',
                        action='store_false',
                        default=True,
                        help='single-end reads (default: paired-end reads)')
    parser.add_argument('-r', '--read-length',
                        dest='read_len',
                        action='store',
                        type=int,
                        default=100,
                        help='read length (default: 100)')
    parser.add_argument('-f', '--fragment-length',
                        dest='frag_len',
                        action='store',
                        type=int,
                        default=250,
                        help='fragment length (default: 250)')
    parser.add_argument('-n', '--num-fragment',
                        dest='num_frag',
                        action='store',
                        type=int,
                        default=1000000,
                        help='number of fragments (default: 1000000)')
    parser.add_argument('-e', '--expr-profile',
                        dest='expr_profile',
                        action='store',
                        type=str,
                        default='flux',
                        help='expression profile: flux or constant (default: flux)')
    parser.add_argument('--error-rate',
                        dest='error_rate',
                        action='store',
                        type=float,
                        default=0.0,
                        help='per-base sequencing error rate (%%) (default: 0.0)')
    parser.add_argument('--max-mismatch',
                        dest='max_mismatch',
                        action='store',
                        type=int,
                        default=3,
                        help='max mismatches due to sequencing errors (default: 3)')
    parser.add_argument('--random-seed',
                        dest='random_seed',
                        action='store',
                        type=int,
                        default=0,
                        help='random seeding value (default: 0)')
    parser.add_argument('--sanity-check',
                        dest='sanity_check',
                        action='store_true',
                        help='sanity check')
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')
    parser.add_argument('--version', 
                        action='version',
                        version='%(prog)s 2.0.0-alpha')
    args = parser.parse_args()
    if not args.index_fname:
        parser.print_help()
        exit(1)
    if not args.dna:
        print >> sys.stderr, "Error: --rna is not implemented."
        exit(1)
    # if args.dna:
    #    args.expr_profile = "constant"
    simulate_reads(args.index_fname, args.base_fname, \
                       args.dna, args.paired_end, args.read_len, args.frag_len, \
                       args.num_frag, args.expr_profile, args.error_rate, args.max_mismatch, \
                       args.random_seed, args.sanity_check, args.verbose)
