#! /Users/karol/miniconda3/bin/python
#
#  SIMPLE MODE:
#  convert vcf to table
#
#  REGENOTYPE MODE:
#  Add coverage information to after merging two VCF's
#  OUTPUT: a new recalculated VCF on stdout
#          a table as named argument
#  Karol

import sys, getopt, re
import pandas as pd
import numpy as np
from numpy import array
#import pysam
import argparse
import vcf


verbose = False
parent_parser = argparse.ArgumentParser()

subparsers = parent_parser.add_subparsers(title = "mode", help = "sub-command help",dest = 'mode')
subparsers.required = True

parser_simple = subparsers.add_parser("simple",help ="run simple vcf to table conversion")
parser_simple.add_argument("-o","--out", help = "file to output table", required = True)
parser_simple.add_argument("-i","--input_vcf", help = "input vcf to be parsed", required = True)
parser_simple.add_argument('-debug',action = "store_true")
parser_simple.add_argument('-tab',action = "store_true")
parser_simple.add_argument('--build', help = "genome build")

parser_reGeno = subparsers.add_parser("reGenotype", help ="convert vcf to table and add missing fields from additional input")
parser_reGeno.add_argument("-i","--input_vcf", help = "input vcf to be parsed", required = True)
parser_reGeno.add_argument("-o","--out", help = "file to output table", required = True)
parser_reGeno.add_argument("--t1_bam", help = "tumor 1 bam file", required = True)
parser_reGeno.add_argument("--t2_bam", help = "tumor 2 bam file"  , required = True)
parser_reGeno.add_argument("--normal", help = "ID of normal sample"  , required = True)
parser_reGeno.add_argument('-debug',action = "store_true")
parser_reGeno.add_argument('-tab',action = "store_true")
parser_reGeno.add_argument('--build', help = "genome build")

# parser_reGeno.add_argument("-m","--mutect_vcf", help = "mutect vcf", required = True)
# parser_reGeno.add_argument("-v","--vardict_vcf",help = "vardict vcf", required = True)
# parser_reGeno.add_argument("-s","--sample", help = "sample name", required = True)


args = parent_parser.parse_args()

build = 'hs37d5'
build = '.'
if args.build :
    build = args.build

def eprint(str):
    print(str, file = sys.stderr)


def load_vcf(file_name, sample, vc):
    vcf_reader = vcf.Reader(filename=file_name)
    var_dtype  = [('CHROM','S50'),('POS',int),('REF','S50'),('ALT','S50'),('tum_ref_ad',int),('tum_alt_ad',int),('norm_ref_ad',int),('norm_alt_ad',int)]
    var_df     = np.empty([0,8],dtype=var_dtype)
    tumor_sample  = sample +"_tumor." + vc
    normal_sample = sample +"_normal." + vc

    for record in vcf_reader:
        try:
            tum_ref_ad,  tum_alt_ad   = record.genotype(tumor_sample)["AD"]
            norm_ref_ad, norm_alt_ad  = record.genotype(normal_sample)["AD"]
            eprint("CHROM:{}, POS:{}, REF:{}, ALT:{}, tum_ref_ad:{}, tum_alt_ad:{}, norm_ref_ad:{}, norm_alt_ad:{}".format(record.CHROM, record.POS, record.REF, record.ALT[0], tum_ref_ad, tum_alt_ad, norm_ref_ad, norm_alt_ad))
            tmp = np.array((record.CHROM, record.POS, record.REF, record.ALT[0], tum_ref_ad, tum_alt_ad, norm_ref_ad, norm_alt_ad), dtype=var_dtype)
            var_df = np.append(var_df,tmp)
        except KeyError:
            print("Error while parsing genotype fields: {}".format(KeyError))

    return pd.DataFrame(var_df)

def extract_gt(gth, field, df, chrom, pos, alt, sm):
    gtf = re.split(r':',field)
    depth = gtf[gth.index('DP')]
    freq  = gtf[gth.index('AF')]
    if ('AD' in gth):
        ad    = re.split(r',',gtf[gth.index('AD')])
    else:
        ad    = df.loc[(    df["CHROM"] == chrom.encode() ) & (    df["POS"] ==        int(pos)) & (    df["ALT"] == alt.encode() )][[sm + '_ref_ad', sm + '_alt_ad']]
        if(ad.empty):
            ad = ['.','.']
        else:
            ad = ad.values.tolist()[0]
    return(depth,freq,ad)

def extract_gt_strelka(gth, field, ref, alt):
    gtf = re.split(r':',field)

    if(len(ref) > 1 or len(alt) > 1):
        ad = ['.',re.split(r',',gtf[gth.index('TIR')])[0]]
    else:
        ad = [re.split(r',',gtf[gth.index(ref + 'U')])[0],re.split(r',',gtf[gth.index(alt +'U')])[0]]
    depth = gtf[gth.index('DP')]
    freq  = gtf[gth.index('AF')]
    return(depth,freq,ad)

def extract_gt_vardict(gth, field):
    gtf = re.split(r':',field)
    #eprint('Chrom: {}, Pos: {}, Alt: {}'.format(chrom, pos, alt))
    depth = gtf[gth.index('DP')]
    freq  = gtf[gth.index('AF')]
    ad    = re.split(r',',gtf[gth.index('AD')])
    return(depth,freq,ad)

consequence_order = \
      { "transcript_ablation":1,
        "splice_acceptor_variant":2,
        "splice_donor_variant":3,
        "stop_gained":4,
        "frameshift_variant":5,
        "stop_lost":6,
        "start_lost":7,
        "transcript_amplification":8,
        "inframe_insertion":9,
        "inframe_deletion":10,
        "missense_variant":11,
        "protein_altering_variant":12,
        "splice_region_variant":13,
        "incomplete_terminal_codon_variant":14,
        "start_retained_variant":15,
        "stop_retained_variant":16,
        "synonymous_variant":17,
        "coding_sequence_variant":18,
        "mature_miRNA_variant":19,
        "5_prime_UTR_variant":20,
        "3_prime_UTR_variant":21,
        "non_coding_transcript_exon_variant":22,
        "intron_variant":23,
        "NMD_transcript_variant":24,
        "non_coding_transcript_variant":25,
        "upstream_gene_variant":26,
        "downstream_gene_variant":27,
        "TFBS_ablation":28,
        "TFBS_amplification":29,
        "TF_binding_site_variant":30,
        "regulatory_region_ablation":31,
        "regulatory_region_amplification":32,
        "feature_elongation":33,
        "regulatory_region_variant":34,
        "feature_truncation":35,
        "intergenic_variant":36 }

def source_value(source):
    ret = 3
    if source == "RefSeq":
        ret = 2
    if source == "Ensembl":
        ret = 1
    return(ret)

def sort_annotation(ann, found_anno):
    # assuming the consequences are ordered when concatenated with &
    ann['Variant_Classification'] = ann['Consequence'].apply(lambda x: re.split('&',x)[0])
    ann["ord_cs"] = ann['Variant_Classification'].apply(lambda x : consequence_order[x])
    if "SOURCE" in found_anno:
        ann["ord_source"] = ann["SOURCE"].apply(lambda x : source_value(x))
    else:
        ann["ord_source"] = 0
    if "CANONICAL" in found_anno:
        ann["ord_canon"]  = ann["CANONICAL"].apply(lambda x: 1 if (x=="YES") else 2)
    else:
        ann["ord_canon"] = 0
    ann = ann.sort_values(["ord_cs","ord_source","ord_canon"],ascending = [True,True,True])

    #with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    #    print(ann)
    return(ann)


def squish(ann):
    #ann['HGVSp']     = ann["HGVSp"].apply(lambda x: "" if (re.split(":",x) == ['']) else re.split(":",x)[1])
    #ann['HGVSc'] = ann["HGVSc"].apply(lambda x: "" if (re.split(":",x) == [''] ) else re.split(":",x)[1])

    x = ann[['Feature','HGVSc','HGVSp']].to_string(header=False, index=False, index_names=False).split("\n")
    vals = [';'.join(ele.split()) for ele in x]
    #print("===============",file = sys.stderr)
    #print("|".join(vals), file=sys.stderr)
    #print("===============",file = sys.stderr)
    return('|'.join(vals))


f = open(args.input_vcf, 'r')
lines = f.readlines()


if args.mode not in {'simple','reGenotype'}:
    sample = args.sample

    mut_df = load_vcf(args.mutect_vcf, sample, 'mutect2')
    var_df = load_vcf(args.vardict_vcf, sample, 'vardict')

ann_header_s = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|REFSEQ_MATCH|SOURCE|GIVEN_REF|USED_REF|BAM_EDIT|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE"

ann_header = ""
found_anno = ""

if args.mode == 'simple':
    # !!! new column name filter was added
    new_variant = { 'comment':'', 'gene_symbol':'.', 'genome_build':build,
                'Chromosome':'.',    'Start_Position':'.', 'vcf_filter':'.',
                'Variant_Classification':'.', 'Variant_Type':'.',
                'Reference_Allele':'.',       'Tumor_Seq_Allele2':'.',
                'HGVSc':'.', 'HGVSp':'.', 'Mutation_description':'.', 'NM_transcriptID':'.' ,'feature_type':'.',
                'Transcript_ID':'.', 'Exon_Number':'.',
                'tumor_DP':'.', 'tumor_AD_ref':'.', 'tumor_AD_alt':'.','tumor_AF':'.', 'strand_bias':'.',
                'alt_transcripts':'.', 'Codons':'.', 'Existing_variation':'.', 'BIOTYPE':'.',
                'CANONICAL':'.',
                'SIFT':'.', 'PolyPhen':'.', 'clinvar':'.',
                'IMPACT':'.',
                'af':'.',
                'eur_af':'.',    'gnomad_af':'.',               'gnomad_nfe_af':'.',
                'max_af':'.',    'max_af_pops':'.',        'pubmed':'.',
                }

if args.mode == 'reGenotype':

    t1_samfile = pysam.Samfile(args.t1_bam, "rb")
    t1_sample_id = t1_samfile.header.get("RG")[0]["SM"]

    t2_samfile = pysam.Samfile(args.t2_bam, "rb")
    t2_sample_id = t2_samfile.header.get("RG")[0]["SM"]

    print("t1 Sample id:{} ".format(t1_sample_id),file=sys.stderr)
    print("t2 Sample id:{} ".format(t2_sample_id),file=sys.stderr)

    n_id = args.normal
    new_variant = { 'comment':'', 'gene_symbol':'.', 'genome_build':build,
                    'Chromosome':'.',    'Start_Position':'.',
                    'Variant_Classification':'.', 'Variant_Type':'.',
                    'Reference_Allele':'.',       'Tumor_Seq_Allele2':'.',
                    'HGVSc':'.', 'HGVSp':'.', 'feature_type':'.',
                    'Transcript_ID':'.', 'Exon_Number':'.',
                    '{}_DP'.format(n_id):'.', '{}_AD_ref'.format(n_id):'.', '{}_AD_alt'.format(n_id):'.', '{}_AF'.format(n_id):'.',
                    '{}_DP'.format(t1_sample_id):'.', '{}_AD_ref'.format(t1_sample_id):'.', '{}_AD_alt'.format(t1_sample_id):'.', '{}_AF'.format(t1_sample_id):'.',
                    '{}_DP'.format(t2_sample_id):'.', '{}_AD_ref'.format(t2_sample_id):'.', '{}_AD_alt'.format(t2_sample_id):'.', '{}_AF'.format(t2_sample_id):'.',
                    'TP2-TP1.AF':'.', 'strand_bias':'.',
                    'alt_transcripts':'.', 'Codons':'.', 'Existing_variation':'.', 'BIOTYPE':'.',
                    'CANONICAL':'.',
                    'SIFT':'.', 'PolyPhen':'.', 'clinvar':'.',
                    'IMPACT':'.',
                    'af':'.',
                    'eur_af':'.',    'gnomad_af':'.',               'gnomad_nfe_af':'.',
                    'max_af':'.',    'max_af_pops':'.',        'pubmed':'.',      'set':'.',
                    }

idx =  list(new_variant.keys())
df = pd.DataFrame(columns=idx)

add_variant = new_variant.copy()

for line in lines:
    line = line.strip()
    if(line.startswith("#")):
        if( "ID=CSQ" in line ):
            m = re.search("Format: (.+?)\">",line)
            if m:
                ann_header_s = m.group(1);
            else:
                eprint("WARNING! Unable to extract VEP annotatin found!");
                ann_header_s = ""

            ann_header = re.split("\|",ann_header_s)
            expected_anno = {"SYMBOL", "HGVSc", "HGVSp", "IMPACT", "BIOTYPE", "CANONICAL", "Feature_type", "Feature", \
                             "EXON", "Codons", "Consequence", "SIFT", "PolyPhen", "PUBMED", "Existing_variation", 'MAX_AF_POPS', \
                             'VARIANT_CLASS', 'AF', 'EUR_AF', 'MAX_AF', 'gnomAD_AF', 'gnomAD_NFE_AF', "SOURCE", "CANONICAL", "CLIN_SIG" }
            found_anno = expected_anno.copy()
            for anno in expected_anno:
                if anno not in ann_header:
                    found_anno.remove(anno)
                    eprint("WARNING! '{}' missing from VEP annotation!".format(anno))

        if(line.startswith("#CHROM")):
            if(ann_header == ""):
                eprint("WARNING! No VEP annotatin found!");

            fields = re.split(r'\t',line)
            fields_len = len(fields)

            #reGenotype mode
            if(args.mode == 'reGenotype'):
                n_pos   = fields.index(n_id)
                t1_pos  = fields.index(t1_sample_id)
                t2_pos  = fields.index(t2_sample_id)

            if(args.mode != 'simple'):
                print(line)
            continue
        if(args.mode != 'simple'):
            print(line)
        continue

    #empty new variant dictionary
    for key in add_variant:
        add_variant[key] = '.'
    add_variant['genome_build'] = build
    add_variant['comment'] = ''
    fields = re.split(r'\t+',line)

    INFO_field = fields[7] + ";"


    chr      = fields[0]
    pos      = fields[1]
    #print("chr: {}, pos: {}.".format(chr,pos), file = sys.stderr)
    ref_base = fields[3]
    alt_base = fields[4]
    # !!! add vcf filter info
    vcf_filt = fields[6]
    add_variant['Chromosome'] = chr
    add_variant['Start_Position']   = pos
    add_variant['Reference_Allele']   = ref_base
    add_variant['Tumor_Seq_Allele2']   = alt_base
    add_variant['vcf_filter'] = vcf_filt

    # reGenotype
    if args.mode != 'simple':
        set_from = re.search("set=(.+?);",INFO_field).group(1)
        add_variant['set'] = set_from

    gth = re.split(r':',fields[8])

    missing = 0
    present = 0

    dp_pos = gth.index("DP")
    af_pos = -1
    ao_pos = -1
    if "AF" in gth:
        af_pos = gth.index("AF")
    ad_pos = -1
    if "AD" in gth:
        ad_pos = gth.index("AD")

    if (missing > 0):
        gtf = fields[missing]
    if args.mode == 'reGenotype':

        if set_from in t1_sample_id:
            missing = t2_pos
            samfile = t2_samfile
            present = t1_pos
        if set_from in t2_sample_id:
            missing = t1_pos
            samfile = t1_samfile
            present = t2_pos


        if (missing > 0) :
            gtf = ["."]*len(gth)
            ref = 0
            alt = 0
            depth = 0
            insertion = 0
            deletion = 0
            if len(ref_base) > 1:
                deletion = 1
            if len(alt_base) > 1:
                insertion = 1

            for pileupcolumn in samfile.pileup(chr,int(pos)-1,int(pos)+1):
                for pileupread in pileupcolumn.pileups:
                    if(insertion == 0 & deletion == 0):
                        if pileupcolumn.pos == (int(pos) -1): # 0 based coordinates ??
                            if pileupread.query_position != None:
                                depth = depth + 1
                                base = pileupread.alignment.query_sequence[pileupread.query_position]
                                if base.upper() == ref_base:
                                    ref = ref + 1
                                elif base.upper() == alt_base:
                                    alt =  alt + 1
                    else:
                        if pileupcolumn.pos == int(pos)-1: # 0 based coordinates ??
                            if pileupread.query_position != None:
                                depth = depth + 1
                                if pileupread.indel != 0:
                                    alt = alt + 1
                                else:
                                    ref = ref + 1

            # if int(pos) == 4794873:
            #     print("ref:{}, alt:{}, depth:{}, insertion:{}, deletion:{}".format(ref,alt,depth,insertion,deletion),file = sys.stderr)
            #     print("ref_base:{}, alt_base:{}".format(ref_base,alt_base),file = sys.stderr)
            #     print("len(ref_base):{}".format(len(ref_base)),file = sys.stderr )
            #     print("len(alt_base):{}".format(len(alt_base)),file = sys.stderr )

            gtf[dp_pos] = str(depth)
            gtf[ad_pos] = "{},{}".format(ref,alt)

            if depth > 0:
                af = str("%f"%(float(alt)/float(depth)))
            else:
                af = '0'
            gtf[af_pos] = str(af)
            fields[missing] = ":".join(gtf)

            print('\t'.join(fields[0:fields_len]))

        gtf = re.split(r':',fields[n_pos])
        add_variant['{}_AD_ref'.format(n_id)] = re.split(',',gtf[ad_pos])[0]
        add_variant['{}_AD_alt'.format(n_id)] = re.split(',',gtf[ad_pos])[1]
        add_variant['{}_DP'.format(n_id)] = gtf[dp_pos]
        add_variant['{}_AF'.format(n_id)] = gtf[af_pos]

        gtf = re.split(r':',fields[t1_pos])
        add_variant['{}_AD_ref'.format(t1_sample_id)] = re.split(',',gtf[ad_pos])[0]
        add_variant['{}_AD_alt'.format(t1_sample_id)] = re.split(',',gtf[ad_pos])[1]
        add_variant['{}_DP'.format(t1_sample_id)] = gtf[dp_pos]
        add_variant['{}_AF'.format(t1_sample_id)] = gtf[af_pos]
        tp1_af = gtf[af_pos]

        gtf = re.split(r':',fields[t2_pos])
        add_variant['{}_AD_ref'.format(t2_sample_id)] = re.split(',',gtf[ad_pos])[0]
        add_variant['{}_AD_alt'.format(t2_sample_id)] = re.split(',',gtf[ad_pos])[1]
        add_variant['{}_DP'.format(t2_sample_id)] = gtf[dp_pos]
        add_variant['{}_AF'.format(t2_sample_id)] = gtf[af_pos]
        tp2_af = gtf[af_pos]

        diff = float(tp2_af) - float(tp1_af)
        add_variant["TP2-TP1.AF"] = "{0:.3f}".format(diff)

    # if args.debug:
    #     eprint("{}:{}:{}".format(chr,pos,alt_base))

    if(args.mode == 'simple'):

        gtf = re.split(r':',fields[9])

        ad_pos = gth.index("AD")
        dp_pos = gth.index("DP")
        if "AF" in gth:
            af_pos = gth.index("AF")
            af_val = gtf[af_pos]
        else:
            af_val = "."
            #eprint("gtf ad_pos before splitting: {}".format(gtf[ad_pos]))
            #eprint("whole gth: {}".format(gth))
            #eprint("whole gtf: {}".format(gtf))
            if gtf[ad_pos] != "." and  len(re.split(',',gtf[ad_pos])) > 1 and gtf[dp_pos] != "." and gtf[dp_pos] != "0":
                af_val = "{0:.3f}".format(float(re.split(',',gtf[ad_pos])[1]) / float(gtf[dp_pos]))

        add_variant['tumor_AD_ref'] = re.split(',',gtf[ad_pos])[0]
        add_variant['tumor_AD_alt'] = re.split(',',gtf[ad_pos])[1]
        add_variant['tumor_DP'] = gtf[dp_pos]
        add_variant['tumor_AF'] = af_val


    ### Annotation
    annotation = ""
    search = re.search("CSQ=(.+?);",INFO_field)
    if search != None :
        annotation = search.group(1)
        annotations = re.split(",", annotation)
        ann_list   = []

        for a in annotations:
            spl = re.split("\|",a)
            ann_list = ann_list + [spl]

        if ann_header != "" :
            ann_pd = pd.DataFrame(ann_list, columns = ann_header)
            sel = sort_annotation(ann_pd,found_anno).iloc[0]


        for anno in found_anno:
            #eprint(anno)
            if anno == "HGVSc" or anno == "HGVSp":
                coding = sel[anno]
                if coding != "":
                    coding = re.split(":",coding)[1]
                    add_variant[anno] = coding.replace('%3D','=')
                continue

            if anno in {"MAX_AF", "gnomAD_AF", "gnomAD_NFE_AF"}:
                af = sel[anno]
                if af != "":
                    add_variant[anno.lower()] = str("%f"%float(af))
                continue

            if anno in {"PUBMED", "MAX_AF_POPS","AF", "EUR_AF", "Feature_type"}:
                add_variant[anno.lower()] = sel[anno]
                continue

            if anno == "SYMBOL":
                add_variant["gene_symbol"] =  sel["SYMBOL"]
            if anno == "Feature":
                add_variant["Transcript_ID"] = sel["Feature"]
            if anno == "EXON":
                add_variant["Exon_Number"] = sel["EXON"]
            if anno == "VARIANT_CLASS":
                add_variant['Variant_Type'] = sel['VARIANT_CLASS']
            if anno == "Consequence":
                add_variant["Variant_Classification"] = sel[anno]
            if anno =='CLIN_SIG':
                add_variant["clinvar"] = sel[anno]

            add_variant[anno] = sel[anno]

            add_variant['alt_transcripts'] = squish(ann_pd).replace('%3D','=')

            # extract refseq NM ID
            transcriptID_list = add_variant['alt_transcripts'].split("|")
            add_variant['NM_transcriptID'] = ','.join([sub.split(';')[0] for sub in [transcriptID_list[i] for i, s in enumerate(transcriptID_list) if 'NM' in s]])

    # add Mutation_description, if there are two AF, use first one a round to 2 decimal
    af_tmp = float(af_val.split(",")[0])*100
    add_variant['Mutation_description'] = add_variant['HGVSc'] + " " + add_variant['HGVSp'] + " " + str("%.2f"%af_tmp) + "%"

    clinvar = re.search("clinvar_sig=(.+?);",INFO_field)
    if clinvar != None:
        add_variant["clinvar"] = clinvar.group(1)

    #strand bias
    saaf = re.search("SAAF=(.+?);",INFO_field)
    if saaf != None:
        sapp = re.search("SAPP=(.+?);",INFO_field)
        if sapp != None:
            add_variant['strand_bias'] = "{};{}".format(sapp.group(1),saaf.group(1))

    df = df.append(pd.DataFrame(add_variant, index=[0]),sort=True)

    # remove all multiallelic mutation except 20:31022441
    #df = df[(df['Start_Position'] == '31022441') | (~df['vcf_filter'].str.contains("multiallelic"))]


if(args.tab):
    df[idx].to_csv(args.out,index=False,sep="\t")
else:
    df[idx].to_csv(args.out,index=False)


        #Deadcode

        #reAnnotate mode
        # if("mutect" in set_from):
        #     #eprint("extracting mutect from: {}".format(set_from))
        #     ret = extract_gt(gth, fields[mtct_tum_pos], mut_df, chr, pos, alt_base,'tum')
        #     add_variant['Mutect2_tumor_DP'] = ret[0]
        #     add_variant['Mutect2_tumor_AF'] = ret[1]
        #     add_variant['Mutect2_tumor_AD_ref'] = ret[2][0]
        #     add_variant['Mutect2_tumor_AD_alt'] = ret[2][1]
        #
        #     ret = extract_gt(gth, fields[mtct_nor_pos], mut_df, chr, pos, alt_base,'norm')
        #     add_variant['Mutect2_normal_DP'] = ret[0]
        #     add_variant['Mutect2_normal_AF'] = ret[1]
        #     add_variant['Mutect2_normal_AD_ref'] = ret[2][0]
        #     add_variant['Mutect2_normal_AD_alt'] = ret[2][1]
        #
        #
        # if("strelka" in set_from):
        #     #eprint("extracting strelka from: {}".format(set_from));
        #     ret = extract_gt_strelka(gth, fields[strl_nor_pos], ref_base, alt_base)
        #     add_variant['Strelka2_tumor_DP'] = ret[0]
        #     add_variant['Strelka2_tumor_AF'] = ret[1]
        #     add_variant['Strelka2_tumor_AD_ref'] = ret[2][0]
        #     add_variant['Strelka2_tumor_AD_alt'] = ret[2][1]
        #
        #     ret = extract_gt_strelka(gth, fields[strl_tum_pos], ref_base, alt_base)
        #     add_variant['Strelka2_normal_DP'] = ret[0]
        #     add_variant['Strelka2_normal_AF'] = ret[1]
        #     add_variant['Strelka2_normal_AD_ref'] = ret[2][0]
        #     add_variant['Strelka2_normal_AD_alt'] = ret[2][1]
        #
        #
        # if('vardict' in set_from):
        #
        #     ret = extract_gt(gth, fields[vdct_tum_pos], var_df, chr, pos, alt_base,'tum')
        #     add_variant['Vardict_tumor_DP'] = ret[0]
        #     add_variant['Vardict_tumor_AF'] = ret[1]
        #     add_variant['Vardict_tumor_AD_ref'] = ret[2][0]
        #     add_variant['Vardict_tumor_AD_alt'] = ret[2][1]
        #
        #     ret = extract_gt(gth, fields[vdct_nor_pos], var_df, chr, pos, alt_base,'norm')
        #     add_variant['Vardict_normal_DP'] = ret[0]
        #     add_variant['Vardict_normal_AF'] = ret[1]
        #     add_variant['Vardict_normal_AD_ref'] = ret[2][0]
        #     add_variant['Vardict_normal_AD_alt'] = ret[2][1]
