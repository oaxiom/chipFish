"""

format specifiers. Part of glbase

This is now the approved way to get at the format specifiers:

format=format.bed

TODO:
Implement all the fileformats from::

http://genome.ucsc.edu/FAQ/FAQformat#format5.1

See below for the catalogue of file formats

"""

import csv

from format_container import fc

# ------------------- standard formats:
bed = fc(name="bed",
    format={"loc": "location(chr=column[0], left=column[1], right=column[2])",  
        "name": 3, "score": 4, "strand": 5,
        "force_tsv": True, "skiplines": -1},
    description= "The 'normal' 6 column definition of a BED file, containing location, name, score and strand")
    
minimal_bed = fc(name="minimal_bed",
    format={"loc": "location(chr=column[0], left=column[1], right=column[2])", "force_tsv": True,
        "skiplines": -1},
    description= "A minimal BED file contains just the location data, in columns 0, 1 and 2.")

full_bed = fc(name="full_bed", # The full formal definition of a BED file.
    format={"loc": "location(chr=column[0], left=column[1], right=column[2])",  
        "name": 3, "score": 4, "strand": 5, "thickStart": 6, "thickEnd": 7,
        "itemRgb": 8, "blockCount": 9, "blockSizes": 10, "blockStarts": 11, 
        "force_tsv": True, "skiplines": -1},
    description= "This is the full 12-column definition of a BED file")
    
# Would be useful to add a "optional" and "required" for bed files? 
    
bed_no_strand = minimal_bed # Old name, kept for compatability

psl = fc(name="psl",
    format={"skiplines": -1, "force_tsv": True,
        "matches": 0, "misMatches": 1, "repMatches": 2, "nCount": 3,
        "qNumInsert": 4, "qBaseInsert": 5, "tNumInsert": 6, "tBaseInsert": 7,
        "strand": 8, "qName": 9, "qSize": 10, "qStart": 11, "qEnd": 12,
        "tName": 13, "tSize": 14, "tStart": 15, "tEnd": 16,"blockCount": 17, 
        "blockSizes": 18, "qStarts": 19, "tStarts": 20},
    description="PSL files, as produced by BLAT")
    
fasta = fc(name="fasta",
    description="FASTA sequence file",
    format={"special": "fasta"})

gtf = fc(name="gtf",
    description="GTF, gene transfer file format.",
    format={"feature_type": 1, "feature": 2, "gtf_decorators": 8,
	    "loc": "location(chr=column[0], left=column[3], right=column[4])", 
	    "strand": 6, "skiplines": -1, "force_tsv": True})

snp = fc(name="snp",
    format=dict(bin=0, name=4, score=5, strand=6, refNCBI=7, refUCSC=8, observed=9,
        molType=10, cclass=11, valid=12, avHet=13, avHetSE=14, func=15, locType=16, weight=17,
        loc={"code": "location(chr=column[1], left=column[2], right=column[3])"}, force_tsv=True),
    description="snp.txt file format reader")

pgsnp = fc(name="pgsnp",
    format={"loc": "location(chrom=column[0], left=column[1], right=column[2]",
        "name": 3, "alleleCount": 4, "alleleFreq": 5,
        "alleleScores": 6, "skiplines": -1, "forse_tsv": True},
    description="Personal Genome SNP file format")

encode_rna_expn = fc("encode_rna_expn",
    format={"loc": "location(chr=column[0], left=column[1], right=column[2])",  
        "name": 3, "score": 4, "strand": 5,
        "level": 6, "signif": 7, "score2": 8,
        "force_tsv": True, "skiplines": -1},
    description="ENCODE RNA elements: BED6 + 3 scores format")

encode_broadpeak = fc("encode_broadpeak",
    format={"loc": "location(chr=column[0], left=column[1], right=column[2])",  
        "name": 3, "score": 4, "strand": 5,
        "signalValue": 6, "pValue": 7, "qValue": 8,
        "force_tsv": True, "skiplines": -1},
    description="ENCODE broadPeak: Broad Peaks (or Regions) format")

# --------------------- motif discovery formats
fimo_out = fc(name="fimo_out", 
    description="Load in the fimo.txt file output by FIMO, part of the MEME suite.",
    format={"force_tsv": True, "pattern-name": 0, "sequence": 1, "left": 2, "right": 3, "score": 4,
        "p-value": 5, "q-value": 6, "sequence": 7, 
        "loc": "location(chr=column[1], left=min(column[2], column[3]), right=max(column[2], column[3]))"})

# --------------------- ChIP-seq file-types
macs_output = {"loc": {"code": "location(chr=column[0], left=column[1], right=column[2])"}, "tag_height": 5,
    "skiptill": "chr", "force_tsv": True, "fold_change": 7}

macs_summit = {"loc": {"code": "location(chr=column[0], left=int(column[1])+int(column[4]), right=int(column[1])+int(column[4]))"},
	"fold_change": 7, "tag_count": 5, 
	"skiptill": "chr", 
	"force_tsv": True}

ccat_output = {"loc": {"code": "location(chr=column[0], left=column[1], right=column[1])"}, "force_tsv": True,
    "tag_height": 4, "fold": 6, "skiplines": -1}

ccat_output_csv = {"loc": {"code": "location(chr=column[0], left=column[1], right=column[1])"},
    "tag_height": 4, "fold": 6, "fold_change": {"code": "barSplitter(3)[2]"}, "skiplines": -1}

# load in SISSRS file.
sissrs_output = {"loc": {"code": "location(chr=column[0], left=column[1], right=column[2])"}, "tag_height": 3,
    "fold": 4, "p-value": 5, "force_tsv": True, "skiplines": 57}

# --------------------- next-gen sequencing formats:
# This is not a full implementation of the sam file specification.
# but it will accept a sam file as outputted by tophat.
sam_tophat = fc(name="sam_tophat",
    format={"name": 0, "loc": {"code": "location(chr=column[2], left=column[3], right=column[3])"},
        "mapq": 3, "seq": 9, "force_tsv": True, "skiplines": -1},
    description="SAM file. Note this only implements SAM files as output by tophat")

exporttxt_loc_only = fc(name="exporttxt_loc_only",
    format={"loc": {"code": "location(chr=column[10].strip(\".fa\"), left=column[12], right=int(column[12])+25)"},
        "strand": 13,
        "force_tsv": True}, # export.txt file (output from the illumina pipeline), but only loads the location and strand.
    description="Load in the export.txt file as output by the Illumina HTS pipeline.\n\t\tThis variant loads only the read location")

exporttxt_all = fc(name="exporttxt_all",
    format={"loc": {"code": "location(chr=column[10].strip(\".fa\"), left=column[12], right=int(column[12])+25)"},
        "strand": 13, "seq": 6, "quality_score": 7,
        "force_tsv": True}, # export.txt file (output from the illumina pipeline), but only loads the location and strand.
    description="Load in the export.txt file as output by the Illumina HTS pipeline.\n\t\tThis variant loads the location, sequence and quality score")

# --------------------- miscellaneous file formats
ann_list = {"loc": 4, "strand": 6, "name": 9, "refseq": 11, "entrez": 12, "tag_height": 1}# output format from my annotation script
peak_loc = {"loc": 4, "name": 9, "refseq": 11, "entrez": 12, "tag_height": 1} # A list of peak locations, annotated?
genespring_out = {"refseq": 4, "entrez": 8, "fold_change": 1, "array_systematic_name": 0, "force_tsv": True}# A list out of GeneSpring
array_simplified = {"refseq": 10, "array_systematic_name": 0, "entrez": 8}
illumina_anotations = {"array_systematic_name":0, "refseq": 3, "entrez": 7}

# --------------------- snpXXX.txt file
snp_txt = dict(bin=0, name=4, score=5, strand=6, refNCBI=7, refUCSC=8, observed=9,
    molType=10, cclass=11, valid=12, avHet=13, avHetSE=14, func=15, locType=16, weight=17,
    loc={"code": "location(chr=column[1], left=column[2], right=column[3])"}, force_tsv=True)

# --------------------- microarray-like file formats
array_tsv = {"refseq": 0, "entrez": 1, "symbol": 2,
        "conditions": {"code": "column[4:]"}, "array_systematic_name": 1, "duplicates_key": False,
        "force_tsv": True}

array_csv = {"refseq": 0, "entrez": 1, "symbol": 2,
        "conditions": {"code": "column[4:]"}, "array_systematic_name": 1, "duplicates_key": False}

mm8_refgene = fc(name="mm8_refgene",
    format={"loc": {"code": "location(chr=column[0], left=column[2], right=column[3])"},
        "strand": 1, "name": 4, "description": 5, "dialect" : csv.excel_tab,
        "refseq": 6, "tss_loc": {"code": "strandSorter(column[0], column[2], column[3], column[1])"},
        },
    description="The mm8 refGene table downloaded from UCSC Genome Browser")

mm9_refgene = fc(name="mm9_refgene",
    format={"loc": "location(chr=column[2], left=column[4], right=column[5])",
        "strand": 3, "name": 12, "force_tsv": True,
        "refseq": 1, "tss_loc": {"code": "strandSorter(column[2], column[4], column[5], column[3])"},
        "cds_loc": {"code": "location(chr=column[2], left=column[6], right=column[7])"},
        "exons_count": 8
        }, # description is lost from the mm9 table?
    description="The mm9 refGene table downloaded from UCSC Genome Browser")

ensembl = {"loc": {"code": "location(chr=column[2], left=column[4], right=column[5])"},
	"tss_loc": {"code": "strandSorter(column[2], column[4], column[5], column[3])"},
	"ensmbl": 1, "name": 12, "exon_count": 8,
	"force_tsv": True, "skiplines": 0} 
	
# hg18 default refseq export.
hg18_refseq = fc(name="hg18_refseq",
    format={"loc": {"code": "location(chr=column[0], left=column[1], right=column[2])"},
        "strand": 5, "force_tsv": True,
        "refseq": 3, "tss_loc": {"code": "strandSorter(column[0], column[2], column[3], column[1])"}},
    description="The Hg18 RefSeq gene table as downloaded from the UCSC Genome Browser")

homer_annotated= {"loc": "location(chr=column[0], left=column[1], right=column[2])",
    "peak_id":3,"motif_score":4,"strand": 5,"motif_seq":6,"summit_dist":7,
    "summit_height":8,"neares_gene":9,"TSS_dist":10,"annotation":11, 
    "force_tsv": True,"skiplines": -1} 

MACS_combined={"loc": "location(chr=column[0], left=column[1], right=column[2])","peak_id":3,
    "summit_height":4,"p-value_score": 5,"number_tags":6,"fold_enrichment":7,"FDR":8, "force_tsv": True,
    "skiplines": -1}

# --------------------- class container

class fccatalogue():
    def __init__(self, formats):
        # put in dict by name
        
        self.formats = {}
        for item in formats:
            self.formats[item.name] = item
        
    def __str__(self):
        print "Found %s formats" % len(self.formats)
        a = []
        for key in self.formats:
            a.append("format.{name:24} - {desc:5}".format(name=key, desc=self.formats[key].description))
        return("\n".join(a))
    
    def __iter__(self):
        for k in self.formats:
            yield self.formats[k]
            
    def __len__(self):
        return(len(self.formats))
        
    def find(self, value):
        """
        **Purpose**
            find a particular motif based on value
            
        **Arguments**
            value
                searches through the list of formats and returns possible matches, searches the name
                and description to find relevant formats.
        """

        a = []
        for key in self.formats:
            if value.lower() in self.formats[key].name.lower() or value in self.formats[key].description.lower():
                a.append("format.{name:22} - {desc:6}".format(name=key, desc=self.formats[key].description))
        if a:
            print "\n".join(a)
        else:
            print "None found"

catalogue = fccatalogue([fimo_out, fasta, gtf, bed, full_bed, minimal_bed,
    exporttxt_loc_only, exporttxt_all, mm8_refgene, mm9_refgene, hg18_refseq,
    snp, pgsnp, psl, encode_rna_expn, encode_broadpeak, sam_tophat])
         
if __name__ == "__main__":
    print catalogue
    #print "Find:"
    #catalogue.find("bed")
        