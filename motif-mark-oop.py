#!/usr/bin/env python

import argparse
import re
import cairo

#argparse expressions to take in user input
def get_args():
    parser = argparse.ArgumentParser(description = " A program to find motifs within a fasta file")
    parser.add_argument("-f", "--file", help ="SAM file", required = True)
    parser.add_argument("-m", "--motifs", help = "file containing motifs that you are searching for", required = True)
    return parser.parse_args()
args = get_args()
fasta_file = args.file
motif_file = args.motifs 



#this line determines output file name based on input file name
output_file = "".join(fasta_file.split('.')[:-1]) + '.svg'
# setting available rbg colors for motifs
color_list = [(0.2,0.4,1), (1,0.6,0), (1,0,0.1), (1,0.3,1), (1,1,0.2),(0.1,1,0),(0.6,0.2,0.7)]
                # Blue       orange      red       pink      yellow     green       purple
                # I encoded for more than required...just in case!

GENE_HEIGHT = 65
Y_OFFSET = 30
xax = 65

with open(fasta_file, "r") as fasta:
    fasta_dict = {}
    gene_names = []
    for line in fasta:
        line = line.strip()
        if line.startswith('>'):
            header = line
            gene_names.append(header)
            read = ''
        else:
            read += line
        i = {header:read}
        fasta_dict.update(i)


with open(motif_file, 'r') as motif_file:
    motif_list = []
    for motif in motif_file:
        motif = motif.strip()
        motif_list.append(motif)


yax = 50
width, height = 1000,1500
surface = cairo.SVGSurface(output_file, width, height)
context = cairo.Context(surface)

#defining all dictionaries that will be used and liast that will be used
motif_dict = {}
color_dict_trans_motif = {}     #color dictionaryfor IUPAC translated motifs
color_dict = {}                 #color dictionary for original motif provided - colors will correspond to translated motifs



class Exon:
    '''
    An Exon object should be given enough information that it can figure out how to
    draw itself. Exons should be black.
    '''
    def __init__(self, line, context, fasta_dict, gene_number):
        self.line = line
        self.context = context
        self.fasta_dict = fasta_dict
        self.gene_number = gene_number
        self.xax = xax
        self.exon_coord = self.find_exon(line)
    def find_exon(self, line):
        '''
        A function to find the exon coordinates in a fasta file
        '''
        caps = "([AUTCG]+)"
        line = line.strip()
    #print(line)
        exon_indicies = [(m.start(0), m.end(1)-1) for m in re.finditer(caps, line)]
        return(exon_indicies)
    def draw_exon(self):
        '''
        This function is responsible for drawing the svg output of the input fasta and motif files
        '''
        yax = GENE_HEIGHT * self.gene_number + Y_OFFSET
        x = self.exon_coord[0][0]
        y = self.exon_coord[0][1]
        context.set_source_rgb(0, 0, 0)
        context.set_line_width(30)
        context.move_to(x+xax, yax)
        context.line_to(y+xax, yax)
        context.stroke()


class Motif:
    '''
    A Motif object should be given enough information that it can figure out how to
    draw itself. Each distinct Motif should have a different color 
    '''
    def __init__(self, motif_list, gene_number, color_dict, header):
        self.color_dict = color_dict
        self.gene_number = gene_number
        self.header = header
        self.iupac_dict = {
                        #this defines a dictionary full of IUPAC degeneragte bases to be mathed to giiven motifs
                        # keys are IUPAC Bases, values are Regex expressions
                        "A":"[Aa]",  #adenine
                        "C":"[Cc]",  #cytosine
                        "G":"[Gg]",  #Guanine
                        "T":'[TtUu]',  #thymine/uracil
                        "U":'[UuTt]',      #uracil/thymine
                        "W":'[AaTtUu]',  #weak
                        "S":'[CcGg]',  #strong
                        "M":'[AaCc]',  #amino
                        "K":'[GgTtUu]',  #keto
                        "R":'[AaGg]',  #Purine
                        "Y":'[CcTtUu]',  #pyrimidine
                        "B":'[CcGgTtUu]',  #not A
                        "D":'[AaGgTtUu]',  #not C
                        "H":'[AaCcTtUu]',  #not G
                        "V":'[AaCcGg]',  #not T
                        "N":'[AaCcGgTtUu]',  #Any one base
                        "Z":'[]'}   #zero
        self.iupac = self.iupac_finder(motif_list)
        self.motif_list = motif_list
        self.indicies_motif = self.motif_indicies()
        self.xax = 65
    def iupac_finder(self, motif_list):
        '''
        A function to define the IUPAC nucleotide base pairs and return the regex expression for a given motif 
        '''
        conv_dict = {}
        for motif in motif_list:
            conv_motif = ''
            motif = motif.upper()
            for mote in motif:
                if mote in self.iupac_dict:
                    conv_motif += self.iupac_dict[mote]
                else:
                    print("ERROR, mystery iupac detected: is that really what you wanted to do?")
            conv_dict[conv_motif] = len(motif)
        return(conv_dict)
    def motif_indicies(self):
        mo_ind = []
        for motif in self.iupac.keys():
            indicies_motif = [(m.start(0)) for m in re.finditer(rf'(?=({motif}))', fasta_dict[self.header], re.IGNORECASE)]
            mo_ind.append(indicies_motif)
            #print(indicies_motif)
        #print(mo_ind)
        return mo_ind
    def draw_motifs(self):
        '''
        This function is responsible for drawing the svg output of the input fasta and motif files
        '''
        for i, motif in enumerate(self.indicies_motif):
            yax = GENE_HEIGHT * self.gene_number + Y_OFFSET
            r,b,g = self.color_dict[motif_list[i]]
            context.set_source_rgb(r,b,g)
            line_width = 15
            for index in motif:
                context.rectangle(index + xax, yax-7, len(motif_list[i]), line_width)
                context.fill()
                context.stroke()




class Fastaheader:
    '''
    A FastaHeader object should be given enough information that it can figure out
    how to draw itself. FastaHeaders should be black.
    '''
    def __init__(self, header, gene_number):
        self.header = header
        self.gene_number = gene_number
        #this code creates a dictionary with keys as the header of the fasta(lines that start with a > and contain the gene name)
        #and the values of the dictionary as thee entire fasta sequence given, to be used for indexing and motif placement
    def draw_fastaheader(self):
        yax = GENE_HEIGHT * self.gene_number + Y_OFFSET
        context.set_source_rgb(0,0,0)   
        context.move_to(xax , yax - 30)
        context.show_text(self.header)

class Gene:
    def __init__(self, sequence ,gene_number):
        self.gene_number = gene_number
        self.sequence = sequence
        self.gene_length = len(sequence)
    def draw_gene(self):
        yax = GENE_HEIGHT * self.gene_number + Y_OFFSET
        context.set_source_rgba(0, 0, 0)
        context.set_line_width(2)
        context.move_to(xax,yax)
        context.line_to(self.gene_length + xax, yax)
        context.stroke()



class Genegroup:
    '''
    A GeneGroup object is responsible for telling its children (e.g., Genes,
    Motifs) to go draw themselves. MAKE ME DUMB
    '''
    def __init__(self, exon, motifs, gene, header):
        self.exon_indicies = exon
        self.motif_indicies = motifs
        self.gene_indicies = gene
        self.header = header
    def draw_geneobjects(self):
        self.exon_indicies.draw_exon()
        self.motif_indicies.draw_motifs()
        self.gene_indicies.draw_gene()
        self.header.draw_fastaheader()

color_number = 0
for motif in motif_list:
    color_dict[motif] = (color_list[color_number])
    color_number += 1

big_info_list = []
gene_number = 1
for header in fasta_dict.keys():
    exons_object = Exon(fasta_dict[header], context, fasta_dict, gene_number)
    #print("hi")
    motif_object = Motif(motif_list, gene_number, color_dict,header)
    gene_indicies = Gene(fasta_dict[header], gene_number)
    header = Fastaheader(header, gene_number)
    gene_groupies = Genegroup(exons_object, motif_object, gene_indicies, header)
    big_info_list.append(gene_groupies)
    gene_number += 1



for i in big_info_list:
    i.draw_geneobjects()

surface.finish()

# this code block fills lists and calls the defined functions above to fill the lists/dictionaries that are specified above
#it also sets corresponding colors for the svg outut
# gene_number = 1
# for header,sequence in fasta_dict.items():
#     exons_object = exon(sequence, context, fasta_dict, gene_number, xax, exon_dict)
#     exons_object.find_exon(sequence)
#     e = {header:exons_object}
#     exon_dict.update(e)
#     motif_object = motif(motif_list, iupac_dict)
#     motif_object.iupac_finder(motif_list)
#     motif_regex = iupac_finder(motif_list)
#     LN = 0
#     line = 0
#     for motif, length in motif_regex.items():
#         indicies_motif = [(m.start(0)) for m in re.finditer(rf'(?=({motif}))', sequence, re.IGNORECASE)]
#         info_list.append((header, motif, length, indicies_motif))
#         color_dict_trans_motif[motif] = (color_list[LN])
#         LN += 1
#     for motif in motif_list:
#         color_dict[motif] = (color_list[line])
#         line += 1
#         gene_number += 1
# def draw_svg():
#     '''
#     This function is responsible for drawing the svg output of the input fasta and motif files
#     '''
#     for motif in motif_list:
#         x,y,z = color_dict[motif]
#         context.set_source_rgb(0,0,0)
#         context.move_to(65 , yax)
#         context.show_text(motif)
#         context.set_source_rgb(x,y,z)
#         context.set_line_width(10)
#         context.move_to(52, yax-10)
#         context.line_to(60, yax+10)
#         context.stroke()
#         yax += 25
#     context.set_source_rgb(0,0,0)
#     context.move_to(65 , yax)
#     context.show_text("Exon")
#     context.set_source_rgb(0,0,0)
#     context.set_line_width(10)
#     context.move_to(52, yax-10)
#     context.line_to(60, yax+10)
#     context.stroke()
#     ###
#     for header, sequence in fasta_dict.items():
#         sequence = sequence.strip()
#         yax += 100
#         xax = 65
#         context.set_source_rgba(0, 0, 0)
#         context.set_line_width(2)
#         context.move_to(xax,yax)
#         context.line_to(len(sequence)+xax,yax)
#         context.stroke()
#         x = exon_dict[header][0][0]
#         y = exon_dict[header][0][1]
#         context.set_source_rgb(0, 0, 0)
#         context.set_line_width(30)
#         context.move_to(x+xax,yax)
#         context.line_to(y+xax,yax)
#         context.stroke()
#         for item in info_list:
#             if item[0] == header:
#                 for m, length in motif_regex.items():
#                     if m == item[1]:
#                         r,b,g = color_dict_trans_motif[m]
#                         for index in item[3]:
#                             line_width = 15
#                             context.set_source_rgb(r,b,g)
#                             context.rectangle(index+xax, yax-7, length, line_width)
#                             context.fill()
#                             context.set_source_rgb(0,0,0)   
#                             context.move_to(xax , yax -30 )
#                             context.show_text(header)
surface.finish()

#drawing the final output!!
# draw_svg()
