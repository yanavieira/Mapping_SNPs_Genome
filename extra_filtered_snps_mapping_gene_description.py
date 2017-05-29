#!/usr/bin/python

import argparse
import matplotlib.pyplot as plt

from collections import defaultdict

parser = argparse.ArgumentParser(description="Mapping Cg_Nara5 - SNPs "
    "charaterization - coding or non coding region and gene identification")

parser.add_argument("-in_sam", dest="sam_infile", help="The .sam file from the"
                    " Mapping", required=True)

parser.add_argument("-in_vcf", dest="vcf_infile",
                    help=" VCF file with the SNPs conserved in Ck and totally"
                    "differentiated from UG1/2")

parser.add_argument("-in_gff3", dest="gff3_infile", help="GFF3-file with"
                    "scaffold information from Cg_Nara5 (gene location,"
                    "description, intron/exons, etc")

parser.add_argument("-in_loci", dest="loci_infile", nargs="*", help="open all "
                    "loci files to parse the sequence and has the correct codon")

arg = parser.parse_args()


def sam_parsed(sam_file):
     """ Open the .sam file provided from the users and retain the information"
     "about the Cg_scaffold that the locis mapped."""

     sam_file= open(sam_file)

     sam_dic = {}
     read_frame_dic ={}
     count = 0
     counter_1 = 0
     counter_2 = 0
     #.sam file parsed - crucial information was retrited (scaffold information)
     # start - the starting position of the locus_sequence
     #reading_frame - locus in the correct sense [0] or CR [16]
     #sequence_locus - locus sequence information

     for line in sam_file:

         if line.startswith("@"):
             pass

         else:
             line_information = line.strip().split()
             scaffold = line_information[2]
             loci = line_information[0]
             mapping_beginning = line_information[3]
             read_frame = line_information [1]
             locus_sequence = line_information [9]
             cigar = line_information [5]
             if "D" in cigar or "I" in cigar:
                 count += 1
             if "D" in cigar and "I" in cigar:
                  counter_2 +=1
             a = count - counter_2
             if scaffold != '*':
                 sam_dic[loci] = {"scaffold": scaffold,
                                  "start": int(mapping_beginning),
                                  "reading_frame": read_frame,
                                  "sequence_locus": locus_sequence,
                                  "cigar": cigar}
                 counter_1 +=1
     print ("Number of loci mappead on Cg: {}".format(len(sam_dic)))

     print ("Step 1 - Parse the .sam file -- Done")

     #The sam_dic return a dictionary where the key is the locus(read) and the
     #value has the scaffold information, the position of the gene beginin,
     #the correct read frame of the gene, and finally the sequence of locus, in
     #the same reading frame of the Cg

     # 
    #  print ("Number of locus with insertion or deletion " + str(count))
    #  print ("Number of locus with insertion and deletion " + str(counter_2))
    #  print ("Number of locus with problems " + str(a))
     return sam_dic


def vcf_parsed (vcf_file, sam_dic):
     """ Open the .vcf file provided from the users and retain the information"
         "about the real SNPs localization """
     #real SNPs localiation - snp localiation in Cg_Nara5 genome
     #snps variation - the real modification between Ck and UG1 and UG2

     vcf_file = open(arg.vcf_infile)
     vcf_dic = {}
     snps_variation_dic ={}
     counter = 0
     counter_1 = 0
     for line in vcf_file:
         if line.startswith("##"):
             pass
         elif line.startswith("#CHROM"):
             pass
         else:
             line_information = line.strip().split()

             loci_chrom = line_information[0].split()
             loci_pos = line_information[1]
             snps_variation = line_information [3]
             snps_variation_2 = line_information [4]

             for i in loci_chrom:
                if i not in vcf_dic:
                    vcf_dic[i] = [int(loci_pos)]
                    snps_variation_dic[i] = [(snps_variation, snps_variation_2)]

                else:
                     vcf_dic[i].append(int(loci_pos))
                     snps_variation_dic[i].append((snps_variation, snps_variation_2))

     for locus, position_list in vcf_dic.items():

        if locus in sam_dic:
            sam_dic[locus]["loci_position"] = position_list
            counter += len(position_list)

            real_snp = []
            real_snp_cr = []
            for snp in position_list:
                if sam_dic[locus]["reading_frame"] == "0":
                    real_snp_loc = (sam_dic[locus]["start"] + snp)
                    real_snp.append(real_snp_loc)
                    sam_dic[locus]["real_snp_localization"] = real_snp
                else:
                    sequence_length = len(sam_dic[locus]["sequence_locus"])
                    real_snp_loc_cr = ((sequence_length - snp) + sam_dic[locus]["start"]) + 1
                    real_snp.append(real_snp_loc_cr)
                    sam_dic[locus]["real_snp_localization"] = real_snp

     for locus, variation in snps_variation_dic.items():
         if locus in sam_dic:
             sam_dic[locus]["snp_variation"]= snps_variation_dic[locus]




     print ("Number of SNPs mapped on Cg: {}".format(counter))

     print ("Step 2 - Parse the .vcf file -- Done")

     # The sam_dic return: Key - loci name; value - the same information previously
     # described as well as the real snp localization and the snp variation
     # present in the initial dataset (Ck vs UG1/2)

     return sam_dic

def gff3_parsed (gff3_file, sam_dic):
    """ Open the .gff3 file provided from the users and retain the information"
        about the Cg_Nara5 scaffolds (genes, exons, pseudogene).
        Add to the main dictionary the information related with the location
        of the SNPs in Cg genome as well as the gene description."""

    #A special type of dictionary in which the values were saved in a list
    gff_dic = defaultdict(list)

    gff3_file = open(arg.gff3_infile)
    gff3_dic = {}

    gene_dic = {}
    exon_list = []
    gene_idx = 1

    counter_1 = 0
    counter_2 = 0
    counter_3 = 0
    counter_4 = 0
    counter_5 = 0
    counter_6 = 0
    counter_7 = 0
    idx_pseudogene = 0

    #A dictionary
    gene_idexes = {"gene": gene_idx, "exon": gene_idx,
                   "pseudogene": "pseudogene"}


    for line in gff3_file:
        if line.startswith("##"):
            pass
        elif line.startswith("#!"):
            pass
        else:
            line_information = line.strip().split()

            # Make a dic with the genes present on Gg genome and its anotattion
            if line_information[2] == ("gene"):
                # deal with the PREVIOUS gene
                #This peace of code add to the gff3_dic(the main dic of gff3 file)
                #the information of which are the exons of one particular gene
                #Note: this happends at the same time that the gene information
                #were parsed
                if exon_list:
                    gff3_dic[gene_idx]["exon_list"] = exon_list
                    gene_idx += 1

                    exon_list = []
                #parse the gene information and add this information to a new dic (gff3_dic)
                #with all the information related to the genes present in gff3 file (Cg_Nara5)
                # deal with CURRENT gene
                scaffold = line_information [0]
                gene_beg = line_information[3]
                gene_end = line_information [4]
                gene_loc = [gene_beg, gene_end]
                gene_strand = line_information[6]
                gene_information = line_information [8]
                gene_information = line.strip().split(";")
                gene_description = [gene_information[2]]
                gff3_dic[gene_idx] = {"scaffold": scaffold,
                                      "gene_range": gene_loc,
                                      "description": gene_description,
                                      "exon_list": None,
                                      "strand": gene_strand}

            # Make a list with the exons-genes present on Gg genome and its anotattion
            # If in this line the "gene" keyword is not present but the "exon"
            #keyword are append the range information to the exon list which
            # will be added to main gff3 dic
            elif line_information[2] == ("exon"):
                exon_beg = line_information[3]
                exon_end = line_information [4]
                exon_loc = (exon_beg, exon_end)
                exon_list.append(exon_loc)

                exon_information = line_information [8]
                exon_information = line.strip().split()[8].split(";")[0]
                gff3_dic[gene_idx]["exon_reference"] = exon_information
            #At the same time - regardless the previous code if the line has
            #any of this keywords the information of the gene_range were added
            # to the gff_dic.
            if line_information[2] in ["gene", "exon", "pseudogene"]:

                gene_range = (line_information[3], line_information[4])

                #Note: this peace of code happends because the gene description
                #of the gene is not the same as the exon description. Therefore,
                #the gene description has to be recovered

                if line_information[2] == "gene":
                    gene_information = line_information [8]
                    gene_information = line.strip().split(";")
                    gene_description = [gene_information[2]]

                # Example:
                # gff_dic[scaffold1] = [[1, "gene", (82, 1159), description],
                #                        1, "exon", (82, 603), description],
                #                        2, "gene", (1440, 4998), description
                #                        pseudogene_idx, pseudogene, (1999, 3000)]]

                #To keep only the information regardless gene_idx (gene index)
                #to the gene or the exons present in this gene. When I have
                #pseudogenes, the gene index is replaced for pseudogene
                if line_information[2] in ["exon", "gene"]:
                    idx = gene_idx
                else:
                    idx_pseudogene += 1
                    idx = "pseudogene_"+ str(idx_pseudogene)

                #add the previous information in a different format in which
                #the key is the sacffold and the values are the index (to easly
                #acess the information present in gff3 dictionary), the keyword
                #(gene, exon, pseudogene), the range, and the description.
                #All these informations will be used to perfome the SNP range
                # discover only within the true scaffold and not in all the scaffolds
                #present in the gff3 file. Making the code mor efficient and realibel
                gff_dic[line_information[0]].append([idx,
                                                     line_information[2],
                                                     gene_range,
                                                     gene_description])

    # Add last exon list to last gene index\
    else:
        if exon_list:
            gff3_dic[gene_idx]["exon_list"] = exon_list

    print ("Step 3a - Parse the .gff3 file -- Done")


    for locus, info_dict in sam_dic.items():

        # Get all info from current scaffold
        # scaffold_info is a list containing all genes, exons and pseudogenes
        # of the scaffold in sam_dic

        scaffold_info = gff_dic[info_dict["scaffold"]]
        #we create two different "values" in the sam_dic dictionary with the len
        #of the real snp location in which all the "values" begin with "intergenic" or None
        #and as we make the check codes this values will be replaced for new
        # values or will be remain like this

        info_dict["element_type"] = ["intergenic"] * len(info_dict["real_snp_localization"])
        info_dict["element_range"] = [None] * len(info_dict["real_snp_localization"])
        info_dict["gene_index"] = "intergenic"

        # Check if locus is in any range
        # The enumerate function give the value of the "value" as well as the
        #position of the value. Example: l = ["a", "b", "c"]
        #enumerate (l) --- (0, "a"); (1, "b"); (2, "c")
        #pos - the position of the snp in the list
        #snp - is the real snp localization under analyse

        # Get the position of the snp in the list. This position will
        # be used to create a key for the gene_inf_dic.
        for pos, snp in enumerate(info_dict["real_snp_localization"]):
            # The "element" is the several lists present in the gff_dic.
            #Note: all the lists regardless the type has exactly the same length.
            # Example : [10459, 'gene', ('18930', '23805'), ['description=LysM domain-containing protein']
            #So for each list we will check if the SNP is in the range
            for element in scaffold_info:
                element_beg = int(element[2][0])
                element_end = int(element[2][1])
                element_range= range(element_beg, element_end)


                # YAY, one of the SNP matches one element of the scaffold
                if snp in element_range:

                    info_dict["gene_index"] = element[0]

                    # ELEMENT KEY:
                    # "exon": The SNP is in a coding region
                    # "gene": The SNP is in an intron
                    # "pseudogene": The SNP is in a pseudogene
                    info_dict["element_type"][pos] = element[1]

                    info_dict["element_range"][pos] = element[2]

                    info_dict["description"] = element[3]



    #Get the main statistics from our dataset

    for locus, locus_info in sam_dic.items():

        element_type = locus_info["element_type"]

        # Adding information for loci in a intergenic region
        #The set return an object with only 1 "element" in that case "intergenic"
        #So if the locus has 2 snps 1 in a intergenic region and other in a gene
        # this locus will not count as a intergenic locus, because the set will
        #have two elenets {"intergenic", "gene"} and not only 1 {"intergenic"}.
        #Note: The set works for each element_type present in sam_dic (loop)
        if set(element_type) == {"intergenic"}:
            counter_1 += 1

        # Adding information for SNPs in intergenic region
        #This counter gives the number of times the intergenic word appears
        counter_2 += element_type.count("intergenic")

        # Adding information for loci in pseudogenes
        if "pseudogene" in element_type:
            counter_3 += 1

        #Adding information for SNPs in pseudogene
        counter_4 += element_type.count("pseudogene")

        #Adding information for loci in genes
        #As previously refered the gene information were recorded in two different formats
        #gene- when the SNP were in a gene but not in a exon (aka intron)
        #exon - when the SNP were in a gene and in a specific exon
        #So in order to have the statistics for the gene we need to search
        #booth keywords on the element_type . Not in this particular case the set
        #doesn't work because the set don't has an order (gene, exon) or (exon, gene)

        if "gene" in element_type or "exon" in element_type:
            counter_5 += 1

        #Adding information for SNPs in gene

        counter_6 += element_type.count("exon") + element_type.count("gene")

        #Adding information for SNPs in exons

        counter_7 += element_type.count("exon")



    print("Data resume:")
    print("Number of loci in a non coding region: {}".format(counter_1))
    print("Number of SNPs in a non coding region: {}".format(counter_2))

    print("Number of loci located in pseudogenes:{}".format(counter_3))
    print("Number of SNPs located in pseudogenes:{}".format(counter_4))

    print("Number of loci located in genes: {}".format(counter_5))
    print("Number of SNPs located in genes: {}".format(counter_6))
    print("Number of SNPs located in exons: {}".format(counter_7))

#    print(gff3_dic[6207])
    #print(sam_dic["gene_idx"])
    return (sam_dic, gff3_dic)
#

def loci_parsed(loci_file):
    """ Open the directory where the locis files are and parse all .loci files
    having a dictionary with all the real sequences for the Ck, Ug1 and Ug2 """
    #
    ck_list = ["Rua_1", "Ang_6", "Cam_5", "Eti_20", "Cam_1", "Zim_14", "Cam_2","Que_2", "Ang_30", "Cam_8", "Que_48","Zim_12", "Zim_1", "Tan_13", "Que_72","UAG_5", "UAG_1", "UAG_3", "Ang_67","Ang_29"]

    ug1_list = ["Cg12063", "Cg125212", "Cg126212", "Cg12758", "Cg_432"]

    ug2_list = ["Cg128233", "Cg12824", "Cg12881", "Cg1291", "Cg8801"]

    loci_dic = {}

    loci_list = {"ck": None, "ug1": None, "ug2": None}

    for files in loci_file:
        name= files.strip().split ("/")
        name_loci = name[12].split("_")
        name_loci_1 = name_loci[1].split(".")
        real_name_loci = name_loci_1[0]


        loci_file = open(files)

        for line in loci_file:

            if line[:1] in '0123456789':
                pass
            else:

                line_information = line.strip().split()
                isolate = line_information[0]
                sequence = line_information [1]
                # if "-" in sequence:
                #     sequence = sequence.replace ("-", "")

                if isolate in ck_list and loci_list["ck"] == None:
                    loci_list["ck"] = sequence
                if isolate in ug1_list and loci_list["ug1"] == None:
                    loci_list["ug1"] = sequence
                if isolate in ug2_list and loci_list["ug2"] == None:
                    loci_list["ug2"] = sequence
                    loci_dic[real_name_loci] = loci_list



        loci_list = {"ck": None, "ug1": None, "ug2": None}

    return loci_dic


def get_codon_position_positivestrand(exon_list, snp_position):

    offset = 0
    prev_end = 0

    # Exon start
    start = int(exon_list[0][0])


    for rng in exon_list:

        if prev_end:
            offset += (int(rng[0]) - prev_end) - 1


        if snp_position in range(int(rng[0]), int(rng[1])):
            # print(snp_position, offset, start, rng)
            #
            loci_snp_position = ((snp_position - offset) - start)
            codon_pos = int((snp_position - offset) - start) % 3
            return (codon_pos, loci_snp_position)

        prev_end = int(rng[1])
    return ("error", "error")


def get_codon_position_negat_16(exon_list, start, snp_neg, rad_length):
    # print (exon_list)
    exon_list= [exon_list_inv [::-1] for exon_list_inv in exon_list]
    exon_list_inv = exon_list[::-1]
    # print(exon_list_inv)

    # beg_teste = exon_list[0][1]

    offset = 0

    prev_end = 0

    # Exon start
    start_gene = int(exon_list_inv[0][0])

    real_snp = start + (rad_length - snp_neg)

    for rng in exon_list_inv:

        if prev_end:
            offset += (prev_end - int(rng[0])) - 1

#
# [('134382', '133633'), ('134734', '134494'), ('135009', '134789')]
# [('135009', '134789'), ('134734', '134494'), ('134382', '133633')]

        if real_snp in range(int(rng[1]), int(rng[0])):
            rad_beg = ((start_gene) - (start + offset)) + 1
            # print("rad_begining " + str(rad_beg))
            loci_snp_position = rad_beg - (rad_length - snp_neg)
            codon_pos = int(rad_beg - (rad_length - snp_neg)) % 3
            #
            # print("offset " + str(offset))
            return (codon_pos, loci_snp_position)


        prev_end = int(rng[1])

        # if int(beg_teste) >= int(start):
        #     print ("coco")
    return ("error", "error")

    # print("offset" + str(offset))


def get_codon_position_negat_0(exon_list, start, snp_neg, rad_length):
    exon_list= [exon_list_inv [::-1] for exon_list_inv in exon_list]
    exon_list_inv = exon_list[::-1]
    # print(exon_list_inv)

    offset = 0

    prev_end = 0

    # Exon start
    start_gene = int(exon_list_inv[0][0])

    real_snp = start + snp_neg

    for rng in exon_list_inv:

        if prev_end:
            offset += (prev_end - int(rng[0])) - 1

#
# [('134382', '133633'), ('134734', '134494'), ('135009', '134789')]
# [('135009', '134789'), ('134734', '134494'), ('134382', '133633')]

        if real_snp in range(int(rng[1]), int(rng[0])):
            rad_beg = ((start_gene) - (start + offset)) + 1
            # print("rad_begining " + str(rad_beg))
            loci_snp_position = (rad_beg - snp_neg) + 1
            codon_pos = (int(rad_beg - snp_neg)+1) % 3

            # print("offset " + str(offset))
            return (codon_pos, loci_snp_position)


        prev_end = int(rng[1])

        # if int(beg_teste) >= int(start):
        #     print ("coco")

    # print("offset" + str(offset))
    return ("error", "error")

    # return ("hello world")


def get_rc_loci_dic (loci_dic):

    alt_map = {'ins':'0'}

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', "-": "-"}

    loci_dic_rc = {}

    loci_list_rc = {"ck": None, "ug1": None, "ug2": None}

    for loci, values in loci_dic.items():
        for k,v in alt_map.iteritems():
            ck_rc = str(values["ck"]).replace(k,v)
            bases = list(ck_rc)
            bases = reversed([complement.get(base,base) for base in bases])
            bases = ''.join(bases)
            for k,v in alt_map.iteritems():
                bases = bases.replace(v,k)
            loci_list_rc["ck"] = bases
        for k,v in alt_map.iteritems():
            ug1_rc = str(values["ug1"]).replace(k,v)
            bases = list(ug1_rc)
            bases = reversed([complement.get(base,base) for base in bases])
            bases = ''.join(bases)
            for k,v in alt_map.iteritems():
                bases = bases.replace(v,k)
            loci_list_rc["ug1"] = bases

        for k,v in alt_map.iteritems():
            ug2_rc = str(values["ug2"]).replace(k,v)
            bases = list(ug2_rc)
            bases = reversed([complement.get(base,base) for base in bases])
            bases = ''.join(bases)
            for k,v in alt_map.iteritems():
                bases = bases.replace(v,k)
            loci_list_rc["ug2"] = bases
        loci_dic_rc[loci] = loci_list_rc
        loci_list_rc = {"ck": None, "ug1": None, "ug2": None}



    return loci_dic_rc

def mutation_type(sam_dic, gff3_dic, loci_dic, loci_dic_rc):
    "This fucntion use the information present on the read_frame_1 to infer"
    "the type of mutation present(synonymous or non-synonymous). "
    "Note: 16- The read is on reverse complemented mode so the sequence has to"
    "reversed complemented before the analysis; "
    "0 - The read is in the same sense of the Cg sequence."

    start_position = []
    end_position = []
    percentage =[]
    count = 0
    count1 = 0
    count2 = 0
    count3 = 0
    bad_locus_list_final = []
    count_cigar_problem = 0
    count_syn_final = 0
    count_nsyn_final = 0
    count_error_final = 0
    count_len_exon_error = 0
    count_codon_position = 0

    #  for idx, gene_inf in gff3_dic.items():
    # locus_list= []
    for locus, snp_values in sam_dic.items():
        snp_values["aa_information"] = [None] * len(snp_values["real_snp_localization"])


        if "exon" in snp_values["element_type"]:

            if "I" in snp_values["cigar"]:
                count_cigar_problem += 1
                pass
            else:

                try:
                    exon_list = gff3_dic[snp_values["gene_index"]]["exon_list"]

                    exon_length = sum([(int(end) - int(start)) + 1 for start, end in exon_list])

                    if exon_length % 3 != 0:
                        count_len_exon_error += 1
                        # print(locus)

                        pass

                    else:
                        if gff3_dic[snp_values["gene_index"]]["strand"] == "+":

                            for p, snp in enumerate(snp_values["real_snp_localization"]):
                                if snp_values["element_type"][p] == "exon":

                                    codon_pos, loci_snp_position = get_codon_position_positivestrand(exon_list, snp)
                                    if codon_pos == "error":
                                        count_codon_position +=1

                                    if snp_values["reading_frame"] == "0":
                                        # print ("lalal")

                                        real_codon, count_error = find_real_codon(loci_dic[locus], codon_pos, sam_dic[locus]["loci_position"][p], locus)

                                        bad_locus_list, aa_dic_final, count_syn, count_nsyn, mutation = translate_dna(real_codon, locus)

                                        count_syn_final += count_syn
                                        count_nsyn_final += count_nsyn
                                        count_error_final += count_error


                                        if bad_locus_list != None:
                                            bad_locus_list_final.append(bad_locus_list)
                                        snp_values["aa_information"][p] = mutation

                                    # real_codon =  translate_dna(sam_dic, codon_pos, loci_snp_position, gff3_dic)
                                        # print ("--------------------------------------")
                                        # print ("snp position " + str(loci_snp_position))
                                        # print("codon position "+ str(codon_pos))
                                        # print("locus " + str(locus))
                                        # # print(exon_list, exon_length, snp)
                                        # print("gene identification " + str(gff3_dic[snp_values["gene_index"]]["exon_reference"]))
                                        # print("snp variation " + str(snp_values["snp_variation"]))
                                        # print("reading frame " + str(snp_values["reading_frame"]))
                                        # # rad_pos = snp_values["loci_position"][p]
                                        # # print(rad_pos)
                                        # print("sequence locus " + str(snp_values["sequence_locus"]))
                                        # print("strand " + str(gff3_dic[snp_values["gene_index"]]["strand"]))
                                        # print("real codon " + str(real_codon))
                                        # print("loci_dic " + str(loci_dic[locus]))
                                        # print("aa_information" + str(aa_dic_final))
                                        # count += 1
                                        # print("0+" +str(count))

                                    if snp_values["reading_frame"] == "16":
                                        # print ("lalal")
                                        new_loci_position = (len(sam_dic[locus]["sequence_locus"]) - int(sam_dic[locus]["loci_position"][p]))+1
                                        real_codon, count_error = find_real_codon(loci_dic_rc[locus], codon_pos, new_loci_position, locus)

                                        bad_locus_list, aa_dic_final, count_syn, count_nsyn, mutation = translate_dna(real_codon, locus)
                                        count_syn_final += count_syn
                                        count_nsyn_final += count_nsyn
                                        count_error_final += count_error

                                        if bad_locus_list != None:
                                            bad_locus_list_final.append(bad_locus_list)
                                        snp_values["aa_information"][p] = mutation
                                        # real_codon =  translate_dna(sam_dic, codon_pos, loci_snp_position, gff3_dic)
                                        # print ("--------------------------------------")
                                        # print ("snp position " + str(loci_snp_position))
                                        # print("codon position "+ str(codon_pos))
                                        # print("locus " + str(locus))
                                        # # print(exon_list, exon_length, snp)
                                        # print("gene identification " + str(gff3_dic[snp_values["gene_index"]]["exon_reference"]))
                                        # print("snp variation " + str(snp_values["snp_variation"]))
                                        # print("reading frame " + str(snp_values["reading_frame"]))
                                        # # rad_pos = snp_values["loci_position"][p]
                                        # # print(rad_pos)
                                        # print("sequence locus " + str(snp_values["sequence_locus"]))
                                        # print("strand " + str(gff3_dic[snp_values["gene_index"]]["strand"]))
                                        # print("real codon " + str(real_codon))
                                        # print("loci_dic_rc " + str(loci_dic_rc[locus]))
                                        # print("aa_information" + str(aa_dic_final))
                                        # count1 += 1
                                        # print("16+" +str(count1))



                        else:
                            for p, snp_neg in enumerate(snp_values["loci_position"]):
                                if snp_values["element_type"][p] == "exon" and snp_values["reading_frame"] == "16":
                                    start = snp_values["start"]
                                    rad_length = len(sam_dic[locus]["sequence_locus"])

                                    codon_pos, loci_snp_position = get_codon_position_negat_16(exon_list, start, snp_neg, rad_length)
                                    if codon_pos == "error":
                                        count_codon_position +=1
                                    real_codon, count_error = find_real_codon(loci_dic[locus], codon_pos, sam_dic[locus]["loci_position"][p], locus)

                                    bad_locus_list, aa_dic_final, count_syn, count_nsyn, mutation = translate_dna(real_codon, locus)

                                    count_syn_final += count_syn
                                    count_nsyn_final += count_nsyn
                                    count_error_final += count_error

                                    if bad_locus_list != None:
                                        bad_locus_list_final.append(bad_locus_list)
                                    snp_values["aa_information"][p] = mutation
                                    # print ("--------------------------------------")
                                    # print ("snp position " + str(loci_snp_position))
                                    # print("codon position " + str(codon_pos))
                                    # print("locus " + str(locus))
                                    # print("gene identification " + str(gff3_dic[snp_values["gene_index"]]["exon_reference"]))
                                    # print("snp variation "+ str(snp_values["snp_variation"]))
                                    # print("reading frame " + str(snp_values["reading_frame"]))
                                    # print("sequence locus " + str(snp_values["sequence_locus"]))
                                    # print("strand " + str(gff3_dic[snp_values["gene_index"]]["strand"]))
                                    #
                                    # print("real codon " + str(real_codon))
                                    # print("loci_dic" + str(loci_dic[locus]))
                                    # print("aa_information" + str(aa_dic_final))
                                    # count2 += 1
                                    # print("16-" +str(count2))

                                else:
                                    start = snp_values["start"]
                                    rad_length = len(sam_dic[locus]["sequence_locus"])

                                    codon_pos, loci_snp_position = get_codon_position_negat_0(exon_list, start, snp_neg, rad_length)
                                    if codon_pos == "error":
                                        count_codon_position +=1
                                    new_loci_position = (len(sam_dic[locus]["sequence_locus"]) - int(sam_dic[locus]["loci_position"][p]))+1

                                    real_codon, count_error = find_real_codon(loci_dic_rc[locus], codon_pos, new_loci_position, locus)

                                    bad_locus_list, aa_dic_final, count_syn, count_nsyn, mutation = translate_dna(real_codon, locus)

                                    count_syn_final += count_syn
                                    count_nsyn_final += count_nsyn
                                    count_error_final += count_error
                                    if bad_locus_list != None:
                                        bad_locus_list_final.append(bad_locus_list)

                                    snp_values["aa_information"][p] = mutation
                                    # print ("--------------------------------------")
                                    # print ("snp position " + str(loci_snp_position))
                                    # print("codon position " + str(codon_pos))
                                    # print("locus " + str(locus))
                                    # print("gene identification " + str(gff3_dic[snp_values["gene_index"]]["exon_reference"]))
                                    # print("snp variation "+ str(snp_values["snp_variation"]))
                                    # print("reading frame " + str(snp_values["reading_frame"]))
                                    # print("sequence locus " + str(snp_values["sequence_locus"]))
                                    # print("strand " + str(gff3_dic[snp_values["gene_index"]]["strand"]))
                                    # print("real codon " + str(real_codon))
                                    # print("loci_dic_rc " + str(loci_dic[locus]))
                                    # print("aa_information" + str(aa_dic_final))
                                # count3 += 1
                                # print("0-" +str(count3))



                except KeyError:
                    pass
    # print (bad_locus_list_final)
    print("Synonymous " + str(count_syn_final))
    print("Non-synonymous " + str(count_nsyn_final))
    #Insertion and Delition in .sam files
    print("Number of snps with cigar problems " + str(count_cigar_problem))
    #Not all the codons has 3 letters. Somethimes the sequences ends before the codon ends. 
    #Sometimes the codon only has 1 letter. This make the codon traduction to a.a impossible.
    print("Number of snps with real codon problems " + str(len(bad_locus_list_final)))
    #Not all the groups ("Ug1", "Ug2" and "Ck") has an codon - Missing data
    print("Number of snps without information for all type of groups " + str(count_error_final))
    #The exon division by 3 in .gff3 file is different from 0. Something is wrong on .gff3 file
    print("Len of exons with problems " + str(count_len_exon_error))
    #The function of codon position sometimes give us an error.
    print("codon position error " + str(count_codon_position))

    return(codon_pos, loci_snp_position, sam_dic)



def find_real_codon (loci_value, codon_pos, snp_position, locus):
    real_codon = {"ck": None, "ug1": None, "ug2": None}
    real_codon_final = {}
    count_error = 0
    if loci_value["ug1"] == None or loci_value["ck"] == None or loci_value["ug2"] == None:
        count_error += 1
        pass

    else:

        if codon_pos == 0:
            real_codon["ck"] = loci_value["ck"][snp_position-3:snp_position]
            real_codon["ug1"] = loci_value["ug1"][snp_position-3:snp_position]
            real_codon["ug2"] = loci_value["ug2"][snp_position-3:snp_position]

        elif codon_pos == 2:
            real_codon["ck"] = loci_value["ck"][snp_position-2: snp_position+1]
            real_codon["ug1"] = loci_value["ug1"][snp_position-2: snp_position+1]
            real_codon["ug2"] = loci_value["ug2"][snp_position-2: snp_position+1]

        elif codon_pos == 1:
            real_codon["ck"] = loci_value["ck"][snp_position-1:snp_position+2]
            real_codon["ug1"] = loci_value["ug1"][snp_position-1:snp_position+2]
            real_codon["ug2"] = loci_value["ug2"][snp_position-1:snp_position+2]

        real_codon_final[locus] = real_codon

    return (real_codon_final, count_error)




def translate_dna(real_codon, locus):
    count = 0
    bad_locus_list = None
    aa_dic = {}
    aa_dic_final = {}
    count_syn = 0
    count_nsyn = 0
    mutation = None
    codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    if len(str(real_codon[locus]["ck"])) < 3:
        bad_locus_list = locus

    else:
        aa_dic["aa_ck"] = codontable[real_codon[locus]["ck"]]
        aa_dic["aa_ug1"] = codontable[real_codon[locus]["ug1"]]
        aa_dic["aa_ug2"] = codontable[real_codon[locus]["ug2"]]
        if len(set([aa_dic["aa_ug2"], aa_dic["aa_ck"], aa_dic["aa_ug1"]])) == 1:
            count_syn += 1
            mutation = "synonymous"
        elif len(set([aa_dic["aa_ug2"], aa_dic["aa_ck"], aa_dic["aa_ug1"]])) and (aa_dic["aa_ck"] != aa_dic["aa_ug2"] or aa_dic["aa_ck"] != aa_dic["aa_ug1"]):
            count_nsyn += 1
            mutation ="non-synonymous"
    aa_dic_final[locus] = aa_dic

    return(bad_locus_list, aa_dic_final, count_syn, count_nsyn, mutation)


def statistics(sam_dic, gff3_dic):

    gene_dic = defaultdict(int)
    gene_dic_exon = defaultdict(int)
    syn_per_gene = defaultdict(int)
    nsys_per_gene = defaultdict(int)

    for locus, snp_values in sam_dic.items():

        gene_dic[sam_dic[locus]["gene_index"]] += sam_dic[locus]["element_type"].count("exon")
        syn_per_gene[sam_dic[locus]["gene_index"]] += sam_dic[locus]["aa_information"].count("synonymous")
        nsys_per_gene[sam_dic[locus]["gene_index"]] += sam_dic[locus]["aa_information"].count("non-synonymous")
        # gene_dic[sam_dic[locus]["gene_index"]] += sam_dic[locus]["element_type"].count("exon") + sam_dic[locus]["element_type"].count("gene")

        if nsys_per_gene[sam_dic[locus]["gene_index"]] >= 2 or gene_dic[sam_dic[locus]["gene_index"]] >=3:
            print(sam_dic[locus]["gene_index"],locus, gene_dic[sam_dic[locus]["gene_index"]], nsys_per_gene[sam_dic[locus]["gene_index"]], sam_dic[locus]["description"])

        # if locus in ['8027','5876','32847','33281', '24324','2364','39782','9642']:
        #     print(locus, gene_dic[sam_dic[locus]["gene_index"]], nsys_per_gene[sam_dic[locus]["gene_index"]], sam_dic[locus]["description"])


             #graphics - number of SNPs located on genes per scaffold

    mpl_fig5 = plt.figure()
    ax = mpl_fig5.add_subplot(111)
    plt.bar(range(len(nsys_per_gene)), nsys_per_gene.values(), align='center')
    plt.xticks(range(len(nsys_per_gene)), nsys_per_gene.keys())
    # ax.set_ylim(0, 35)
    ax.set_ylabel('Nr of non-synonymous mutation')
    ax.set_xlabel('gene')
    # ax.set_title("Number of SNPs per gene")
    plt.savefig("nr_non_syn_mutation.png")

 
def main():

    sam_dic_1 = sam_parsed(arg.sam_infile)
    sam_dic_2= vcf_parsed(arg.vcf_infile, sam_dic_1)
    sam_dic_3, gff3_dic_1= gff3_parsed(arg.gff3_infile, sam_dic_2)
    loci_dic = loci_parsed(arg.loci_infile)
    loci_dic_rc = get_rc_loci_dic (loci_dic)
    codon_pos_1, loci_snp_position_1, sam_dic_4 = mutation_type(sam_dic_3, gff3_dic_1, loci_dic, loci_dic_rc)
    statistics_1 = statistics(sam_dic_4, gff3_dic_1)

    #
    # scaffold_dic = scaffold_information(sam_dic_3, gff3_dic_scaffold)

main()
