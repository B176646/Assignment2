#!/usr/bin/python3
import subprocess,os,re

#Obtain protein sequence data
os.mkdir("raw_protein_sequences")
print("Begin search")

#Define a function to search relevant seqences and returns the size of path
def search_proteinseqs():
    toxonomics = input("Please enter the taxonomic group:")
    protein_family = input("Please enter the protein family:")
    myfile = "./raw_protein_sequences/proteinseqs.fa"
    subprocess.call("esearch -db protein -query \"%s[organism] AND %s[PROT] NOT Partial\" | efetch -format fasta >./raw_protein_sequences/proteinseqs.fa" %(toxonomics,protein_family),shell=True)
    size = os.path.getsize(myfile)
    return size

#Check whether the sequences are obtained
size_of_path = search_proteinseqs()
while size_of_path == 0:
    print("Input error,please try again")
    size_of_path = search_proteinseqs()

#calculate number of raw sequences
open_proteinseqs = open("./raw_protein_sequences/proteinseqs.fa")
seqs_number = 0
for line in open_proteinseqs:
    if line.startswith(">"):
       seqs_number = seqs_number + 1
count_seqs = open("./raw_protein_sequences/number_of_proteinseqs" + "txt","w")
count_seqs.write("The number of raw protein sequences is " + str(seqs_number) + "." + "\n")
count_seqs.close()

#Determine the level of conservation
#Align all protein sequences
os.mkdir("conservation_level")
subprocess.call("clustalo -i ./raw_protein_sequences/proteinseqs.fa -o ./conservation_level/clustalo_seqs.fa -v",shell = True)

#Make dictionary for all target sequences
sequences = {}
f = open("./conservation_level/clustalo_seqs.fa")
line = f.readline()
while line:
    if line.startswith(">"):
        name = line.rstrip("\n")
        sequences[name] = ""
    else:
        sequences[name] = sequences[name] + line.rstrip("\n")
    line = f.readline()
f.close()

#Define a function to calculate the "-" number of each sequence 
def calculate():
    f = open("./conservation_level/similarity_calculate.txt","a")
    for key,value in sequences.items():
        print(str(value.count("-")) + "\t" + key)
        f.writelines(str(value.count("-")) + "\t" + key + "\n")
    f.close()
    print("Calculation Finished")
calculate()

#Select 250 sequences with least number of "-" and put it in a txt 
subprocess.call("sort -n -k1 ./conservation_level/similarity_calculate.txt | head -250 > ./conservation_level/most_similar.txt",shell = True)
subprocess.call("rm -f ./conservation_level/similarity_calculate.txt",shell = True)

#Remove the number and partial information to get a txt only contained the names
subprocess.call("cut -f 2 ./conservation_level/most_similar.txt > ./conservation_level/prepare_for_blastp.txt",shell = True)
subprocess.call("rm -f ./conservation_level/most_similar.txt",shell = True)
a = open("./conservation_level/prepare_for_blastp.txt").read().split()
b = open("./conservation_level/only_name" + ".txt","w")
for m in a:
    if re.search(r'^>',m):
        b.write(m.strip(">") + "\n")
b.close()

#Use the name file to extract protein sequences with pullseq
m = open("./conservation_level/only_name.txt").read()
subprocess.call("/localdisk/data/BPSM/Assignment2/pullseq -i ./raw_protein_sequences/proteinseqs.fa -n ./conservation_level/only_name.txt > ./conservation_level/pul_seqs.fa",shell = True)
subprocess.call("rm -f ./conservation_level/only_name.txt",shell = True)

#Build index to prepare for blastp
subprocess.call("makeblastdb -in ./raw_protein_sequences/proteinseqs.fa -dbtype prot -out ./conservation_level/prodb", shell=True)

#Calculate the sequences that are most similar to the target sequences 
subprocess.call("blastp -db ./conservation_level/prodb -query ./conservation_level/pul_seqs.fa -outfmt 7 -out ./conservation_level/blastp.txt",shell = True)

#Extract information of sequences with protein names
subprocess.call("sed -e '/#/d' ./conservation_level/blastp.txt | cut -f 2 > ./conservation_level/protein_names",shell = True)
n = open("./conservation_level/protein_names").read()
subprocess.call("/localdisk/data/BPSM/Assignment2/pullseq -i ./conservation_level/clustalo_seqs.fa -n ./conservation_level/protein_names > ./conservation_level/pulresult.fa",shell = True)

#Plot the conservation between sequences with plotcon
subprocess.call("plotcon -sequences ./conservation_level/pulresult.fa -winsize 10 -graph svg",shell = True)
os.system("display plotcon.svg")
subprocess.call("mv ./plotcon.svg ./conservation_level/",shell = True)

#Determine the known motifs
#Scan protein sequences with patmatmotifs and extract information of the report
os.mkdir("motifs")
subprocess.call("patmatmotifs -sequence ./conservation_level/pulresult.fa -outfile ./motifs/report",shell = True)

#Define a function to get index of items from a list
def get_index(list=None, item=""):
    return [i for i in range(len(list)) if list[i] == item]

#Put overall hitcount number and motif names into hc_motif.txt
r = open("./motifs/report").read().split()
h = open("./motifs/hc_motif" + ".txt","w")
i = r.index("HitCount:")
h.write("The overall HitCount number is: " +  str(r[i+1]) + "." + "\n")
index_m = get_index(r, "Motif")
number = 1
for m in index_m:
    h.write("Name of motif" + str(number) + " is " + str(r[m+2]) + "." + "\n")
    number += 1
h.close()
h.close()

#Calculates statistics of protein properties
os.mkdir("protein_properties")
open("./conservation_level/pul_seqs.fa").read()
subprocess.call("pepstats -sequence ./conservation_level/pul_seqs.fa -outfile ./protein_properties/pepstats_of_seq",shell = True)

subprocess.call("rm -f ./conservation_level/prodb*",shell = True)
subprocess.call("rm -f ./conservation_level/prepare_for_blastp.txt",shell = True)
subprocess.call("rm -f ./conservation_level/protein_names",shell = True)
subprocess.call("rm -f ./conservation_level/clustalo_seqs.fa",shell = True)
subprocess.call("rm -f ./conservation_level/pul*",shell = True)
print("Congratulations! All tasks have been done.")
