import pandas as pd
import tkinter as tk
import csv
import os
from tkinter import Entry, Button, Label, ttk, filedialog

my_data = []

#fxn for selecting the RE of interest in gui
def select_enzyme():
    global re_of_interest  # Declare re_of_interest as a global variable
    enzyme_name = enzyme_entry.get()
    re_of_interest = enzyme_name

def get_table(list1, list2, list3, list4, list5):
    global my_data
    for i in range(len(list1)):  
        Sequence = list1[i]
        Index = list2[i]
        Mismatch = list3[i]
        AA_Seq = list4[i]
        Silent_Mut_Seq = list5[i]
        data = (Sequence, Index, Mismatch, AA_Seq, Silent_Mut_Seq)
        pandas_to_tree.insert(parent="", index=0, values=data)
        my_data.append(data)
        
def export():
    fln = filedialog.asksaveasfilename(initialdir = os.getcwd(), title= "Save CSV", filetypes=(("CSV File", "*.csv"), ("All Files", "*")))
    with open(fln, mode= 'w') as my_file:
        exp_writer = csv.writer(my_file, delimiter = ",")
        for i in my_data:
            exp_writer.writerow(i)

def clean_dna_sequence(dna_sequence):
    # Remove non-alphabetic characters (e.g., spaces, tabs, etc.)
    global cleaned_sequence
    cleaned_sequence = ''.join(filter(str.isalpha, dna_sequence))
    potential_cutsite(cleaned_sequence)
    

def potential_cutsite(sequence):
    re_dict = {"hindIII": "AGCTT", "ecorI": "GAATTC", "pstI": "CTGCAG", "notI" : "GCGGCCGC", "bamhI": "GGATCC", "pmeI" : "GTTTAAAC", "xbaI" : "TCTAGA", "ecorV" : "GATATC", "ndeI" : "CATATG", "fspI" : "TGGCA", "draI" : "TTTAAA", "bstbI" : "TTCGAA", "bgIII" : "AGATCT"}
    potential_sites_i = []  # List to store potential cut site indices
    potential_sites_seq = []  # List to store the matching sequences
    potential_sites_mismatch = [] #list to store location of mismatched bp
    print ()
    #re_of_interest = input("What restriction enzyme are you interested in using? Note: provide text with NO capitalization OR spacing. If applicable, use roman numerals for numbers (i.e. hindIII is proper vs Hind3) ")
    print()

    recognition_site = re_dict.get(re_of_interest, None)  # Get the recognition site 
    if recognition_site is None: 
            print()
            print("Recognition site not found for the specified enzyme. Add enzyme to library.")
            print()
            return [], []

    window_size = len(recognition_site)  # Set the window size based on the recognition site
     

    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]  # Create a window of variable length
        count = 0


        for j in range(window_size):  # Iterate through the indices of the window and the recognition site simultaneously
            if window[j].upper() == recognition_site[j]:
                count += 1
                
        
        percent_match = (count / (window_size)) * 100
        if percent_match >= 60 and len(window) % 3 == 0:
            potential_sites_i.append(i)  # Append the index where the match is found
            potential_sites_seq.append(window)  # Append the matching sequence

        if percent_match >= 60 and len(window) % 3 != 0:
            if len(window) == 4:
                window = sequence[i:i + window_size + 2]
                potential_sites_i.append(i)  # Append the index where the match is found
                potential_sites_seq.append(window)  # Append the matching sequence
            elif len(window) == 5:
                window = sequence[i:i + window_size + 1]
                potential_sites_i.append(i)  # Append the index where the match is found
                potential_sites_seq.append(window)  # Append the matching sequence
            elif len(window) == 7:
                window = sequence[i:i + window_size + 2]
                potential_sites_i.append(i)  # Append the index where the match is found
                potential_sites_seq.append(window)  # Append the matching sequence
            elif len(window) == 8:
                window = sequence[i:i + window_size + 1]
                potential_sites_i.append(i)  # Append the index where the match is found
                potential_sites_seq.append(window)  # Append the matching sequence

    results_dict = {seq: index for seq, index in zip(potential_sites_seq, potential_sites_i)}
    filtered_results_dict = {}
    for seq, index in results_dict.items():
        if index % 3 == 0:  # Check if the index is evenly divisible by 3
            filtered_results_dict[seq] = index
    
    print()
    for key in filtered_results_dict:
        mismatched_seq = ""
        recognition_site_adjusted = ""
        if len(key) == len(recognition_site_adjusted):
            for j in range(len(key)):
                if key[j].upper() != recognition_site_adjusted[j]:
                    mismatched_seq.append('*')
                else:
                    mismatched_seq += key[j]
        else:
            recognition_site_adjusted = list(recognition_site)  # Convert recognition_site to a list

            while len(key) > len(recognition_site_adjusted):
                recognition_site_adjusted.append('-')

            for y in range(len(recognition_site_adjusted)):
                if key[y].upper() != recognition_site_adjusted[y]:
                    mismatched_seq += "*"
                else:
                    mismatched_seq += key[y]

        potential_sites_mismatch.append(mismatched_seq)
    
    list_index = []
    list_seq = []
    for key, value in filtered_results_dict.items():
        list_seq.append(key)
        list_index.append(value)
    
    codon_dict = {'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'C': ['TGT', 'TGC'], 'W': ['TGG'], 'E': ['GAA', 'GAG'], 'D': ['GAT', 'GAC'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'N': ['AAT', 'AAC'], 'M': ['ATG'], 'K': ['AAA', 'AAG'], 'Y': ['TAT', 'TAC'], 'I': ['ATT', 'ATC', 'ATA'], 'Q': ['CAA', 'CAG'], 'F': ['TTT', 'TTC'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'], '*': ['TAA', 'TAG', 'TGA'], 'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'G': ['GGT', 'GGC', 'GGA', 'GGG'], 'H': ['CAT', 'CAC']}
    translated_aa = []
    
    for i in list_seq:
        codon_1 = i[0:3]
        codon_2 = i[3:6]
        codon_3 = i[6:9] if len(i) > 6 else None
        aa_1 = None
        aa_2 = None
        aa_3 = None
        
        for amino_acid, codons in codon_dict.items():
            if codon_1 in codons:
                aa_1 = amino_acid
            if codon_2 in codons:
                aa_2 = amino_acid
            if codon_3 is not None and codon_3 in codons:
                aa_3 = amino_acid
            if codon_3 == None:
                aa_3 = 'No AA'
        
        # Append the translated amino acids to the list
        translated_aa.append((aa_1, aa_2, aa_3))
   
    #still not translating the second aa correctly. only giving the DNA seq, no aa
    recognition_site_aa = []
    codon_1_rs = recognition_site[0:3]
    codon_2_rs = recognition_site[3:6]

    codon_3_rs = recognition_site[6:9] if len(recognition_site) > 6 else None
    aa_1_rs = None
    aa_2_rs = codon_2_rs
    aa_3_rs = None
       
    for amino_acid, codons in codon_dict.items():
        if codon_1_rs in codons:
            aa_1_rs = amino_acid
        if codon_2_rs in codons:
            aa_2_rs = amino_acid
        if codon_3_rs is not None and codon_3_rs in codons:
            aa_3_rs = amino_acid
        if codon_3_rs is not None and codon_3_rs is not codons:
            aa_3_rs = codon_3_rs
        
        # Append the translated amino acids to the list
    recognition_site_aa.append((aa_1_rs, aa_2_rs, aa_3_rs))
    recog_site_aa = recognition_site_aa[0] #list of all aa in recog. seq
    recog_site_aa_2 = recog_site_aa[1] #only the second aa in recog. seq
    recog_site_aa_3 = recog_site_aa[2]
    #print("AA Residues, Recog Site:", recog_site_aa)
    
    #creates pandas table with all the information from above
    pandas_list = []
    for i in range(len(list_seq)):
            pandas_list.append([list_seq[i], list_index[i], potential_sites_mismatch[i], (translated_aa[i])])
      
    df = pd.DataFrame(pandas_list, columns=["Sequence:", "Index:", "Mismatch:", "AA Seq:"])


    #Filters through df and only keeps the sequences where aa1 seq == aa1 RE recog. site    
    filtered_df = df[df['AA Seq:'].apply(lambda x: x[0]) == recog_site_aa[0]]
    
    codon_dict = {'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'C': ['TGT', 'TGC'], 'W': ['TGG'], 'E': ['GAA', 'GAG'], 'D': ['GAT', 'GAC'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'N': ['AAT', 'AAC'], 'M': ['ATG'], 'K': ['AAA', 'AAG'], 'Y': ['TAT', 'TAC'], 'I': ['ATT', 'ATC', 'ATA'], 'Q': ['CAA', 'CAG'], 'F': ['TTT', 'TTC'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'], '*': ['TAA', 'TAG', 'TGA'], 'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'G': ['GGT', 'GGC', 'GGA', 'GGG'], 'H': ['CAT', 'CAC']}
    
    #Silent mutation 
    final_pandas_list = [] 
    for index, row in filtered_df.iterrows():  # Iterate through DataFrame rows
        second_aa = row['AA Seq:'][1]
        second_codon = row["Sequence:"][3:6]
        third_aa = row['AA Seq:'][2]

        #case w hindIII when second codon is not full
        if third_aa == "No AA" :              
            if second_aa in codon_dict:  # Check if second_aa exists in codon_dict
                for codon in codon_dict[second_aa]:# Iterate through the codons for second_aa
                    if codon.startswith(recog_site_aa_2):  # Check if codon starts with recog_site_aa_2
                        codon_2 = codon
                        codon_1 = recognition_site[0:3]
                        combined_codon = codon_1 + codon_2
                        final_pandas_list.append([row['Sequence:'], row["Index:"], row['Mismatch:'], row["AA Seq:"], combined_codon])
                        break
            if second_aa in codon_dict and len(recog_site_aa_2) < 2:  # Check if second_aa exists in codon_dict
                for codon in codon_dict[recog_site_aa_2]:# Iterate through the codons for second_aa
                        if second_codon == codon:
                            codon_2 = recognition_site[3:6]
                            codon_1 = recognition_site[0:3]
                            combined_codon = codon_1 + codon_2
                            final_pandas_list.append([row['Sequence:'], row["Index:"], row['Mismatch:'], row["AA Seq:"], combined_codon])
                            break

        #case when there is partial 3rd codon
        if third_aa != "No AA":                    #trying to set up a filter
            if third_aa in codon_dict and second_aa == recog_site_aa_2:  
                for codon in codon_dict[third_aa]:
                    if codon.startswith(recog_site_aa_3):
                        codon_3 = codon                       
                        codon_2 = recognition_site[3:6]
                        codon_1 = recognition_site[0:3]
                        combined_codon = codon_1 + codon_2 + codon_3
                        final_pandas_list.append([row['Sequence:'], row["Index:"], row['Mismatch:'], row["AA Seq:"], combined_codon])
                        break
        
        #case when there is a full 3rd codon, needs more work           
        elif third_aa != "No AA":                    #trying to set up a filter
            if third_aa in codon_dict and second_aa == recog_site_aa_2:  
                for codon in codon_dict[third_aa]:
                    if codon == codon_dict[recog_site_aa_3]:
                        codon_3 = codon                       
                        codon_2 = recognition_site[3:6]
                        codon_1 = recognition_site[0:3]
                        combined_codon = codon_1 + codon_2 + codon_3
                        final_pandas_list.append([row['Sequence:'], row["Index:"], row['Mismatch:'], row["AA Seq:"], combined_codon])
                        break

        
# Create a new DataFrame from the list of matching rows
    final_df = pd.DataFrame(final_pandas_list, columns=["Sequence:", "Index:", "Mismatch:", "AA Seq:", "Silent Mut Seq:"])
    
    sequence_list = final_df["Sequence:"].tolist()
    index_list = final_df["Index:"].tolist()
    mismatch_list = final_df["Mismatch:"].tolist()
    aa_seq_list = final_df["AA Seq:"].tolist()
    silent_mut_seq_list = final_df["Silent Mut Seq:"].tolist()
    get_table(sequence_list, index_list, mismatch_list, aa_seq_list,silent_mut_seq_list)
   
    
#___________________________________________________________________________________________________________

root = tk.Tk()
root.title("CodonInsight")
root.geometry("600x600")

#function for passing DNA Seq to clean_dna_sequence(dna_sequence) via button1
def seq_input_button():
    dna_sequence = entry.get()
    clean_dna_sequence(dna_sequence)
    potential_cutsite(cleaned_sequence)



#DNA Seq Frame
dna_frame = tk.Frame(root)
dna_frame.pack()

label = Label(dna_frame, text="Provide DNA Sequence Here")
label.pack()

entry = Entry(dna_frame, width = 100, borderwidth = 2 ) #text field to pass through to program
entry.pack(padx= 30)

button1 = Button(dna_frame, text = "Submit", command=seq_input_button)
button1.pack()

        #__________

# Enzyme Frame
enzyme_frame = tk.Frame(root)
enzyme_frame.pack()

enzyme_label = Label(enzyme_frame, text="Select Restriction Enzyme:") # Add a label for enzyme selection
enzyme_label.pack()

enzyme_entry = Entry(enzyme_frame, width=30, borderwidth=2) # Entry field for enzyme name
enzyme_entry.pack()


enzyme_button = Button(enzyme_frame, text="Submit Enzyme", command=select_enzyme) # Button to submit enzyme selection
enzyme_button.pack()

        #_________

expt_buton = Button(root, text = "Export CSV", command= export)
expt_buton.pack()
                    
#setting up Treeview for pandas df    
pandas_to_tree = ttk.Treeview(root, columns= ("Sequence:", "Index:", "Mismatch:", "AA Seq:", "Silent Mut Seq:") , show="headings")
pandas_to_tree.heading("Sequence:", text= 'Sequence:')
pandas_to_tree.heading("Index:", text = "Index:")
pandas_to_tree.heading("Mismatch:", text = "Mismatch:")
pandas_to_tree.heading("AA Seq:", text = "AA Seq:")
pandas_to_tree.heading("Silent Mut Seq:", text= "Silent Mut Seq:")
pandas_to_tree.pack()



root.mainloop()



