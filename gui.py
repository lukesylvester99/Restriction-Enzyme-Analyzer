import tkinter as tk
from tkinter import Entry, Button, Label, ttk

#fxn for selecting the RE of interest in gui
def select_enzyme():
    global re_of_interest  # Declare re_of_interest as a global variable
    enzyme_name = enzyme_entry.get()
    re_of_interest = enzyme_name


root = tk.Tk()
root.title("CodonInsight")
root.geometry("600x600")

#function for passing DNA Seq to clean_dna_sequence(dna_sequence) via button1
def seq_input_button():
    dna_sequence = entry.get()
    clean_dna_sequence(dna_sequence)




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

#setting up Treeview for pandas df    
#pandas_to_tree = ttk.Treeview(root, columns=list(final_df.columns), show="headings")
#result_tree.pack()

#or col in pandas_to_tree.columns:
       # pandas_to_tree.heading(col, text=col)
        #pandas_to_tree.column(col, width=100)

#for i, row in final_df.iterrows():
        #pandas_to_tree.insert("", "end", values=row.tolist())    

root.mainloop()

