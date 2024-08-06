#Importing relevant libraries and packages
import tkinter as tk
from tkinter import ttk
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from matplotlib import pyplot as plt
import io
import os
from Chem_database_Functions import DatabaseManager, TreeManager, FilterParameters
import argparse

#Setting up argparse to parse optional arguments for creating, loading and quering database using flags, and a positional argument to specify the database filename for all functions.
parser=argparse.ArgumentParser(description='Creating chemical databases, and querying data in userfriendly interface.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--database', help="Name of database file.")
parser.add_argument('--ddl', help="Path of SQL DDL statement file.")
parser.add_argument('--output_folder', help="Path of new directory to store output files in.")
parser.add_argument('--sdf', help="Path of molecule SDF file.")
args=parser.parse_args()

#Assigning args to variables
SQL_create = args.ddl
output_folder=args.output_folder
database_name = os.path.join(output_folder, args.database)
sdf=args.sdf

#Creating output directory
if not os.path.isdir(output_folder):
    os.mkdir(output_folder)

#Making database manager object
db=DatabaseManager(database_name)

#Creating database
db.create_db(SQL_create)

#Pulling molecule properties from SDF
parameters= []
molecules= Chem.ForwardSDMolSupplier(sdf)
for mol in molecules:
    mol_name = mol.GetProp("Name")
    atom_number = int(mol.GetNumAtoms())
    mol_weight = Descriptors.ExactMolWt(mol)
    cd_id = mol.GetProp("CdId")
    log_d = float(mol.GetProp("LogD"))
    mol_id = mol.GetProp("Mol_ID")
    mol_formula = mol.GetProp("Formula")
    rotate_bonds = int(Chem.rdMolDescriptors.CalcNumRotatableBonds(mol))
    number_rings = int(mol.GetRingInfo().NumRings())
    smiles = Chem.MolToSmiles(mol)
    log_p = float(Chem.Crippen.MolLogP(mol))
    h_donors = int(Chem.rdMolDescriptors.CalcNumHBD(mol))
    h_acceptors= int(Chem.rdMolDescriptors.CalcNumHBA(mol))
    psa= Chem.rdMolDescriptors.CalcTPSA(mol)            

    #Following 3 lines was suggested on various patforms before looping over the rings
    Chem.AssignStereochemistry(mol, force=True, flagPossibleStereoCenters=True)
    Chem.SanitizeMol(mol)
    Chem.Kekulize(mol)
    
    #Finding aromatic fused rings
    fused_rings = 0
    if number_rings >=2:
        #Looping over ring indexes as the following function don't take the same input
        for index in range(0, number_rings):
            if mol.GetRingInfo().IsRingFused(index):
                if all(mol.GetBondWithIdx(bond).GetIsAromatic() for bond in mol.GetRingInfo().BondRings()[index]):        
                    fused_rings += 1
                
    #Counting fullfilled criteria for bioavailability
    biobool=0
    if mol_weight <= 500:
        biobool += 1
    if log_p <= 5:
        biobool += 1
    if h_donors <= 5:
        biobool += 1
    if h_acceptors <= 10:
        biobool += 1
    if rotate_bonds <= 10:
        biobool += 1
    if psa <= 200:
        biobool += 1
    if fused_rings <=5:
        biobool += 1
    if biobool >= 6:
        bioavailability = "Yes"
    else:
        bioavailability = "No"

    #Checking lipinsky's rule of 5
    if mol_weight <= 500 and log_p <= 5 and h_donors <= 5 and h_acceptors <= 10:
        lipinsky = "Yes"
    else:
        lipinsky = "No"

    #Checking lead-likeness
    if mol_weight <= 450 and log_d >= -4 and log_d <= 4 and number_rings <= 4 and rotate_bonds <= 10 and h_donors <= 5 and h_acceptors <= 8:
        lead_like = "Yes"
    else:
        lead_like = "No"

    #Drawing 2D structures and parsing them into bytes for the SQL database
    image = Draw.MolToImage(mol, size=(150, 150), mols_per_row=2)
    image_binary= io.BytesIO()
    image.save(image_binary, format = 'PNG')
    binary=image_binary.getvalue()

    #Saving the 2D structure in PNG format as well
    image.save(os.path.join(output_folder, f'{cd_id}.png'))
    
    #Making list of everything to go into sql database
    parameters.append([cd_id, smiles, round(mol_weight, 2), round(log_p,2), atom_number, h_donors, rotate_bonds, number_rings, bioavailability, lipinsky, lead_like, binary, h_acceptors, mol_name])

#Inserting data
sql = (f'INSERT INTO Molecule VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)')
db.insert_db(sql, parameters)


## GUI ####
#Creating main frame
window = tk.Tk()
window.configure(bg="#AAD7D9")  
window.minsize(1250, 650)
window.maxsize(1250, 650)
tk.Label(window, text="Database of chemical compounds", font=("Copperplate Gothic Bold", 25), bg="#AAD7D9").pack()
#Creating left frame
left_frame = tk.Frame(window, width=200, height=650, bg="#92C7CF")
left_frame.pack(side="left")
tk.Label(left_frame, text="Querying tools", bg="#92C7CF").pack()

#Creating right frame
right_frame = tk.Frame(window, width=1050, height=650, bg="#92C7CF")
right_frame.pack(side="right")

tk.Label(right_frame, text="Hits", bg="#92C7CF").pack()

database_frame=tk.Frame(right_frame, width=180, height=650, bg="#92C7CF")
database_frame.pack()

# Create The Treeview
my_tree = ttk.Treeview(database_frame, selectmode="extended", height=3)
tree_mang= TreeManager(my_tree)
s= ttk.Style()
s.configure('Treeview', background="#92C7CF", rowheight=150)
# Create a Treeview Scrollbar
tree_scrollx = tk.Scrollbar(database_frame, orient="horizontal", command=my_tree.xview)
tree_scrollx.pack(side="bottom", fill="x")
my_tree.configure(xscrollcommand=tree_scrollx.set)

#Defining treeview columns
my_tree['columns'] = ("CdID", "Name", "Smiles", "Molecular weight", "LogP", "Atom", "H-bond donors", "H-bond acceptors", "Rotatable bonds", "Rings", "Bioavailable", "Lipinsky", "Lead-like")

tree_mang.tree_column("#0", 180)
tree_mang.tree_column("CdID", 40)
tree_mang.tree_column("Name", 100)
tree_mang.tree_column("Smiles", 150)
tree_mang.tree_column("Molecular weight", 100)
tree_mang.tree_column("LogP", 70)
tree_mang.tree_column("Atom", width=60)
tree_mang.tree_column("H-bond donors", 100)
tree_mang.tree_column("H-bond acceptors", 100)
tree_mang.tree_column("Rotatable bonds", 100)
tree_mang.tree_column("Rings", 70)
tree_mang.tree_column("Bioavailable", 100)
tree_mang.tree_column("Lipinsky", 100)
tree_mang.tree_column("Lead-like", 100)

my_tree.heading("#0", text="2D structure", anchor= "center")
for header in my_tree['columns']:
    tree_mang.tree_header(header)

#Content of left frame
tool_bar = tk.Frame(left_frame, width=180, height=650, bg="#AAD7D9")
tool_bar.pack()

tk.Label(tool_bar, text="Parameter", bg="#AAD7D9").grid(column=0, row=1)
tk.Label(tool_bar, text="Minimum", bg="#AAD7D9").grid(column=1, row=1)
tk.Label(tool_bar, text="Maximum", bg="#AAD7D9").grid(column=2, row=1)

#Making a class object to set up filtering entries
filter = FilterParameters(my_tree, tool_bar)
#CdID
cdid_filter = filter.param_label("CdID", 2)
cdid_minentry= filter.param_min(2)
cdid_maxentry= filter.param_max(2)
#Mol_weight
weight_filter = filter.param_label("Molecular weight", 3)
mol_weight_minentry= filter.param_min(3)
mol_weight_maxentry= filter.param_max(3)
#LogP
logp_filter= filter.param_label("Log P", 4)
logp_minentry= filter.param_min(4)
logp_maxentry= filter.param_max(4)
#Atom number
atom_filter=filter.param_label("Atoms", 5)
atom_minentry= filter.param_min(5)
atom_maxentry= filter.param_max(5)
#h bond donors
donors_filter=filter.param_label("H-bond donors", 6)
donors_minentry= filter.param_min(6)
donors_maxentry= filter.param_max(6)
#h bond acceptors
acceptors_filter=filter.param_label("H-bond acceptors", 7)
acceptors_minentry= filter.param_min(7)
acceptors_maxentry= filter.param_max(7)
#rotatable bonds
rotatable_filter=filter.param_label("Rotatable bonds", 8)
rotate_minentry= filter.param_min(8)
rotate_maxentry= filter.param_max(8)
#Ring count
ring_filter=filter.param_label("Rings", 9)
rings_minentry= filter.param_min(9)
rings_maxentry= filter.param_max(9)

entries_list= (cdid_minentry, cdid_maxentry, mol_weight_minentry, mol_weight_maxentry, logp_minentry, logp_maxentry, atom_minentry, atom_maxentry, donors_minentry, donors_maxentry, acceptors_minentry, acceptors_maxentry, rotate_minentry, rotate_maxentry, rings_minentry, rings_maxentry)
filters_list= [[cdid_filter, cdid_minentry, cdid_maxentry], [weight_filter, mol_weight_minentry, mol_weight_maxentry], [logp_filter, logp_minentry, logp_maxentry], [atom_filter, atom_minentry, atom_maxentry], [donors_filter, donors_minentry, donors_maxentry], [acceptors_filter, acceptors_minentry, acceptors_maxentry], [rotatable_filter, rotate_minentry, rotate_maxentry], [ring_filter, rings_minentry, rings_maxentry]]

#Creating cheackboxes
tk.Label(tool_bar, text="Fullfills", bg="#AAD7D9").grid(row=10, column=1, padx=5, pady=5)

#Bioavailability, lipinsky, lead
biocheck=tk.IntVar()
tk.Label(tool_bar, text="Bioavailable", bg="#AAD7D9").grid(row=11, column=0, padx=5, pady=5)
bioavailability_button= tk.Checkbutton(tool_bar, variable=biocheck, onvalue=1, offvalue=0, bg="#AAD7D9").grid(row=12, column=0)

lipinskycheck= tk.IntVar()
tk.Label(tool_bar, text="Lipinsky", bg="#AAD7D9").grid(row=11, column=1, padx=5, pady=5)
lipinski_button= tk.Checkbutton(tool_bar, variable=lipinskycheck, onvalue=1, offvalue=0, bg="#AAD7D9").grid(row=12, column=1)

leadcheck=tk.IntVar()
tk.Label(tool_bar, text="Lead-like", bg="#AAD7D9").grid(row=11, column=2, padx=5, pady=5)
lead_button= tk.Checkbutton(tool_bar, variable=leadcheck, onvalue=1, offvalue=0, bg="#AAD7D9").grid(row=12, column=2)

check_list=(biocheck, lipinskycheck, leadcheck)

#Creating sorting dropdown menus
tk.Label(tool_bar, text="Sort by", bg="#AAD7D9").grid(row=13, column=0, columnspan=3, padx=5, pady=5)

tk.Label(tool_bar, text="Parameter", bg="#AAD7D9").grid(row=14, column=0, padx=5, pady=5)
drop_param_string = tk.StringVar()
drop_param_string.set("CdID")
drop_param = tk.OptionMenu(tool_bar, drop_param_string, *['CdID', 'Molecular weight', 'Log P', 'Atoms', 'H-bond donors', 'H-bond acceptors', 'Rotatable bonds', 'Rings'] )
drop_param.grid(row=14, column=1, padx=20)
drop_param["menu"].configure(bg= "white")

tk.Label(tool_bar, text="Order", bg="#AAD7D9").grid(row=15, column=0, padx=5, pady=5)
drop_order_string = tk.StringVar()
drop_order_string.set("Ascending")
drop_order = tk.OptionMenu(tool_bar, drop_order_string, *["Ascending", "Descending"])
drop_order.grid(row=15, column=1, padx=20)
drop_order["menu"].configure(bg= "white")

#Creating buttons for querying, clearing filtering parameters and resetting database to default
tk.Label(tool_bar, text="Update query", bg="#AAD7D9").grid(row=16, column=0, columnspan=3, padx=5, pady=5)

search_button= tk.Button(tool_bar, text="Query", cursor="heart", command= lambda: tree_mang.query_treeview(entries_list, check_list, drop_param_string, drop_order_string, hits_label, window, database_name), background="white")
search_button.grid(row=17, column=0)
reset_button= tk.Button(tool_bar, text="Reset", cursor="pirate", command= lambda: tree_mang.reset_treeview(database_name, total_hits, entries_list, tsv_entry, check_list, drop_param_string, drop_order_string, hits_label, window), background="white")
reset_button.grid(row=17, column=2)
clear_button= tk.Button(tool_bar, text="Clear", cursor="pirate", command= lambda: tree_mang.clear_entries(entries_list, tsv_entry, check_list, drop_param_string, drop_order_string), background="white")
clear_button.grid(row=17, column=1)

tk.Label(tool_bar, text="Filename", bg="#AAD7D9").grid(row=18, column=0, padx=5, pady=5)
tsv_entry=tk.Entry(tool_bar)
tsv_entry.grid(row=18, column=1, padx=5, pady=5)
tsv_button=tk.Button(tool_bar, text="Export to tsv", command= lambda: tree_mang.export_tsv(my_tree['columns'], output_folder, tsv_entry))
tsv_button.grid(row=18, column=2, padx=5, pady=5)

#Adding labels of total molemules and molecules within query
total_hits = tk.Label(tool_bar, bg="#AAD7D9")
total_hits.grid(row=20, column=1, padx=5, pady=5)

hits_label = tk.Label(tool_bar, bg="#AAD7D9")
hits_label.grid(row=20, column=2, padx=5, pady=5)

#Fetching all the molecules and adding them to the GUI
tree_mang.fetch_whole_database(database_name, total_hits)

#Runs GUI
window.mainloop()
