#import seaborn as sns

# Import data from shared.py
#from shared import df
import pandas as pd
from pathlib import Path
from shiny.express import input, render, ui
import shinyswatch

#RDKit imports
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import base64
from io import BytesIO



#def abxdata():
 #   infile = Path(__file__).parent / "anitbioticsdata.csv"
  #  return pd.read_csv(infile)

abxdata = pd.read_csv(Path(__file__).parent / "antibiotics_data.csv") #imports csv



ui.page_opts(title="PenicillinX!",
             window_title="PynicillinX",
             fillable=True,
             theme=shinyswatch.theme.superhero,
             style="text-align: center;"
             )



with ui.tags.head(): # Edits the title to centre it
    ui.tags.style("""
        .navbar-brand {
            text-align: center;
            font-size: 38px;
            flex-grow: 1;
            
        }
        .navbar-header {
            display: flex;
            justify-content: center;
            width: 100%;
        }
    """)




with (ui.sidebar(open="always")):
    ui.input_selectize(
        id="abx1",
        label="Select Antibiotic 1:",
        selected="",
        choices=[""] + abxdata["Antibiotic"].tolist(),
        multiple=False,
        options={"Placeholder": "CHECK12", "allowEmptyOption": True})

    ui.input_selectize(
        id="abx2",
        label="Select Antibiotic 2:",
        selected="Select Antibiotic",
        choices=[""] + abxdata["Antibiotic"].tolist(),
        options={"Placeholder": "Select Please", "allowEmptyOption": False})

    ui.input_selectize(
        id="abx3",
        label="Select Antibiotic 3:",
        selected="Select Antibiotic",
        choices=[""] + abxdata["Antibiotic"].tolist(),
        options={"Placeholder": "Select Please", "allowEmptyOption": False})


#Output panel
@render.ui
def selectedlist1():
    # Get the selected antibiotic from the input
    selected_abx = input.abx1()
    if not selected_abx:
        return ui.HTML("<p>Please select an antibiotic.</p>")  # Handle empty selection

    # Fetch the SMILES string corresponding to the selected antibiotic
    smiles = abxdata.loc[abxdata["Antibiotic"] == selected_abx, "SMILESFull"].values
    if len(smiles) == 0:
        return ui.HTML("<p>SMILES not found for the selected antibiotic.</p>")

    # Generate the molecule from the retrieved SMILES
    mol = Chem.MolFromSmiles(smiles[0])
    if not mol:
        return ui.HTML("<p>Invalid SMILES format for the selected antibiotic.</p>")
    # Generate the molecule image
    svg = Draw.MolsToGridImage(
        [mol],
        molsPerRow=1,
        subImgSize=(300, 300),
        legends=[selected_abx],
        useSVG=True,
        )
    return ui.HTML(svg)

@render.ui
def selectedlist2():
    # Get the selected antibiotic from the input
    selected_abx = input.abx2()
    if not selected_abx:
        return ui.HTML("<p>Please select an antibiotic.</p>")  # Handle empty selection

    # Fetch the SMILES string corresponding to the selected antibiotic
    smiles = abxdata.loc[abxdata["Antibiotic"] == selected_abx, "SMILESFull"].values
    if len(smiles) == 0:
        return ui.HTML("<p>SMILES not found for the selected antibiotic.</p>")

    # Generate the molecule from the retrieved SMILES
    mol = Chem.MolFromSmiles(smiles[0])
    if not mol:
        return ui.HTML("<p>Invalid SMILES format for the selected antibiotic.</p>")
    # Generate the molecule image
    svg = Draw.MolsToGridImage(
        [mol],
        molsPerRow=1,
        subImgSize=(300, 300),
        legends=[selected_abx],
        useSVG=True,
        )
    return ui.HTML(svg)

@render.ui
def selectedlist3():
    # Get the selected antibiotic from the input
    selected_abx = input.abx3()
    if not selected_abx:
        return ui.HTML("<p>Please select an antibiotic.</p>")  # Handle empty selection

    # Fetch the SMILES string corresponding to the selected antibiotic
    smiles = abxdata.loc[abxdata["Antibiotic"] == selected_abx, "SMILESFull"].values
    if len(smiles) == 0:
        return ui.HTML("<p>SMILES not found for the selected antibiotic.</p>")

    # Generate the molecule from the retrieved SMILES
    mol = Chem.MolFromSmiles(smiles[0])
    if not mol:
        return ui.HTML("<p>Invalid SMILES format for the selected antibiotic.</p>")
    # Generate the molecule image
    svg = Draw.MolsToGridImage(
        [mol],
        molsPerRow=1,
        subImgSize=(300, 300),
        legends=[selected_abx],
        useSVG=True,
        )
    return ui.HTML(svg)