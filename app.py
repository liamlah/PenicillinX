# import seaborn as sns

# Import data from shared.py
# from shared import df
import pandas as pd
from pathlib import Path

from rdkit.Chem.Draw.IPythonConsole import drawOptions
from shiny.express import input, render, ui
import shinyswatch

# RDKit imports
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D



# def abxdata():
#   infile = Path(__file__).parent / "anitbioticsdata.csv"
#  return pd.read_csv(infile)

abxdata = pd.read_csv(Path(__file__).parent / "antibiotics_data.csv")  # imports csv

ui.page_opts(title="PenicillinX!",
             window_title="PynicillinX",
             fillable=True,
             theme=shinyswatch.theme.superhero,
             style="text-align: center;"
             )

with ui.tags.head():  # Edits the title to centre it
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

with ui.layout_column_wrap(fixed_width=False, width="330px"):

    with ui.card(fill=False, min_height="300px"):
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
            # RDKit Mol drawing with transparent background
            drawer = rdMolDraw2D.MolDraw2DSVG(300, 300, )  # Set canvas size
            opts = drawer.drawOptions()
            opts.setBackgroundColour((1, 1, 1, 0))  # Transparent background (RGBA)
            atom_palette = {atom.GetIdx(): (1.0, 1.0, 1.0) for atom in mol.GetAtoms()}
            opts.setAtomPalette(atom_palette)
            opts.setLegendColour((1,1,1))
            opts.legendFontSize = 20
            drawer.DrawMolecule(mol, legend=selected_abx)
            drawer.FinishDrawing()

            # Extract and clean up the SVG
            svg = drawer.GetDrawingText()
            return ui.HTML(svg)




    with ui.card(fill=False):
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
            # RDKit Mol drawing with transparent background
            drawer = rdMolDraw2D.MolDraw2DSVG(300, 300, )  # Set canvas size
            opts = drawer.drawOptions()
            opts.setBackgroundColour((1, 1, 1, 0))  # Transparent background (RGBA)
            atom_palette = {atom.GetIdx(): (1.0, 1.0, 1.0) for atom in mol.GetAtoms()}
            opts.setAtomPalette(atom_palette)
            opts.setLegendColour((1, 1, 1))
            opts.legendFontSize = 20
            drawer.DrawMolecule(mol, legend=selected_abx)
            drawer.FinishDrawing()

            # Extract and clean up the SVG
            svg = drawer.GetDrawingText()
            return ui.HTML(svg)


    with ui.card(fill=False):
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
            # RDKit Mol drawing with transparent background
            drawer = rdMolDraw2D.MolDraw2DSVG(300, 300, )  # Set canvas size
            opts = drawer.drawOptions()
            opts.setBackgroundColour((1, 1, 1, 0))  # Transparent background (RGBA)
            atom_palette = {atom.GetIdx(): (1.0, 1.0, 1.0) for atom in mol.GetAtoms()}
            opts.setAtomPalette(atom_palette)
            opts.setLegendColour((1, 1, 1))
            opts.legendFontSize = 20
            drawer.DrawMolecule(mol, legend=selected_abx)
            drawer.FinishDrawing()

            # Extract and clean up the SVG
            svg = drawer.GetDrawingText()
            return ui.HTML(svg)

# todo
# kekerize
# transparent background SVG
# SVG alignment


