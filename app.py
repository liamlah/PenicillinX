import pandas as pd
from pathlib import Path

# Shiny imports
from shiny.express import input, render, ui
import shinyswatch
from shiny import reactive

# RDKit imports
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from favicon_utils import generate_favicon

# Favicon generation
favicon_data_uri = generate_favicon()

# Import CSV
abxData = pd.read_csv(Path(__file__).parent / "antibiotics_data.csv")

# Page level options
ui.page_opts(title="PenicillinX",
             window_title="PenicillinX",
             #fillable=True,
             theme=shinyswatch.theme.superhero,
             style="text-align: center;"
             )




with ui.tags.head():
    if favicon_data_uri:
        ui.tags.link(rel="icon", href=favicon_data_uri, type="image/png")
    ui.tags.link(rel="stylesheet",
                 href="https://fonts.googleapis.com/css2?family=Ubuntu:ital,wght@0,300;0,400;0,500;0,700&display=swap")
    ui.tags.link(rel="stylesheet", href="style.css")



#================================================== SIDEBAR ==================================================
with (ui.sidebar(open="always")):
    ui.input_selectize(
        id="abx1",
        label="Select Antibiotic 1:",
        selected="",
        choices=[""] + abxData["Antibiotic"].tolist(),
        multiple=False,
        options={"Placeholder": "CHECK12", "allowEmptyOption": True, "plugins":["clear_button"]})

    ui.input_selectize(
        id="abx2",
        label="Select Antibiotic 2:",
        selected="Select Antibiotic",
        choices=[""] + abxData["Antibiotic"].tolist(),
        options={"Placeholder": "Select Please", "allowEmptyOption": True, "plugins":["clear_button"]})

    ui.input_selectize(
        id="abx3",
        label="Select Antibiotic 3:",
        selected="Select Antibiotic",
        choices=[""] + abxData["Antibiotic"].tolist(),
        options={"Placeholder": "Select Please", "allowEmptyOption": True, "plugins":["clear_button"]})


    with ui.tooltip(id="btn_tooltip", placement="bottom"):
        ui.input_action_button("btn", "Info")
        ui.HTML(
                "This application allows you to compare the chemical structures of various antibiotics to . "
                "Select up to three antibiotics from the dropdown. The application will then analyze the selected antibiotics and display any matching features they share. Eventually, this app will allow clinicians "
                "to identify the common structures behind a patient's antibiotic allergy, and allow them to safely prescribe alternatives. "
                "<strong>NOTE: App currently in development. This should absolutely not be used in any clinical setting.</strong> <a href='https://github.com/liamlah/penicillinX' target='_blank' style='color: #ffffff;'>More information can be found on GitHub</a>")

@reactive.effect # Ensures antibiotic can't be selected twice - remove during RDKIT debugging ********************************************
def update_dropdown_choices():
    # Grabs the antibiotics selected in each dropdown
    selected_abx1 = input.abx1()
    selected_abx2 = input.abx2()
    selected_abx3 = input.abx3()

    # Get the full list of antibiotics
    all_choices = abxData["Antibiotic"].tolist()

    # Exclude selected antibiotics for each dropdown, ensures that current selection is not removed
    choices_for_abx1 = (
        [selected_abx1] + [abx for abx in all_choices if abx != selected_abx2 and abx != selected_abx3]
        if selected_abx1 else
        [""] + [abx for abx in all_choices if abx != selected_abx2 and abx != selected_abx3]
    )
    choices_for_abx2 = (
        [selected_abx2] + [abx for abx in all_choices if abx != selected_abx1 and abx != selected_abx3]
        if selected_abx2 else
        [""] + [abx for abx in all_choices if abx != selected_abx1 and abx != selected_abx3]
    )
    choices_for_abx3 = (
        [selected_abx3] + [abx for abx in all_choices if abx != selected_abx1 and abx != selected_abx2]
        if selected_abx3 else
        [""] + [abx for abx in all_choices if abx != selected_abx1 and abx != selected_abx2]
    )

    # Dynamically update the dropdown choices while maintaining the current selection
    ui.update_selectize("abx1", choices=choices_for_abx1, selected=selected_abx1)
    ui.update_selectize("abx2", choices=choices_for_abx2, selected=selected_abx2)
    ui.update_selectize("abx3", choices=choices_for_abx3, selected=selected_abx3)


    #================================================== MAIN PANEL ==================================================
with ui.tags.div(class_="main-content"):
    with ui.layout_column_wrap(fixed_width=True, width="330px", class_="ubuntu-regular"):
        with ui.card(fill=True, min_height="250px"):
            @render.ui
            def selectedlist1():
                # Get the selected antibiotic from the input
                selected_abx = input.abx1()
                selected_abx2 = input.abx2()
                selected_abx3 = input.abx3()

                if not selected_abx:
                    return ui.HTML("<p>Please select an antibiotic.</p>")


                # Get the SMILES string corresponding to the selected antibiotic
                smiles = abxData.loc[abxData["Antibiotic"] == selected_abx, "SMILESFull"].values
                if len(smiles) == 0:
                    return ui.HTML("<p>SMILES not found for the selected antibiotic.</p>")

                # Turns the SMILES into a MOL with coordinates
                mol = Chem.MolFromSmiles(smiles[0])
                if not mol:
                    return ui.HTML("<p>Invalid SMILES format for the selected antibiotic.</p>")

                # Beta-Lactam structure, not explicitly defined in antibiotics_data.csv
                lactamRing = Chem.MolFromSmarts("O=[#6]-1-[#6]-[#6]-[#7]-1")
                lactam_atoms = list(mol.GetSubstructMatch(lactamRing))  # Lactam atom indices

                # RDKit commands to generate image
                drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)  # Set canvas size
                opts = drawer.drawOptions()
                opts.setBackgroundColour((1, 1, 1, 0))  # Transparent background (RGBA)

                # Makes the molecule white to match the dark theme
                atom_palette = {atom.GetIdx(): (1.0, 1.0, 1.0) for atom in mol.GetAtoms()}
                opts.setAtomPalette(atom_palette)
                opts.fontFile = "www/fonts/Ubuntu-Regular.ttf"
                opts.setLegendColour((1, 1, 1))  # White legend color
                opts.legendFontSize = 20

                # Highlight dictionaries
                highlight_colors = {}  # Atom highlight colors
                highlight_atoms = []  # List of atoms to highlight

                messages = []  # Collects the messages

                # For blue Beta-Lactam ring
                if lactam_atoms:
                    highlight_atoms.extend(lactam_atoms)
                    for atom_idx in lactam_atoms:
                        highlight_colors[atom_idx] = (0.4, 0.4, 1, 0.5)  # Blue RGBA
                    messages.append(f"• {selected_abx} contains a <mark style='background-color:#454eb4; color: white;'>Beta-Lactam ring</mark>.")

                # 2. Substructure comparison with other abx
                for other_abx in [selected_abx2, selected_abx3]:
                    if other_abx:
                        substructure_smiles = abxData.loc[
                            abxData["Antibiotic"] == other_abx,
                            ["SMILESR1", "SMILESR2", "SMILESR3", "SMILESRing"]
                        ]

                        for col_name, sub_smiles in substructure_smiles.iloc[0].items():
                            # Skip invalid entries where side chains don't exist
                            if pd.isna(sub_smiles):
                                continue

                            substructure = Chem.MolFromSmiles(sub_smiles)
                            if substructure:
                                matching_atoms = list(mol.GetSubstructMatch(substructure))
                                if matching_atoms:
                                    highlight_atoms.extend(matching_atoms)

                                    if col_name == "SMILESRing":  # Specifically ring match
                                        for atom_idx in matching_atoms:
                                            highlight_colors[atom_idx] = (0.0, 0.8, 0.0, 0.5)  # Green RGBA
                                        messages.append(f"• {selected_abx} shares a similar <mark style='background-color:#08791c; color: white;'>core</mark> with {other_abx}.")
                                    else:  # R1, R2, or R3 match
                                        for atom_idx in matching_atoms:
                                            highlight_colors[atom_idx] = (0.9, 0.0, 0.0, 0.5)  # Red RGBA
                                        messages.append(f"• {selected_abx} shares a similar <mark style='background-color:#940e15; color: white;'>side chain</mark> with {other_abx}.")

                # Puts the highlights on the molecule
                drawer.DrawMolecule(
                    mol,
                    highlightAtoms=highlight_atoms,
                    highlightAtomColors=highlight_colors,
                    legend=selected_abx
                )
                drawer.FinishDrawing()

                # Extract and clean up the SVG
                svg = drawer.GetDrawingText()

                # Combine the SVG and messages
                message_html = "".join(f"<p>{msg}</p>" for msg in messages)
                return ui.HTML(svg + message_html)


        with ui.card(fill=True, min_height="250px"):
            @render.ui
            def selectedlist2():
                # Get the selected antibiotic from the input
                selected_abx = input.abx2()
                selected_abx1 = input.abx1()
                selected_abx3 = input.abx3()

                if not selected_abx:
                    return ui.HTML("<p>""</p>")  # No output if selection empty

                # Get the SMILES string corresponding to the selected antibiotic
                smiles = abxData.loc[abxData["Antibiotic"] == selected_abx, "SMILESFull"].values
                if len(smiles) == 0:
                    return ui.HTML("<p>SMILES not found for the selected antibiotic.</p>")

                # Turns the SMILES into a MOL with coordinates
                mol = Chem.MolFromSmiles(smiles[0])
                if not mol:
                    return ui.HTML("<p>Invalid SMILES format for the selected antibiotic.</p>")

                # Beta-lactam structure, not explicitly defined in antibiotics_data.csv
                lactamRing = Chem.MolFromSmarts("O=[#6]-1-[#6]-[#6]-[#7]-1")
                lactam_atoms = list(mol.GetSubstructMatch(lactamRing))  # Lactam atom indices

                # RDKit commands to generate image
                drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)  # Set canvas size
                opts = drawer.drawOptions()
                opts.setBackgroundColour((1, 1, 1, 0))  # Transparent background (RGBA)

                # Makes the molecule white to match the dark theme
                atom_palette = {atom.GetIdx(): (1.0, 1.0, 1.0) for atom in mol.GetAtoms()}
                opts.setAtomPalette(atom_palette)
                opts.fontFile = "www/fonts/Ubuntu-Regular.ttf"
                opts.setLegendColour((1, 1, 1))  # White legend color
                opts.legendFontSize = 20

                # Highlight dictionaries
                highlight_colors = {}  # Atom highlight colors
                highlight_atoms = []  # List of atoms to highlight

                messages = []  # Collects the messages

                # For blue Beta-Lactam ring
                if lactam_atoms:
                    highlight_atoms.extend(lactam_atoms)
                    for atom_idx in lactam_atoms:
                        highlight_colors[atom_idx] = (0.4, 0.4, 1, 0.5)  # Blue RGBA
                    messages.append(f"• {selected_abx} contains a <mark style='background-color:#454eb4; color: white;'>Beta-Lactam ring</mark>.")

                # 2. Substructure comparison with other abx
                for other_abx in [selected_abx1, selected_abx3]:
                    if other_abx:
                        substructure_smiles = abxData.loc[
                            abxData["Antibiotic"] == other_abx,
                            ["SMILESR1", "SMILESR2", "SMILESR3", "SMILESRing"]
                        ]

                        for col_name, sub_smiles in substructure_smiles.iloc[0].items():
                            # Skip invalid entries where side chains don't exist
                            if pd.isna(sub_smiles):
                                continue

                            substructure = Chem.MolFromSmiles(sub_smiles)
                            if substructure:
                                matching_atoms = list(mol.GetSubstructMatch(substructure))
                                if matching_atoms:
                                    highlight_atoms.extend(matching_atoms)

                                    if col_name == "SMILESRing":  # Specifically ring match
                                        for atom_idx in matching_atoms:
                                            highlight_colors[atom_idx] = (0.0, 0.8, 0.0, 0.5)  # Green RGBA
                                        messages.append(f"• {selected_abx} shares a similar <mark style='background-color:#08791c; color: white;'>core</mark> with {other_abx}.")
                                    else:  # R1, R2, or R3 match
                                        for atom_idx in matching_atoms:
                                            highlight_colors[atom_idx] = (0.9, 0.0, 0.0, 0.5)  # Red RGBA
                                        messages.append(f"• {selected_abx} shares a similar <mark style='background-color:#940e15; color: white;'>side chain</mark> with {other_abx}.")

                # Puts the highlights on the molecule
                drawer.DrawMolecule(
                    mol,
                    highlightAtoms=highlight_atoms,
                    highlightAtomColors=highlight_colors,
                    legend=selected_abx
                )
                drawer.FinishDrawing()

                # Extract and clean up the SVG
                svg = drawer.GetDrawingText()

                # Combine the SVG and messages
                message_html = "".join(f"<p>{msg}</p>" for msg in messages)
                return ui.HTML(svg + message_html)

        with ui.card(fill=True, min_height="250px"):
            @render.ui
            def selectedlist3():
                # Get the selected antibiotic from the input
                selected_abx = input.abx3()
                selected_abx1 = input.abx1()
                selected_abx2 = input.abx2()

                if not selected_abx:
                    return ui.HTML("<p>""</p>")  # No output if selection empty

                # Get the SMILES string corresponding to the selected antibiotic
                smiles = abxData.loc[abxData["Antibiotic"] == selected_abx, "SMILESFull"].values
                if len(smiles) == 0:
                    return ui.HTML("<p>SMILES not found for the selected antibiotic.</p>")

                # Turns the SMILES into a MOL with coordinates
                mol = Chem.MolFromSmiles(smiles[0])
                if not mol:
                    return ui.HTML("<p>Invalid SMILES format for the selected antibiotic.</p>")

                # Beta lactam structure, not explicitly defined in antibiotics_data.csv
                lactamRing = Chem.MolFromSmarts("O=[#6]-1-[#6]-[#6]-[#7]-1")
                lactam_atoms = list(mol.GetSubstructMatch(lactamRing))  # Lactam atom indices

                # RDKit commands to generate image
                drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)  # Set canvas size
                opts = drawer.drawOptions()
                opts.setBackgroundColour((1, 1, 1, 0))  # Transparent background (RGBA)
                opts.fontFile= "www/fonts/Ubuntu-Regular.ttf"
                # Makes the molecule white to match the dark theme
                atom_palette = {atom.GetIdx(): (1.0, 1.0, 1.0) for atom in mol.GetAtoms()}
                opts.setAtomPalette(atom_palette)
                opts.setLegendColour((1, 1, 1))  # White legend color
                opts.legendFontSize = 20

                # Highlight dictionaries
                highlight_colors = {}  # Atom highlight colors
                highlight_atoms = []  # List of atoms to highlight

                messages = []  # Collects the messages

                # For blue Beta-Lactam ring
                if lactam_atoms:
                    highlight_atoms.extend(lactam_atoms)
                    for atom_idx in lactam_atoms:
                        highlight_colors[atom_idx] = (0.4, 0.4, 1, 0.5)  # Blue RGBA
                    messages.append(f"• {selected_abx} contains a <mark style='background-color:#454eb4; color: white;'>Beta-Lactam ring</mark>.")

                # 2. Substructure comparison with other abx
                for other_abx in [selected_abx1, selected_abx2]:
                    if other_abx:
                        substructure_smiles = abxData.loc[
                            abxData["Antibiotic"] == other_abx,
                            ["SMILESR1", "SMILESR2", "SMILESR3", "SMILESRing"]
                        ]

                        for col_name, sub_smiles in substructure_smiles.iloc[0].items():
                            # Skip invalid entries where side chains don't exist
                            if pd.isna(sub_smiles):
                                continue

                            substructure = Chem.MolFromSmiles(sub_smiles)
                            if substructure:
                                matching_atoms = list(mol.GetSubstructMatch(substructure))
                                if matching_atoms:
                                    highlight_atoms.extend(matching_atoms)

                                    if col_name == "SMILESRing":  # Specifically ring match
                                        for atom_idx in matching_atoms:
                                            highlight_colors[atom_idx] = (0.0, 0.8, 0.0, 0.5)  # Green RGBA
                                        messages.append(f"• {selected_abx} shares a similar <mark style='background-color:#08791c; color: white;'>core</mark> with {other_abx}.")
                                    else:  # R1, R2, or R3 match
                                        for atom_idx in matching_atoms:
                                            highlight_colors[atom_idx] = (0.9, 0.0, 0.0, 0.5)  # Red RGBA
                                        messages.append(f"• {selected_abx} shares a similar <mark style='background-color:#940e15; color: white;'>side chain</mark> with {other_abx}.")

                # Puts the highlights on the molecule
                drawer.DrawMolecule(
                    mol,
                    highlightAtoms=highlight_atoms,
                    highlightAtomColors=highlight_colors,
                    legend=selected_abx
                )
                drawer.FinishDrawing()

                # Extract and cleans up the SVG
                svg = drawer.GetDrawingText()

                # Combines SVG and text matches
                message_html = "".join(f"<p>{msg}</p>" for msg in messages)
                return ui.HTML(svg + message_html)









# todo
# Improve small screen and mobile layount: sidebar layout - on top in mobile
# New text box at bottom with summary of matches
# refactor 3 functions to one if possible

# Bugs to fix
# Fix one way matches. e.g ampicillin penicillinG
# Fix repeating match output for different molecules e.g amox + cefalex + ampicillin
# Vancomycin rendering is borked
