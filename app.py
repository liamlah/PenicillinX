#import seaborn as sns

# Import data from shared.py
#from shared import df
import pandas as pd
from pathlib import Path
from shiny.express import input, render, ui
import shinyswatch


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
ui.markdown("Selected Antibiotics:")
@render.code
def selectedlist1():
    return input.abx1()
@render.code
def selectedlist2():
    return input.abx2()
@render.code
def selectedlist3():
    return input.abx3()