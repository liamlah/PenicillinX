import os

import pandas as pd
from htmltools import HTML
from shiny import render, reactive, ui
import shinyswatch

# Define the UI
app_ui = ui.page_fluid(
    ui.tags.head(
        # Inserted JavaScript
        ui.tags.title("PenicillinX"),
        ui.tags.script(
            HTML(
                """
            $(document).ready(function() {
                $('#reset_button').click(function() {
                    location.reload();
                });
                $('.instructions-toggle').click(function() {
                    $(this).next('.instructions-text').slideToggle();
                    // Toggle the button text between '> Instructions' and '< Instructions'
                    $(this).text(function(i, text){
                        return text === "> Instructions" ? "< Instructions" : "> Instructions";
                    });
                });
            });
            """
            )
        ),
    ),

    ui.panel_title(
        ui.h1("PenicillinX", align="center"),
    ),
    shinyswatch.theme.superhero(),

    # Instructions section
    ui.div(
        "> Instructions",
        class_="instructions-toggle",
        style="text-align: left; font-weight: bold; cursor: pointer; margin-bottom: 10px;",
    ),
    ui.div(
        "This application allows you to compare different antibiotics based on their chemical structure. Select up to three antibiotics from the dropdown. The application will then analyze the selected antibiotics and display any matching features they share. Eventually this app will allow clinicians to identify the common structures behind a patient's antibiotic allergy, and allow them to safely prescribe alternatives.",
        class_="instructions-text",
        style="display: none; text-align: justify;",
    ),
    ui.layout_sidebar(
        ui.sidebar(
            # Dropdown menu for selecting the first antibiotic with placeholder
            ui.input_select(
                id="antibiotic1",
                label="Select Antibiotic 1:",
                choices=[],  # Will be populated by the csv
                selected="",  # Default no selection
            ),

            # Dropdown menu for selecting the second antibiotic with placeholder
            ui.input_select(
                id="antibiotic2",
                label="Select Antibiotic 2:",
                choices=[],  # Will be populated by the csv
                selected="",  # Default no selection
            ),

            # Dropdown menu for selecting the third antibiotic with placeholder
            ui.input_select(
                id="antibiotic3",
                label="Select Antibiotic 3:",
                choices=[],  # Will be populated by the csv
                selected="",  # Default no selection
            ),
            ui.input_action_button("reset_button", "Reset", style="color: #4e5d6c; background-color: #fff; border: 1px solid #4e5d6c; border-radius: 3px;"),
        ),
        ui.panel_main(
            # Where the antibiotic pictures will be shown
            ui.h3("Selected Antibiotics:"),
            ui.output_ui("antibiotic_images"),

            # Display the matching attributes
            ui.h3("Matching Features:"),
            ui.output_text_verbatim("matches"),

            # If there are common features among 3 antibiotics, then this will show up
            ui.output_ui("common_matches_output")
        )
    ),

    ui.tags.footer(  # Disclaimer and link to github
        ui.p(
            ui.HTML(
                "App currently in development. This should absolutely not be used in any clinical setting. <a href='https://github.com/liamlah/PenicillinX' target='_blank'>More information can be found on the GitHub page</a>"
            )
        ),
        style="position:center; bottom: ; left: 0%; transform: translateX(-0%); width: 100%; height: auto; color: white; padding: 0px; background-color: ; z-index: 100;"
    )
)

# Define the server
def server(input, output, session):
    # Reads from the CSV ---- semi-dummy data at present
    antibiotics_data = pd.read_csv("antibiotics_data.csv")

    # Dropdown with antibiotic choices
    @reactive.Effect
    def _():
        antibiotic_choices = ["Select antibiotic"] + antibiotics_data["Antibiotic"].tolist()
        ui.update_select(session, "antibiotic1", choices=antibiotic_choices, selected="")
        ui.update_select(session, "antibiotic2", choices=antibiotic_choices, selected="")
        ui.update_select(session, "antibiotic3", choices=antibiotic_choices, selected="")

    # User can select 2 or more antibiotics
    @reactive.Calc
    def selected_antibiotics():
        selections = [input.antibiotic1(), input.antibiotic2(), input.antibiotic3()]
        selections = [s for s in selections if s != ""]  # Remove empty selections
        return selections

    # Finds the common features among two antibiotics
    def find_common(antibiotic_a, antibiotic_b):
        row_a = antibiotics_data[antibiotics_data["Antibiotic"] == antibiotic_a].iloc[0]
        row_b = antibiotics_data[antibiotics_data["Antibiotic"] == antibiotic_b].iloc[0]

        side_chain = row_a["SideChain"] if row_a["SideChain"] == row_b["SideChain"] else None
        core_ring = row_a["CoreRing"] if row_a["CoreRing"] == row_b["CoreRing"] else None

        return {"side_chain": side_chain, "core_ring": core_ring}

    # Finds common features among the selected antibiotics
    def find_common_all(selected):
        selected_data = antibiotics_data[antibiotics_data["Antibiotic"].isin(selected)]

        # Find common features using a loop
        common_side_chain = selected_data["SideChain"].iloc[0]
        common_core_ring = selected_data["CoreRing"].iloc[0]
        for i in range(1, len(selected_data)):
            common_side_chain = common_side_chain if common_side_chain in selected_data["SideChain"].iloc[i] else None
            common_core_ring = common_core_ring if common_core_ring in selected_data["CoreRing"].iloc[i] else None

        return {"side_chain": common_side_chain, "core_ring": common_core_ring}

    # Formats the feature matches
    def format_matches(features):
        output_strings = []

        if features["side_chain"] is not None:
            output_strings.append(f"Side Chain: {features['side_chain']}")

        if features["core_ring"] is not None:
            output_strings.append(f"Core Ring: {features['core_ring']}")

        if len(output_strings) > 0:
            return " | ".join(output_strings)
        else:
            return "No matching features."

    # Text for the output of matching features
    @output
    @render.text
    def matches():
        selected = selected_antibiotics()
        num_selected = len(selected)

        if num_selected < 2:
            return "Please select at least two antibiotics for comparison."

        # Ensures someone hasnt picked the same antibiotic twice
        if len(set(selected)) < num_selected:
            return "Please select distinct antibiotics for comparison."

        output_text = ""

        # Compares in pairs
        for i in range(num_selected):
            for j in range(i + 1, num_selected):
                features = find_common(selected[i], selected[j])
                match_info = format_matches(features)
                output_text += f"🔹 Matching features between {selected[i]} and {selected[j]}:\n{match_info}\n\n"

        # output of the result
        return output_text

    # Common features of all three. in separate box for readability and distinctness
    @output
    @render.ui
    def common_matches_output():
        selected = selected_antibiotics()
        num_selected = len(selected)

        # Ensures three antibiotics are selected before doing this for aesthetics
        if num_selected == 3:
            common_features = find_common_all(selected)
            common_features_formatted = format_matches(common_features)

            # Prepares the output for the common features
            if common_features["side_chain"] is not None or common_features["core_ring"] is not None:
                return ui.div(
                    style="margin-top: 20px; padding: 10px; border: 1px solid #ccc; background-color: #f0f8ff; color: #000;",
                    children=[
                        ui.strong(
                            f"🔹 Matching features among all three antibiotics ({', '.join(selected)}):"
                        ),
                        common_features_formatted,
                    ],
                )
            else:
                return ui.div(
                    style="margin-top: 20px; padding: 10px; border: 1px solid #ccc; background-color: #f0f8ff; color: #000;",
                    children=[ui.strong("❗ No matching features among all three antibiotics.")],
                )
        else:
            return None  # No output if not all three antibiotics are selected

    # Images and labels for the antibiotics
    @output
    @render.ui
    def antibiotic_images():
        selected = selected_antibiotics()

        if len(selected) == 0:
            return ui.p("No antibiotics selected.")

        # Create a list to hold image-tag pairs
        image_tags = []
        for antibiotic in selected:
            img_src = f"ImagesInv/{antibiotic}.svg"
            # Placeholder image if I decide to add more antibiotics
            img_path = f"www/Images/{antibiotic}.svg"
            if not os.path.exists(img_path):
                img_src = "Images/placeholder.png"

            image_tags.append(
                ui.div(
                    style="text-align: center; margin: 10px;",
                    children=[
                        ui.img(src=img_src, height="150px", alt=f"Image of {antibiotic}"),
                        ui.tags.br(),
                        ui.strong(antibiotic),
                    ],
                )
            )

        # Arrange images in a fluidRow with columns
        return ui.row_fluid(
            *[ui.column(width=4, children=[tag]) for tag in image_tags]
        )

app = shiny.App(app_ui, server)