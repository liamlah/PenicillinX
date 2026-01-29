# PenicillinX

A Python Shiny web app that compares the molecular components of common $\beta$-lactam antibiotics with the goal to assist clinicians with antibiotic selection for patients with previous history of reactions.

## Usage

The clinician selects two or more antibiotics based on previous reactions reported by a patient. PenicillinX assesses and highlights common structures between the antibiotics. With this method it may be possible to
determine the structure responsible for a particular patient's previous reactions, and therefore allow the clinician to determine antibiotics that could be administered with a lower risk of reactions.  

## Warning:

This program is currently in development and ***MUST NOT be used in any clinical decision-making under any circumstance***. Validation of the program's outputs against existing empirical data is ongoing, and a number of
bugs are currently present during molecular substructure comparisons.

## Testing

If you wish to give feedback, you can test a live version of the app [here.](https://penicillinX.liamjones.science)

## Dependencies and sources

- Molecular structures were found using [Reaxys](https://www.reaxys.com/). Their structures and substructures are stored using the [SMILES](https://archive.epa.gov/med/med_archive_03/web/html/smiles.html) format.
- PenicillinX relies heavily on the [RDKit](https://www.rdkit.org/) open source Cheminformatics package for Python.
- The [Shiny for Python](https://github.com/posit-dev/py-shiny/) framework was used to build PenicillinX, including [Shiny R](https://github.com/rstudio/shiny) for the early iterations of the program.

## Licence

PenicillinX is licenced under the GNU Affero General Public Licence. See the [LICENCE](LICENSE) file for more details.


