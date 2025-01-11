from rdkit import Chem
from rdkit.Chem import Draw
from io import BytesIO
import base64


def generate_favicon(filepath: str = 'www/PX.mol', img_size: tuple = (192, 192)) -> str:
    """
    Generate a favicon from a MOL file using RDKit with a white background
    and black molecules (default rendering).

    Args:
        filepath (str): Path to the MOL file to generate the molecule structure.
        img_size (tuple): Image size for the favicon (width, height).

    Returns:
        str: Base64 encoded image string for the favicon or an empty string if generation fails. Does not fail gracefully though
    """
    try:
        # Load the molecule from the file
        mol = Chem.MolFromMolFile(filepath)
        if not mol:
            print(f"Failed to load molecule from file: {filepath}")
            return ""

        # Set up MolDraw2DCairo (with custom size)
        width, height = img_size
        drawer = Draw.MolDraw2DCairo(width, height)  # PNG generation
        opts = drawer.drawOptions()
        drawer.drawOptions().useBWAtomPalette()
        # Set white background color
        opts.setBackgroundColour((1.0, 1.0, 1.0, 0.8))  # White background

        # Disable atom labels
        opts.noAtomLabels = False
        opts.bondLineWidth = 5


        # Render the molecule
        drawer.DrawMolecule(mol)  # Default coloring (black atoms/bonds on white)
        drawer.FinishDrawing()

        # Retrieve the PNG image and encode as base64
        png_data = drawer.GetDrawingText()
        base64_data = base64.b64encode(png_data).decode("utf-8")

        return f"data:image/png;base64,{base64_data}"

    except Exception as e:
        print(f"Error generating favicon: {e}")
        return ""
