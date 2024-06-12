from flask import Flask, render_template, request, redirect, url_for, session, jsonify
from flask_sqlalchemy import SQLAlchemy
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolTransforms
import pubchempy as pcp
import pandas as pd

app = Flask(__name__)

# Set a secret key for session
app.secret_key = b'?\xaf\x13\xca+`\x82\xd4\xc7\x97\x0ek\x1f:\x99Q1\x9e\x88?\x98\xf4x\t'
# start of mol predict
def get_smiles_from_formula(formula):
    compounds = pcp.get_compounds(formula, 'formula')
    if compounds:
        return compounds[0].canonical_smiles
    else:
        return None

def process_molecular_formula(formula):
    smiles = get_smiles_from_formula(formula)
    if smiles is None:
        return {"Error": "Failed to find a compound with the given molecular formula"}
    return process_smiles(smiles)

def process_smiles(smiles):
    if isinstance(smiles, float):
        return {
            "SMILES": "",
            "Molecule Name": "",
            "Molecular Formula": "",
            "Bond Information": "",
            "Functional Groups": "",
            "Atom Information": "",
            "Symmetry": "",
            "Planarity": "",
            "Isomers": "",
            "Molecular weight": "",
            "PiPi Stacking Potential": "",
            "ESP": "",
            "oxygen_balance": "",
            "NO2 Count": "",
            "NH2 Count": "",
            "H Count": "",
            "Graphite-like Structure": "",
            "Specific bond": "",
            "Explosive groups": "",
            "Energetics": ""
        }

    molecule_name = get_molecule_name(smiles)

    molecule = Chem.MolFromSmiles(smiles)
    if molecule:
        bonds = molecule.GetBonds()
        bond_info = []
        if bonds:
            for bond in bonds:
                bond_info.append(f"{bond.GetBeginAtomIdx()}-{bond.GetEndAtomIdx()} ({bond.GetBondType()})")
        else:
            bond_info.append("No bonds")

        functional_groups = get_functional_groups(molecule)

        atoms = molecule.GetAtoms()
        atom_info = [f"{atom.GetIdx()}: {atom.GetSymbol()}" for atom in atoms]

        molecular_formula = Chem.rdMolDescriptors.CalcMolFormula(molecule)

        symmetry = check_symmetry(smiles)
        planarity = check_planarity(smiles)
        isomers = generate_isomers(smiles)
        homo_energy = homo_e(smiles)
        pipi_stacking_potential = pipi_stacking(smiles)
        esp_charge = calculate_esp(smiles)
        oxygen_balance = calculate_oxygen_balance(smiles)
        substituent_counts = count_specific_substituents(smiles)
        graphite_like_structure = has_graphite_like_structure(smiles)
        specific_bond = has_specific_bonds(smiles)
        exp_groups = check_for_explosive_groups(smiles)

        energetics = (symmetry == 1) & (planarity == 1) & (oxygen_balance >= -80) & (oxygen_balance <= -40) & (exp_groups == 1) & (specific_bond == 1)

        return {
            "SMILES": smiles,
            "Molecule Name": molecule_name,
            "Molecular Formula": molecular_formula,
            "Bond Information": ", ".join(bond_info),
            "Functional Groups": ", ".join(functional_groups),
            "Atom Information": ", ".join(atom_info),
            "Symmetry": symmetry,
            "Planarity": planarity,
            "Isomers": isomers,
            "Molecular weight": Descriptors.MolWt(molecule),
            "PiPi Stacking Potential": pipi_stacking_potential,
            "ESP": esp_charge,
            "oxygen_balance": oxygen_balance,
            "NO2 Count": substituent_counts["[N+](=O)[O-]"],
            "NH2 Count": substituent_counts["[NH2]"],
            "H Count": substituent_counts["[H]"],
            "Graphite-like Structure": graphite_like_structure,
            "Specific bond": specific_bond,
            "Explosive groups": exp_groups,
            "Energetics": int(energetics)
        }
    else:
        return {
            "SMILES": smiles,
            "Molecule Name": "Failed to create molecule",
            "Molecular Formula": "",
            "Bond Information": "",
            "Functional Groups": "",
            "Atom Information": "",
            "Symmetry": "",
            "Planarity": "",
            "Isomers": "",
            "Molecular weight": "",
            "PiPi Stacking Potential": "",
            "ESP": "",
            "oxygen_balance": "",
            "NO2 Count": "",
            "NH2 Count": "",
            "H Count": "",
            "Graphite-like Structure": False,
            "Specific bond": "",
            "Explosive groups": "",
            "Energetics": ""
        }


# end of mol predict2


# Database configuration
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///app.db'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False

# Initialize SQLAlchemy
db = SQLAlchemy(app)

# Define User model
class User(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(100), unique=True, nullable=False)
    email = db.Column(db.String(100), unique=True, nullable=False)
    password = db.Column(db.String(100), nullable=False)

# Ensure db.create_all() is executed within the application context
with app.app_context():
    db.create_all()

# Route for login
@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        username = request.form['username']
        session['username'] = username
        return redirect(url_for('home'))
    return render_template('index.html')

# Route for registration
@app.route('/register', methods=['GET', 'POST'])
def register():
    if request.method == 'POST':
        username = request.form['username']
        email = request.form['email']
        password = request.form['password']
        confirmPassword = request.form['confirmPassword']

        if password != confirmPassword:
            return "Passwords do not match. Please try again."

        existing_user = User.query.filter_by(username=username).first()
        existing_email = User.query.filter_by(email=email).first()
        if existing_user or existing_email:
            return "Username or email already exists. Please choose another one."

        new_user = User(username=username, email=email, password=password)
        db.session.add(new_user)
        db.session.commit()

        return redirect(url_for('index'))

    return render_template('register.html')

# Route to display all registered users
@app.route('/users')
def users():
    users = User.query.all()
    return render_template('users.html', users=users)

@app.route('/home')
def home():
    username = session.get('username')
    return render_template('home.html', username=username)

@app.route('/cards')
def cards():
    return render_template('cards.html')

@app.route('/tab2')
def tab2():
    return render_template('tab2.html')

@app.route('/tab3')
def tab3():
    return render_template('tab3.html')

@app.route('/predict1', methods=['POST'])
def single_smile():
    if request.method == 'POST':
        data = request.json
        smiles_input = data.get('smiles', '')

        result = process_smiles(smiles_input)
        return jsonify({'result': result})

def process_smiles(smiles):
    # Example processing: reversing the SMILES string and adding mock data
    processed_result = {
        'Original SMILES': smiles,
        'Processed SMILES': smiles[::-1],  # Reverse the SMILES string as an example
        'Additional Info': 'Sample info'  # Replace with actual processing details
    }
    return processed_result


@app.route('/predict2', methods=['POST'])
def predict_molecular():
    data = request.json
    molecular_formula = data.get('formula', '')
    result = process_molecular_formula(molecular_formula)
    return jsonify({'result': result})




@app.route('/pubchem')
def pubchem():
    return render_template('pubchem.html')

@app.route('/tecompounds')
def tecompounds():
    return render_template('tecompounds.html')

@app.route('/ffcompounds')
def ffcompounds():
    return render_template('ffcompounds.html')

def get_molecule_name(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule:
        return Chem.MolToSmiles(molecule)
    else:
        return "Invalid SMILES"

def get_functional_groups(molecule):
    functional_groups = []

    if molecule.HasSubstructMatch(Chem.MolFromSmarts('[C,c]-[C,c]')):
        functional_groups.append("Alkane")

    if molecule.HasSubstructMatch(Chem.MolFromSmarts('[C,c]=C')):
        functional_groups.append("Alkene")

    if molecule.HasSubstructMatch(Chem.MolFromSmarts('[C,c]#C')):
        functional_groups.append("Alkyne")

    if molecule.HasSubstructMatch(Chem.MolFromSmarts('a')):
        functional_groups.append("Aromatic")

    return functional_groups

def check_planarity(smiles):
    molecule = Chem.MolFromSmiles(smiles)

    if molecule is None:
        print("Invalid Smiles input.")
        return

    molecule = Chem.AddHs(molecule)
    AllChem.EmbedMolecule(molecule, randomSeed=42)

    has_non_planar_rings = False
    for ring in molecule.GetRingInfo().AtomRings():
        if len(ring) > 3:
            dihedral_angles = []
            for i in range(len(ring) - 3):
                dihedral = rdMolTransforms.GetDihedralDeg(
                    molecule.GetConformer(), ring[i], ring[i + 1], ring[i + 2], ring[i + 3]
                )
                dihedral_angles.append(dihedral)
            planarity_threshold = 10.0
            has_non_planar_rings |= any(
                abs(angle) > planarity_threshold for angle in dihedral_angles
            )

    if has_non_planar_rings:
        return 0
    else:
        return 1

def check_symmetry(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "Invalid SMILES string"

    AllChem.Compute2DCoords(mol)
    bond = mol.GetBondBetweenAtoms(0, 1)
    if bond is None:
        return "Symmetric"

    bond.SetStereo(Chem.BondStereo.STEREOE)
    rotated_smiles = Chem.MolToSmiles(mol)

    if rotated_smiles == smiles:
        return 1
    else:
        return 0

def generate_isomers(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")

        isomers = []

        isomers.extend(generate_structural_isomers(mol))
        isomers.extend(generate_stereoisomers(mol))

        return isomers
    except Exception as e:
        raise ValueError(f"Error generating isomers: {e}")

def generate_structural_isomers(mol):
    isomers = []
    isomeric_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    isomers.append((isomeric_smiles, "Structural"))
    return isomers

def generate_stereoisomers(mol):
    isomers = []
    stereo_mols = AllChem.EnumerateStereoisomers(mol)
    for stereo_mol in stereo_mols:
        isomers.append((Chem.MolToSmiles(stereo_mol, isomericSmiles=True), "Stereoisomer"))
    return isomers

def homo_e(smiles):
    mol = Chem.MolFromSmiles(smiles)
    homo_energy = Descriptors.MolWt(mol)
    return homo_energy

def pipi_stacking(smiles):
    mol = Chem.MolFromSmiles(smiles)
    aromatic_rings = Chem.GetSSSR(mol)
    homo_energy = Descriptors.MolWt(mol)

    pipi_stacking_potential = len(aromatic_rings) * homo_energy
    return pipi_stacking_potential

def count_specific_substituents(smiles):
    substituents = ["[N+](=O)[O-]", "[NH2]", "[H]"]
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")

    counts = {substituent: 0 for substituent in substituents}
    for substituent in substituents:
        matches = mol.GetSubstructMatches(Chem.MolFromSmarts(substituent))
        counts[substituent] = len(matches)

    return counts

def calculate_oxygen_balance(smiles):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    oxygen_atoms = sum([atom.GetAtomicNum() == 8 for atom in mol.GetAtoms()])
    molecular_weight = Descriptors.MolWt(mol)

    oxygen_balance = ((oxygen_atoms * 16) - molecular_weight) / molecular_weight * 100

    return oxygen_balance

def calculate_esp(smiles):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        print("Invalid SMILES input.")
        return None

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.ComputeGasteigerCharges(mol)

    esp_charge = sum([atom.GetDoubleProp('_GasteigerCharge') for atom in mol.GetAtoms()])

    return esp_charge

def has_graphite_like_structure(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0

    rings = Chem.GetSymmSSSR(mol)
    benzene_like = [ring for ring in rings if len(ring) == 6 and all(mol.GetAtomWithIdx(atom).GetAtomicNum() == 6 for atom in ring)]
    return int(len(benzene_like) >= 2)

def has_specific_bonds(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0

    substructure1 = Chem.MolFromSmarts("[N+](=O)[O-]")
    substructure2 = Chem.MolFromSmarts("O[N+](=O)[O-]")
    contains_ono2 = mol.HasSubstructMatch(substructure1)
    contains_nno2 = mol.HasSubstructMatch(substructure2)

    return int(contains_ono2 or contains_nno2)

def check_for_explosive_groups(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0

    explosive_functional_groups = ['[N+](=O)[O-]', 'N(=O)=O', 'N=[N+]=[N-]', 'N=[N+]=N']
    for group in explosive_functional_groups:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(group)):
            return 1

    return 0

def process_smiles(smiles):
    if isinstance(smiles, float):
        return {
            "SMILES": "",
            "Molecule Name": "",
            "Molecular Formula": "",
            "Bond Information": "",
            "Functional Groups": "",
            "Atom Information": "",
            "Symmetry": "",
            "Planarity": "",
            "Isomers":"",
            "Molecular weight": "",
            "PiPi Stacking Potential": "",
            "ESP": "",
            "oxygen_balance": "",
            "NO2 Count": "",
            "NH2 Count": "",
            "H Count": "",
            "Graphite-like Structure":"",
            "Specific bond":"",
            "Explosive groups":" ",
            "Energetics": ""
        }

    molecule_name = get_molecule_name(smiles)

    molecule = Chem.MolFromSmiles(smiles)
    if molecule:
        bonds = molecule.GetBonds()
        bond_info = []
        if bonds:
            for bond in bonds:
                bond_info.append(f"{bond.GetBeginAtomIdx()}-{bond.GetEndAtomIdx()} ({bond.GetBondType()})")
        else:
            bond_info.append("No bonds")

        functional_groups = get_functional_groups(molecule)

        atoms = molecule.GetAtoms()
        atom_info = [f"{atom.GetIdx()}: {atom.GetSymbol()}" for atom in atoms]

        molecular_formula = Chem.rdMolDescriptors.CalcMolFormula(molecule)

        symmetry = check_symmetry(smiles)

        planarity = check_planarity(smiles)

        isomers =  generate_isomers(smiles)

        homo_energy = homo_e(smiles)

        pipi_stacking_potential = pipi_stacking(smiles)

        esp_charge = calculate_esp(smiles)

        oxygen_balance = calculate_oxygen_balance(smiles)

        substituent_counts = count_specific_substituents(smiles)

        graphite_like_structure = has_graphite_like_structure(smiles)

        specific_bond = has_specific_bonds(smiles)

        exp_groups = check_for_explosive_groups(smiles)

        energetics = (symmetry == 1) & (planarity == 1) & (oxygen_balance >= -80) & (oxygen_balance <= -40) & (exp_groups == 1) & (specific_bond == 1)

        return {
            "SMILES": smiles,
            "Molecule Name": molecule_name,
            "Molecular Formula": molecular_formula,
            "Bond Information": ", ".join(bond_info),
            "Functional Groups": ", ".join(functional_groups),
            "Atom Information": ", ".join(atom_info),
            "Symmetry": symmetry,
            "Planarity": planarity,
            "Isomers": isomers,
            "Molecular weight": Descriptors.MolWt(molecule),
            "PiPi Stacking Potential": pipi_stacking_potential,
            "ESP": esp_charge,
            "oxygen_balance": oxygen_balance,
            "NO2 Count": substituent_counts["[N+](=O)[O-]"],
            "NH2 Count": substituent_counts["[NH2]"],
            "H Count": substituent_counts["[H]"],
            "Graphite-like Structure": graphite_like_structure,
            "Specific bond": specific_bond,
            "Explosive groups": exp_groups,
            "Energetics": int(energetics)
        }
    else:
        return {
            "SMILES": smiles,
            "Molecule Name": "Failed to create molecule",
            "Molecular Formula": "",
            "Bond Information": "",
            "Functional Groups": "",
            "Atom Information": "",
            "Symmetry": "",
            "Planarity": "",
            "Isomers" :"",
            "Molecular weight": "",
            "PiPi Stacking Potential": "",
            "ESP": "",
            "oxygen_balance": "",
            "NO2 Count": "",
            "NH2 Count": "",
            "H Count": "",
            "Graphite-like Structure": False,
            "Specific bond":"",
            "Explosive groups":"",
            "Energetics": ""
        }

@app.route('/upload', methods=['POST'])
def upload():
    file = request.files['file']
    if file:
        try:
            contents = []
            for line in file:
                contents.append(line.strip().split(','))  # Assuming CSV format with comma-separated values
            return jsonify({'contents': contents})
        except Exception as e:
            return jsonify({'error': str(e)}), 400
    return jsonify({'error': 'No file uploaded'}), 400

if __name__ == '__main__':
    app.run(debug=True, port=5002)
