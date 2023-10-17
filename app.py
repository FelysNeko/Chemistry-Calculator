from flask import Flask, render_template, request, flash
import api.chemistry as chem


app = Flask(__name__)
app.secret_key = "ELYSIA IS MY WAIFU"


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/equation", methods=['POST', 'GET'])
def equation():
    action = str(request.form["action"])
    equation = chem.Equation([request.form["molecule1"], request.form["molecule2"]]).rectify(True)

    if True in [i==chem.Molecule.null('1') for i in equation.molecule]:
        flash("At least one molecule is incorrect")
        return render_template("index.html")
    
    if (action == "predict"):
        result = ' + '.join(equation.short) + ' -> ' + ' + '.join(equation.prediction.short)
    elif (action == "balance"):
        re, pr = equation.balance(manual=False)
        result = ' + '.join(re.short) + ' -> ' + ' + '.join(pr.short)
    else:
        result = "error"

    flash(result)
    return render_template("index.html")


@app.route("/property", methods=['POST', 'GET'])
def property():
    action = str(request.form["action"])
    molecule = chem.Molecule(request.form["molecule"])

    if molecule == chem.Molecule.null('1').rectify(True):
        flash("Invalid molecule")
        return render_template("index.html")
    
    if (action == "mass"):
        result = molecule.mass
    elif (action == "solubility"):
        result = molecule.solubility
    elif (action == "bond"):
        result = molecule.bond
    else:
        result = "error"

    flash(molecule.short + ': ' + str(result))
    return render_template("index.html")
