from flask import Flask, render_template, request, flash
import api.chemistry as chem

app = Flask(__name__)
app.secret_key = "manbearpig_MUDMAN888"

@app.route("/")
def index():
    
    return render_template("index.html")

@app.route("/equation", methods=['POST', 'GET'])
def equation():
    
    action = str(request.form["action"])
    m1 = str(request.form["molecule1"])
    m2 = str(request.form["molecule2"])
    if (not (bool(chem.expand(m1)) and bool(chem.expand(m2)))):
        flash("At least one molecule is incorrect.")
        return render_template("index.html")
    equation = chem.Equation([m1, m2])
    if (action == "predict"):
        result = equation.predict().short[0]
        for i in range(1, len(equation.predict().short)):
            result = result + '+' + equation.predict().short[i]

    elif (action == "balance"):
        re, pr = equation.balance()
        re = re.short
        pr = pr.short
        result = re[0]
        for i in range(1, len(re)):
            result = result + '+' + re[i]
        result = result + ' → '
        result = result + pr[0]
        for i in range(1, len(pr)):
            result = result + '+' + pr[i]
    else:
        result = "error"
    flash(result)

    return render_template("index.html")

@app.route("/property", methods=['POST', 'GET'])
def property():
    action = str(request.form["action"])
    molecule = chem.Molecule(str(request.form["molecule"]))
    if (not bool(chem.expand(str(request.form["molecule"])))):
        flash("Invalid molecule. ")
        return render_template("index.html")
    elif (action == "mass"):
        result = molecule.mass
    elif (action == "solubility"):
        result = molecule.solubility
    elif (action == "bond"):
        result = molecule.bond
    else:
        result = "error"
    flash(str(result))
    return render_template("index.html")
