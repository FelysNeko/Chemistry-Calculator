from flask import Flask, render_template, request, flash, jsonify
from flask_cors import CORS
import api.chemistry as chem
from api.form import res


app = Flask(__name__)

CORS(app)

#/api/equation
@app.route("/api/equation", methods=['POST'])
def equation():
    if not request.is_json: 
        return "Request type not supported", 400
    data =  request.get_json()

    action = str(data["action"])
    equation = chem.Equation([data["molecule1"], data["molecule2"]]).rectify(True)

    if True in [i==chem.Molecule.null('1') for i in equation.molecule]:
        return jsonify(res('At least one molecule is incorrect', False))
    
    if (action == "predict"):
        result = ' + '.join(equation.short) + ' -> ' + ' + '.join(equation.prediction.short)
    elif (action == "balance"):
        re, pr = equation.balance(manual=False)
        result = ' + '.join(re.short) + ' -> ' + ' + '.join(pr.short)
    else:
        return jsonify(res('error', False))

    return jsonify(res(result, True))


#/api/property
@app.route("/api/property", methods=['POST'])
def property():
    if not request.is_json: 
        return "Request type not supported", 400
    data =  request.get_json()

    action = str(data["action"])
    molecule = chem.Molecule(data["molecule"])

    if molecule == chem.Molecule.null('1').rectify(True):
        return jsonify(res('Invalid molecule', False))
    
    if (action == "mass"):
        result = molecule.mass
    elif (action == "solubility"):
        result = molecule.solubility
    elif (action == "bond"):
        result = molecule.bond
    else:
        return jsonify(res('error', False))

    return jsonify(res(str(result), True))

if __name__ == "__main__":
    app.run(debug=True, port=8080)
