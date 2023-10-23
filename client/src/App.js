import './App.css';
import { useState } from 'react';
import { TextField, Button } from '@mui/material';


function App() {
  return (
    <>
      <EquationForm />
      <MoleculeForm />
    </>
  );
}

export function EquationForm() {
  const[result, setResult] = useState("");
  function handleSubmit(e) {
    // Prevent the browser from reloading the page
    e.preventDefault();

    // Read the form data
    const form = e.target;
    const formData = new FormData(form);

    // Pass formData as a fetch body directly:
    fetch('http://127.0.0.1:8080/api/equation', { method: form.method, body: formData })
      .then((Response)=>Response.json()).then((data) => {
        
        console.log(data);
        if (data.success) {
          setResult(data.result);
          
        } else {
          alert(data.msg);
        }
      });

    // plain object:
    const formJson = Object.fromEntries(formData.entries());
    console.log(formJson);
  }

  return (
    <>
    <form method="post" onSubmit={handleSubmit}>
      <h1>Equation</h1>
      <TextField 
         name="molecule1" label="Molecule 1" />
      <>+</>
      <TextField name="molecule2" label="Molecule 2" />
      <hr />
      <label htmlFor="actution">Calculation: </label>
      <select name="action">
        <option value="balance">Balance</option>
        <option value="predict">Predict</option>
      </select>
      <hr />
      <Button type="submit">Calculate!</Button>
    </form>
    <p>{result}</p>
    </>
  );
}

export function MoleculeForm() {
  const[result, setResult] = useState("");
  function handleSubmit(e) {
    // Prevent the browser from reloading the page
    e.preventDefault();

    // Read the form data
    const form = e.target;
    const formData = new FormData(form);

    // Pass formData as a fetch body directly:
    fetch('http://127.0.0.1:8080/api/property', { method: form.method, body: formData })
      .then((Response)=>Response.json()).then((data) => {
        if (data.success) {
          setResult(data.result);
          
        } else {
          alert(data.msg);
        }
        console.log(data);
      });

    // plain object:
    const formJson = Object.fromEntries(formData.entries());
    console.log(formJson);
  }

  return (
    <>
    <form method="post" onSubmit={handleSubmit}>
      <h1>Substance Property</h1>
      <TextField 
         name="molecule" label="Substance (Molecule/Element)" />
      <hr />
      <label for="action">Determin: </label>
        <select name="action">
          <option value="mass">Mass</option>
          <option value="solubility">Solubility</option>
          <option value="bond">Bond</option>
        </select>
      <hr />
      <Button type="submit">Reveal!</Button>
    </form>
    <p>{result}</p>
    </>
  );
}



export default App;
