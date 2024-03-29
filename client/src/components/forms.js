import { useState } from 'react';
import { TextField, Button } from '@mui/material';

let fetchPostAPI = (url, data) => fetch(url, { 
  method: "POST",
  mode: "cors",
  headers: {
    "Content-Type": "application/json",
  },
  body: JSON.stringify(data),
})

function getFormJson(e) {
  // Read the form data
  const form = e.target;
  const formData = new FormData(form);
  // plain object:
  const formJson = Object.fromEntries(formData.entries());
  console.log(formJson);
  return formJson;
}

export function EquationForm() {
  const[result, setResult] = useState("");
  function handleSubmit(e) {
    // Prevent the browser from reloading the page
    e.preventDefault();

    // Pass formData as a fetch body directly:
    fetchPostAPI('http://127.0.0.1:8080/api/equation', getFormJson(e))
    .then((Response)=>Response.json()).then((data) => {
      console.log(data);
      if (data.success) {
        setResult(data.result);
      } else {
        alert(data.msg);
      }
    });
  }
  
    return (
      <>
      <form method="post" onSubmit={handleSubmit}>
        <h1>Equation</h1>
        <TextField 
           name="molecule1" label="Molecule 1" />
        <>+</>
        <TextField name="molecule2" label="Molecule 2" />
        <br />
        <label htmlFor="actution">Calculation: </label>
        <select name="action">
          <option value="balance">Balance</option>
          <option value="predict">Predict</option>
        </select>
        <br />
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

    fetchPostAPI('http://127.0.0.1:8080/api/property', getFormJson(e))
    .then((Response)=>Response.json()).then((data) => {
      if (data.success) {
        setResult(data.result);
      } else {
        alert(data.msg);
      }
      console.log(data);
    });      
  }

  return (
    <>
    <form method="post" onSubmit={handleSubmit}>
      <h1>Substance Property</h1>
      <TextField 
          name="molecule" label="Substance (Molecule/Element)" />
      <br />
      <label for="action">Determin: </label>
        <select name="action">
          <option value="mass">Mass</option>
          <option value="solubility">Solubility</option>
          <option value="bond">Bond</option>
        </select>
      <br />
      <Button type="submit">Reveal!</Button>
    </form>
    <p>{result}</p>
    </>
  );
}