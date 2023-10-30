# Chemistry-Calculator

A Calculator that Solves Basic Chemistry Problems

## Web Application Version(v0.3)

This app now has separate React front-end and Flask back-end, connected via REST API.

To run locally, run local server of back-end and front-end separately.

### Back-end
Navigate to ./server and run the following code to install dependencies
```shell
python3 -m venv .venv
. .venv/bin/activate
pip install Flask sympy
```
and then run
```shell
python app.py
```

### Front-end
Navigate to ./client and run
```shell
npm start
```

<br>

## API Version(v1.2)

You may import this module and use it in jupyter notebook, see all available exmaples in example.py. The website is based on it. Please only import the class or function you need.

```python3
from api.chemistry import Molecule, Equation, lookup
```
<br>

## Command Line Version(v4.1)

Execute the command line version in terminal. Make sure you run the chemistry.py under 'colnsole' folder. This was the oldest code, and I made the very first version one year ago which is the time that I don't even know OOP. However, it provided me with lots of fundamental support for any further development like the API version.

```shell
python3 chemistry.py
```

**Note: The molarcular calculator seems to be borken, but I don't have time to fix it since the code is nasty and old**
