# Chemistry-Calculator

A Calculator that Solves Basic Chemistry Problems

<br>

## Web Application(v0.3)

This app now has separate React front-end and Flask back-end, connected via REST API.

To run locally, run local server of back-end and front-end separately.

### Back-end
Navigate to [server](server) and run the following code to install dependencies
```shell
python3 -m venv .venv
source .venv/bin/activate
pip3 install Flask sympy
```

and then run
```shell
python3 app.py
```

### Front-end
Navigate to [client](client) and run
```shell
npm start
```

<br>

## Package(v1.2)

You may import this module and use it in jupyter notebook, see all available exmaples in [guide.py](server/guide.py). The website is based on it. Please only import the class or function you need.

```python3
from lib.chemistry import Molecule, Equation, lookup
```

You may also find the future package(v2.0) under folder [future](future) and simply import it by `from chemistry import *`. Note that a lot of functions are removed in this version since I am thinking of transfering this project to a university level. Hence, all the implementation of high school chemistry are not available.

<br>

## Command Line Tool(v4.1)

Execute the command line version in terminal. Make sure you run the chemistry.py under 'console' folder. This was the oldest code, and I made the very first version one year ago which is the time that I don't even know OOP. However, it provided me with lots of fundamental support for any further development like the package version.

```shell
python3 chemistry.py
```

**Note: The molarcular calculator seems to be borken, but I don't have time to fix it since the code is nasty and old**
