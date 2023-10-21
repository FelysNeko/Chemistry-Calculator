# Chemistry-Calculator

A Calculator that Solves Grade 11 Chemistry Problems

## Web Application Version(v0.3)

This app now has separate React front-end and Flask back-end, connected via REST API.

To run locally, run local server of back-end and front-end separately.

### Backend
Navigate to server/ and run
```shell
python3 -m venv .venv
. .venv/bin/activate
pip install Flask
```

### Frontend
Navigate to client/ and run
```shell
npm start
```

<br>

## API Version(v1.1)

You may import this module and use it in jupyter notebook, see all available exmaples in example.py. The website is based on it.

```python3
import api.chemistry as chem
```
<br>

## Command Line Version(v4.1)

Execute the command line version in terminal. Make sure you run the chemistry.py under 'colnsole' folder. This was the oldest code, and I made the very first version one year ago which is the time that I don't even know OOP. However, it provided me with lots of fundamental support for any further development like the API version.

```shell
python3 chemistry.py
```

**Note: The molarcular calculator seems to be borken, but I don't have time to fix it since the code is nasty and old**
