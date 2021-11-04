#!/usr/bin/env python3
"""
    Main support for the website
"""
from flask import Flask
app = Flask(__name__)

@app.route('/')
def hello_world():
    return 'Hello, World!'
