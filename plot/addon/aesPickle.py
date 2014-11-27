#!/usr/bin/python3
# -*- coding: utf-8 -*-
# 
# Author:  Mario S. Könz <mskoenz@gmx.net>
# Date:    30.07.2014 14:18:52 CEST
# File:    aesPickle.py

try:
    import cPickle as pickle
except ImportError:
    import pickle

from Crypto import Random
from Crypto.Cipher import AES

from .helper import *

validate = b"  you may pass  "

def dumps(data, key):
    # serialize data
    encrypt = pickle.dumps(data)
    encrypt = padding(encrypt, 16) + validate
    
    # create initial vector
    rndfile = Random.new()
    iv = rndfile.read(16)
    
    # encrypt data
    obj = AES.new(padding(key, 16), AES.MODE_CBC, iv)
    ciphertext = obj.encrypt(encrypt)
    
    return iv + ciphertext

def dump(data, filename, key):
    f = open(filename, 'wb')
    c = dumps(data, key)
    f.write(c)
    f.close()

def loads(text, key):
    # get initial vector
    iv = text[:16]
    ciphertext = text[16:]
    
    # decrypt data
    obj2 = AES.new(padding(key, 16), AES.MODE_CBC, iv)
    decrypt = obj2.decrypt(ciphertext)
    
    val = decrypt[-16:]
    decrypt = decrypt[:-16]
    
    # check if it worked / deserialize data
    if val == validate:
        data = pickle.loads(decrypt)
        return data
    else:
        ERROR("decryption failed")
        return 0

def load(filename, key):
    f = open(filename, 'rb')
    c = f.read()
    f.close()
    return loads(c, key)
