#!/usr/bin/env python

import logging
import subprocess
import os

__author__ = 'jiri1kolar'
__email__ = "jiri1kolar@gmail.com"
#load logger

def execute_evolver(file_name):
    output_file = file_name[:-4] + '.json'
    logging.info('Executing evolver ...')
    thread = subprocess.Popen(
        ['se_api', '-i', file_name,
         '--all-union', nvolpercell.__str__(),
         '-a', output_file])
    thread.wait()
    logging.info('Executing evolver ... done')
    return output_file

