#!/usr/local/bin/python3

import subprocess
import os
from io import StringIO

def get_full_path(path):
	return os.getcwd() + "/" + path

class HYPHY:
	def __init__(self, hyphy_path = 'HYPHYMP', file = None):
		self.hyphy_path = hyphy_path
		self.file = file
		self.proc = None
	
	def run(self, input):
		if self.file:
			cmd = [self.hyphy_path, '-f', self.file]
		else:
			cmd = [self.hyphy_path]
		self.proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, universal_newlines=True)
		self.proc.communicate(input)
