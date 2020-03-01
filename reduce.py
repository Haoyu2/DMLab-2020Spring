import pandas as pd
import numpy as np




class Reduce:
	def __init__(self, mat):
		self.mat = mat
		self.row_sorted = np.argsort(np.sum(mat,axis=1))
		self.col_sorted = np.argsort(np.sum(mat,axis=0))
	def reduce(self,row,col):
		rows = np.sort(self.row_sorted[-row:])
		cols = np.sort(self.col_sorted[-col:])
		return self.mat[rows][:cols]
