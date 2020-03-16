import pandas as pd
import numpy as np
from IPython.display import display, HTML



class OrBi():
	def __init__(self, mat):
		self.mat = mat
		self.df = pd.DataFrame(mat)  
		
		
	def gini_mu(self, rows,cols, mode=1,base=2.0):
		_, r_mu = self.r_par(cols,mode=mode,base=base)
		_, c_mu = self.c_par(rows,mode=mode,base=base)
		return r_mu+c_mu-r_mu*c_mu
		
		
	def r_par(self,cols, mode=0,base=2.0):    
		
		mat_pat,index, counts = np.unique(self.mat[:,cols],axis=0,return_index=1, return_counts=1)		
		gini = 1 - np.sum(np.square(counts/np.sum(counts)))
		
		if (mode==0):
			pen = np.power(base,np.sum(mat_pat,axis=1))
		if (mode==1):
			pen = np.power(base,2 * np.sum(mat_pat,axis=1) - mat_pat.shape[1])
		gini_m = 1 - np.dot(np.square(counts/np.sum(counts)), pen)

		return gini, gini_m
		
		
		
	def c_par(self,rows,mode=0,base=2.0):    	
		
		mat_pat,index, counts = np.unique(self.mat[rows,:],axis=1,return_index=1, return_counts=1)
		
		gini = 1 - np.sum(np.square(counts/np.sum(counts))) 
		if (mode==0):
			pen = np.power(base,np.sum(mat_pat,axis=0))
		if (mode==1):
			pen = np.power(base,2 * np.sum(mat_pat,axis=0) - mat_pat.shape[0])
		
		
		gini_m = 1 - np.dot(np.square(counts/np.sum(counts)),pen)

		return gini, gini_m
		
	def proj(self, row, cols):
		
		return self.df.style.applymap(self._bg_map,
							   subset=pd.IndexSlice[row,cols])
	def rstr(self, col, rows):
				return self.df.style.applymap(self._bg_map,
							   subset=pd.IndexSlice[rows,col])
	
		
	def _bg_map(self,val, color='yellow'):
		return f'background-color: {color}'
	
	def _powerset(self,iterable):
		"powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
		from itertools import chain,combinations

		s = list(iterable)
		return list(chain.from_iterable(combinations(s, r) for r in range(len(s)+1)))[1:]
	def _allBi(self):
		
		s_r = list(self._powerset( range(self.mat.shape[0])))
		s_c = list(self._powerset( range(self.mat.shape[1])))
		
		from itertools import product
		
		return product(s_r,s_c)
	
	def allBi_giniu(self, mode=1,base=2.0):
		
		all_bi = list(self._allBi())
		
		self.all_bi_str = {}
		zero_one = []
		
		gini_u_ratio = []
		one_ratio = []
		for bi in all_bi:

			self.all_bi_str[self.tuple2str(bi)] = bi
			bic = self.mat[list(bi[0])][:,list(bi[1])]
			ones = np.sum(bic)

			zero_one.append((bic.size-ones,ones,bic.size))
			one_ratio.append(float(ones)/bic.size)
       
		gini_u_all = [self.gini_mu(i[0],i[1],mode=mode,base=base) for i in all_bi]
		gini_u_ratio = [g*ii[1]/ii[2] for g,ii in zip(gini_u_all, zero_one)]

		self.df_gini_all = pd.DataFrame(gini_u_all,index=all_bi,columns=['Gini_u'])
		self.df_gini_all['Zeros_Ones'] = zero_one
		self.df_gini_all['One_ratio'] = one_ratio

		return self.df_gini_all
	def show_bi(self,rows,cols):
		return self.df.style.applymap(self._bg_map,
							   subset=pd.IndexSlice[rows,cols])
	
	def tuple2str(self,bi):
		bi_str1 = [''.join(str(d)) for d in bi]
		
		return ','.join(bi_str1)
		

	
	
	
	




mat = np.array([[1,1,0,1,0],
[1,1,1,1,1],
[1,1,0,0,0],
[0,0,0,1,1],
[1,0,0,0,1],
])
mat

orbi = OrBi(mat)



import matplotlib.pyplot as plt

df_giniu_all = orbi.allBi_giniu(base=2.0)
display(df_giniu_all.describe())
df_giniu_all.plot()
display(df_giniu_all.sort_values(by=['Gini_u']))
print()
print("biclusters with lowest gini_u:")
index_min = df_giniu_all.loc[df_giniu_all['Gini_u']==df_giniu_all.min()[0]]
display(index_min)
in_min = [orbi.all_bi_str[bi] for bi in index_min.index]

[display(orbi.show_bi(list(i[0]),list(i[1]))) for i in list(in_min)]
print()
print()
print("biclusters with highest gini_u:")
index_max = df_giniu_all.loc[df_giniu_all['Gini_u']==df_giniu_all.max()[0]]
display(index_max)

in_max = [orbi.all_bi_str[bi] for bi in index_max.index]
[display(orbi.show_bi(list(i[0]),list(i[1]))) for i in list(in_max)]
print()




print()


display(df_giniu_all.sort_values(by=['Gini_u_ratio']))
print()
print("biclusters with lowest Gini_u_ratio:")
index_min = df_giniu_all.loc[df_giniu_all['Gini_u_ratio']==df_giniu_all['Gini_u_ratio'].min()]
display(index_min)
in_min = [orbi.all_bi_str[bi] for bi in index_min.index]

[display(orbi.show_bi(list(i[0]),list(i[1]))) for i in list(in_min)]
print()
print()
print("biclusters with highest Gini_u_ratio:")
index_max = df_giniu_all.loc[df_giniu_all['Gini_u_ratio']==df_giniu_all['Gini_u_ratio'].max()]
display(index_max)

in_max = [orbi.all_bi_str[bi] for bi in index_max.index]
[display(orbi.show_bi(list(i[0]),list(i[1]))) for i in list(in_max)]
print()

