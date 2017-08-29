import numpy as np
from sklearn.decomposition import PCA
#========================================================================================#
## These are standard sub-routines involved with knn graph construction

def get_vscores(E, min_mean=0, nBins=50, fit_percentile=0.1, error_wt=1):
	''' Sub-routine for choosing highly variable genes'''
	mu_gene = np.mean(E, axis=0)
	gene_ix = np.nonzero(mu_gene > min_mean)[0]
	mu_gene = mu_gene[gene_ix]
	FF_gene = np.var(E[:,gene_ix], axis=0) / mu_gene

	data_x = np.log(mu_gene)
	data_y = np.log(FF_gene / mu_gene)

	x, y = runningquantile(data_x, data_y, fit_percentile, nBins)
	x = x[~np.isnan(y)]
	y = y[~np.isnan(y)]

	gLog = lambda input: np.log(input[1] * np.exp(-input[0]) + input[2])
	h,b = np.histogram(np.log(FF_gene[mu_gene>0]), bins=200)
	b = b[:-1] + np.diff(b)/2
	max_ix = np.argmax(h)
	c = np.max((np.exp(b[max_ix]), 1))
	errFun = lambda b2: np.sum(abs(gLog([x,c,b2])-y) ** error_wt)
	b0 = 0.1
	b = scipy.optimize.fmin(func = errFun, x0=[b0])
	a = c / (1+b) - 1
	v_scores = FF_gene / ((1+a)*(1+b) + b * mu_gene);
	CV_eff = np.sqrt((1+a)*(1+b) - 1);
	CV_input = np.sqrt(b);

	return v_scores, CV_eff, CV_input, gene_ix, mu_gene, FF_gene, a, b


def filter_genes_fitting(E,min_count,min_cells,min_vscore_pctl):
	''' Find highly variable genes using noise-above poisson statistic'''	
	Vscores, CV_eff, CV_input, gene_ix, mu_gene, FF_gene, a, b = get_vscores(E)	
	min_log_vscore = np.percentile(np.log(Vscores),min_vscore_pctl)
	x_min = 0.5*np.min(mu_gene);
	x_max = 2*np.max(mu_gene);
	xTh = x_min * np.exp(np.log(x_max/x_min)*np.linspace(0,1,100));
	yTh = (1 + a)*(1+b) + b * xTh;
	ix = ((np.sum(E[:,gene_ix] >= min_count,axis=0) >= min_cells) & (np.log(Vscores) >= min_log_vscore))
	gene_filter = gene_ix[ix]
	return gene_filter


def row_normalize(X):
	''' Normalize a matrix so that the rows sum to 1 '''
	d = np.sum(X,axis=1)
	return X / np.tile(d[:,None],(1,X.shape[1]))

def Zscore(X):
	''' Xscore normalize the columns of a matrix'''
	means = np.tile(np.mean(X,axis=0)[None,:],(X.shape[0],1))
	stds = np.tile(np.std(X,axis=0)[None,:],(X.shape[0],1))
	return (X - means) / (stds + .0001)

def get_distance_matrix(M):
	''' Compute the distance matrix of M using rows as observations and columns as features'''
	D = np.zeros((M.shape[0],M.shape[0]))
	for i in range(M.shape[0]):
		Mtiled = np.tile(M[i,:][None,:],(M.shape[0],1))
		D[i,:] = np.sqrt(np.sum((Mtiled - M)**2, axis=1))
	return D


def get_adjacency_matrix(D, k):
	''' Get knn adjacency matrix from distances D and num nearest neighbors k'''
	adjmat = np.zeros((D.shape[0],D.shape[0]))
	for i in range(D.shape[0]):
		sorted_nodes = np.argsort(D[i,:])[1:k+1]
		for j in sorted_nodes:
			adjmat[i,j] = 1
			adjmat[j,i] = 1
	return adjmat

#========================================================================================#
# These are higher level functions specific to the Wolball algorithm

def coembed(A,B,k,p):
	'''
	Input:
		A [array; required]: Gene expression matrix
		B [array; required]: Gene expression matrix
		
	Output:
		adjacency matrix of knn graph built from A and B together. 
		A is used to define the PC space for graph construction
	'''
	# row normalize A and B
	A = row_normalize(A)
	B = row_normalize(B)
	
	# Find highly variable genes
	gene_filter = filter_genes_fitting(A*3000,3,3,80)
	
	# Get PC space
	pca = PCA(n_components=p)
	pca.fit(Zscore(A[:,gene_filter]))
	
	# Project A and B into PC space
	Xpca = pca.transform(Zscore(np.vstack((A,B))[:,gene_filter]))

	# Calculate distance matrix and then adjacency matrix
	D = get_distance_matrix(Xpca)	
	return get_adjacency_matrix(D,k)
	
def diffuse(score,adjacency_matrix,s):
	'''
	Input:
		score [1D array; required]: a vector of scores to be diffused
		adjacency_matrix [2D array; required]: adjacency matrix of nearest neighbor graph.
											   must be square matrix with same size as 'score'
		s [float; required]: mean path length for the diffusion process
	
	Output:
		A vector of smoothed scores (same shape as score)
	'''
	adjacency_matrix = row_normalize(adjacency_matrix)
	beta = 1./s; n = adjacency_matrix.shape[0]
	mm = beta * np.linalg.inv(np.identity(n) - (1-beta) * adjacency_matrix)
	return np.dot(mm,score)


def get_wolball_score(A,B,k=5,s=20,p=30):
	'''
	Input:
		A [array; required]: a matrix of good/ambiguous cells (rows=cells, columns=genes)
		B [array; required]: a matrix of 'cells' that are background with high confidence
		k [int; optional]: number of nearest neighbors in the knn graph
		s [int; optional]: average path length of the diffusion operator
		p [int; optional]: number of PCs to use for knn graph construction
	
	Output:
		Array of Wolball scores with length equal to the number of rows in A
	'''
	# Step 1: Verify input
	if A.shape[1] != B.shape[1]: print 'Error: A and B must have the same number of columns'; return
	
	# Step 2: Co-embed A and B in a knn graph; use A to define the PC space for embedding
	print 'Co-embedding cells in knn graph'
	adjacency_matrix = coembed(A,B,k,p)
	
	# Step 3: Perform graph diffusion to propagate background labels from B to A
	print 'Diffusion background labels on the graph'
	raw_score = np.hstack((np.zeros(A.shape[0]),np.ones(B.shape[0])))
	smooth_score = diffuse(raw_score,adjacency_matrix,s)[:A.shape[0]]
	return rescale(smooth_score)

















