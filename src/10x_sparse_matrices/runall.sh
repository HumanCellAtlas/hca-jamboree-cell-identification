# PBMC data
wget http://cf.10xgenomics.com/samples/cell-exp/2.0.1/pbmc8k/pbmc8k_raw_gene_bc
_matrices.tar.gz

tar -xz -f pbmc8k_raw_gene_bc_matrices.tar.gz
mv raw_gene_bc_matrices pbmc8k_raw_gene_bc_matrices

Rscript tenx_matrix_dir_to_sparse_rds.R pbmc8k_raw_gene_bc_matrices 10Xpbmc8k_matrix.rds

# Mouse brain data
wget http://cf.10xgenomics.com/samples/cell-exp/2.0.1/neuron_9k/neuron_9k_raw_gene_bc_matrices.tar.gz
tar -xz -f neuron_9k_raw_gene_bc_matrices.tar.gz
mv raw_gene_bc_matrices neuron9k_raw_gene_bc_matrices

Rscript tenx_matrix_dir_to_sparse_rds.R neuron9k_raw_gene_bc_matrices 10Xmouseneuron9k_matrix.rds --genome mm10

# Pan t-cells
wget http://cf.10xgenomics.com/samples/cell-exp/2.0.1/t_4k/t_4k_raw_gene_bc_matrices.tar.gz
tar -xz -f t_4k_raw_gene_bc_matrices.tar.gz
mv raw_gene_bc_matrices pant4k_raw_gene_bc_matrices

Rscript tenx_matrix_dir_to_sparse_rds.R pant4k_raw_gene_bc_matrices 10Xpant4k_matrix.rds

rm *.tar.gz
