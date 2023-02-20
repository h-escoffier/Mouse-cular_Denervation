import csv
from tqdm import tqdm


# The annotations for the genes comes from the Human Protein Atlas (https://www.proteinatlas.org/)
# TODO: It would be interesting to add an information about the tissues in addition to the cell type.


def extract_from_seurat(path):
    only_genes = []
    with open(path) as genes_per_cluster:
        content = [x.strip('\n').split(',') for x in genes_per_cluster.readlines()][1:]
    for line in tqdm(iterable=content, desc='extract_genes'):
        line = [elm[1:-1] if elm.startswith('"') and elm.endswith('"') else elm for elm in line]
        only_genes.append(line[7])
    return only_genes, content


def gene_annotation(path_annotations, gene_list):
    with open(path_annotations) as all_annotations:
        content = [x.strip('\n').split('\t') for x in all_annotations.readlines()][1:]
    gene_ct_ntpm = []
    for gene in tqdm(iterable=gene_list, desc='annotating'):
        gene_upper = gene.upper()
        for line in content:
            if gene_upper == line[1] and float(line[3]) >= 1:
                gene_ct_ntpm.append([gene, line[2], line[3]])
    return gene_ct_ntpm


def i_m_annotating(genes_tissue_list, genes_clusters):
    final_list = []
    final_list_excel_version = []
    for line in tqdm(iterable=genes_clusters, desc='merge'):
        line = [elm[1:-1] if elm.startswith('"') and elm.endswith('"') else elm for elm in line]
        for gene in genes_tissue_list:
            if line[7] == gene[0]:
                final_list.append([line[6], line[7], str(line[5]), gene[1], float(gene[2])])
                # Excel interpret X.X as a data that why I replace it by X,X
                final_list_excel_version.append([line[6], line[7], str(line[5]), gene[1], gene[2].replace('.', ',')])
    # columns=['cluster', 'gene', 'pvalue_a', 'cell_type', 'nTPM']
    with open('genes_per_cluster_w_cell_type.tsv', 'w', newline='') as f_output:
        content = csv.writer(f_output, delimiter='\t')
        content.writerows(final_list_excel_version)


if __name__ == '__main__':
    print('start')
    # Extract all the genes
    genes, all_content = extract_from_seurat('data/annotation/integrated_top10_genes_per_cluster.csv')
    # Get all the cell type with a nTPM > 1 for each genes
    genes_tissues = gene_annotation('data/annotation/rna_single_cell_type.tsv', genes)
    # Create a tsv with all the genes
    i_m_annotating(genes_tissues, all_content)
    print('end')

