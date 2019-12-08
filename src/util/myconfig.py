"""
This module is to set the common variable for project
Left to update.
"""
import os 
class config(object):
    def __init__(self, root_path, log_dir="", data_type="p", data_path=None):
        # path
        self.log_dir = log_dir
        self.root_path = root_path
        if data_path:
            self.data_path = data_path
        else:
            self.data_path = os.path.join(self.root_path, "data/" + data_type + "-data/")
        self.assay_data_path = os.path.join(self.data_path, "assay-data/begGraph/")
        self.pilot_region_path = os.path.join(self.data_path, "pilot-data")
        self.pilot_region_file = os.path.join(self.data_path, "pilot-data/encodePilotRegions.hg19.bed")
        self.processed_pilot_data_path = os.path.join(self.data_path, "processed-pilot-data/")
        self.blacklist_region_file = os.path.join(self.data_path, "blacklist-data/ENCFF419RSJ.bed")
        self.blacklist_rm_data_path = os.path.join(self.data_path, "blacklist-rm-data/")
        self.experiment_path = os.path.join(self.root_path, "experiment")
        self.experiment_path_new = os.path.join(self.root_path, "experiment-new")
        self.assay_region_file = os.path.join(self.data_path, "assay-data/assay-region.bed")
    
        # others
        self.chromosome_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 
                                'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                                'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 
                                'chr19', 'chr20', 'chr21', 'chr22','chrX']
        self.e003_assays = ['E003-DNase', 
                'E003-H2A.Z', 'E003-H2AK5ac', 
                'E003-H2BK5ac', 'E003-H2BK12ac', 'E003-H2BK15ac', 'E003-H2BK20ac', 'E003-H2BK120ac', 
                'E003-H3K4ac', 'E003-H3K4me1', 'E003-H3K4me2', 'E003-H3K4me3', 
                'E003-H3K9ac', 'E003-H3K9me3', 'E003-H3K14ac', 'E003-H3K18ac', 'E003-H3K23ac', 
                'E003-H3K23me2', 'E003-H3K27ac', 'E003-H3K27me3', 'E003-H3K36me3', 'E003-H3K56ac', 
                'E003-H3K79me1', 'E003-H3K79me2', 'E003-H4K5ac', 'E003-H4K8ac', 'E003-H4K20me1', 
                'E003-H4K91ac']
        self.assay_list = ['E116-H3K27ac', 'E116-H3K4me1', 
                            'E116-H3K36me3', 
                            'E116-H4K20me1', 'E116-H3K9ac', 
                            'E116-H3K4me2', 'E116-H3K4me3', 
                            'E116-H3K27me3']
        self.e003_assays = ['E003-DNase', 
                'E003-H2A.Z', 'E003-H2AK5ac', 
                'E003-H2BK5ac', 'E003-H2BK12ac', 'E003-H2BK15ac', 'E003-H2BK20ac', 'E003-H2BK120ac', 
                'E003-H3K4ac', 'E003-H3K4me1', 'E003-H3K4me2', 'E003-H3K4me3', 
                'E003-H3K9ac', 'E003-H3K9me3', 'E003-H3K14ac', 'E003-H3K18ac', 'E003-H3K23ac', 
                'E003-H3K23me2', 'E003-H3K27ac', 'E003-H3K27me3', 'E003-H3K36me3', 'E003-H3K56ac', 
                'E003-H3K79me1', 'E003-H3K79me2', 'E003-H4K5ac', 'E003-H4K8ac', 'E003-H4K20me1', 
                'E003-H4K91ac']
        self.test_assay = ['E116-H3K4me1', 'E116-H3K27ac']
        self.chromHMM_assay = ['E116-H3K4me3', 'E116-H3K4me1', 'E116-H3K36me3', 'E116-H3K27me3', 'E116-H3K9me3']
        self.gexp_signal_file = self.data_path + "tss-data/Ensembl_v65.Gencode_v10.ENSG.gene_info"
        self.gexp_region_file = self.data_path + "tss-data/57epigenomes.RPKM.pc"
        self.gexp_file = self.data_path + "tss-data/gexp.bed"
        self.chromHMM_annotation_result_dir = os.path.join(self.data_path, "annotation-result/exp26-1/")
        self.ssm_annotation_result_dir = os.path.join(self.data_path, "annotation-result")
        self.annotation_result_dir = os.path.join(self.data_path, "annotation-result")   
        self.enhancer_file = os.path.join(self.data_path, "enhancer-data/enhancer.bed")
        self.filtered_gexp_file = os.path.join(self.experiment_path, "35-0728-18-gexpEnhancerFilter-p", "exp1", "filtered_gexp.bed")
        self.filtered_enhancer_file = os.path.join(self.experiment_path, "35-0728-18-gexpEnhancerFilter-p", "exp2", "filtered_enhancer.bed")
    
    def pilot_index_region(self, resolution_size):
        return os.path.join(self.pilot_region_path, "index_resol" + str(resolution_size) + "bp.bed")

    def openLog(self, log_name = "log"):
        self.log_file = os.path.join(self.log_dir, log_name + ".txt")
        self.log_f = open(self.log_file, 'w')

    def toLog(self, msg):
        # using print just for convenience, not a good implementation anyway
        print(msg, file=self.log_f)

    def closeLog(self):
        self.log_f.close()
