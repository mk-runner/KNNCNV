# -*-coding:utf-8-*-

from preprocessing import preprocessing

from sklearn.mixture import BayesianGaussianMixture
from sklearn.utils.validation import check_array
from pyod.models.knn import KNN
from scipy import stats
import pandas as pd
import numpy as np


def vbgmm(scores):
    """
    The variational Bayesian Gaussian mixture model (VBGMM)
    :param scores: ndarray (n_samples, 1), the outlier scores for all genome segments
    :return: binary labels, where 1 is outlier, 0 stands for inlier
    """
    scores = check_array(scores)
    clf = BayesianGaussianMixture(n_components=2)
    labels = clf.fit_predict(scores)
    outlier_label = np.argmax(clf.means_)
    labels = labels == outlier_label
    return labels.astype(int)


def combiningCNV(seg_chr, seg_start, seg_end, seg_count, labels, mode):
    def _func(x):
        if x == 2:
            return "duplication"
        else:
            return "deletion"

    index = labels == 1
    CNV_chr = seg_chr[index]
    CNVstart = seg_start[index]
    CNVend = seg_end[index]
    CNVRD = seg_count[index, 0]

    type = np.full(len(CNVRD), 1)
    for i in range(len(CNVRD)):
        if CNVRD[i] > mode:  # "duplication"
            type[i] = 2

    for i in range(len(CNVRD) - 1):
        if CNVend[i] + 1 == CNVstart[i + 1] and type[i] == type[i + 1]:
            CNVstart[i + 1] = CNVstart[i]
            type[i] = 0

    index = type != 0
    CNVRD = CNVRD[index]
    CNV_chr = CNV_chr[index]
    CNVstart = CNVstart[index]
    CNVend = CNVend[index]
    CNVtype = type[index]
    CNVtype = [_func(i) for i in CNVtype]
    return CNV_chr, CNVstart, CNVend, CNVRD, CNVtype


def save_result(CNV_chr, CNVstart, CNVend, CNVRD, CNVtype, save_name=None):
    df = pd.DataFrame()
    df['chr'] = CNV_chr
    df["start"] = CNVstart
    df["end"] = CNVend
    df["type"] = CNVtype
    df["RD"] = CNVRD
    if save_name is None:
        df.to_csv("detail_result.txt", index=False)
    else:
        df.to_csv("%s.txt" % save_name, index=False)


def sta_performance(groudtruth_path, result_start, result_end, result_type, is_simulation=False):
    ground_truth = pd.read_table(groudtruth_path)
    if is_simulation:
        truth_type = ground_truth["state"].apply(lambda x: 'duplication' if x == 'gain' else 'deletion').tolist()
        truth_start = ground_truth['start'].tolist()
        truth_end = ground_truth['end'].tolist()
    else:
        truth_type = ground_truth["variant type"].tolist()
        truth_start = ground_truth['start'].tolist()
        truth_end = ground_truth['stop'].tolist()

    count = 0
    for i in range(len(result_type)):
        for j in range(len(truth_type)):
            if truth_start[j] <= result_start[i] <= truth_end[j] and truth_type[j] == result_type[i]:
                if result_end[i] <= truth_end[j]:
                    count += (result_end[i] - result_start[i] + 1)

                elif result_end[i] > truth_end[j]:
                    count += (truth_end[j] - result_start[i] + 1)

            elif truth_start[j] >= result_start[i] and truth_type[j] == result_type[i]:
                if truth_start[j] <= result_end[i] <= truth_end[j]:
                    count += (result_end[i] - truth_start[j] + 1)

                elif result_end[i] >= truth_end[j]:
                    count += (truth_end[j] - truth_start[j] + 1)

    result_count = 0
    for i in range(len(result_start)):
        result_count += (result_end[i] - result_start[i] + 1)

    truth_count = 0
    for i in range(len(truth_start)):
        truth_count += (truth_end[i] - truth_start[i] + 1)

    if result_count == 0:
        precision = 0
    else:
        precision = count / result_count
    sensitivity = count / truth_count

    return [precision, sensitivity, stats.hmean((precision, sensitivity))]


def knncnv(bam_path, fa_path, bin_size=1000, gt_path=None, cbs_imp='python', ncol=50,
         is_simulation=False, seg_path=None, iter_num=5):
    """
    Parameters
    ----------
    bam_path : str
        Local path of the *.bam file (i.e., sequenced sample).

    fa_path : str
        Local path of the *.fasta file or the *.fa file (i.e., reference genome).

    bin_size : int, optional (default=1000)
        The bin size.

    gt_path : str, optional (default=None)
        Local path of the ground truth of the sequenced sample.

    cbs_imp: str, optional (default='python')
        The implementation of CBS algorithm. In addition to "python", cbs_imp can also be "R".

    ncol : int, optional (default=50)
        The number of  partitions to CBS in R.

    is_simulation: bool, optional (default=False)
        is the simulation dataset executed

    seg_path: str, optional (default is None)
        If the CBS program is executed by R language, seg_path must be specified

    iter_num: int, optional (default = 5)
        Repeated times for each sample.

    Returns
    ----------

    """

    for epoch in range(iter_num):

        # Preprocessing
        all_chr, all_start, all_end, all_rd, mode = preprocessing(bam_path, fa_path, bin_size=bin_size, cbs_imp=cbs_imp,
                                                                  ncol=ncol, is_simulation=is_simulation, seg_path=seg_path)

        n_neighbors = np.random.choice(range(int(len(all_rd) * 0.2), int(len(all_rd) * 0.35)))
        clf = KNN(n_neighbors=n_neighbors, method='mean')
        clf.fit(all_rd)
        scores = clf.decision_function(all_rd)
        labels = vbgmm(scores.reshape(-1, 1))

        # Statistics of experimental results
        CNV_chr, CNVstart, CNVend, CNVRD, CNVtype = combiningCNV(all_chr, all_start, all_end, all_rd, labels, mode)

        # save_result(CNV_chr, CNVstart, CNVend, CNVRD, CNVtype)

        if gt_path:
            # statistic the performance(i.e., precision and sensitivity)
            temp_ans = sta_performance(gt_path, CNVstart, CNVend, CNVtype, is_simulation=is_simulation)
            print(temp_ans)


if __name__ == '__main__':

    fa_path = './chr21.fa'
    bam_path = "./NA19238.chrom21.SLX.maq.SRP000032.2009_07.bam"
    gt_path = './NA19238.gt'
    knncnv(bam_path, fa_path, gt_path=gt_path, iter_num=20)
