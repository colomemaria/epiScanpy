import multiprocessing

import numpy as np
import pandas as pd

from Bio import SeqIO
from pyjaspar import jaspardb



def get_sequence(chrom, start, end, genome):

    seq_record = genome[chrom]
    seq = seq_record.seq[start:end]

    return(str(seq).upper())



def get_chrom_sequence(chrom, genome):

    seq_record = genome[chrom]

    return(str(seq_record.seq).upper())



def get_peak_sequences(adata, genome):

    genome = SeqIO.to_dict(SeqIO.parse(genome, format="fasta"))

    tmp = []
    for row in adata.var.to_dict("records"):

        try:
            seq = get_sequence(chrom=row["chr"], start=row["start"], end=row["stop"], genome=genome)
        except KeyError:
            seq = np.nan

        tmp.append(seq)

    adata.var["seq"] = tmp



def fetch_motifs(tf_names, release="JASPAR2022"):

    # connect to jaspar
    jd_obj = jaspardb(release=release)

    # fetch motifs data
    motifs = jd_obj.fetch_motifs(
        collection="CORE",
        tax_group="vertebrates",
        tf_name=tf_names
    )
    
    names = sorted({m.name for m in motifs})
    idx_dict = {name: i for i, name in enumerate(names)}

    motifs = sorted(motifs, key=lambda m: idx_dict[m.name])

    return motifs



def background_distribution(adata):

    seq_concat = "".join([seq for seq in adata.var.seq if not pd.isna(seq)])

    counts = np.zeros(4)

    idx = {
        "A": 0,
        "a": 0,
        "C": 1,
        "c": 1,
        "G": 2,
        "g": 2,
        "T": 3,
        "t": 3
    }

    for base in seq_concat:
        if base in idx:
            counts[idx[base]] += 1

    bg = counts / counts.sum()
    bg = {base: frac for base, frac in zip(["A", "C", "G", "T"], bg)}

    return bg



def search_motifs(adata, tf_names, release="JASPAR2022", pseudocounts=1, threshold_fpr=1e-6, n_jobs=1):

    # fetch motif data
    motifs = fetch_motifs(tf_names=tf_names, release=release)

    # estimate background distribution
    bg = background_distribution(adata)

    n_jobs = multiprocessing.cpu_count() if n_jobs == -1 else n_jobs

    sem = multiprocessing.Semaphore(n_jobs)

    manager = multiprocessing.Manager()
    q = manager.Queue()

    producers = [MotifSearcher(adata, motif, motif_idx, bg, pseudocounts, threshold_fpr, q, sem) for motif_idx, motif in enumerate(motifs)]

    for p in producers:
        sem.acquire()
        p.start()

    for p in producers:
        p.join()

    results = []
    while not q.empty():
        for res in q.get():
            results.append(res)

    print(f"Number of binding sites: {len(results)}")

    adata.uns["motif_search"] = {}
    adata.uns["motif_search"]["tf_motifs"] = motifs
    adata.uns["motif_search"]["results"] = results



class MotifSearcher(multiprocessing.Process):

    def __init__(self, adata, motif, motif_idx, bg, pseudocounts, threshold_fpr, q, sem):
        super().__init__()

        self.adata = adata

        self.motif = motif
        self.motif_idx = motif_idx

        self.bg = bg

        self.pseudocounts = pseudocounts
        self.threshold_fpr = threshold_fpr

        self.q = q
        self.sem = sem


    def run(self) -> None:
        
        # add pseudocount
        self.motif.pseudocounts = self.pseudocounts

        # highest possible score for a motif
        for pos, score in self.motif.pssm.search(self.motif.consensus):
            max_score = score

        # determine threshold
        score_distribution = self.motif.pssm.distribution(background=self.bg)
        threshold = score_distribution.threshold_fpr(self.threshold_fpr)

        results = []
        # check sequences for motif
        for feature_idx, seq in enumerate(self.adata.var.seq.tolist()):

            if pd.isna(seq):
                continue
            
            for pos, score in self.motif.pssm.search(seq, threshold=threshold):
                
                results.append(
                    {
                        "motif_idx": self.motif_idx,
                        "feature_idx": feature_idx,
                        "pos": pos,
                        "score": score/max_score,
                        "motif_length": self.motif.length
                    }
                )

        self.q.put(results)

        self.sem.release()


def binary_motif_mtx(adata):

    motif_presence_mtx = np.zeros((len(adata.uns["motif_search"]["tf_motifs"]), adata.n_vars))

    for tfbs in adata.uns["motif_search"]["results"]:
        motif_presence_mtx[tfbs["motif_idx"], tfbs["feature_idx"]] = 1

    adata.uns["motif_search"]["motif_presence_mtx"] = motif_presence_mtx