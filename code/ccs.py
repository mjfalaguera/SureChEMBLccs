#!/usr/bin/env python

"""
ccs.py: Identify the Core Chemical Structure (CCS) for a collection of distinct molecules (see Falaguera & Mestres 2021).

>> import ccs
>> idx_smiles = {'CHEMBL293613': 'O=C(O)Cn1c(=O)c(=O)[nH]c2cc([N+](=O)[O-])c(-n3ccc(C(=O)O)c3)cc21',
 'CHEMBL429540': 'O=C(O)Cn1c(=O)c(=O)[nH]c2cc([N+](=O)[O-])c(-n3ccc(CNCCCc4ccccc4)c3)cc21',
 'CHEMBL63861': 'O=C(O)Cn1c(=O)c(=O)[nH]c2cc([N+](=O)[O-])c(-n3ccc(CN4CCN(Cc5ccccc5)CC4)c3)cc21',
 'CHEMBL12465': 'NCc1ccn(-c2cc3c(cc2[N+](=O)[O-])[nH]c(=O)c(=O)n3CC(=O)O)c1',
 'CHEMBL67828': 'O=C(O)Cn1c(=O)c(=O)[nH]c2cc([N+](=O)[O-])c(-n3ccc(CN4CCN(CCc5ccccc5)CC4)c3)cc21',
 'CHEMBL416685': 'O=C(O)Cn1c(=O)c(=O)[nH]c2cc([N+](=O)[O-])c(-n3ccc(CN4CCC(c5ccccc5)CC4)c3)cc21',
 'CHEMBL268284': 'O=C(O)Cn1c(=O)c(=O)[nH]c2cc([N+](=O)[O-])c(-n3cccc3)cc21'}
>> CCS = getCCSforMolecules(idx_smiles=idx_smiles)
"""

__author__ = "Maria J. Falaguera"
__date__ = "25 Jan 2024"

import gzip
import itertools

from matplotlib import pyplot as plt
import numpy as np
import random
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import rdFMCS
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

# change input paring function to access group DB


def getInfoForMolecules(idx_smiles, minWeight=150, maxWeight=1000, minRings=1):
    """
    Convert dictionary of smiles into molecules.

    Arg:
        idx_smiles (dict):  Dict mapping molecule identifier to SMILES
        minWeight (int):    Min. weight to admit the molecule
        maxWeight (int):    Max. weight to admit the molecule
        minRings (int):     Min. number of rings to admit the molecule

    Returns:
        idx_mol (dict):     Dict mapping molecule identifier to:
                             - mol: RDKit mol object
                             - fp: RDKit Morgan fingerprint with radius=2
                             - smiles: Canonized RDKit SMILES
                             - mw: Molecular weight (Da)
                             - rings: Number of rings
    """

    idx_mol = {}
    for idx, smiles in idx_smiles.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            idx_mol[idx] = {"mol": None, "idx": idx, "fp": None, "smiles": smiles}
            continue
        Chem.GetSSSR(mol)  # initialize RingInfo (to avoid error)
        mw = Descriptors.MolWt(mol)
        if mw < float(minWeight):
            continue  # discard tiny molecules
        if mw > float(maxWeight):
            continue  # discard huge molecules (for memory reasons)
        rings = len(Chem.GetSSSR(mol))
        if rings < float(minRings):
            continue  # discard non-cyclic molecules
        smiles = Chem.MolToSmiles(mol)  # canonize SMILES
        fp = AllChem.GetMorganFingerprint(mol, 2)
        idx_mol[idx] = {
            "mol": mol,
            "idx": idx,
            "fp": fp,
            "smiles": smiles,
            "mw": mw,
            "rings": rings,
        }

    return idx_mol


def getSampleForMolecules(idx_mol, sample=1000):

    sampleMols = set([idx for idx in idx_mol if idx_mol[idx]["mol"] is not None])

    if len(sampleMols) > sample:
        sampleMols = set(random.sample(list(sampleMols), sample))

    sampleMols = {idx: idx_mol[idx] for idx in sampleMols}

    return sampleMols


def getSimilarityDistForMolecules(idx_mol, plot=False):
    """
    Get molecular pairwise similarity distribution for a collection of molecules (Morgan FP with radius=2 and Dice similarity).

    Args:
        idx_mol (dict):          Dict of molecule mapped to RDKit mol object and other data
        plot (bool):             Plot similarity values distribution

    Returns:
        similarities (dict):     Dict mapping molecular pairs to similarity value

    """

    # Calculate similarities
    pair_sim = {}
    for idx1, idx2 in itertools.combinations(list(idx_mol), 2):
        fp1 = idx_mol[idx1]["fp"]
        fp2 = idx_mol[idx2]["fp"]
        if (fp1 is not None) and (fp2 is not None):
            sim = DataStructs.DiceSimilarity(fp1, fp2)
        else:
            sim = np.nan
        pair_sim[tuple(sorted([idx1, idx2]))] = round(sim, 3)

    # Plot
    if plot:
        plt.hist(list(pair_sim.values()), bins=20)
        plt.xlabel("Pairwise similarity")
        plt.ylabel("#Molecular pairs")
        plt.xlim((0, 1))

    return pair_sim


def getMCSforMolecules(idx_mol, pair_sim, sim_cutoff=0.4, timeout=1):
    """
    Get Maximum Common Substructures (MCS) for molecular pairs.

    Args:
        idx_mol (dict):         Dict of molecules mapped to RDKit mol objects and other data
        pair_sim (dict):        Dict mapping each molecular pair to their molecular similarity
        sim_cutoff (float):     Min. similarity to allow MCS calculaion. This makes the MCS calculation faster

    Returns:
        MCS (dict):             Dict mapping calculated MCSs to pairs of molecules originating them
    """

    MCS = {}
    for idx1, idx2 in pair_sim:
        if pair_sim[(idx1, idx2)] >= sim_cutoff:
            mol1 = idx_mol[idx1]["mol"]
            mol2 = idx_mol[idx2]["mol"]
            mcs = rdFMCS.FindMCS(
                [mol1, mol2],
                ringMatchesRingOnly=False,
                completeRingsOnly=True,
                timeout=timeout,
            ).smartsString
            if mcs == "":
                continue  # no MCS found

            try:
                _ = MCS[mcs]
            except KeyError:
                MCS[mcs] = {"pairs": set()}
            MCS[mcs]["pairs"].update([(idx1, idx2)])

    return MCS


def getMCScoverage(MCS, idx_mol, plot=False):
    """
    Calculate MCS coverage defined as the no. molecules that have this MCS as a substructure divided by the total amount of molecules.

    Args:
        MCS (dict):     Dict mapping calculated MCSs to data
        idx_mol (dict): Dict of molecules mapped to RDKit mol object and other data
        plot (bool):    Plot distribution of coverage values

    Reurns:
        MCS (dict):     Dict mapping calculated MCSs to:
                         - pairs of molecules originating them
                         - distinct molecules covered by them
                         - coverage value
    """

    # Get molecules covered by MCS
    for mcs in MCS:

        # get distinct mols originating (and so covered) by the MCS
        MCS[mcs]["mols"] = set()
        for idx1, idx2 in MCS[mcs]["pairs"]:
            MCS[mcs]["mols"].update([idx1, idx2])

        # generated MCS mol object
        mcsMol = Chem.MolFromSmarts(mcs)
        MCS[mcs]["mol"] = mcsMol
        if mcsMol is not None:

            # get distinct mols not originating the MCS but covered by it
            for idx in set(idx_mol).difference(MCS[mcs]["mols"]):
                if idx_mol[idx]["mol"] is not None:
                    if idx_mol[idx]["mol"].HasSubstructMatch(mcsMol):
                        MCS[mcs]["mols"].update([idx])

            # number of mols covered by MCS
            MCS[mcs]["n_mols"] = len(MCS[mcs]["mols"])

            # coverage = mols covered by MCS / total mols
            MCS[mcs]["cov"] = float(MCS[mcs]["n_mols"]) / float(len(idx_mol))

    # Plot distribution of coverage values
    if plot:
        coverages = [MCS[mcs]["cov"] for mcs in MCS]
        plt.hist(coverages, bins=20)
        plt.xlabel("Coverage")
        plt.ylabel("#MCS")
        plt.xlim((0, 1))

    return MCS


def getMCShomogeneity(MCS, pair_sim, plot=False):
    """
    Calculate MCS internal homogeneity defined as the median of the similarity distribution of the molecules covered by the MCS.

    Args:
        MCS (dict):         Dict mapping calculated MCSs to data
        pair_sim (dict):    Dict mapping each molecular pair to their molecular similarity
        plot (bool):        Plot homogeneity values distribution

    Reurns:
        MCS (dict):         Dict mapping calculated MCSs to:
                             - pairs of molecules originating them
                             - distinct molecules covered by them
                             - homogeneity value
    """

    # Calculate MCS internal homogeneity
    for mcs in MCS:
        hom = []
        for idx1, idx2 in itertools.combinations(MCS[mcs]["mols"], 2):
            sim = pair_sim[tuple(sorted([idx1, idx2]))]
            hom.append(sim)
        MCS[mcs]["hom"] = np.median(hom)

    # Plot
    if plot:
        homologies = [MCS[mcs]["hom"] for mcs in MCS]
        plt.hist(homologies, bins=20)
        plt.xlabel("Homology")
        plt.ylabel("#MCS")
        plt.xlim((0, 1))

    return MCS


def getMCSinclusion(MCS, plot=False):
    """
    Calculate MCS inclusion and congeneric values.
    Inclusion = MCSs included as substructures of the MCS.
    Congeneric = MCSs with the MCS as a substructure.

    Args:
        MCS (dict):         Dict mapping calculated MCSs to data
        plot (bool):        Plot distribution of homogeneity values

    Reurns:
        MCS (dict):         Dict mapping calculated MCSs to:
                             - pairs of molecules originating them
                             - distinct molecules covered by them
                             - inclusion value
                             - congeneric value
    """

    # Assess MCS substructure inclusion and congeneric degree
    for mcs in MCS:
        MCS[mcs]["congeneric"] = set()
        MCS[mcs]["included"] = set()

    for mcs1, mcs2 in itertools.combinations(MCS, 2):
        mol1 = MCS[mcs1]["mol"]
        mol2 = MCS[mcs2]["mol"]

        if (mol1 is not None) and (mol2 is not None):
            if mol1.HasSubstructMatch(mol2, useQueryQueryMatches=True):
                MCS[mcs1]["included"].update(
                    [mcs2]
                )  # mol1 contains mol2 as a substructure (mol1 includes mol2)
                MCS[mcs2]["congeneric"].update(
                    [mcs1]
                )  # mol2 contains mol1 as a substructure (mol1 covers mol2)
            if mol2.HasSubstructMatch(mol1, useQueryQueryMatches=True):
                MCS[mcs2]["included"].update([mcs1])
                MCS[mcs1]["congeneric"].update([mcs2])

    n_mcs = len(MCS)
    for mcs in MCS:
        MCS[mcs]["inc"] = float(len(MCS[mcs]["included"])) / float(n_mcs)

    # Plot
    if plot:
        inclusions = [MCS[mcs]["inc"] for mcs in MCS]
        plt.hist(inclusions, bins=20)
        plt.xlabel("Inclusion")
        plt.ylabel("#MCS")
        plt.xlim(0, 1)

    return MCS


def getMCSscore(MCS, idx_mol, pair_sim, plot=False):
    """
    Calculate MCS score = coverage * inclusion * homogeneity.

    Args:
        MCS (dict):         Dict mapping calculated MCSs to data
        plot (bool):        Plot distribution of homogeneity values

    Reurns:
        MCS (dict):         Dict mapping calculated MCSs to:
                             - pairs of molecules originating them
                             - distinct molecules covered by them
                             - score value
    """

    # Score MCS
    cov_weight = 0.8
    hom_weight = 1
    inc_weight = 1

    MCS = getMCScoverage(MCS=MCS, idx_mol=idx_mol)
    MCS = getMCShomogeneity(MCS=MCS, pair_sim=pair_sim)
    MCS = getMCSinclusion(MCS=MCS)

    for mcs in MCS:
        score = round(
            (
                MCS[mcs]["cov"] ** cov_weight
                * MCS[mcs]["hom"] ** hom_weight
                * MCS[mcs]["inc"] ** inc_weight
            )
            ** (1 / 3),
            4,
        )
        MCS[mcs]["score"] = score

    # Plot
    if plot:
        scores = [MCS[mcs]["score"] for mcs in MCS]
        plt.hist(scores)
        plt.xlabel("Score")
        plt.ylabel("#MCS")
        plt.xlim(0, 1)

    return MCS


def getTopScoringMCS(MCS, quantile=0.95, non_redundant=True):
    """
    Get top scoring MCSs from a list.

    Args:
        MCS (dict):                 Dict mapping calculated MCSs data
        quantile (float):           Quantile representing the score cutoff for the top MCSs selection
        remove_redundant (bool):    If set to True, it collapses low scoring MCSs onto higher scoring ones that contain them as substructures

    Reurns:
        topMCS (dict):      Dict mapping calculated MCSs to data
    """

    # Cutoff
    scores = [MCS[mcs]["score"] for mcs in MCS]

    if not len(scores):  # no MCS found
        return {}

    cutoff = round(np.quantile(scores, quantile), 4)

    # MCSs with score >= cutoff
    topMCS = {mcs: data for mcs, data in MCS.items() if data["score"] >= cutoff}

    # Remove redundant MCS
    if non_redundant:

        # unify MCS
        redundant = set()
        for mcs1 in topMCS:
            if mcs1 not in redundant:
                for mcs2 in set(topMCS).difference(set([mcs1])):
                    # if they are the same MCS, select one
                    if (mcs1 in topMCS[mcs2]["congeneric"]) and (
                        mcs2 in topMCS[mcs1]["congeneric"]
                    ):
                        redundant.update([mcs2])
        topMCS = {mcs: topMCS[mcs] for mcs in set(topMCS).difference(redundant)}

        # remove MCS contained within another MCS
        for mcs in topMCS:
            topMCS[mcs]["congeneric"] = set(topMCS).intersection(
                topMCS[mcs]["congeneric"]
            )

        redundant = set()
        for mcs in sorted(
            topMCS, key=lambda x: len(topMCS[x]["congeneric"]), reverse=True
        ):
            redundant.update(topMCS[mcs]["congeneric"])

        topMCS = {mcs: topMCS[mcs] for mcs in set(topMCS).difference(redundant)}

    return topMCS

    # Remove redundant MCS
    # if non_redundant:
    #    tmp = set()
    #    for mcs1 in sorted(topMCS, key=lambda x: topMCS[x]["score"], reverse=True):
    #        new = True
    #        for mcs2 in tmp:
    #            if (
    #                mcs1 in topMCS[mcs2]["congeneric"]
    #            ):  # mcs1 contains mcs2 as a substructure
    #                new = False
    #        if new:
    #            tmp.update([mcs1])

    #    topMCS = {mcs: topMCS[mcs] for mcs in tmp}

    return topMCS


def drawMCS(MCS, molsPerRow=6):
    """
    Draw grid of MCSs.

    Args:
        MCS (dict):         Dict mapping calculated MCSs to data

    Returns:
        img (rdkit img):    RDKit grid of MCSs labeled with data
    """

    mols = []
    legends = []
    for mcs in sorted(MCS, key=lambda x: MCS[x]["score"], reverse=True):
        mols.append(MCS[mcs]["mol"])
        legends.append(
            "C:{:.1f} | H:{:.1f} | I:{:.1f} | S:{:.2f}".format(
                MCS[mcs]["cov"], MCS[mcs]["hom"], MCS[mcs]["inc"], MCS[mcs]["score"]
            )
        )

    img = Draw.MolsToGridImage(
        mols, legends=legends, molsPerRow=molsPerRow, maxMols=100
    )

    return img


def recoverSimilarMolecules(MCS, idx_mol, pair_sim, sim_cutoff=0.8):
    """
    Recover molecules not covered by a given MCS but highly similar to covered ones.

    Args:
        MCS (dict):         Dict mapping calculated MCSs to data
        idx_mol (dict):     Dict of molecule mapped to SMILES and other data
        pair_sim (dict):    Dict mapping each molecular pair to their molecular similarity
        sim_cutoff (float): Min. similarity to allow molecular recovery

    Reurns:
        MCS (dict):         Dict mapping calculated MCSs to data
    """

    for mcs in MCS:
        coveredMols = set([idx for idx in MCS[mcs]["mols"]])
        toRecoverMols = set(idx_mol).difference(coveredMols)
        recoveredMols = set([])
        for coveredMol in coveredMols:
            for toRecoverMol in toRecoverMols:
                sim = pair_sim[tuple(sorted([coveredMol, toRecoverMol]))]
                if sim >= sim_cutoff:
                    recoveredMols.update([toRecoverMol])
            toRecoverMols = toRecoverMols.difference(recoveredMols)
        coveredMols = coveredMols.union(recoveredMols)
        MCS[mcs]["mols+recov"] = coveredMols

    return MCS


def addFilteredMolecules(MCS, idx_mol):
    for mcs in MCS:
        MCS[mcs]["filteredMolecules"] = {
            idx: idx_mol[idx] for idx in MCS[mcs]["mols+recov"]
        }

    return MCS


def getFilteredMolecules(MCS, idx_mol):
    """
    Get molecules filtered by Core Chemical Structure (CCS).

    Args:
        MCS (dict):     Dict mapping calculated MCSs to data
        idx_mol (dict): Dict of molecule mapped to RDKit mol object and other data
        plot (bool):    Plot grid of filtered molecules

    Returns:
        RDKit grid of filtered mols
    """

    # Get molecules filtered by top MCSs (CCS)
    filtered_idx_mol = {
        idx: idx_mol[idx] for mcs in MCS for idx in MCS[mcs]["mols+recov"]
    }

    return filtered_idx_mol


def getPatentHomogeneity(filteredMolecules, pair_sim, plot=False):
    sims = [
        pair_sim[tuple(sorted([idx1, idx2]))]
        for idx1, idx2 in itertools.combinations(filteredMolecules, 2)
    ]

    if not len(sims):  # no molecules filtered
        return np.nan

    homogeneity = round(np.median(sims), 3)

    # Plot
    if plot:
        plt.hist(sims, bins=20)
        plt.xlabel("Pairwise similarity")
        plt.ylabel("#Molecular pairs")
        plt.xlim((0, 1))

    return homogeneity


def drawMolecules(idx_mol, molsPerRow=10, maxMols=100):
    mols = []
    legends = []
    for idx, mol in idx_mol.items():
        mols.append(mol["mol"])
        legends.append(idx)

    print("Showing only {} molecules:".format(maxMols))
    img = Draw.MolsToGridImage(
        mols, legends=legends, molsPerRow=molsPerRow, maxMols=maxMols
    )
    return img


def getCCSforMolecules(idx_smiles, label=None, plotMCS=False, plotMolecules=False):
    """
    Identify Core Chemical Structure (CCS) of a collection of molecules.

    Args:
        idx_smiles (dcit):  Dict mapping molecule identifier to SMILES
        plot (bool):        Plot grid of MCSs composing the CCS

    Reurns:
        MCS (dict):         Dict mapping CCS to:
                             - pairs:               pairs of molecules originating them
                             - mols:                distinct molecules covered by them
                             - mols+recov:          distinct molecules covered by them and recovered
                             - filteredMolecules:   distinct molecules covered by them and recovered, as a dict with molecule properties
                             - cov:                 coverage value
                             - hom:                 homogeneity value
                             - inc:                 inclusion value
                             - homogeneity:         patent homogeneity value
                             - score:               score value
                             - mol:                 RDKit mol object of the MCS
    """

    idx_mol = getInfoForMolecules(idx_smiles=idx_smiles)
    pair_sim = getSimilarityDistForMolecules(idx_mol=idx_mol)
    MCS = getMCSforMolecules(idx_mol=idx_mol, pair_sim=pair_sim)
    MCS = getMCSscore(MCS=MCS, idx_mol=idx_mol, pair_sim=pair_sim)
    MCS = getTopScoringMCS(MCS=MCS)
    MCS = recoverSimilarMolecules(MCS=MCS, idx_mol=idx_mol, pair_sim=pair_sim)
    MCS = addFilteredMolecules(MCS=MCS, idx_mol=idx_mol)
    filteredMolecules = getFilteredMolecules(MCS=MCS, idx_mol=idx_mol)
    homogeneity = getPatentHomogeneity(
        filteredMolecules=filteredMolecules, pair_sim=pair_sim
    )

    if not len(idx_mol):
        pctFiltered = np.nan
    else:
        pctFiltered = round(float(len(filteredMolecules)) / float(len(idx_mol)), 3)

    if plotMCS:
        drawMCS(MCS=MCS)

    if plotMolecules:
        drawMolecules(idx_mol=filteredMolecules)

    return {
        "label": label,
        "filteredMolecules": filteredMolecules,
        "pctFiltered": pctFiltered,
        "homogeneity": homogeneity,
        "coreChemicalStructure": MCS,
    }


def writeResults(
    results, fileName, coreChemicalStructure=False, molecules=False, molecule_ccs=False
):
    label = results["label"]
    filteredMolecules = results["filteredMolecules"]
    homogeneity = results["homogeneity"]
    pctFiltered = results["pctFiltered"]
    coreChemicalStructure = results["coreChemicalStructure"]

    if not molecule_ccs:
        if coreChemicalStructure:
            with gzip.open(fileName, "wt") as fh:
                for mcs in coreChemicalStructure:
                    fh.write(
                        "{}\t{}\t{}\t{}\n".format(label, homogeneity, pctFiltered, mcs)
                    )

        if molecules:
            with gzip.open(fileName, "wt") as fh:
                for filteredMolecule in filteredMolecules:
                    fh.write(
                        "{}\t{}\t{}\t{}\t{}\n".format(
                            label,
                            homogeneity,
                            pctFiltered,
                            filteredMolecule,
                            filteredMolecules[filteredMolecule],
                        )
                    )

    else:
        with gzip.open(fileName, "wt") as fh:
            for mcs in coreChemicalStructure:
                for idx, mol in coreChemicalStructure[mcs]["filteredMolecules"].items():
                    fh.write(
                        "{}\t{}\t{}\t{}\t{}\t{}\n".format(
                            label, homogeneity, pctFiltered, mcs, idx, mol["smiles"]
                        )
                    )

        print(fileName)
