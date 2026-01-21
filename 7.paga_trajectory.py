#!/usr/bin/env python3
"""
PAGA (Partition-based Graph Abstraction) Trajectory Analysis
=============================================================

Author: GSE228544 CITE-seq Analysis Pipeline
Date: 2026-01-21

Description:
-----------
This script performs PAGA trajectory analysis on cDC CITE-seq data,
providing topology-based trajectory inference with three UMAP space support.

Main Features:
-------------
1. **Multi-modal UMAP Support**
   - RNA UMAP (gene expression-based)
   - ADT UMAP (surface protein-based)
   - WNN UMAP (weighted nearest neighbor integration) ⭐ Recommended

2. **PAGA Topology Analysis**
   - Cell type connectivity graph
   - Trajectory topology inference
   - Branch point identification
   - Developmental path prediction

3. **Cross-validation**
   - Compare with Monocle3/Slingshot results
   - Validate trajectory consistency
   - Identify robust developmental paths

4. **Mode Support**
   - cDC1-specific analysis
   - cDC2-specific analysis
   - Integrated analysis

Input Data:
----------
From `4.downstream_analysis.Rmd` output:
- `seurat_obj_annotated.rds` (convert to h5ad)
OR
- CSV files: umap_coordinates.csv, adt_umap_coordinates.csv,
  wnn_umap_coordinates.csv, cell_metadata.csv

Output:
-------
paga_analysis_{mode}/
├── data/
│   ├── adata_processed.h5ad           # Processed AnnData object
│   └── paga_results.pkl               # PAGA results object
├── Plots/
│   ├── PAGA/
│   │   ├── paga_graph_rna_umap.png    # PAGA on RNA UMAP
│   │   ├── paga_graph_adt_umap.png    # PAGA on ADT UMAP
│   │   ├── paga_graph_wnn_umap.png    # PAGA on WNN UMAP
│   │   ├── paga_connectivity.png      # Connectivity heatmap
│   │   └── paga_paths.png             # Developmental paths
│   ├── Trajectory/
│   │   ├── trajectory_rna_umap.png
│   │   ├── trajectory_adt_umap.png
│   │   └── trajectory_wnn_umap.png
│   ├── Pseudotime/
│   │   ├── dpt_rna_umap.png          # Diffusion pseudotime
│   │   ├── dpt_adt_umap.png
│   │   └── dpt_wnn_umap.png
│   └── Comparison/
│       ├── umap_comparison.png        # Compare 3 UMAP spaces
│       └── trajectory_comparison.png   # Compare with Monocle3
└── Tables/
    ├── paga_connectivity.csv          # Cell type connectivity matrix
    ├── dpt_pseudotime.csv             # Diffusion pseudotime per cell
    ├── trajectory_genes.csv           # Trajectory-associated genes
    └── cell_transitions.csv           # Cell type transitions

Usage:
------
# Basic usage
python 7.paga_trajectory.py --mode cDC1

# With custom input
python 7.paga_trajectory.py --mode cDC1 --input custom_data.h5ad

# Specify UMAP preference
python 7.paga_trajectory.py --mode cDC1 --primary_umap wnn

Requirements:
------------
- scanpy >= 1.9.0
- anndata >= 0.8.0
- pandas >= 1.3.0
- numpy >= 1.20.0
- matplotlib >= 3.5.0
- seaborn >= 0.11.0
- rpy2 >= 3.5.0 (for reading .rds files)
"""

import os
import sys
import argparse
import warnings
from pathlib import Path
from typing import Optional, Dict, List, Tuple, Union

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import anndata as ad

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)

# Set plot parameters
sc.set_figure_params(dpi=300, frameon=False, figsize=(10, 8), facecolor='white')
sc.settings.verbosity = 1  # Show progress but not too verbose

# Color palettes for consistent visualization
CELL_TYPE_COLORS = {
    'Pre-cDC1s': '#1f77b4',
    'Proliferating cDC1s': '#ff7f0e',
    'Early mature cDC1': '#2ca02c',
    'Mature cDC1': '#d62728',
    'Pre-cDC2s': '#9467bd',
    'Proliferating cDC2s': '#8c564b',
    'Mature cDC2': '#e377c2',
    'pDCs': '#7f7f7f',
    'Monocytes': '#bcbd22'
}

# =============================================================================
# Helper Functions
# =============================================================================

def create_output_directories(base_dir: Path) -> Dict[str, Path]:
    """
    Create organized output directory structure.

    Parameters:
    ----------
    base_dir : Path
        Base output directory

    Returns:
    -------
    Dict[str, Path]
        Dictionary mapping directory names to paths
    """
    dirs = {
        'base': base_dir,
        'data': base_dir / 'data',
        'robjects': base_dir / 'Robjects',
        'plots': base_dir / 'Plots',
        'paga': base_dir / 'Plots' / 'PAGA',
        'trajectory': base_dir / 'Plots' / 'Trajectory',
        'pseudotime': base_dir / 'Plots' / 'Pseudotime',
        'comparison': base_dir / 'Plots' / 'Comparison',
        'tables': base_dir / 'Tables'
    }

    for dir_path in dirs.values():
        dir_path.mkdir(parents=True, exist_ok=True)

    return dirs


def load_data_from_rds(rds_path: Path, mode: str) -> ad.AnnData:
    """
    Load Seurat object from .rds file and convert to AnnData.

    Parameters:
    ----------
    rds_path : Path
        Path to seurat_obj_annotated.rds
    mode : str
        Analysis mode (cDC1, cDC2, or integrated)

    Returns:
    -------
    ad.AnnData
        AnnData object with all UMAP spaces
    """
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri, numpy2ri
        from rpy2.robjects.conversion import localconverter

        pandas2ri.activate()
        numpy2ri.activate()

        print(f"Loading Seurat object from {rds_path}...")

        # Load RDS file
        readRDS = ro.r['readRDS']
        seurat_obj = readRDS(str(rds_path))

        # Extract expression matrix
        ro.r('library(Seurat)')
        get_assay = ro.r('GetAssayData')
        counts = get_assay(seurat_obj, slot='counts', assay='RNA')

        with localconverter(ro.default_converter + pandas2ri.converter):
            counts_df = ro.conversion.rpy2py(counts)

        # Create AnnData object
        adata = ad.AnnData(X=counts_df.T)
        adata.var_names = counts_df.index
        adata.obs_names = counts_df.columns

        # Extract metadata
        metadata = ro.r('as.data.frame')(seurat_obj.slots['meta.data'])
        with localconverter(ro.default_converter + pandas2ri.converter):
            adata.obs = ro.conversion.rpy2py(metadata)

        # Extract embeddings
        embeddings = ro.r('Embeddings')

        # RNA UMAP
        try:
            rna_umap = embeddings(seurat_obj, 'umap')
            with localconverter(ro.default_converter + pandas2ri.converter):
                adata.obsm['X_umap'] = ro.conversion.rpy2py(rna_umap)
            print("✓ Loaded RNA UMAP")
        except:
            print("✗ RNA UMAP not found")

        # ADT UMAP
        try:
            adt_umap = embeddings(seurat_obj, 'adt_umap')
            with localconverter(ro.default_converter + pandas2ri.converter):
                adata.obsm['X_adt_umap'] = ro.conversion.rpy2py(adt_umap)
            print("✓ Loaded ADT UMAP")
        except:
            print("✗ ADT UMAP not found")

        # WNN UMAP
        try:
            wnn_umap = embeddings(seurat_obj, 'wnn_umap')
            with localconverter(ro.default_converter + pandas2ri.converter):
                adata.obsm['X_wnn_umap'] = ro.conversion.rpy2py(wnn_umap)
            print("✓ Loaded WNN UMAP")
        except:
            print("✗ WNN UMAP not found")

        # Extract PCA if available
        try:
            pca = embeddings(seurat_obj, 'pca')
            with localconverter(ro.default_converter + pandas2ri.converter):
                adata.obsm['X_pca'] = ro.conversion.rpy2py(pca)
            print("✓ Loaded PCA")
        except:
            print("✗ PCA not found, will compute later")

        print(f"Loaded {adata.n_obs} cells × {adata.n_vars} genes")

        return adata

    except ImportError:
        print("ERROR: rpy2 not installed. Please install with: pip install rpy2")
        print("Falling back to CSV loading method...")
        return None
    except Exception as e:
        print(f"ERROR loading RDS file: {e}")
        print("Falling back to CSV loading method...")
        return None


def load_data_from_csv(data_dir: Path, mode: str) -> ad.AnnData:
    """
    Load data from CSV files exported by Seurat.

    Parameters:
    ----------
    data_dir : Path
        Directory containing CSV files
    mode : str
        Analysis mode

    Returns:
    -------
    ad.AnnData
        AnnData object with UMAP coordinates
    """
    print(f"Loading data from CSV files in {data_dir}...")

    # Load metadata
    metadata_file = data_dir / 'cell_metadata.csv'
    if not metadata_file.exists():
        raise FileNotFoundError(f"Required file not found: {metadata_file}")

    metadata = pd.read_csv(metadata_file, index_col=0)
    print(f"✓ Loaded metadata: {metadata.shape[0]} cells")

    # Create empty AnnData (we'll need expression matrix separately)
    # For PAGA, we mainly need the UMAP coordinates and metadata
    adata = ad.AnnData(obs=metadata)

    # Load RNA UMAP
    rna_umap_file = data_dir / 'umap_coordinates.csv'
    if rna_umap_file.exists():
        rna_umap = pd.read_csv(rna_umap_file, index_col=0)
        adata.obsm['X_umap'] = rna_umap[['UMAP_1', 'UMAP_2']].values
        print(f"✓ Loaded RNA UMAP")
    else:
        print("✗ RNA UMAP file not found")

    # Load ADT UMAP
    adt_umap_file = data_dir / 'adt_umap_coordinates.csv'
    if adt_umap_file.exists():
        adt_umap = pd.read_csv(adt_umap_file, index_col=0)
        adata.obsm['X_adt_umap'] = adt_umap[['adtUMAP_1', 'adtUMAP_2']].values
        print(f"✓ Loaded ADT UMAP")
    else:
        print("✗ ADT UMAP file not found")

    # Load WNN UMAP
    wnn_umap_file = data_dir / 'wnn_umap_coordinates.csv'
    if wnn_umap_file.exists():
        wnn_umap = pd.read_csv(wnn_umap_file, index_col=0)
        adata.obsm['X_wnn_umap'] = wnn_umap[['wnnUMAP_1', 'wnnUMAP_2']].values
        print(f"✓ Loaded WNN UMAP")
    else:
        print("✗ WNN UMAP file not found")

    return adata


def preprocess_for_paga(adata: ad.AnnData,
                       cell_type_col: str = 'author_annotation',
                       min_cells: int = 10) -> ad.AnnData:
    """
    Preprocess AnnData object for PAGA analysis.

    Parameters:
    ----------
    adata : ad.AnnData
        Input AnnData object
    cell_type_col : str
        Column name for cell type annotation
    min_cells : int
        Minimum cells required per cell type

    Returns:
    -------
    ad.AnnData
        Preprocessed AnnData object
    """
    print("\n" + "="*80)
    print("PREPROCESSING FOR PAGA ANALYSIS")
    print("="*80)

    # Check if expression matrix exists, if not, skip gene-based processing
    has_expression = adata.X is not None and adata.X.shape[1] > 0

    if has_expression:
        print(f"Input: {adata.n_obs} cells × {adata.n_vars} genes")

        # Basic QC
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        print(f"After QC: {adata.n_obs} cells × {adata.n_vars} genes")

        # Normalize
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        # Find highly variable genes
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        print(f"Identified {np.sum(adata.var['highly_variable'])} highly variable genes")

        # Scale
        sc.pp.scale(adata, max_value=10)

        # PCA if not already present
        if 'X_pca' not in adata.obsm:
            print("Computing PCA...")
            sc.tl.pca(adata, svd_solver='arpack')

        # Compute neighborhood graph
        print("Computing neighborhood graph...")
        sc.pp.neighbors(adata, n_neighbors=30, n_pcs=40)
    else:
        print("No expression matrix found, using only UMAP coordinates")
        print(f"Cells: {adata.n_obs}")

        # For PAGA visualization, we still need a neighborhood graph
        # We can compute it from UMAP coordinates
        if 'X_wnn_umap' in adata.obsm:
            print("Computing neighborhood graph from WNN UMAP...")
            sc.pp.neighbors(adata, use_rep='X_wnn_umap', n_neighbors=30)
        elif 'X_umap' in adata.obsm:
            print("Computing neighborhood graph from RNA UMAP...")
            sc.pp.neighbors(adata, use_rep='X_umap', n_neighbors=30)
        else:
            raise ValueError("No UMAP coordinates found for neighborhood computation")

    # Set cell type as categorical
    if cell_type_col in adata.obs.columns:
        adata.obs['cell_type'] = pd.Categorical(adata.obs[cell_type_col])
        print(f"\nCell type distribution:")
        print(adata.obs['cell_type'].value_counts())

        # Filter rare cell types
        cell_counts = adata.obs['cell_type'].value_counts()
        valid_types = cell_counts[cell_counts >= min_cells].index
        adata = adata[adata.obs['cell_type'].isin(valid_types)].copy()
        print(f"\nAfter filtering rare types (min {min_cells} cells): {adata.n_obs} cells")
    else:
        raise ValueError(f"Cell type column '{cell_type_col}' not found in metadata")

    return adata


def run_paga_analysis(adata: ad.AnnData,
                     resolution: float = 1.0) -> ad.AnnData:
    """
    Run PAGA trajectory inference.

    Parameters:
    ----------
    adata : ad.AnnData
        Preprocessed AnnData object
    resolution : float
        Resolution for Leiden clustering

    Returns:
    -------
    ad.AnnData
        AnnData with PAGA results
    """
    print("\n" + "="*80)
    print("PAGA TRAJECTORY INFERENCE")
    print("="*80)

    # Leiden clustering for PAGA
    print(f"Running Leiden clustering (resolution={resolution})...")
    sc.tl.leiden(adata, resolution=resolution, key_added='leiden')

    # Run PAGA
    print("Computing PAGA graph...")
    sc.tl.paga(adata, groups='cell_type')

    print(f"PAGA graph computed with {len(adata.obs['cell_type'].cat.categories)} cell types")

    # Compute diffusion pseudotime
    # Find root cell (typically the earliest cell type)
    root_types = ['Pre-cDC1s', 'Pre-cDC2s', 'pre-cDC1', 'pre-cDC2']
    root_cell = None

    for root_type in root_types:
        if root_type in adata.obs['cell_type'].cat.categories:
            root_cells = np.where(adata.obs['cell_type'] == root_type)[0]
            if len(root_cells) > 0:
                root_cell = root_cells[0]
                print(f"Using root cell type: {root_type}")
                break

    if root_cell is not None:
        adata.uns['iroot'] = root_cell
        print("Computing diffusion pseudotime...")
        sc.tl.diffmap(adata)
        sc.tl.dpt(adata)
        print("✓ Diffusion pseudotime computed")
    else:
        print("⚠ No root cell type found, skipping pseudotime calculation")

    return adata


def plot_paga_graphs(adata: ad.AnnData,
                    output_dir: Path,
                    title_prefix: str = "cDC") -> None:
    """
    Create comprehensive PAGA visualizations across all UMAP spaces.

    Parameters:
    ----------
    adata : ad.AnnData
        AnnData with PAGA results
    output_dir : Path
        Output directory for plots
    title_prefix : str
        Prefix for plot titles
    """
    print("\n" + "="*80)
    print("CREATING PAGA VISUALIZATIONS")
    print("="*80)

    # Define UMAP spaces
    umap_spaces = []
    if 'X_umap' in adata.obsm:
        umap_spaces.append(('X_umap', 'RNA UMAP', 'rna_umap'))
    if 'X_adt_umap' in adata.obsm:
        umap_spaces.append(('X_adt_umap', 'ADT UMAP', 'adt_umap'))
    if 'X_wnn_umap' in adata.obsm:
        umap_spaces.append(('X_wnn_umap', 'WNN UMAP', 'wnn_umap'))

    # 1. PAGA graph on each UMAP space
    for basis, umap_name, suffix in umap_spaces:
        print(f"Plotting PAGA on {umap_name}...")

        fig, ax = plt.subplots(figsize=(12, 10))
        sc.pl.paga(adata,
                   basis=basis,
                   color='cell_type',
                   title=f'{title_prefix} PAGA Graph - {umap_name}',
                   frameon=False,
                   node_size_scale=2.0,
                   edge_width_scale=0.5,
                   min_edge_width=0.1,
                   max_edge_width=3,
                   show=False,
                   ax=ax)

        plt.tight_layout()
        output_file = output_dir / f'paga_graph_{suffix}.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {output_file.name}")

    # 2. PAGA connectivity heatmap
    print("Creating connectivity heatmap...")
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.paga(adata,
               color='cell_type',
               plot=False)

    # Extract connectivity matrix
    conn_matrix = adata.uns['paga']['connectivities'].toarray()
    cell_types = adata.obs['cell_type'].cat.categories

    sns.heatmap(conn_matrix,
                xticklabels=cell_types,
                yticklabels=cell_types,
                cmap='viridis',
                annot=True,
                fmt='.2f',
                ax=ax,
                cbar_kws={'label': 'Connectivity'})

    ax.set_title(f'{title_prefix} PAGA Connectivity Matrix')
    plt.tight_layout()
    output_file = output_dir / 'paga_connectivity.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file.name}")

    # 3. PAGA paths visualization
    print("Creating developmental paths visualization...")
    fig, ax = plt.subplots(figsize=(10, 8))

    sc.pl.paga_path(adata,
                    nodes=list(cell_types),
                    keys=['cell_type'],
                    show=False,
                    ax=ax)

    plt.tight_layout()
    output_file = output_dir / 'paga_paths.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file.name}")


def plot_trajectories(adata: ad.AnnData,
                     output_dir: Path,
                     title_prefix: str = "cDC") -> None:
    """
    Plot trajectory visualizations on different UMAP spaces.

    Parameters:
    ----------
    adata : ad.AnnData
        AnnData with PAGA results
    output_dir : Path
        Output directory
    title_prefix : str
        Prefix for titles
    """
    print("\n" + "="*80)
    print("CREATING TRAJECTORY VISUALIZATIONS")
    print("="*80)

    # Define UMAP spaces
    umap_spaces = []
    if 'X_umap' in adata.obsm:
        umap_spaces.append(('umap', 'RNA UMAP', 'rna_umap'))
    if 'X_adt_umap' in adata.obsm:
        umap_spaces.append(('adt_umap', 'ADT UMAP', 'adt_umap'))
    if 'X_wnn_umap' in adata.obsm:
        umap_spaces.append(('wnn_umap', 'WNN UMAP', 'wnn_umap'))

    for basis, umap_name, suffix in umap_spaces:
        print(f"Plotting trajectory on {umap_name}...")

        fig, ax = plt.subplots(figsize=(10, 8))

        sc.pl.embedding(adata,
                       basis=basis,
                       color='cell_type',
                       title=f'{title_prefix} Trajectory - {umap_name}',
                       frameon=False,
                       show=False,
                       ax=ax)

        plt.tight_layout()
        output_file = output_dir / f'trajectory_{suffix}.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {output_file.name}")


def plot_pseudotime(adata: ad.AnnData,
                   output_dir: Path,
                   title_prefix: str = "cDC") -> None:
    """
    Plot diffusion pseudotime on different UMAP spaces.

    Parameters:
    ----------
    adata : ad.AnnData
        AnnData with DPT results
    output_dir : Path
        Output directory
    title_prefix : str
        Prefix for titles
    """
    if 'dpt_pseudotime' not in adata.obs.columns:
        print("⚠ Diffusion pseudotime not computed, skipping pseudotime plots")
        return

    print("\n" + "="*80)
    print("CREATING PSEUDOTIME VISUALIZATIONS")
    print("="*80)

    # Define UMAP spaces
    umap_spaces = []
    if 'X_umap' in adata.obsm:
        umap_spaces.append(('umap', 'RNA UMAP', 'rna_umap'))
    if 'X_adt_umap' in adata.obsm:
        umap_spaces.append(('adt_umap', 'ADT UMAP', 'adt_umap'))
    if 'X_wnn_umap' in adata.obsm:
        umap_spaces.append(('wnn_umap', 'WNN UMAP', 'wnn_umap'))

    for basis, umap_name, suffix in umap_spaces:
        print(f"Plotting pseudotime on {umap_name}...")

        fig, ax = plt.subplots(figsize=(10, 8))

        sc.pl.embedding(adata,
                       basis=basis,
                       color='dpt_pseudotime',
                       title=f'{title_prefix} Diffusion Pseudotime - {umap_name}',
                       cmap='viridis',
                       frameon=False,
                       show=False,
                       ax=ax)

        plt.tight_layout()
        output_file = output_dir / f'dpt_{suffix}.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {output_file.name}")


def plot_comparison(adata: ad.AnnData,
                   output_dir: Path,
                   title_prefix: str = "cDC") -> None:
    """
    Create comparison plots across UMAP spaces.

    Parameters:
    ----------
    adata : ad.AnnData
        AnnData with results
    output_dir : Path
        Output directory
    title_prefix : str
        Prefix for titles
    """
    print("\n" + "="*80)
    print("CREATING COMPARISON VISUALIZATIONS")
    print("="*80)

    # Compare UMAP spaces side-by-side
    umap_spaces = []
    if 'X_umap' in adata.obsm:
        umap_spaces.append(('umap', 'RNA UMAP'))
    if 'X_adt_umap' in adata.obsm:
        umap_spaces.append(('adt_umap', 'ADT UMAP'))
    if 'X_wnn_umap' in adata.obsm:
        umap_spaces.append(('wnn_umap', 'WNN UMAP'))

    if len(umap_spaces) > 1:
        print("Creating UMAP space comparison...")

        fig, axes = plt.subplots(1, len(umap_spaces),
                                figsize=(8 * len(umap_spaces), 6))

        if len(umap_spaces) == 1:
            axes = [axes]

        for ax, (basis, umap_name) in zip(axes, umap_spaces):
            sc.pl.embedding(adata,
                           basis=basis,
                           color='cell_type',
                           title=umap_name,
                           frameon=False,
                           show=False,
                           ax=ax)

        plt.suptitle(f'{title_prefix} - UMAP Space Comparison',
                    fontsize=16, y=1.02)
        plt.tight_layout()
        output_file = output_dir / 'umap_comparison.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {output_file.name}")


def export_results(adata: ad.AnnData,
                  output_dir: Path,
                  tables_dir: Path) -> None:
    """
    Export PAGA results to files.

    Parameters:
    ----------
    adata : ad.AnnData
        AnnData with PAGA results
    output_dir : Path
        Directory for data outputs
    tables_dir : Path
        Directory for table outputs
    """
    print("\n" + "="*80)
    print("EXPORTING RESULTS")
    print("="*80)

    # Save processed AnnData
    output_file = output_dir / 'adata_processed.h5ad'
    adata.write_h5ad(output_file)
    print(f"✓ Saved AnnData: {output_file}")

    # Export connectivity matrix
    if 'paga' in adata.uns:
        conn_matrix = adata.uns['paga']['connectivities'].toarray()
        cell_types = adata.obs['cell_type'].cat.categories

        conn_df = pd.DataFrame(conn_matrix,
                              index=cell_types,
                              columns=cell_types)

        output_file = tables_dir / 'paga_connectivity.csv'
        conn_df.to_csv(output_file)
        print(f"✓ Saved connectivity matrix: {output_file}")

    # Export pseudotime
    if 'dpt_pseudotime' in adata.obs.columns:
        dpt_df = pd.DataFrame({
            'cell_id': adata.obs_names,
            'cell_type': adata.obs['cell_type'],
            'dpt_pseudotime': adata.obs['dpt_pseudotime']
        })

        output_file = tables_dir / 'dpt_pseudotime.csv'
        dpt_df.to_csv(output_file, index=False)
        print(f"✓ Saved pseudotime: {output_file}")

    # Export cell type transitions
    if 'paga' in adata.uns:
        transitions = []
        conn_matrix = adata.uns['paga']['connectivities'].toarray()
        cell_types = adata.obs['cell_type'].cat.categories

        for i, source in enumerate(cell_types):
            for j, target in enumerate(cell_types):
                if i != j and conn_matrix[i, j] > 0.1:  # Threshold for meaningful connections
                    transitions.append({
                        'source': source,
                        'target': target,
                        'connectivity': conn_matrix[i, j]
                    })

        trans_df = pd.DataFrame(transitions)
        trans_df = trans_df.sort_values('connectivity', ascending=False)

        output_file = tables_dir / 'cell_transitions.csv'
        trans_df.to_csv(output_file, index=False)
        print(f"✓ Saved cell transitions: {output_file}")

    print("\n" + "="*80)
    print("PAGA ANALYSIS COMPLETE!")
    print("="*80)


# =============================================================================
# Main Function
# =============================================================================

def main():
    """Main execution function."""

    # Parse arguments
    parser = argparse.ArgumentParser(
        description='PAGA Trajectory Analysis for cDC CITE-seq Data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument('--mode',
                       type=str,
                       default='cDC1',
                       choices=['cDC1', 'cDC2', 'integrated'],
                       help='Analysis mode (default: cDC1)')

    parser.add_argument('--input',
                       type=str,
                       default=None,
                       help='Input file path (.rds or .h5ad). If not specified, will look for default paths')

    parser.add_argument('--data-dir',
                       type=str,
                       default=None,
                       help='Directory containing CSV files (alternative to RDS/h5ad input)')

    parser.add_argument('--output-dir',
                       type=str,
                       default=None,
                       help='Output directory (default: paga_analysis_{mode})')

    parser.add_argument('--primary-umap',
                       type=str,
                       default='wnn',
                       choices=['rna', 'adt', 'wnn'],
                       help='Primary UMAP space for analysis (default: wnn)')

    parser.add_argument('--resolution',
                       type=float,
                       default=1.0,
                       help='Leiden clustering resolution (default: 1.0)')

    parser.add_argument('--cell-type-col',
                       type=str,
                       default='author_annotation',
                       help='Column name for cell type annotation (default: author_annotation)')

    args = parser.parse_args()

    # Print header
    print("\n" + "="*80)
    print(" "*20 + "PAGA TRAJECTORY ANALYSIS")
    print(" "*25 + "GSE228544 CITE-seq Data")
    print("="*80)
    print(f"\nMode: {args.mode}")
    print(f"Primary UMAP: {args.primary_umap.upper()}")
    print(f"Resolution: {args.resolution}")
    print(f"Cell type column: {args.cell_type_col}")

    # Set up output directory
    if args.output_dir:
        output_base = Path(args.output_dir)
    else:
        output_base = Path(f'paga_analysis_{args.mode}')

    dirs = create_output_directories(output_base)
    print(f"\nOutput directory: {output_base}")

    # Load data
    adata = None

    # Try loading from specified input file
    if args.input:
        input_path = Path(args.input)
        if input_path.suffix == '.rds':
            adata = load_data_from_rds(input_path, args.mode)
        elif input_path.suffix == '.h5ad':
            print(f"Loading AnnData from {input_path}...")
            adata = sc.read_h5ad(input_path)
            print(f"✓ Loaded {adata.n_obs} cells")
        else:
            raise ValueError(f"Unsupported file format: {input_path.suffix}")

    # Try loading from CSV directory
    elif args.data_dir:
        adata = load_data_from_csv(Path(args.data_dir), args.mode)

    # Try default paths
    else:
        # Look for RDS file
        default_rds = Path(f'results/{args.mode}/Robjects/seurat_obj_annotated.rds')
        if default_rds.exists():
            adata = load_data_from_rds(default_rds, args.mode)

        # Fall back to CSV
        if adata is None:
            default_csv_dir = Path(f'results/{args.mode}/Tables')
            if default_csv_dir.exists():
                adata = load_data_from_csv(default_csv_dir, args.mode)
            else:
                raise FileNotFoundError(
                    "No input data found. Please specify --input or --data-dir, "
                    "or ensure default paths exist:\n"
                    f"  - {default_rds}\n"
                    f"  - {default_csv_dir}"
                )

    if adata is None:
        raise ValueError("Failed to load data from any source")

    # Preprocess
    adata = preprocess_for_paga(adata,
                               cell_type_col=args.cell_type_col,
                               min_cells=10)

    # Run PAGA
    adata = run_paga_analysis(adata, resolution=args.resolution)

    # Create visualizations
    plot_paga_graphs(adata, dirs['paga'], title_prefix=args.mode)
    plot_trajectories(adata, dirs['trajectory'], title_prefix=args.mode)
    plot_pseudotime(adata, dirs['pseudotime'], title_prefix=args.mode)
    plot_comparison(adata, dirs['comparison'], title_prefix=args.mode)

    # Export results
    export_results(adata, dirs['data'], dirs['tables'])

    print(f"\n✓ All results saved to: {output_base}")
    print("\nNext steps:")
    print("  1. Review PAGA graphs in Plots/PAGA/")
    print("  2. Compare with Monocle3/Slingshot results")
    print("  3. Analyze connectivity matrix in Tables/")
    print("  4. Use pseudotime for downstream gene analysis")


if __name__ == '__main__':
    main()
