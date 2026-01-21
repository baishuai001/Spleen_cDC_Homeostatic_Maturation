#!/usr/bin/env python3
"""
RNA Velocity Analysis for cDC Development
=========================================

Author: GSE228544 CITE-seq Analysis Pipeline
Date: 2026-01-21

Description:
-----------
This script performs RNA velocity analysis to infer the directionality and
speed of cellular differentiation in cDC development.

RNA velocity uses the ratio of unspliced to spliced mRNA to predict the
future state of individual cells, providing insights into:
- Direction of cell state transitions
- Speed of differentiation
- Stability of cell states
- Regulatory dynamics

Main Features:
-------------
1. **Velocity Computation**
   - Spliced/unspliced mRNA quantification
   - Velocity vector calculation
   - Phase portrait analysis

2. **Multi-UMAP Support**
   - RNA UMAP visualization
   - ADT UMAP visualization
   - WNN UMAP visualization ⭐ Recommended

3. **Integration with Trajectory Analysis**
   - Compare with Monocle3/Slingshot/PAGA results
   - Validate developmental directions
   - Identify stable vs transitional states

4. **Comprehensive Visualization**
   - Velocity stream plots
   - Phase portraits
   - Speed and coherence analysis
   - Gene-specific dynamics

Input Data:
----------
Required:
- Velocyto .loom file (from velocyto CLI preprocessing)
- UMAP coordinates (from Seurat analysis)
- Cell type annotations

Optional:
- Monocle3/Slingshot pseudotime for validation

Velocyto Preprocessing:
----------------------
Run velocyto on BAM files before using this script:

```bash
# For 10x data
velocyto run10x -m repeat_msk.gtf \\
  /path/to/sample_folder \\
  /path/to/genes.gtf

# Output: sample_folder/velocyto/sample.loom
```

Output:
-------
rna_velocity_{mode}/
├── data/
│   ├── adata_with_velocity.h5ad       # AnnData with velocity
│   └── velocity_genes.txt             # Velocity genes list
├── Plots/
│   ├── Velocity/
│   │   ├── velocity_stream_rna_umap.png    # Velocity on RNA UMAP
│   │   ├── velocity_stream_adt_umap.png    # Velocity on ADT UMAP
│   │   ├── velocity_stream_wnn_umap.png    # Velocity on WNN UMAP ⭐
│   │   ├── velocity_grid_rna_umap.png      # Grid velocity
│   │   └── velocity_embedding_*.png        # Embedding plots
│   ├── Phase/
│   │   ├── phase_portraits_top_genes.png   # Top dynamic genes
│   │   └── phase_portraits_markers_*.png   # Marker genes
│   ├── Dynamics/
│   │   ├── velocity_length.png             # Velocity magnitude
│   │   ├── velocity_confidence.png         # Prediction confidence
│   │   └── latent_time.png                 # Latent time
│   └── Validation/
│       ├── velocity_vs_pseudotime.png      # Compare with trajectory
│       └── speed_by_celltype.png           # Speed distribution
└── Tables/
    ├── velocity_scores.csv                 # Per-cell velocity scores
    ├── latent_time.csv                     # RNA velocity time
    ├── velocity_genes.csv                  # Dynamic genes
    └── validation_metrics.csv              # Validation results

Usage:
------
# Basic usage (requires velocyto .loom file)
python 8.rna_velocity.py --mode cDC1 \\
  --loom data/cDC1.loom \\
  --umap results/cDC1/Tables/wnn_umap_coordinates.csv

# With full integration
python 8.rna_velocity.py --mode cDC1 \\
  --loom data/cDC1.loom \\
  --seurat results/cDC1/Robjects/seurat_obj_annotated.rds \\
  --pseudotime monocle3_analysis_cDC1/Tables/pseudotime_by_cell.csv

# Custom parameters
python 8.rna_velocity.py --mode cDC1 \\
  --loom data/cDC1.loom \\
  --umap results/cDC1/Tables/wnn_umap_coordinates.csv \\
  --n-top-genes 3000 \\
  --min-shared-counts 30

Requirements:
------------
- scvelo >= 0.2.5
- scanpy >= 1.9.0
- anndata >= 0.8.0
- loompy >= 3.0.7
- pandas, numpy, matplotlib, seaborn
"""

import os
import sys
import argparse
import warnings
from pathlib import Path
from typing import Optional, Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import scvelo as scv
import anndata as ad

# Suppress warnings
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)

# Set parameters
scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo')
sc.settings.verbosity = 1

# Color palette
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
    """Create organized output directory structure."""
    dirs = {
        'base': base_dir,
        'data': base_dir / 'data',
        'plots': base_dir / 'Plots',
        'velocity': base_dir / 'Plots' / 'Velocity',
        'phase': base_dir / 'Plots' / 'Phase',
        'dynamics': base_dir / 'Plots' / 'Dynamics',
        'validation': base_dir / 'Plots' / 'Validation',
        'tables': base_dir / 'Tables'
    }

    for dir_path in dirs.values():
        dir_path.mkdir(parents=True, exist_ok=True)

    return dirs


def load_velocyto_loom(loom_path: Path) -> ad.AnnData:
    """
    Load velocyto .loom file.

    Parameters:
    ----------
    loom_path : Path
        Path to .loom file

    Returns:
    -------
    ad.AnnData
        AnnData with spliced/unspliced counts
    """
    print(f"\nLoading velocyto data from {loom_path}...")

    if not loom_path.exists():
        raise FileNotFoundError(
            f"Loom file not found: {loom_path}\n\n"
            "Please run velocyto preprocessing first:\n"
            "velocyto run10x -m repeat_msk.gtf sample_folder genes.gtf"
        )

    # Load loom file
    adata = scv.read(loom_path, cache=True)

    print(f"✓ Loaded {adata.n_obs} cells × {adata.n_vars} genes")
    print(f"  Layers: {list(adata.layers.keys())}")

    return adata


def add_umap_coordinates(adata: ad.AnnData,
                        umap_file: Optional[Path] = None,
                        seurat_rds: Optional[Path] = None,
                        umap_type: str = 'wnn') -> ad.AnnData:
    """
    Add UMAP coordinates from external sources.

    Parameters:
    ----------
    adata : ad.AnnData
        AnnData object
    umap_file : Path, optional
        CSV file with UMAP coordinates
    seurat_rds : Path, optional
        Seurat RDS file
    umap_type : str
        UMAP type (rna, adt, wnn)

    Returns:
    -------
    ad.AnnData
        AnnData with UMAP coordinates
    """
    print("\nAdding UMAP coordinates...")

    # Try loading from Seurat RDS
    if seurat_rds and seurat_rds.exists():
        try:
            import rpy2.robjects as ro
            from rpy2.robjects import pandas2ri
            from rpy2.robjects.conversion import localconverter

            pandas2ri.activate()

            print(f"Loading from Seurat RDS: {seurat_rds}")
            ro.r('library(Seurat)')
            seurat_obj = ro.r['readRDS'](str(seurat_rds))
            embeddings = ro.r['Embeddings']

            # Load all available UMAPs
            umap_map = {
                'rna': 'umap',
                'adt': 'adt_umap',
                'wnn': 'wnn_umap'
            }

            for key, reduction_name in umap_map.items():
                try:
                    umap = embeddings(seurat_obj, reduction_name)
                    with localconverter(ro.default_converter + pandas2ri.converter):
                        umap_coords = ro.conversion.rpy2py(umap)
                        adata.obsm[f'X_{key}_umap'] = umap_coords
                        print(f"  ✓ Added {key.upper()} UMAP")
                except:
                    print(f"  ✗ {key.upper()} UMAP not found")

            # Load metadata
            try:
                metadata = ro.r('as.data.frame')(seurat_obj.slots['meta.data'])
                with localconverter(ro.default_converter + pandas2ri.converter):
                    meta_df = ro.conversion.rpy2py(metadata)

                # Match cells by barcode
                common_cells = set(adata.obs_names).intersection(set(meta_df.index))
                if len(common_cells) > 0:
                    adata = adata[list(common_cells)].copy()
                    adata.obs = adata.obs.join(meta_df.loc[common_cells])
                    print(f"  ✓ Matched {len(common_cells)} cells with metadata")

            except Exception as e:
                print(f"  ⚠ Could not load metadata: {e}")

        except ImportError:
            print("  ✗ rpy2 not available, falling back to CSV")

    # Load from CSV file
    if umap_file and umap_file.exists():
        print(f"Loading UMAP from CSV: {umap_file}")
        umap_df = pd.read_csv(umap_file, index_col=0)

        # Detect UMAP column names
        umap_cols = [col for col in umap_df.columns if 'UMAP' in col or 'umap' in col]

        if len(umap_cols) >= 2:
            # Match cell barcodes
            common_cells = set(adata.obs_names).intersection(set(umap_df.index))

            if len(common_cells) > 0:
                adata = adata[list(common_cells)].copy()
                adata.obsm[f'X_{umap_type}_umap'] = umap_df.loc[common_cells, umap_cols[:2]].values
                print(f"  ✓ Matched {len(common_cells)} cells")
                print(f"  ✓ Added {umap_type.upper()} UMAP: {umap_cols[:2]}")

                # Add cell type if available
                if 'cell_type' in umap_df.columns:
                    adata.obs['cell_type'] = umap_df.loc[common_cells, 'cell_type']
                    print(f"  ✓ Added cell type annotations")

            else:
                print("  ✗ No matching cells found between loom and UMAP file")
        else:
            print(f"  ✗ Could not find UMAP columns in {umap_file}")

    # Check if we have at least one UMAP
    umap_found = any(key in adata.obsm for key in ['X_rna_umap', 'X_adt_umap', 'X_wnn_umap'])

    if not umap_found:
        print("\n⚠ WARNING: No UMAP coordinates found!")
        print("Velocity visualization will be limited.")

    return adata


def preprocess_velocity_data(adata: ad.AnnData,
                            n_top_genes: int = 2000,
                            min_shared_counts: int = 30,
                            n_pcs: int = 30,
                            n_neighbors: int = 30) -> ad.AnnData:
    """
    Preprocess data for velocity analysis.

    Parameters:
    ----------
    adata : ad.AnnData
        Input AnnData
    n_top_genes : int
        Number of top variable genes
    min_shared_counts : int
        Minimum shared counts for velocity genes
    n_pcs : int
        Number of PCs
    n_neighbors : int
        Number of neighbors

    Returns:
    -------
    ad.AnnData
        Preprocessed AnnData
    """
    print("\n" + "="*80)
    print("PREPROCESSING FOR VELOCITY ANALYSIS")
    print("="*80)

    print(f"Input: {adata.n_obs} cells × {adata.n_vars} genes")

    # Basic filtering
    print("\n1. Filtering genes and cells...")
    scv.pp.filter_and_normalize(adata,
                                min_shared_counts=min_shared_counts,
                                n_top_genes=n_top_genes)

    print(f"After filtering: {adata.n_obs} cells × {adata.n_vars} genes")

    # Find variable genes if not already computed
    if 'highly_variable' not in adata.var.columns:
        print("\n2. Identifying highly variable genes...")
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
        n_hvg = adata.var['highly_variable'].sum()
        print(f"Found {n_hvg} highly variable genes")

    # Compute moments
    print("\n3. Computing moments for velocity estimation...")

    # Check if we need to compute PCA and neighbors
    if 'X_pca' not in adata.obsm:
        print("  Computing PCA...")
        sc.tl.pca(adata, n_comps=n_pcs)

    if 'neighbors' not in adata.uns:
        print("  Computing neighbor graph...")
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    # Compute first and second order moments
    scv.pp.moments(adata, n_neighbors=n_neighbors)

    print("✓ Preprocessing complete")

    return adata


def compute_velocity(adata: ad.AnnData,
                    mode: str = 'stochastic') -> ad.AnnData:
    """
    Compute RNA velocity.

    Parameters:
    ----------
    adata : ad.AnnData
        Preprocessed AnnData
    mode : str
        Velocity mode ('deterministic', 'stochastic', or 'dynamical')

    Returns:
    -------
    ad.AnnData
        AnnData with velocity
    """
    print("\n" + "="*80)
    print(f"COMPUTING RNA VELOCITY ({mode.upper()} MODE)")
    print("="*80)

    if mode == 'deterministic':
        print("Using deterministic model (fast, assumes steady-state)...")
        scv.tl.velocity(adata, mode='deterministic')

    elif mode == 'stochastic':
        print("Using stochastic model (balanced speed and accuracy)...")
        scv.tl.velocity(adata, mode='stochastic')

    elif mode == 'dynamical':
        print("Using dynamical model (most accurate, slower)...")
        scv.tl.recover_dynamics(adata)
        scv.tl.velocity(adata, mode='dynamical')

    else:
        raise ValueError(f"Unknown mode: {mode}")

    # Compute velocity graph
    print("\nComputing velocity graph...")
    scv.tl.velocity_graph(adata)

    # Compute additional metrics
    print("Computing velocity metrics...")

    # Velocity length (speed)
    scv.tl.velocity_length(adata)

    # Velocity confidence
    scv.tl.velocity_confidence(adata)

    # Latent time (RNA velocity pseudotime)
    try:
        scv.tl.latent_time(adata)
        print("✓ Latent time computed")
    except:
        print("⚠ Could not compute latent time")

    print("✓ Velocity computation complete")

    return adata


def plot_velocity_embeddings(adata: ad.AnnData,
                            output_dir: Path,
                            basis: str = 'wnn_umap',
                            title_prefix: str = "cDC") -> None:
    """
    Create velocity stream plots on UMAP embeddings.

    Parameters:
    ----------
    adata : ad.AnnData
        AnnData with velocity
    output_dir : Path
        Output directory
    basis : str
        Embedding basis (e.g., 'wnn_umap')
    title_prefix : str
        Prefix for titles
    """
    print("\n" + "="*80)
    print("CREATING VELOCITY VISUALIZATIONS")
    print("="*80)

    # Determine available UMAP spaces
    available_bases = []
    for umap_type in ['rna_umap', 'adt_umap', 'wnn_umap']:
        if f'X_{umap_type}' in adata.obsm:
            available_bases.append(umap_type)

    if not available_bases:
        print("⚠ No UMAP coordinates found, skipping embedding plots")
        return

    print(f"Available UMAP spaces: {', '.join(available_bases)}")

    # Plot velocity stream for each UMAP space
    for umap_basis in available_bases:
        umap_name = umap_basis.replace('_', ' ').upper()
        print(f"\nPlotting velocity on {umap_name}...")

        # Stream plot
        fig, ax = plt.subplots(figsize=(10, 8))
        scv.pl.velocity_embedding_stream(
            adata,
            basis=umap_basis,
            color='cell_type' if 'cell_type' in adata.obs else None,
            title=f'{title_prefix} Velocity - {umap_name}',
            size=50,
            alpha=0.8,
            arrow_size=2,
            arrow_length=3,
            legend_loc='right margin',
            show=False,
            ax=ax
        )
        plt.tight_layout()
        output_file = output_dir / f'velocity_stream_{umap_basis}.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Saved: {output_file.name}")

        # Grid plot
        fig, ax = plt.subplots(figsize=(10, 8))
        scv.pl.velocity_embedding_grid(
            adata,
            basis=umap_basis,
            color='cell_type' if 'cell_type' in adata.obs else None,
            title=f'{title_prefix} Velocity Grid - {umap_name}',
            size=50,
            alpha=0.8,
            arrow_size=2,
            arrow_length=3,
            show=False,
            ax=ax
        )
        plt.tight_layout()
        output_file = output_dir / f'velocity_grid_{umap_basis}.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Saved: {output_file.name}")

        # Velocity embedding (colored by velocity)
        fig, ax = plt.subplots(figsize=(10, 8))
        scv.pl.velocity_embedding(
            adata,
            basis=umap_basis,
            arrow_length=3,
            arrow_size=2,
            title=f'{title_prefix} Velocity Magnitude - {umap_name}',
            show=False,
            ax=ax
        )
        plt.tight_layout()
        output_file = output_dir / f'velocity_embedding_{umap_basis}.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Saved: {output_file.name}")


def plot_phase_portraits(adata: ad.AnnData,
                        output_dir: Path,
                        marker_genes: Optional[List[str]] = None,
                        n_top_genes: int = 20) -> None:
    """
    Create phase portraits for dynamic genes.

    Parameters:
    ----------
    adata : ad.AnnData
        AnnData with velocity
    output_dir : Path
        Output directory
    marker_genes : List[str], optional
        Specific marker genes to plot
    n_top_genes : int
        Number of top dynamic genes to plot
    """
    print("\n" + "="*80)
    print("CREATING PHASE PORTRAITS")
    print("="*80)

    # Identify top likelihood genes
    scv.tl.rank_velocity_genes(adata, groupby='cell_type' if 'cell_type' in adata.obs else None)

    # Plot top dynamic genes
    print(f"Plotting top {n_top_genes} dynamic genes...")

    top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).head(n_top_genes).index

    fig = scv.pl.scatter(
        adata,
        basis=top_genes[:12],  # Show top 12 in grid
        ncols=4,
        frameon=False,
        title='Top Dynamic Genes - Phase Portraits',
        show=False,
        figsize=(16, 12)
    )

    output_file = output_dir / 'phase_portraits_top_genes.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved: {output_file.name}")

    # Plot marker genes if provided
    if marker_genes:
        available_markers = [g for g in marker_genes if g in adata.var_names]

        if available_markers:
            print(f"Plotting {len(available_markers)} marker genes...")

            fig = scv.pl.scatter(
                adata,
                basis=available_markers[:12],
                ncols=4,
                frameon=False,
                title='Marker Genes - Phase Portraits',
                show=False,
                figsize=(16, 12)
            )

            output_file = output_dir / 'phase_portraits_markers.png'
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"  ✓ Saved: {output_file.name}")


def plot_velocity_dynamics(adata: ad.AnnData,
                          output_dir: Path,
                          basis: str = 'wnn_umap') -> None:
    """
    Plot velocity dynamics metrics.

    Parameters:
    ----------
    adata : ad.AnnData
        AnnData with velocity
    output_dir : Path
        Output directory
    basis : str
        UMAP basis for plotting
    """
    print("\n" + "="*80)
    print("CREATING DYNAMICS VISUALIZATIONS")
    print("="*80)

    # Check available basis
    if f'X_{basis}' not in adata.obsm:
        # Try to find any available UMAP
        for alt_basis in ['wnn_umap', 'rna_umap', 'adt_umap']:
            if f'X_{alt_basis}' in adata.obsm:
                basis = alt_basis
                break

    # Velocity length (speed)
    if 'velocity_length' in adata.obs:
        print("Plotting velocity magnitude...")

        fig, ax = plt.subplots(figsize=(10, 8))
        scv.pl.scatter(
            adata,
            basis=basis,
            c='velocity_length',
            cmap='coolwarm',
            title='RNA Velocity Magnitude',
            show=False,
            ax=ax
        )
        plt.tight_layout()
        output_file = output_dir / 'velocity_length.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Saved: {output_file.name}")

    # Velocity confidence
    if 'velocity_confidence' in adata.obs:
        print("Plotting velocity confidence...")

        fig, ax = plt.subplots(figsize=(10, 8))
        scv.pl.scatter(
            adata,
            basis=basis,
            c='velocity_confidence',
            cmap='coolwarm',
            title='RNA Velocity Confidence',
            show=False,
            ax=ax
        )
        plt.tight_layout()
        output_file = output_dir / 'velocity_confidence.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Saved: {output_file.name}")

    # Latent time
    if 'latent_time' in adata.obs:
        print("Plotting latent time...")

        fig, ax = plt.subplots(figsize=(10, 8))
        scv.pl.scatter(
            adata,
            basis=basis,
            c='latent_time',
            cmap='viridis',
            title='RNA Velocity Latent Time',
            show=False,
            ax=ax
        )
        plt.tight_layout()
        output_file = output_dir / 'latent_time.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Saved: {output_file.name}")


def validate_with_pseudotime(adata: ad.AnnData,
                            output_dir: Path,
                            pseudotime_file: Optional[Path] = None) -> None:
    """
    Validate velocity with trajectory pseudotime.

    Parameters:
    ----------
    adata : ad.AnnData
        AnnData with velocity
    output_dir : Path
        Output directory
    pseudotime_file : Path, optional
        CSV file with pseudotime values
    """
    if pseudotime_file is None or not pseudotime_file.exists():
        print("\n⚠ No pseudotime file provided, skipping validation")
        return

    print("\n" + "="*80)
    print("VALIDATING VELOCITY WITH PSEUDOTIME")
    print("="*80)

    # Load pseudotime
    pt_df = pd.read_csv(pseudotime_file, index_col=0)

    # Match cells
    common_cells = set(adata.obs_names).intersection(set(pt_df.index))

    if len(common_cells) < 10:
        print(f"⚠ Too few matching cells ({len(common_cells)}), skipping validation")
        return

    print(f"Matched {len(common_cells)} cells")

    # Add pseudotime to adata
    adata.obs['trajectory_pseudotime'] = pt_df.loc[common_cells].iloc[:, 0]

    # Compare with latent time
    if 'latent_time' in adata.obs:
        from scipy.stats import pearsonr

        valid_mask = ~(adata.obs['latent_time'].isna() | adata.obs['trajectory_pseudotime'].isna())
        valid_cells = adata.obs_names[valid_mask]

        if len(valid_cells) > 10:
            corr, pval = pearsonr(
                adata.obs.loc[valid_cells, 'latent_time'],
                adata.obs.loc[valid_cells, 'trajectory_pseudotime']
            )

            print(f"\nCorrelation between latent time and pseudotime:")
            print(f"  r = {corr:.3f}, p = {pval:.2e}")

            # Plot comparison
            fig, ax = plt.subplots(figsize=(8, 6))
            sns.scatterplot(
                data=adata.obs.loc[valid_cells],
                x='latent_time',
                y='trajectory_pseudotime',
                hue='cell_type' if 'cell_type' in adata.obs else None,
                alpha=0.6,
                ax=ax
            )
            ax.set_xlabel('RNA Velocity Latent Time')
            ax.set_ylabel('Trajectory Pseudotime')
            ax.set_title(f'Velocity vs Pseudotime (r = {corr:.3f})')
            plt.tight_layout()

            output_file = output_dir / 'velocity_vs_pseudotime.png'
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"  ✓ Saved: {output_file.name}")

    # Plot speed by cell type
    if 'velocity_length' in adata.obs and 'cell_type' in adata.obs:
        print("\nPlotting velocity speed by cell type...")

        fig, ax = plt.subplots(figsize=(10, 6))
        sns.violinplot(
            data=adata.obs,
            x='cell_type',
            y='velocity_length',
            ax=ax
        )
        ax.set_xlabel('Cell Type')
        ax.set_ylabel('Velocity Magnitude')
        ax.set_title('RNA Velocity Speed by Cell Type')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()

        output_file = output_dir / 'speed_by_celltype.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Saved: {output_file.name}")


def export_results(adata: ad.AnnData,
                  data_dir: Path,
                  tables_dir: Path) -> None:
    """
    Export velocity results.

    Parameters:
    ----------
    adata : ad.AnnData
        AnnData with velocity
    data_dir : Path
        Data output directory
    tables_dir : Path
        Tables output directory
    """
    print("\n" + "="*80)
    print("EXPORTING RESULTS")
    print("="*80)

    # Save AnnData
    output_file = data_dir / 'adata_with_velocity.h5ad'
    adata.write_h5ad(output_file)
    print(f"✓ Saved AnnData: {output_file}")

    # Export velocity scores
    if 'velocity_length' in adata.obs or 'latent_time' in adata.obs:
        score_cols = ['cell_type'] if 'cell_type' in adata.obs else []
        if 'velocity_length' in adata.obs:
            score_cols.append('velocity_length')
        if 'velocity_confidence' in adata.obs:
            score_cols.append('velocity_confidence')
        if 'latent_time' in adata.obs:
            score_cols.append('latent_time')

        scores_df = adata.obs[score_cols].copy()
        scores_df.index.name = 'cell_id'

        output_file = tables_dir / 'velocity_scores.csv'
        scores_df.to_csv(output_file)
        print(f"✓ Saved velocity scores: {output_file}")

    # Export velocity genes
    if 'fit_likelihood' in adata.var:
        genes_df = adata.var[['fit_likelihood']].copy()
        genes_df = genes_df.sort_values('fit_likelihood', ascending=False)
        genes_df.index.name = 'gene'

        output_file = tables_dir / 'velocity_genes.csv'
        genes_df.to_csv(output_file)
        print(f"✓ Saved velocity genes: {output_file}")

        # Save top genes list
        top_genes = genes_df.head(100).index.tolist()
        output_file = data_dir / 'velocity_genes.txt'
        with open(output_file, 'w') as f:
            f.write('\n'.join(top_genes))
        print(f"✓ Saved top 100 velocity genes: {output_file}")

    print("\n" + "="*80)
    print("RNA VELOCITY ANALYSIS COMPLETE!")
    print("="*80)


# =============================================================================
# Main Function
# =============================================================================

def main():
    """Main execution function."""

    parser = argparse.ArgumentParser(
        description='RNA Velocity Analysis for cDC CITE-seq Data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument('--mode',
                       type=str,
                       default='cDC1',
                       choices=['cDC1', 'cDC2', 'integrated'],
                       help='Analysis mode (default: cDC1)')

    parser.add_argument('--loom',
                       type=str,
                       required=True,
                       help='Path to velocyto .loom file (required)')

    parser.add_argument('--umap',
                       type=str,
                       default=None,
                       help='CSV file with UMAP coordinates')

    parser.add_argument('--seurat',
                       type=str,
                       default=None,
                       help='Seurat RDS file (alternative to --umap)')

    parser.add_argument('--pseudotime',
                       type=str,
                       default=None,
                       help='CSV file with pseudotime for validation')

    parser.add_argument('--output-dir',
                       type=str,
                       default=None,
                       help='Output directory (default: rna_velocity_{mode})')

    parser.add_argument('--velocity-mode',
                       type=str,
                       default='stochastic',
                       choices=['deterministic', 'stochastic', 'dynamical'],
                       help='Velocity computation mode (default: stochastic)')

    parser.add_argument('--n-top-genes',
                       type=int,
                       default=2000,
                       help='Number of top variable genes (default: 2000)')

    parser.add_argument('--min-shared-counts',
                       type=int,
                       default=30,
                       help='Minimum shared counts (default: 30)')

    args = parser.parse_args()

    # Print header
    print("\n" + "="*80)
    print(" "*25 + "RNA VELOCITY ANALYSIS")
    print(" "*25 + "GSE228544 CITE-seq Data")
    print("="*80)
    print(f"\nMode: {args.mode}")
    print(f"Velocity mode: {args.velocity_mode}")
    print(f"Loom file: {args.loom}")

    # Set up output directory
    if args.output_dir:
        output_base = Path(args.output_dir)
    else:
        output_base = Path(f'rna_velocity_{args.mode}')

    dirs = create_output_directories(output_base)
    print(f"Output directory: {output_base}")

    # Load velocyto data
    adata = load_velocyto_loom(Path(args.loom))

    # Add UMAP coordinates
    adata = add_umap_coordinates(
        adata,
        umap_file=Path(args.umap) if args.umap else None,
        seurat_rds=Path(args.seurat) if args.seurat else None
    )

    # Preprocess
    adata = preprocess_velocity_data(
        adata,
        n_top_genes=args.n_top_genes,
        min_shared_counts=args.min_shared_counts
    )

    # Compute velocity
    adata = compute_velocity(adata, mode=args.velocity_mode)

    # Create visualizations
    plot_velocity_embeddings(adata, dirs['velocity'], title_prefix=args.mode)
    plot_phase_portraits(adata, dirs['phase'])
    plot_velocity_dynamics(adata, dirs['dynamics'])

    # Validation
    if args.pseudotime:
        validate_with_pseudotime(
            adata,
            dirs['validation'],
            Path(args.pseudotime)
        )

    # Export results
    export_results(adata, dirs['data'], dirs['tables'])

    print(f"\n✓ All results saved to: {output_base}")
    print("\nNext steps:")
    print("  1. Review velocity stream plots in Plots/Velocity/")
    print("  2. Check phase portraits for dynamic genes")
    print("  3. Validate with Monocle3/Slingshot pseudotime")
    print("  4. Identify velocity-driven regulatory genes")


if __name__ == '__main__':
    main()
