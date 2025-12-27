"""
FlashDeconv wrapper for spacedeconv integration.

FlashDeconv is a high-performance spatial transcriptomics deconvolution method
designed for atlas-scale data. It uses structure-preserving randomized sketching
and spatial graph regularization for efficient deconvolution.

GitHub: https://github.com/cafferychen777/flashdeconv
Paper: https://doi.org/10.64898/2025.12.22.696108
"""

import numpy as np
import pandas as pd


def py_deconvolute_flashdeconv(
    sp_obj,
    ref_obj,
    cell_type_col="cell_type",
    sketch_dim=512,
    lambda_spatial=5000.0,
    n_hvg=2000,
    n_markers_per_type=50,
    layer_st=None,
    layer_ref=None,
    spatial_key="spatial",
):
    """
    Deconvolute spatial transcriptomics data using FlashDeconv.

    Parameters
    ----------
    sp_obj : AnnData
        Spatial transcriptomics AnnData object with coordinates in obsm[spatial_key].
    ref_obj : AnnData
        Single-cell reference AnnData object with cell type labels in obs[cell_type_col].
    cell_type_col : str, default="cell_type"
        Column name in ref_obj.obs containing cell type annotations.
    sketch_dim : int, default=512
        Dimension of the sketched space. Higher values preserve more information.
    lambda_spatial : float, default=5000.0
        Spatial regularization strength. Higher values encourage spatially smooth
        cell type distributions. Typical range: 1000-10000 for Visium.
    n_hvg : int, default=2000
        Number of highly variable genes to select.
    n_markers_per_type : int, default=50
        Number of marker genes per cell type.
    layer_st : str, optional
        Layer in sp_obj to use for counts. Uses .X if None.
    layer_ref : str, optional
        Layer in ref_obj to use for counts. Uses .X if None.
    spatial_key : str, default="spatial"
        Key in sp_obj.obsm for spatial coordinates.

    Returns
    -------
    pd.DataFrame
        Cell type proportions with spots as rows and cell types as columns.
        Values sum to 1 for each spot.
    """
    try:
        import flashdeconv as fd
    except ImportError:
        raise ImportError(
            "FlashDeconv is not installed. Please install it with:\n"
            "  pip install flashdeconv\n"
            "For more information, see: https://github.com/cafferychen777/flashdeconv"
        )

    # Make a copy to avoid modifying the original object
    adata_st = sp_obj.copy()

    # Run FlashDeconv deconvolution
    fd.tl.deconvolve(
        adata_st,
        ref_obj,
        cell_type_key=cell_type_col,
        sketch_dim=sketch_dim,
        lambda_spatial=lambda_spatial,
        n_hvg=n_hvg,
        n_markers_per_type=n_markers_per_type,
        layer_st=layer_st,
        layer_ref=layer_ref,
        spatial_key=spatial_key,
        key_added="flashdeconv",
        copy=False,
    )

    # Extract results as DataFrame
    # flashdeconv stores results in obsm['flashdeconv'] as a DataFrame
    proportions = adata_st.obsm["flashdeconv"]

    # Ensure it's a pandas DataFrame with proper index
    if not isinstance(proportions, pd.DataFrame):
        proportions = pd.DataFrame(
            proportions,
            index=adata_st.obs_names,
            columns=adata_st.uns["flashdeconv_params"]["cell_type_names"],
        )

    return proportions
