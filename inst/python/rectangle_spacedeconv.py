import rectanglepy as rectangle


def py_build_rectangle_signatures(
    sc_obj,
    bulks=None,
    cell_type_col="cell_type",
    layer=None,
    raw=False,
    optimize_cutoffs=True,
    p=0.015,
    lfc=1.5,
    n_cpus=None,
    gene_expression_threshold=0.5,
):
    return rectangle.pp.build_rectangle_signatures(
        adata=sc_obj,
        bulks=bulks,
        cell_type_col=cell_type_col,
        layer=layer,
        raw=raw,
        optimize_cutoffs=optimize_cutoffs,
        p=p,
        lfc=lfc,
        n_cpus=n_cpus,
        gene_expression_threshold=gene_expression_threshold,
    )


def py_deconvolute_rectangle(
    signature_result,
    bulks,
    correct_mrna_bias=True,
    n_cpus=None,
):
    estimations, bulk_err_df = rectangle.tl.deconvolution(
        signatures=signature_result,
        bulks=bulks,
        correct_mrna_bias=correct_mrna_bias,
        n_cpus=n_cpus,
    )

    return estimations
