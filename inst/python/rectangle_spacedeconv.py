import rectanglepy as rectangle


def py_deconvolute_rectangle(
    sc_obj,
    bulks,
    cell_type_col="cell_type",
    layer=None,
    raw=False,
    correct_mrna_bias=True,
    optimize_cutoffs=True,
    p=0.015,
    lfc=1.5,
    n_cpus=None,
    gene_expression_threshold=0.5,
):
    estimations, signature_result = rectangle.rectangle(
        adata=sc_obj,
        bulks=bulks,
        cell_type_col=cell_type_col,
        layer=layer,
        raw=raw,
        correct_mrna_bias=correct_mrna_bias,
        optimize_cutoffs=optimize_cutoffs,
        p=p,
        lfc=lfc,
        n_cpus=n_cpus,
        gene_expression_threshold=gene_expression_threshold,
    )

    return estimations
