toy_basis <- function(order = NULL) {
  set.seed(1)
  w <- matrix(rnorm(120), nrow = 40, ncol = 3,
              dimnames = list(as.character(seq_len(40)),
                              paste0("P", 1:3)))
  if (!is.null(order)) w <- w[order, , drop = FALSE]
  w <- validate_program_basis(w, min_genes = 20)
  structure(list(weights = w, basis_type = "signed_template",
                 audit = attr(w, "proret_audit")),
            class = "proret_basis")
}

toy_drugs <- function(basis, disease_gene) {
  x <- cbind(
    same_1 = disease_gene,
    same_2 = disease_gene * 0.95,
    reverse_1 = -disease_gene,
    reverse_2 = -disease_gene * 0.95,
    other_1 = rev(disease_gene),
    other_2 = rev(disease_gene)
  )
  rownames(x) <- names(disease_gene)
  md <- data.frame(
    sig_id = colnames(x),
    pert_id = rep(c("A", "B", "C"), each = 2),
    pert_iname = rep(c("same", "reverse", "other"), each = 2),
    cell_iname = rep(c("CELL1", "CELL2"), 3),
    pert_type = "trt_cp", pert_itime = "24 h", pert_idose = "1 uM",
    is_hiq = 1
  )
  project_lincs_matrix(x, md, basis, gene_space_ids = rownames(x))
}
